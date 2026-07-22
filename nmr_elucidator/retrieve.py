"""候选检索层: 在 nmrshiftdb2 SDF 中找与实测谱兼容的结构.

匹配策略：
- 硬筛：分子式、溶剂、重原子数
- 软筛：化学位移窗口匹配（每个实测峰在库内对应谱上找最近峰）

性能：
- 默认流式扫描 58187 记录; 内存友好
- 命中后即返回（不返回前 N 个），由 score 层排序

调用模式：
    for entry in stream_nmrshiftdb2_sdf(sdf_path):
        hit = retrieve_one(entry, problem, constraints)
        if hit: yield hit
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, List, Optional

from .constraints import Constraints, formula_match, solvent_compatible
from .io_sdf import NmrShiftDb2Entry


# 默认容差（13C 化学位移在 100 MHz 下典型 0.5 ppm, 高场下 1-2 ppm）
DEFAULT_13C_TOL = 2.0
DEFAULT_1H_TOL = 0.3


@dataclass
class ProblemPeaks:
    """观测峰集合（从 BrukerSpectrum 提取）."""

    carbon_shifts: List[float] = field(default_factory=list)
    proton_shifts: List[float] = field(default_factory=list)
    proton_integrals: List[tuple] = field(default_factory=list)  # [(shift, integral), ...]
    solvent: Optional[str] = None
    field_mhz: Optional[float] = None

    @property
    def has_13c(self) -> bool:
        return bool(self.carbon_shifts)

    @property
    def has_1h(self) -> bool:
        return bool(self.proton_shifts)


@dataclass
class RetrievalHit:
    """单条候选的命中信息."""

    entry: NmrShiftDb2Entry
    # 13C 匹配
    c_match_count: int = 0
    c_total_obs: int = 0
    c_total_lib: int = 0
    c_coverage: float = 0.0  # matched / max(obs, lib)
    c_rms: float = 0.0  # 已匹配对的位移差均方根
    # 1H 匹配
    h_match_count: int = 0
    h_total_obs: int = 0
    h_total_lib: int = 0
    h_coverage: float = 0.0
    h_rms: float = 0.0
    # 硬约束
    formula_ok: bool = True
    solvent_ok: bool = True
    heavy_atoms: int = 0
    # DEPT 匹配 (供 score 维度使用)
    dept_total: int = 0  # 候选预测 C 原子总数
    dept_coverage: float = 0.0  # DEPT 类型 full match 比例
    # H 总数 (供 score 维度使用)
    h_count_expected: Optional[int] = None  # DEPT 推断的 H 总数
    h_count_lib: Optional[int] = None  # 候选分子 H 总数


# ---- 单条匹配 ----

def _shifts_lib(entry: NmrShiftDb2Entry, nucleus: str) -> List[float]:
    if nucleus == '13C':
        return [s.shift for s in entry.spectrum_13c]
    if nucleus == '1H':
        return [s.shift for s in entry.spectrum_1h]
    return []


def _greedy_match(obs: List[float], lib: List[float], tol: float) -> tuple:
    """贪心匹配: 对每个 obs 峰找最近的 lib 峰, 一对一, 距离 <= tol 才计入.

    返回 (matched_count, rms_of_matched_pairs).
    """
    if not obs or not lib:
        return 0, 0.0
    lib_sorted = sorted(lib)
    matched = 0
    sqsum = 0.0
    used = [False] * len(lib_sorted)
    for o in obs:
        # 二分最近
        lo, hi = 0, len(lib_sorted) - 1
        # 找最近未用
        best_i = -1
        best_d = None
        # 简单线性扫描（库通常 <= 50 条, 线性即可）
        for i, l in enumerate(lib_sorted):
            if used[i]:
                continue
            d = abs(o - l)
            if d <= tol and (best_d is None or d < best_d):
                best_d = d
                best_i = i
        if best_i >= 0:
            used[best_i] = True
            matched += 1
            sqsum += (best_d or 0.0) ** 2
    rms = (sqsum / matched) ** 0.5 if matched else 0.0
    return matched, rms


def retrieve_one(
    entry: NmrShiftDb2Entry,
    problem: ProblemPeaks,
    constraints: Constraints,
    c_tol: float = DEFAULT_13C_TOL,
    h_tol: float = DEFAULT_1H_TOL,
    obs_dept_summary: Optional[List[dict]] = None,
) -> Optional[RetrievalHit]:
    """对单条候选做匹配判定. 通过返回 hit, 不通过返回 None.

    obs_dept_summary: 实测 DEPT 推断 [{ppm, mult, sign}, ...] (可选).
                     提供后会同时算 dept_coverage 和 h_count 字段,
                     供 score_candidates 综合打分使用.
    """
    # 硬约束：分子式
    if constraints.formula:
        try:
            mol = entry.mol()
        except Exception:
            mol = None
        if mol is None:
            return None
        if not formula_match(constraints.formula, mol):
            return None
        formula_ok = True
    else:
        formula_ok = True
        try:
            mol = entry.mol()
        except Exception:
            mol = None

    # 硬约束：溶剂
    if constraints.solvent:
        if not solvent_compatible(constraints.solvent, entry.solvent):
            return None
    solvent_ok = True

    # 重原子数
    if mol is None:
        heavy = entry.atom_count()
    else:
        heavy = mol.GetNumHeavyAtoms()
    if constraints.heavy_atom_min is not None and heavy < constraints.heavy_atom_min:
        return None
    if constraints.heavy_atom_max is not None and heavy > constraints.heavy_atom_max:
        return None

    # 13C 匹配
    c_match, c_rms = _greedy_match(problem.carbon_shifts, _shifts_lib(entry, '13C'), c_tol)
    # 1H 匹配
    h_match, h_rms = _greedy_match(problem.proton_shifts, _shifts_lib(entry, '1H'), h_tol)

    # 至少需要 13C 或 1H 中任一谱有匹配
    obs_has_c = bool(problem.carbon_shifts)
    obs_has_h = bool(problem.proton_shifts)
    if obs_has_c:
        if c_match / len(problem.carbon_shifts) < 0.15:
            return None
    if obs_has_h:
        if h_match / len(problem.proton_shifts) < 0.15:
            return None
    if (c_match + h_match) < 2:
        return None

    c_total = len(problem.carbon_shifts)
    h_total = len(problem.proton_shifts)
    c_lib = len(entry.spectrum_13c)
    h_lib = len(entry.spectrum_1h)
    c_cov = c_match / max(c_total, c_lib) if max(c_total, c_lib) else 0.0
    h_cov = h_match / max(h_total, h_lib) if max(h_total, h_lib) else 0.0

    # DEPT 匹配 + H 总数 (可选, 用于 score 维度)
    dept_total = 0
    dept_coverage = 0.0
    h_count_expected: Optional[int] = None
    h_count_lib: Optional[int] = None
    if obs_dept_summary and mol is not None and c_lib > 0:
        from .evidence import build_dept_evidence, compute_h_count_evidence
        d_evid = build_dept_evidence(
            [s.shift for s in entry.spectrum_13c],
            problem.carbon_shifts, obs_dept_summary, mol, c_tol,
        )
        dept_total = len(d_evid)
        if dept_total:
            full_match = sum(1 for d in d_evid if d.match_type == 'full')
            dept_coverage = full_match / dept_total
        h_ce = compute_h_count_evidence(obs_dept_summary, mol)
        if h_ce is not None:
            h_count_expected = h_ce.expected
            h_count_lib = h_ce.lib_heavy_h

    return RetrievalHit(
        entry=entry,
        c_match_count=c_match,
        c_total_obs=c_total,
        c_total_lib=c_lib,
        c_coverage=c_cov,
        c_rms=c_rms,
        h_match_count=h_match,
        h_total_obs=h_total,
        h_total_lib=h_lib,
        h_coverage=h_cov,
        h_rms=h_rms,
        formula_ok=formula_ok,
        solvent_ok=solvent_ok,
        heavy_atoms=heavy,
        dept_total=dept_total,
        dept_coverage=dept_coverage,
        h_count_expected=h_count_expected,
        h_count_lib=h_count_lib,
    )


# ---- 顶层接口 ----

def retrieve_candidates(
    sdf_path: str,
    problem: ProblemPeaks,
    constraints: Constraints,
    c_tol: float = DEFAULT_13C_TOL,
    h_tol: float = DEFAULT_1H_TOL,
    max_hits: Optional[int] = None,
    progress_every: int = 5000,
    obs_dept_summary: Optional[List[dict]] = None,
) -> Iterator[RetrievalHit]:
    """流式检索 nmrshiftdb2.

    max_hits 限制返回数（达到后停止扫描）.
    progress_every 每 N 条打印一次进度（仅输出到 stderr, 不污染 stdout）.
    obs_dept_summary: 实测 DEPT 推断 (可选, 提供后用于打分维度).
    """
    from .io_sdf import stream_nmrshiftdb2_sdf
    import sys
    seen = 0
    returned = 0
    for entry in stream_nmrshiftdb2_sdf(sdf_path):
        seen += 1
        hit = retrieve_one(entry, problem, constraints, c_tol=c_tol, h_tol=h_tol,
                           obs_dept_summary=obs_dept_summary)
        if hit is not None:
            yield hit
            returned += 1
            if max_hits is not None and returned >= max_hits:
                return
        if progress_every and seen % progress_every == 0:
            print(f"[retrieve] scanned={seen} hits={returned}", file=sys.stderr)