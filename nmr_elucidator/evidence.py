"""打分与证据链.

score_candidates:
    - 综合打分: 覆盖率权重 + RMS 惩罚 + 峰数差异 + DEPT 匹配 + H 总数一致性
    - 排序返回 top-N

build_evidence:
    - 每条候选: matched/unmatched/conflict 证据
    - DEPT 推断表: 候选预测谱的每个 C 位置与实测 DEPT 类型对照
    - H 总数一致性: 候选 H 总数 vs DEPT 推断的 H 总数
    - 重点: 未匹配峰列表 (留待 ML 模型 + 化学家人工解析)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .constraints import formula_match
from .correlation_2d import Correlation2DScore
from .io_sdf import NmrShiftDb2Entry
from .retrieve import RetrievalHit


# ---------- score ----------

@dataclass
class ScoredCandidate:
    """打分后的候选."""

    hit: RetrievalHit
    score: float
    rank: int = 0


def score_candidates(
    hits: List[RetrievalHit],
    c_weight: float = 0.6,
    h_weight: float = 0.4,
    dept_match_bonus: float = 30.0,
    h_count_match_bonus: float = 20.0,
    use_2d_scoring: bool = False,
    use_noe_scoring: bool = False,
    noe_scores: dict = None,
) -> List[ScoredCandidate]:
    """综合打分.

    无 2D 数据时 (use_2d_scoring=False):
    score = (c_coverage * c_weight + h_coverage * h_weight) * 100
            - 5 * c_rms
            - 5 * h_rms
            - 2 * peak_count_diff
            + dept_match_bonus  * dept_coverage
            + h_count_match_bonus * h_count_match_score

    有 2D 数据时 (use_2d_scoring=True):
    score = 1D_coverage * 40
            + hsqc_pair_score * 30
            - 5 * hsqc_rms
            + dept_bonus * 5
            + h_count_bonus * 5
            - 2 * peak_count_diff

    有 NOE 数据时 (use_noe_scoring=True):
    score += noe_match_rate * 10 - 5 * noe_conflict_count
    """
    scored: List[ScoredCandidate] = []
    for h in hits:
        if use_2d_scoring and h.hsqc_total_obs > 0:
            # 2D 模式: HSQC 配对匹配为核心
            cov_1d = (h.c_coverage * c_weight + h.h_coverage * h_weight) * 40.0
            hsqc_score = h.hsqc_coverage * 30.0
            hsqc_pen = 5.0 * h.hsqc_rms
            diff_c = abs(h.c_total_lib - h.c_total_obs)
            diff_h = abs(h.h_total_lib - h.h_total_obs)
            diff_pen = 2.0 * (diff_c + diff_h)
            # DEPT + H 总数 (权重降低)
            dept_bonus = 0.0
            if h.dept_total > 0:
                dept_bonus = 5.0 * h.dept_coverage
            h_count_bonus = 0.0
            if h.h_count_expected is not None and h.h_count_lib is not None and h.h_count_expected > 0:
                ratio = min(h.h_count_lib, h.h_count_expected) / max(h.h_count_lib, h.h_count_expected)
                h_count_bonus = 5.0 * ratio
            s = cov_1d + hsqc_score - hsqc_pen - diff_pen + dept_bonus + h_count_bonus
        else:
            # 1D 模式: 原始公式
            cov_score = (h.c_coverage * c_weight + h.h_coverage * h_weight) * 100.0
            c_pen = 5.0 * h.c_rms
            h_pen = 5.0 * h.h_rms
            diff_c = abs(h.c_total_lib - h.c_total_obs)
            diff_h = abs(h.h_total_lib - h.h_total_obs)
            diff_pen = 2.0 * (diff_c + diff_h)
            dept_bonus = 0.0
            if h.dept_total > 0:
                dept_bonus = dept_match_bonus * h.dept_coverage
            h_count_bonus = 0.0
            if h.h_count_expected is not None and h.h_count_lib is not None and h.h_count_expected > 0:
                ratio = min(h.h_count_lib, h.h_count_expected) / max(h.h_count_lib, h.h_count_expected)
                h_count_bonus = h_count_match_bonus * ratio
            s = cov_score - c_pen - h_pen - diff_pen + dept_bonus + h_count_bonus

        # NOE 加分 (若有 NOE 验证数据)
        if use_noe_scoring and noe_scores:
            smiles = h.entry.canonical_smiles or h.entry.smiles
            noe_info = noe_scores.get(smiles)
            if noe_info and noe_info.get('n_total', 0) > 0:
                match_rate = noe_info['n_matched'] / noe_info['n_total']
                noe_bonus = match_rate * 10 - 5 * noe_info.get('n_conflict', 0)
                s += noe_bonus

        scored.append(ScoredCandidate(hit=h, score=s))
    # 降序
    scored.sort(key=lambda x: x.score, reverse=True)
    for i, s in enumerate(scored, start=1):
        s.rank = i
    return scored


# ---------- evidence ----------

@dataclass
class PeakAssignment:
    """单条峰的指认证据."""

    observed_shift: float
    matched: bool
    closest_lib_shift: Optional[float] = None
    delta_ppm: Optional[float] = None
    matched_atom_index: Optional[int] = None


@dataclass
class DeptEvidence:
    """单条 13C 候选的 DEPT 匹配证据."""

    candidate_carbon_ppm: float  # 候选预测的 13C 化学位移
    observed_carbon_ppm: Optional[float]  # 实测 13C 化学位移 (匹配上的, 否则 None)
    delta_ppm: Optional[float]  # 偏差
    observed_mult: Optional[str]  # 实测 DEPT 推断的碳类型 (CH/CH2/CH3/Cq/?)
    expected_mult: Optional[str]  # 候选该 C 的碳类型 (从 SMILES 算)
    mult_match: bool  # 碳类型是否一致
    match_type: str  # 'full' (位置+类型都匹配) / 'partial' (位置匹配但类型不一致) / 'unmatched'


@dataclass
class HCountEvidence:
    """H 总数交叉验证证据."""

    expected: int  # DEPT 推断的 H 总数 (CH*1 + CH2*2 + CH3*3)
    lib_heavy_h: int  # 候选分子中 H 原子数 (RDKit)
    delta: int  # 差值 (lib - expected, 正表示候选 H 多)
    match_score: float  # 0-1 一致性
    notes: List[str] = field(default_factory=list)


@dataclass
class CandidateEvidence:
    """单条候选的证据链."""

    smiles: str
    nmrshiftdb2_id: Optional[int]
    inchi: Optional[str]
    inchikey: Optional[str]
    solvent_lib: Optional[str]
    heavy_atoms: int
    formula: Optional[str]
    # 候选预测谱 (用于绘图与 ML 二次校验)
    lib_13c_shifts: List[float] = field(default_factory=list)
    lib_1h_shifts: List[float] = field(default_factory=list)
    rank: int = 0
    score: float = 0.0
    carbon_evidence: List[PeakAssignment] = field(default_factory=list)
    proton_evidence: List[PeakAssignment] = field(default_factory=list)
    # 重点关注: 实测谱中未匹配上的峰 (可能是新结构)
    unmatched_carbon_shifts: List[float] = field(default_factory=list)
    unmatched_proton_shifts: List[float] = field(default_factory=list)
    # DEPT 推断表
    dept_evidence: List[DeptEvidence] = field(default_factory=list)
    # H 总数一致性
    h_count_evidence: Optional[HCountEvidence] = None
    # 2D NMR 关联验证分 (HSQC/HMBC/COSY)
    correlation_2d_score: Optional[Correlation2DScore] = None
    # 冲突: 分子式/溶剂/碳类型不一致等
    conflicts: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def _greedy_match_pairs(obs: List[float], lib: List[float], tol: float):
    """返回 (pairs, unmatched_obs)."""
    if not lib:
        return [], list(obs)
    if not obs:
        return [], []
    lib_sorted = sorted(lib)
    used = [False] * len(lib_sorted)
    pairs = []
    unmatched = []
    for o in obs:
        best_i = -1
        best_d = None
        for i, l in enumerate(lib_sorted):
            if used[i]:
                continue
            d = abs(o - l)
            if d <= tol and (best_d is None or d < best_d):
                best_d = d
                best_i = i
        if best_i >= 0:
            used[best_i] = True
            pairs.append((o, lib_sorted[best_i], best_d))
        else:
            nearest = min(lib_sorted, key=lambda x: abs(x - o))
            unmatched.append(o)
    return pairs, unmatched


def _expected_carbon_type(atom, mol) -> str:
    """从 RDKit 原子推算预期碳类型: CH/CH2/CH3/Cq/?"""
    if atom is None:
        return '?'
    if atom.GetSymbol() != 'C':
        return '?'  # 只处理碳
    if atom.GetTotalNumHs() == 0:
        return 'Cq'
    if atom.GetTotalNumHs() == 1:
        return 'CH'
    if atom.GetTotalNumHs() == 2:
        return 'CH2'
    if atom.GetTotalNumHs() == 3:
        return 'CH3'
    return 'CH3+'  # > 3 不会发生


def build_dept_evidence(
    lib_c_shifts: List[float],
    obs_c_shifts: List[float],
    obs_dept_summary: List[dict],
    entry_mol,
    c_tol: float,
) -> List[DeptEvidence]:
    """为候选分子构建 DEPT 匹配证据.

    lib_c_shifts: 候选预测的 13C 化学位移列表
    obs_c_shifts: 实测 13C 主谱峰
    obs_dept_summary: 实测 DEPT 推断 [{ppm, mult, sign}, ...]
    entry_mol: RDKit mol 对象 (从 SDF)
    c_tol: 13C 匹配容差

    返回每条候选 C 的 DeptEvidence.
    """
    out: List[DeptEvidence] = []
    if not lib_c_shifts:
        return out

    # 构造 obs ppm -> (mult, sign) 字典 (按 ppm 找最近)
    obs_dept: List[Tuple[float, str, str]] = []
    for d in obs_dept_summary:
        obs_dept.append((float(d.get('ppm', 0)), d.get('mult', '?'), d.get('sign', 'unknown')))

    # 候选每个 C 原子 (按 SMILES 顺序)
    n_atoms = entry_mol.GetNumAtoms() if entry_mol is not None else 0
    for i, c_ppm in enumerate(lib_c_shifts):
        # 找对应原子 (按 SDF 顺序, i 通常对应 mol.GetAtomWithIdx(i))
        atom = entry_mol.GetAtomWithIdx(i) if (entry_mol is not None and i < n_atoms) else None
        expected_mult = _expected_carbon_type(atom, entry_mol) if atom else '?'
        # 找最近 obs ppm
        if not obs_c_shifts:
            out.append(DeptEvidence(
                candidate_carbon_ppm=c_ppm,
                observed_carbon_ppm=None, delta_ppm=None,
                observed_mult=None, expected_mult=expected_mult,
                mult_match=False, match_type='unmatched',
            ))
            continue
        nearest_obs = min(obs_c_shifts, key=lambda x: abs(x - c_ppm))
        delta = nearest_obs - c_ppm
        if abs(delta) > c_tol:
            out.append(DeptEvidence(
                candidate_carbon_ppm=c_ppm,
                observed_carbon_ppm=None, delta_ppm=None,
                observed_mult=None, expected_mult=expected_mult,
                mult_match=False, match_type='unmatched',
            ))
            continue
        # 在 obs_dept 中找最近的多重性
        if obs_dept:
            nearest_dept = min(obs_dept, key=lambda x: abs(x[0] - nearest_obs))
            obs_mult = nearest_dept[1] if abs(nearest_dept[0] - nearest_obs) < c_tol else '?'
        else:
            obs_mult = '?'
        mult_match = (obs_mult == expected_mult) or (obs_mult == 'CH2/CH3' and expected_mult in ('CH2', 'CH3'))
        match_type = 'full' if mult_match else 'partial'
        out.append(DeptEvidence(
            candidate_carbon_ppm=c_ppm,
            observed_carbon_ppm=round(nearest_obs, 3), delta_ppm=round(delta, 3),
            observed_mult=obs_mult, expected_mult=expected_mult,
            mult_match=mult_match, match_type=match_type,
        ))
    return out


def compute_h_count_evidence(
    obs_dept_summary: List[dict],
    lib_mol,
) -> Optional[HCountEvidence]:
    """从 DEPT 推断的 H 总数 vs RDKit 算的候选 H 总数.

    DEPT 推断: CH*1 + CH2*2 + CH3*3 = 实测 H 总数 (减去活泼 H).
    候选 H 总数: RDKit CalcNumH + GetTotalNumHs 之和 (含 N-H, O-H 等活泼氢).
    """
    if not obs_dept_summary:
        return None
    expected = 0
    for d in obs_dept_summary:
        m = d.get('mult', '?')
        if m == 'CH':
            expected += 1
        elif m == 'CH2':
            expected += 2
        elif m == 'CH3':
            expected += 3
        # Cq / ? / CH2/CH3: 不计入 (Cq=0, ? 跳过, CH2/CH3 保守按 0)
    if expected == 0:
        return None

    # 候选 H 总数与分子式
    if lib_mol is None:
        return None
    from rdkit.Chem import rdMolDescriptors
    formula_str = ''
    try:
        formula_str = rdMolDescriptors.CalcMolFormula(lib_mol)
        # 提取 H 数量
        import re
        m = re.search(r'H(\d+)', formula_str)
        lib_heavy_h = int(m.group(1)) if m else 0
    except Exception:
        lib_heavy_h = 0
    if lib_heavy_h == 0:
        return None

    delta = lib_heavy_h - expected
    # match_score: 1 - |delta|/max(expected, lib_heavy_h) ; 0 起步
    denom = max(expected, lib_heavy_h, 1)
    match_score = max(0.0, 1.0 - abs(delta) / denom)
    notes: List[str] = []
    if delta > 0:
        notes.append(f"候选比实测多 {delta} 个 H (可能含 OH/NH/COOH 等活泼 H)")
    elif delta < 0:
        notes.append(f"候选比实测少 {-delta} 个 H (可能丢失部分结构或实测 DEPT 漏检)")
    else:
        notes.append("H 总数完全一致")

    # SymPy 辅助: 候选不饱和度与 H/C 比 (帮助化学家判断结构类型)
    if formula_str:
        try:
            from .formula_solver import (
                FormulaParseError, formula_dou, formula_hc_ratio, parse_formula,
            )
            parsed = parse_formula(formula_str)
            dou = formula_dou(parsed)
            hc = formula_hc_ratio(parsed)
            if dou is not None:
                notes.append(f"候选不饱和度 DoU = {dou:.2f} (苯环/双键/环数总和)")
            if hc is not None:
                notes.append(
                    f"候选 H/C 比 = {hc:.2f} "
                    "(饱和烷烃 ~2.2, 烯烃 ~1.5-1.8, 芳烃 ~1.0)"
                )
        except FormulaParseError:
            pass  # RDKit 输出的分子式 SymPy 解析失败 (如含同位素标记), 不阻断主流程
        except Exception:
            pass  # 其他 SymPy 异常, 静默跳过, 不影响主打分

    return HCountEvidence(
        expected=expected,
        lib_heavy_h=lib_heavy_h,
        delta=delta,
        match_score=match_score,
        notes=notes,
    )


def build_evidence(
    scored: ScoredCandidate,
    observed_c: List[float],
    observed_h: List[float],
    c_tol: float,
    h_tol: float,
    formula_constraint: Optional[str],
    obs_dept_summary: Optional[List[dict]] = None,
) -> CandidateEvidence:
    """为单个打分后的候选构建证据链."""
    hit = scored.hit
    entry = hit.entry

    # 13C 证据
    lib_c = [s.shift for s in entry.spectrum_13c]
    pairs_c, unmatched_c = _greedy_match_pairs(observed_c, lib_c, c_tol)
    c_evid = []
    for o, l, d in pairs_c:
        c_evid.append(PeakAssignment(
            observed_shift=o, matched=True, closest_lib_shift=l, delta_ppm=d,
        ))
    for o in unmatched_c:
        if lib_c:
            nearest = min(lib_c, key=lambda x: abs(x - o))
            delta = nearest - o
        else:
            nearest = None
            delta = None
        c_evid.append(PeakAssignment(
            observed_shift=o, matched=False,
            closest_lib_shift=nearest, delta_ppm=delta,
        ))

    # 1H 证据
    lib_h = [s.shift for s in entry.spectrum_1h]
    pairs_h, unmatched_h = _greedy_match_pairs(observed_h, lib_h, h_tol)
    h_evid = []
    for o, l, d in pairs_h:
        h_evid.append(PeakAssignment(
            observed_shift=o, matched=True, closest_lib_shift=l, delta_ppm=d,
        ))
    for o in unmatched_h:
        if lib_h:
            nearest = min(lib_h, key=lambda x: abs(x - o))
            delta = nearest - o
        else:
            nearest = None
            delta = None
        h_evid.append(PeakAssignment(
            observed_shift=o, matched=False,
            closest_lib_shift=nearest, delta_ppm=delta,
        ))

    # 冲突
    conflicts: List[str] = []
    if formula_constraint:
        try:
            mol = entry.mol()
        except Exception:
            mol = None
        if mol is not None and not formula_match(formula_constraint, mol):
            conflicts.append(f"分子式不匹配 (实测 {formula_constraint} vs 候选 {hit.entry.extra_tags.get('nmrshiftdb2 ID', '?')})")
    if not hit.solvent_ok:
        conflicts.append(f"溶剂不兼容 (实测 {hit.entry.solvent} vs 候选 {hit.entry.solvent})")

    # 计算分子式
    formula = None
    try:
        mol = entry.mol()
        if mol is not None:
            from rdkit.Chem import rdMolDescriptors
            formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        formula = None

    # DEPT 证据
    dept_evidence: List[DeptEvidence] = []
    if obs_dept_summary:
        try:
            entry_mol = entry.mol()
        except Exception:
            entry_mol = None
        dept_evidence = build_dept_evidence(
            lib_c, observed_c, obs_dept_summary, entry_mol, c_tol,
        )

    # H 总数交叉验证
    h_count_evidence: Optional[HCountEvidence] = None
    if obs_dept_summary:
        try:
            entry_mol = entry.mol()
        except Exception:
            entry_mol = None
        h_count_evidence = compute_h_count_evidence(obs_dept_summary, entry_mol)

    notes: List[str] = []
    if unmatched_c:
        notes.append(
            f"未匹配 13C 峰数: {len(unmatched_c)} ({unmatched_c}) — 可能是新结构片段, 需人工复核"
        )
    if unmatched_h:
        notes.append(
            f"未匹配 1H 峰数: {len(unmatched_h)} ({unmatched_h}) — 可能是新结构片段, 需人工复核"
        )
    if dept_evidence:
        full_match = sum(1 for d in dept_evidence if d.match_type == 'full')
        notes.append(
            f"DEPT 碳类型匹配: {full_match}/{len(dept_evidence)} (full)"
        )
    if h_count_evidence and h_count_evidence.notes:
        notes.extend(h_count_evidence.notes)

    return CandidateEvidence(
        smiles=entry.canonical_smiles or entry.smiles,
        nmrshiftdb2_id=entry.rec_id,
        inchi=entry.inchi,
        inchikey=entry.inchikey,
        solvent_lib=entry.solvent,
        heavy_atoms=hit.heavy_atoms,
        formula=formula,
        lib_13c_shifts=lib_c,
        lib_1h_shifts=lib_h,
        rank=scored.rank,
        score=scored.score,
        carbon_evidence=c_evid,
        proton_evidence=h_evid,
        unmatched_carbon_shifts=unmatched_c,
        unmatched_proton_shifts=unmatched_h,
        dept_evidence=dept_evidence,
        h_count_evidence=h_count_evidence,
        conflicts=conflicts,
        notes=notes,
    )
