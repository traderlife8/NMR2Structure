"""2D NMR 关联提取与候选结构验证.

从 2D 谱 (HSQC/HMBC/COSY) 的交叉峰中提取原子关联信息,
并与候选分子结构进行匹配验证, 产出 2D 验证分.

对标: MestReNova 解谱指南 第五章 (二维谱关联) + 第七章 (Verify 打分维度):
  - HSQC 关联匹配率
  - HMBC 远程相关匹配率
  - COSY 邻接匹配率

匹配策略:
  1. 从 nmrshiftdb2 候选的 per-atom 归属构建预期 2D 关联
  2. 对每个观测 2D 峰, 检查候选是否有对应的原子对能解释该峰
  3. 匹配率 = 被解释的观测峰数 / 观测峰总数

依赖: RDKit (分子拓扑分析), numpy
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .io_bruker_2d import Peak2D, Spectrum2D
from .io_sdf import NmrShiftDb2Entry, SpectrumSignal


# ---- 数据结构 ----

@dataclass
class Correlation2DScore:
    """2D NMR 关联验证分."""

    hsqc_match_rate: float = 0.0    # 观测 HSQC 峰被候选解释的比例
    hmbc_match_rate: float = 0.0    # 观测 HMBC 峰被候选解释的比例
    cosy_match_rate: float = 0.0    # 观测 COSY 峰被候选解释的比例
    combined_2d_score: float = 0.0  # 加权综合 (HSQC 0.5 + HMBC 0.3 + COSY 0.2)
    n_hsqc_obs: int = 0             # 观测 HSQC 峰数 (过滤后)
    n_hmbc_obs: int = 0             # 观测 HMBC 峰数 (过滤后)
    n_cosy_obs: int = 0             # 观测 COSY 峰数 (过滤后)
    n_hsqc_matched: int = 0
    n_hmbc_matched: int = 0
    n_cosy_matched: int = 0
    notes: List[str] = field(default_factory=list)


# ---- 关联提取 (从观测 2D 谱) ----

def extract_hsqc_pairs(
    spectrum: Spectrum2D,
    min_snr: float = 20.0,
    max_peaks: int = 200,
) -> List[Tuple[float, float]]:
    """从 HSQC 谱提取 (1H_ppm, 13C_ppm) 直接连接对.

    HSQC 的 F1=13C, F2=1H. 返回 (1H, 13C) 对.
    过滤: SNR >= min_snr, 最多 max_peaks 个.
    """
    if spectrum.experiment_type != 'HSQC':
        return []
    # 按 SNR 过滤
    peaks = [p for p in spectrum.peaks if (p.snr or 0) >= min_snr]
    # 按强度降序取 top
    peaks.sort(key=lambda p: p.intensity, reverse=True)
    peaks = peaks[:max_peaks]
    # HSQC: F1=13C (间接维), F2=1H (直接维)
    return [(p.f2_ppm, p.f1_ppm) for p in peaks]


def extract_hmbc_pairs(
    spectrum: Spectrum2D,
    min_snr: float = 10.0,
    max_peaks: int = 300,
) -> List[Tuple[float, float]]:
    """从 HMBC 谱提取 (1H_ppm, 13C_ppm) 远程关联对.

    HMBC 的 F1=13C, F2=1H. 返回 (1H, 13C) 对.
    过滤: SNR >= min_snr, 最多 max_peaks 个.
    """
    if spectrum.experiment_type != 'HMBC':
        return []
    peaks = [p for p in spectrum.peaks if (p.snr or 0) >= min_snr]
    peaks.sort(key=lambda p: p.intensity, reverse=True)
    peaks = peaks[:max_peaks]
    return [(p.f2_ppm, p.f1_ppm) for p in peaks]


def extract_cosy_pairs(
    spectrum: Spectrum2D,
    min_snr: float = 15.0,
    diagonal_tol: float = 0.1,
    max_peaks: int = 200,
) -> List[Tuple[float, float]]:
    """从 COSY 谱提取 (1H_ppm_a, 1H_ppm_b) 邻接对.

    COSY 的 F1=1H, F2=1H. 去除对角线峰 (|F1-F2| < diagonal_tol).
    去重复: (a,b) 和 (b,a) 只保留一个 (a < b).
    """
    if spectrum.experiment_type != 'COSY':
        return []
    pairs: Set[Tuple[float, float]] = set()
    # 按 SNR 过滤
    peaks = [p for p in spectrum.peaks if (p.snr or 0) >= min_snr]
    peaks.sort(key=lambda p: p.intensity, reverse=True)
    for p in peaks:
        # 去对角线
        if abs(p.f1_ppm - p.f2_ppm) < diagonal_tol:
            continue
        # 规范化: 小值在前
        a, b = min(p.f1_ppm, p.f2_ppm), max(p.f1_ppm, p.f2_ppm)
        pairs.add((round(a, 3), round(b, 3)))
        if len(pairs) >= max_peaks:
            break
    return sorted(pairs)


# ---- 从候选分子构建预期 2D 关联 ----

def _build_atom_shift_maps(
    entry: NmrShiftDb2Entry,
) -> Tuple[Dict[int, float], Dict[int, float]]:
    """从 nmrshiftdb2 entry 构建原子索引 -> 化学位移映射.

    Returns:
        (c_shifts_by_atom, h_shifts_by_atom)
        c_shifts_by_atom: {atom_idx: 13C_ppm}
        h_shifts_by_atom: {atom_idx: 1H_ppm}
    """
    c_map: Dict[int, float] = {}
    h_map: Dict[int, float] = {}
    for sig in entry.spectrum_13c:
        if sig.atom_index is not None:
            c_map[sig.atom_index] = sig.shift
    for sig in entry.spectrum_1h:
        if sig.atom_index is not None:
            h_map[sig.atom_index] = sig.shift
    return c_map, h_map


def _get_bond_graph(mol) -> Dict[int, Set[int]]:
    """构建分子键图 (邻接表)."""
    graph: Dict[int, Set[int]] = {i: set() for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[i].add(j)
        graph[j].add(i)
    return graph


def _bfs_distance(graph: Dict[int, Set[int]], start: int, max_dist: int = 4) -> Dict[int, int]:
    """BFS 计算从 start 到所有原子的最短距离 (限制最大距离)."""
    dist: Dict[int, int] = {start: 0}
    queue = [start]
    head = 0
    while head < len(queue):
        node = queue[head]
        head += 1
        if dist[node] >= max_dist:
            continue
        for nb in graph.get(node, set()):
            if nb not in dist:
                dist[nb] = dist[node] + 1
                queue.append(nb)
    return dist


def compute_expected_hsqc_from_entry(
    entry: NmrShiftDb2Entry,
) -> List[Tuple[float, float]]:
    """从 nmrshiftdb2 entry 的 per-atom 归属构建预期 HSQC 关联.

    对每个带 H 的 C 原子: 若同时有 13C 和 1H 归属, 产生 (1H, 13C) 对.
    需要 RDKit 解析 mol_block 获取 C-H 键连接关系.
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolBlock(entry.mol_block, sanitize=True, removeHs=False)
        if mol is None:
            return []
    except Exception:
        return []

    c_map, h_map = _build_atom_shift_maps(entry)
    if not c_map or not h_map:
        return []

    pairs: List[Tuple[float, float]] = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        atom1 = mol.GetAtomWithIdx(a1)
        atom2 = mol.GetAtomWithIdx(a2)
        # 找 C-H 键
        c_idx, h_idx = None, None
        if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'H':
            c_idx, h_idx = a1, a2
        elif atom2.GetSymbol() == 'C' and atom1.GetSymbol() == 'H':
            c_idx, h_idx = a2, a1
        if c_idx is not None and c_idx in c_map and h_idx in h_map:
            pairs.append((h_map[h_idx], c_map[c_idx]))
    return pairs


def compute_expected_hmbc_from_entry(
    entry: NmrShiftDb2Entry,
) -> List[Tuple[float, float]]:
    """从 nmrshiftdb2 entry 构建预期 HMBC 关联 (2-3 键 C...H).

    枚举所有 C 和 H 原子对, 若最短路径距离为 2 或 3, 产生 (1H, 13C) 对.
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolBlock(entry.mol_block, sanitize=True, removeHs=False)
        if mol is None:
            return []
    except Exception:
        return []

    c_map, h_map = _build_atom_shift_maps(entry)
    if not c_map or not h_map:
        return []

    graph = _get_bond_graph(mol)

    # 收集所有 C 和 H 原子
    c_atoms = [idx for idx in c_map.keys()]
    h_atoms = [idx for idx in h_map.keys()]

    pairs: Set[Tuple[float, float]] = set()
    for c_idx in c_atoms:
        dist = _bfs_distance(graph, c_idx, max_dist=3)
        for h_idx in h_atoms:
            d = dist.get(h_idx, 99)
            if 2 <= d <= 3:
                pairs.add((round(h_map[h_idx], 3), round(c_map[c_idx], 3)))
    return sorted(pairs)


def compute_expected_cosy_from_entry(
    entry: NmrShiftDb2Entry,
) -> List[Tuple[float, float]]:
    """从 nmrshiftdb2 entry 构建预期 COSY 关联 (3 键 H-H, vicinal).

    枚举所有 H 原子对, 若最短路径距离为 3 (H-C-C-H), 产生 (1H_a, 1H_b) 对.
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolBlock(entry.mol_block, sanitize=True, removeHs=False)
        if mol is None:
            return []
    except Exception:
        return []

    _, h_map = _build_atom_shift_maps(entry)
    if not h_map:
        return []

    graph = _get_bond_graph(mol)
    h_atoms = sorted(h_map.keys())

    pairs: Set[Tuple[float, float]] = set()
    for i, h1 in enumerate(h_atoms):
        dist = _bfs_distance(graph, h1, max_dist=3)
        for h2 in h_atoms[i + 1:]:
            d = dist.get(h2, 99)
            if d == 3:  # vicinal: H-C-C-H
                a = min(h_map[h1], h_map[h2])
                b = max(h_map[h1], h_map[h2])
                pairs.add((round(a, 3), round(b, 3)))
    return sorted(pairs)


# ---- 匹配打分 ----

def _match_hsqc(
    observed: List[Tuple[float, float]],
    expected: List[Tuple[float, float]],
    h_tol: float = 0.3,
    c_tol: float = 2.5,
) -> int:
    """计算观测 HSQC 峰中被预期关联解释的数量.

    observed/expected: [(1H_ppm, 13C_ppm), ...]
    """
    if not observed or not expected:
        return 0
    matched = 0
    for h_obs, c_obs in observed:
        for h_exp, c_exp in expected:
            if abs(h_obs - h_exp) <= h_tol and abs(c_obs - c_exp) <= c_tol:
                matched += 1
                break
    return matched


def _match_hmbc(
    observed: List[Tuple[float, float]],
    expected: List[Tuple[float, float]],
    h_tol: float = 0.4,
    c_tol: float = 3.0,
) -> int:
    """计算观测 HMBC 峰中被预期关联解释的数量.

    HMBC 容差比 HSQC 稍大 (远程关联化学位移不确定性更高).
    """
    if not observed or not expected:
        return 0
    matched = 0
    for h_obs, c_obs in observed:
        for h_exp, c_exp in expected:
            if abs(h_obs - h_exp) <= h_tol and abs(c_obs - c_exp) <= c_tol:
                matched += 1
                break
    return matched


def _match_cosy(
    observed: List[Tuple[float, float]],
    expected: List[Tuple[float, float]],
    h_tol: float = 0.3,
) -> int:
    """计算观测 COSY 峰中被预期关联解释的数量.

    observed/expected: [(1H_a, 1H_b), ...] (已规范化 a < b)
    """
    if not observed or not expected:
        return 0
    matched = 0
    for h_a_obs, h_b_obs in observed:
        for h_a_exp, h_b_exp in expected:
            # 两种顺序都检查
            if (abs(h_a_obs - h_a_exp) <= h_tol and abs(h_b_obs - h_b_exp) <= h_tol):
                matched += 1
                break
            if (abs(h_a_obs - h_b_exp) <= h_tol and abs(h_b_obs - h_a_exp) <= h_tol):
                matched += 1
                break
    return matched


# ---- 顶层接口 ----

def score_candidate_2d(
    entry: NmrShiftDb2Entry,
    spectra_2d: Dict[str, Spectrum2D],
    h_tol: float = 0.3,
    c_tol: float = 2.5,
    hsqc_min_snr: float = 20.0,
    hmbc_min_snr: float = 10.0,
    cosy_min_snr: float = 15.0,
) -> Correlation2DScore:
    """对单个候选分子进行 2D NMR 关联验证.

    Args:
        entry: nmrshiftdb2 候选 (含 per-atom 归属)
        spectra_2d: {experiment_type: Spectrum2D} 字典
        h_tol: 1H 容差 ppm
        c_tol: 13C 容差 ppm
        hsqc_min_snr: HSQC 峰 SNR 过滤阈值
        hmbc_min_snr: HMBC 峰 SNR 过滤阈值
        cosy_min_snr: COSY 峰 SNR 过滤阈值

    Returns:
        Correlation2DScore 验证结果
    """
    result = Correlation2DScore()

    # 提取观测关联
    obs_hsqc: List[Tuple[float, float]] = []
    obs_hmbc: List[Tuple[float, float]] = []
    obs_cosy: List[Tuple[float, float]] = []

    if 'HSQC' in spectra_2d:
        obs_hsqc = extract_hsqc_pairs(spectra_2d['HSQC'], min_snr=hsqc_min_snr)
    if 'HMBC' in spectra_2d:
        obs_hmbc = extract_hmbc_pairs(spectra_2d['HMBC'], min_snr=hmbc_min_snr)
    if 'COSY' in spectra_2d:
        obs_cosy = extract_cosy_pairs(spectra_2d['COSY'], min_snr=cosy_min_snr)

    result.n_hsqc_obs = len(obs_hsqc)
    result.n_hmbc_obs = len(obs_hmbc)
    result.n_cosy_obs = len(obs_cosy)

    # 无观测数据时返回零分
    if not obs_hsqc and not obs_hmbc and not obs_cosy:
        result.notes.append('无可用 2D 观测数据')
        return result

    # 构建候选的预期关联
    exp_hsqc = compute_expected_hsqc_from_entry(entry)
    exp_hmbc = compute_expected_hmbc_from_entry(entry)
    exp_cosy = compute_expected_cosy_from_entry(entry)

    if not exp_hsqc and not exp_hmbc and not exp_cosy:
        result.notes.append('候选缺少 per-atom 归属, 无法构建预期 2D 关联')
        return result

    # 匹配
    if obs_hsqc and exp_hsqc:
        result.n_hsqc_matched = _match_hsqc(obs_hsqc, exp_hsqc, h_tol, c_tol)
        result.hsqc_match_rate = result.n_hsqc_matched / len(obs_hsqc)

    if obs_hmbc and exp_hmbc:
        result.n_hmbc_matched = _match_hmbc(obs_hmbc, exp_hmbc, h_tol + 0.1, c_tol + 0.5)
        result.hmbc_match_rate = result.n_hmbc_matched / len(obs_hmbc)

    if obs_cosy and exp_cosy:
        result.n_cosy_matched = _match_cosy(obs_cosy, exp_cosy, h_tol)
        result.cosy_match_rate = result.n_cosy_matched / len(obs_cosy)

    # 加权综合分 (仅有数据的维度参与加权)
    weights = {}
    scores = {}
    if obs_hsqc:
        weights['hsqc'] = 0.5
        scores['hsqc'] = result.hsqc_match_rate
    if obs_hmbc:
        weights['hmbc'] = 0.3
        scores['hmbc'] = result.hmbc_match_rate
    if obs_cosy:
        weights['cosy'] = 0.2
        scores['cosy'] = result.cosy_match_rate

    if weights:
        total_w = sum(weights.values())
        result.combined_2d_score = sum(
            weights[k] * scores[k] for k in weights
        ) / total_w

    # 注释
    if obs_hsqc:
        result.notes.append(
            f'HSQC: {result.n_hsqc_matched}/{result.n_hsqc_obs} '
            f'({result.hsqc_match_rate:.0%})')
    if obs_hmbc:
        result.notes.append(
            f'HMBC: {result.n_hmbc_matched}/{result.n_hmbc_obs} '
            f'({result.hmbc_match_rate:.0%})')
    if obs_cosy:
        result.notes.append(
            f'COSY: {result.n_cosy_matched}/{result.n_cosy_obs} '
            f'({result.cosy_match_rate:.0%})')

    return result
