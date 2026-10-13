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

import csv
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
    noe_match_rate: Optional[float] = None   # NOE 命中率 (0-1)
    noe_conflict_count: int = 0              # NOE 冲突数
    noe_rms_distance: Optional[float] = None # NOE RMS 距离偏差 (A)
    # Part 2.A: HMBC ¹J 残留过滤统计
    n_hmbc_filtered_oneJ: int = 0            # 滤掉的 ¹J 直接连接残留数
    # Part 2.C: Karplus ³J 二面角验证
    j_match_rate: Optional[float] = None     # J 命中率 (0-1)
    j_conflict_count: int = 0                # J 冲突数
    j_rms_error_hz: Optional[float] = None   # J RMS 偏差 (Hz)
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
    # 按 SNR 过滤. snr=None 表示峰来自显式峰表 (如 MNova 导出, 无噪声估计), 放行.
    peaks = [p for p in spectrum.peaks
             if p.snr is None or p.snr >= min_snr]
    # 按强度降序取 top
    peaks.sort(key=lambda p: p.intensity, reverse=True)
    peaks = peaks[:max_peaks]
    # HSQC: F1=13C (间接维), F2=1H (直接维)
    return [(p.f2_ppm, p.f1_ppm) for p in peaks]


def extract_hmbc_pairs(
    spectrum: Spectrum2D,
    min_snr: float = 10.0,
    max_peaks: int = 300,
    hsqc_pairs: Optional[List[Tuple[float, float]]] = None,
    h_tol: float = 0.05,
    c_tol: float = 1.0,
) -> Tuple[List[Tuple[float, float]], int]:
    """从 HMBC 谱提取 (1H_ppm, 13C_ppm) 远程关联对.

    HMBC 的 F1=13C, F2=1H. 返回 (1H, 13C) 对.
    过滤: SNR >= min_snr, 最多 max_peaks 个.
    """
    if spectrum.experiment_type != 'HMBC':
        return [], 0
    # snr=None 表示峰来自显式峰表 (如 MNova 导出), 放行.
    peaks = [p for p in spectrum.peaks
             if p.snr is None or p.snr >= min_snr]
    peaks.sort(key=lambda p: p.intensity, reverse=True)
    peaks = peaks[:max_peaks]
    pairs = [(p.f2_ppm, p.f1_ppm) for p in peaks]

    # Part 2.A: 过滤已知 ¹J 直接连接残留
    if hsqc_pairs:
        hsqc_normalized = [(round(h, 3), round(c, 3)) for h, c in hsqc_pairs]
        filtered: List[Tuple[float, float]] = []
        n_filtered = 0
        for h_ppm, c_ppm in pairs:
            is_oneJ = False
            for kh, kc in hsqc_normalized:
                if abs(h_ppm - kh) <= h_tol and abs(c_ppm - kc) <= c_tol:
                    is_oneJ = True
                    break
            if is_oneJ:
                n_filtered += 1
            else:
                filtered.append((h_ppm, c_ppm))
        return filtered, n_filtered

    return pairs, 0


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
    if spectrum.experiment_type not in ('COSY', 'ROESY', 'NOESY', 'TOCSY'):
        return []
    pairs: Set[Tuple[float, float]] = set()
    # 按 SNR 过滤. snr=None 表示峰来自显式峰表 (如 MNova 导出), 放行.
    peaks = [p for p in spectrum.peaks
             if p.snr is None or p.snr >= min_snr]
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
        # Part 2.A: 传入 HSQC 配对, 过滤已知 ¹J 残留
        obs_hmbc, n_oneJ_filtered = extract_hmbc_pairs(
            spectra_2d['HMBC'], min_snr=hmbc_min_snr,
            hsqc_pairs=obs_hsqc if obs_hsqc else None,
        )
        result.n_hmbc_filtered_oneJ = n_oneJ_filtered  # 类型需要定义
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


# ---- 全归属表生成 (化学家的最终交付物) ----

@dataclass
class AssignmentRow:
    """归属表的一行."""

    atom_label: str       # 'C-1', 'C-2', ...
    c_ppm: Optional[float] = None
    multiplicity: str = 'unknown'  # 'CH'/'CH2'/'CH3'/'Cq'
    h_ppm: Optional[float] = None
    h_mult: Optional[str] = None  # 's'/'d'/'t'/'q'/'m'
    hsqc_match: bool = False
    cosy_partners: List[str] = field(default_factory=list)
    hmbc_to_c: List[str] = field(default_factory=list)
    # 产物 G: 1H 峰来源标记 (measured/estimated/'')
    h_source: str = ''


@dataclass
class FullAssignmentTable:
    """完整归属表."""

    rows: List[AssignmentRow] = field(default_factory=list)
    structure_smiles: str = ''
    confidence: float = 0.0
    unassigned_obs_c: List[float] = field(default_factory=list)
    unassigned_obs_h: List[float] = field(default_factory=list)


def build_full_assignment_table(
    entry: NmrShiftDb2Entry,
    interpretation: 'Interpretation2D',
    observed_c13: List[float],
    observed_h1: List[float],
    h_tol: float = 0.3,
    c_tol: float = 2.5,
) -> FullAssignmentTable:
    """为确定结构生成完整 NMR 归属表.

    对标论文中的 Table 1 格式:
    | 位置 | δC | 碳级数 | δH | HSQC | COSY | HMBC |

    Args:
        entry: nmrshiftdb2 候选 (含 per-atom 归属)
        interpretation: 2D 解谱中间结果
        observed_c13: 实测 13C 位移
        observed_h1: 实测 1H 位移
        h_tol: 1H 容差
        c_tol: 13C 容差
    """
    table = FullAssignmentTable(
        structure_smiles=entry.canonical_smiles or entry.smiles,
    )

    # 从 entry 构建原子归属
    c_map, h_map = _build_atom_shift_maps(entry)

    # 解析分子获取 C-H 键关系
    try:
        from rdkit import Chem
        mol = Chem.MolFromMolBlock(entry.mol_block, sanitize=True, removeHs=False)
    except Exception:
        mol = None

    if mol is None:
        table.unassigned_obs_c = list(observed_c13)
        table.unassigned_obs_h = list(observed_h1)
        return table

    # 构建 C-H 键映射: c_idx -> [h_idx, ...]
    c_to_h: Dict[int, List[int]] = {}
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        at1, at2 = mol.GetAtomWithIdx(a1), mol.GetAtomWithIdx(a2)
        if at1.GetSymbol() == 'C' and at2.GetSymbol() == 'H':
            c_to_h.setdefault(a1, []).append(a2)
        elif at2.GetSymbol() == 'C' and at1.GetSymbol() == 'H':
            c_to_h.setdefault(a2, []).append(a1)

    # 确定每个 C 的级数
    def _get_mult(c_idx: int) -> str:
        n_h = len(c_to_h.get(c_idx, []))
        return {0: 'Cq', 1: 'CH', 2: 'CH2', 3: 'CH3'}.get(n_h, 'unknown')

    # 为每个有归属的 C 原子生成一行
    assigned_c_shifts: Set[float] = set()
    assigned_h_shifts: Set[float] = set()
    row_idx = 0

    for c_idx in sorted(c_map.keys(), key=lambda i: c_map[i], reverse=True):
        c_ppm = c_map[c_idx]
        mult = _get_mult(c_idx)
        row_idx += 1
        label = f'C-{row_idx}'

        # 找对应的 H
        h_ppm_val: Optional[float] = None
        h_indices = c_to_h.get(c_idx, [])
        if h_indices and h_indices[0] in h_map:
            h_ppm_val = h_map[h_indices[0]]

        # HSQC 匹配: 检查观测 HSQC 中是否有对应
        # 产物 G: 估算峰 (h_source='estimated') 使用 1.5 倍宽松容差
        hsqc_ok = False
        matched_h_source = ''
        for pair in interpretation.ch_pairs:
            pair_h_tol = h_tol * 1.5 if pair.h_source == 'estimated' else h_tol
            if abs(pair.c_ppm - c_ppm) <= c_tol:
                if h_ppm_val is not None and abs(pair.h_ppm - h_ppm_val) <= pair_h_tol:
                    hsqc_ok = True
                    matched_h_source = pair.h_source
                    break
                elif h_ppm_val is None:
                    break

        # COSY partners: 找观测 COSY 中与该 H 相关的其他 H
        cosy_parts: List[str] = []
        if h_ppm_val is not None:
            for ss in interpretation.spin_systems:
                if any(abs(h_ppm_val - sh) <= h_tol for sh in ss.h_ppms):
                    for other_h in ss.h_ppms:
                        if abs(other_h - h_ppm_val) > h_tol:
                            for other_pair in interpretation.ch_pairs:
                                if abs(other_pair.h_ppm - other_h) <= h_tol:
                                    cosy_parts.append(
                                        other_pair.assigned_h_label or f'H({other_h:.2f})')
                                    break
                    break

        # HMBC: 找观测 HMBC 中从其他 H 到该 C 的远程连接
        hmbc_labels: List[str] = []
        for conn in interpretation.hmbc_connections:
            if abs(conn.c_ppm - c_ppm) <= c_tol:
                for pair in interpretation.ch_pairs:
                    if abs(pair.h_ppm - conn.h_ppm) <= h_tol:
                        hmbc_labels.append(
                            pair.assigned_h_label or f'H({conn.h_ppm:.2f})')
                        break

        row = AssignmentRow(
            atom_label=label,
            c_ppm=round(c_ppm, 2),
            multiplicity=mult,
            h_ppm=round(h_ppm_val, 3) if h_ppm_val is not None else None,
            hsqc_match=hsqc_ok,
            cosy_partners=cosy_parts[:5],
            hmbc_to_c=hmbc_labels[:5],
            h_source=matched_h_source,
        )
        table.rows.append(row)
        assigned_c_shifts.add(round(c_ppm, 1))
        if h_ppm_val is not None:
            assigned_h_shifts.add(round(h_ppm_val, 2))

    # 未归属的观测峰
    table.unassigned_obs_c = [
        c for c in observed_c13
        if not any(abs(c - ac) <= c_tol for ac in assigned_c_shifts)
    ]
    table.unassigned_obs_h = [
        h for h in observed_h1
        if not any(abs(h - ah) <= h_tol for ah in assigned_h_shifts)
    ]

    # 置信度: 已归属 C 占总观测 C 的比例 (封顶 1.0)
    if observed_c13:
        table.confidence = min(1.0, len(assigned_c_shifts) / len(observed_c13))

    return table


def write_assignment_table_csv(
    table: FullAssignmentTable,
    output_path: str,
    rank: int = 0,
    score: Optional[float] = None,
) -> None:
    """把 FullAssignmentTable 写成论文 Table 1 格式的 CSV (化学家可直接复核).

    输出列:
      位置 | δC (ppm) | 碳级数 | δH (ppm) | H多重峰 | HSQC ✓/— | COSY partners | HMBC to C

    文件还包含元信息头注释: Rank / 结构 SMILES / 置信度 / 未归属观测峰.
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        # 元信息头
        w.writerow(['# Rank', rank])
        if score is not None:
            w.writerow(['# Score', f'{score:.4f}'])
        w.writerow(['# Structure SMILES', table.structure_smiles or ''])
        w.writerow(['# Confidence', f'{table.confidence:.2%}'])
        if table.unassigned_obs_c:
            w.writerow(['# Unassigned 13C (ppm)',
                        '; '.join(f'{c:.3f}' for c in table.unassigned_obs_c)])
        if table.unassigned_obs_h:
            w.writerow(['# Unassigned 1H (ppm)',
                        '; '.join(f'{h:.3f}' for h in table.unassigned_obs_h)])
        w.writerow([])
        # 表头
        w.writerow(['位置', 'δC (ppm)', '碳级数', 'δH (ppm)', 'H多重峰',
                    'HSQC', 'COSY partners', 'HMBC to C', '1H来源'])
        # 数据行
        for r in table.rows:
            w.writerow([
                r.atom_label,
                f'{r.c_ppm:.2f}' if r.c_ppm is not None else '',
                r.multiplicity,
                f'{r.h_ppm:.3f}' if r.h_ppm is not None else '',
                r.h_mult or '',
                '✓' if r.hsqc_match else '—',
                '; '.join(r.cosy_partners),
                '; '.join(r.hmbc_to_c),
                r.h_source,
            ])


# ---- NOE 3D 构象验证 ----

# Part 2.C: Karplus ³J 三维二面角验证
# Karplus 公式: ³J(H-C-C-H) = A·cos²(θ) + B·cos(θ) + C
# 化学家经典常数:
#   脂肪饱和: A=7, B=-1, C=0 (Karplus 1963)
#   烯烃:      A=7, B=-1.5, C=2
#   芳香:      A=8.5, B=-1, C=0.5
_DEFAULT_KARPLUS = {
    'aliphatic': (7.0, -1.0, 0.0),
    'olefinic': (7.0, -1.5, 2.0),
    'aromatic':  (8.5, -1.0, 0.5),
}


@dataclass
class JVerificationResult:
    """Karplus 3D 二面角 J 耦合验证结果 (Part 2.C).

    区别于 NOE 距离验证: 走二面角 + Karplus 方程路向验证.
    与 NOE 3D 距离一起构成“距离 + 二面角”双维度反向验证.
    """

    n_total: int = 0                # 可验证的 J 约束数 (两个 H 都有归属 + 有 4 键路径)
    n_matched: int = 0              # 预测 J 与实测 J 偏差 ≤ j_tolerance_hz 的
    n_conflict: int = 0             # 偏差 ≥ 2 * j_tolerance_hz 的
    rms_j_error_hz: float = 0.0     # 预测 vs 实测 J 的 RMS (Hz)
    conflicts: List[dict] = field(default_factory=list)  # {h_ppm, j_obs, j_pred, theta_deg}
    best_conformer_idx: int = 0
    verification_status: str = ''   # 'good' / 'partial' / 'conflict' / 'fallback_2d' / 'skipped'
    karplus_constants: dict = field(default_factory=dict)  # 使用的 A/B/C 参数, e.g. {'aliphatic': (7,-1,0)}


@dataclass
class NOEVerificationResult:
    """候选结构 NOE 约束验证结果."""

    n_total: int = 0
    n_matched: int = 0
    n_conflict: int = 0
    rms_distance: float = 0.0
    conflicts: List[dict] = field(default_factory=list)
    best_conformer_idx: int = 0
    verification_status: str = ''  # 'good'/'partial'/'conflict'/'fallback_2d'/'skipped'


def verify_candidate_3d_noe(
    candidate_smiles: str,
    noe_contacts: List,
    h_assignments: dict,
    max_distance: float = 5.0,
    distance_tolerance: float = 1.0,
    num_confs: int = 10,
) -> NOEVerificationResult:
    """对候选结构用 RDKit 嵌入多构象, 验证 NOE 约束.

    算法:
    1. RDKit Chem.AddHs -> AllChem.EmbedMultipleConfs(mol, numConfs)
    2. 对每个构象: MMFFOptimizeMolecule -> 力场优化
    3. 遍历每个 NOEContact:
       - 根据 h_assignments 找到两个 H 原子 index
       - 计算 3D 距离
       - 与 expected_distance 比较 (容差 +/- distance_tolerance)
       - 距离 > max_distance -> 标为冲突
    4. 取 NOE 命中率最高的构象为 best_conformer
    5. 输出: 命中率, RMS 距离偏差, 冲突列表

    Args:
        candidate_smiles: 候选结构 SMILES
        noe_contacts: NOEContact 列表 (从 interpret_2d 提取)
        h_assignments: {dH_ppm: atom_index} 从归属表构建
        max_distance: NOE 最大距离阈值 (A)
        distance_tolerance: 距离容差 (A)
        num_confs: 嵌入构象数

    Returns:
        NOEVerificationResult
    """
    result = NOEVerificationResult()

    if not noe_contacts or not h_assignments:
        result.verification_status = 'skipped'
        return result

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        result.verification_status = 'skipped'
        result.conflicts.append({'error': 'RDKit 未安装'})
        return result

    mol = Chem.MolFromSmiles(candidate_smiles)
    if mol is None:
        result.verification_status = 'skipped'
        return result

    mol = Chem.AddHs(mol)

    # 嵌入多构象
    try:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol, numConfs=num_confs, randomSeed=42,
            useBasicKnowledge=True, pruneRmsThresh=0.5)
        if len(conf_ids) == 0:
            result.verification_status = 'fallback_2d'
            return result
    except Exception:
        result.verification_status = 'fallback_2d'
        return result

    # 力场优化每个构象
    for conf_id in conf_ids:
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=int(conf_id))
        except Exception:
            pass  # 优化失败不阻断

    # 筛选可验证的 NOE 对 (两个 H 都有归属)
    verifiable: List[tuple] = []  # (idx1, idx2, expected_dist, contact)
    for nc in noe_contacts:
        idx1 = h_assignments.get(round(nc.h1_ppm, 3))
        idx2 = h_assignments.get(round(nc.h2_ppm, 3))
        if idx1 is not None and idx2 is not None and idx1 != idx2:
            verifiable.append((idx1, idx2, nc.expected_distance, nc))

    result.n_total = len(verifiable)
    if result.n_total == 0:
        result.verification_status = 'skipped'
        return result

    # 对每个构象计算 NOE 命中率
    best_rate = -1.0
    best_conf_idx = 0
    best_matched = 0
    best_conflicts_list: List[dict] = []
    best_rms = 0.0

    for conf_id in conf_ids:
        conf = mol.GetConformer(int(conf_id))
        n_matched = 0
        n_conflict = 0
        conflicts_list: List[dict] = []
        dist_errors: List[float] = []

        for idx1, idx2, exp_dist, nc in verifiable:
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            dist = pos1.Distance(pos2)

            if dist <= max_distance:
                n_matched += 1
                dist_errors.append(abs(dist - exp_dist))
            else:
                n_conflict += 1
                conflicts_list.append({
                    'h1_ppm': nc.h1_ppm,
                    'h2_ppm': nc.h2_ppm,
                    'distance': round(dist, 2),
                    'expected': exp_dist,
                })

        rate = n_matched / result.n_total if result.n_total > 0 else 0.0
        if rate > best_rate:
            best_rate = rate
            best_conf_idx = int(conf_id)
            best_matched = n_matched
            best_conflicts_list = conflicts_list
            best_rms = (sum(d**2 for d in dist_errors) / len(dist_errors))**0.5 if dist_errors else 0.0

    result.n_matched = best_matched
    result.n_conflict = result.n_total - best_matched
    result.rms_distance = round(best_rms, 3)
    result.conflicts = best_conflicts_list
    result.best_conformer_idx = best_conf_idx

    # 状态判定
    if best_rate >= 0.8:
        result.verification_status = 'good'
    elif best_rate >= 0.5:
        result.verification_status = 'partial'
    else:
        result.verification_status = 'conflict'

    return result


# ---- Part 2.C: Karplus 3D J 二面角验证 ----


def _find_h_c_c_h_path(mol, h1_idx: int, h2_idx: int,
                        max_bond_path: int = 5) -> Optional[List[int]]:
    """BFS 找 h1 -> h2 间长度 ≤ max_bond_path (默认 5 = 4 键) 的最短键路径.

    返回 4 原子路径 [h1, c1, c2, h2]. 若路径含 4 个原子 (即 H-C-C-H), 可用于
    二面角计算. 失败返回 None.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return None

    if h1_idx == h2_idx:
        return None
    # BFS
    from collections import deque
    queue = deque()
    visited = {h1_idx}
    # parent: idx -> (prev_idx, bond_path_atom_indices)
    parent: Dict[int, Tuple[int, List[int]]] = {h1_idx: (-1, [h1_idx])}
    queue.append(h1_idx)

    # 收集每个原子的邻居
    def neighbors(idx: int):
        atom = mol.GetAtomWithIdx(idx)
        out = []
        for bond in atom.GetBonds():
            out.append(bond.GetOtherAtomIdx(idx))
        return out

    while queue:
        cur = queue.popleft()
        path = parent[cur][1]
        if len(path) > max_bond_path:
            continue
        for nxt in neighbors(cur):
            if nxt in visited:
                continue
            new_path = path + [nxt]
            parent[nxt] = (cur, new_path)
            if nxt == h2_idx:
                # 检查路径原子数 = 4 (H-C-C-H), 即 len=4 的路径包含 3 个键
                if len(new_path) == 4:
                    return new_path
                # 若更长, 不在 4 键路径范围内, 继续找
                # (BFS 已保证最短, 不会再短于此)
                if len(new_path) < max_bond_path:
                    visited.add(nxt)
                    queue.append(nxt)
                # 否则放弃
                continue
            visited.add(nxt)
            queue.append(nxt)
    return None


def _calculate_dihedral_rad(conf, atom_indices: List[int]) -> Optional[float]:
    """从 conformer 取 4 原子位置计算二面角 (弧度)."""
    try:
        from rdkit.Chem import rdMolTransforms
        if len(atom_indices) != 4:
            return None
        return rdMolTransforms.GetDihedralRad(conf, *atom_indices)
    except (ImportError, ValueError, RuntimeError):
        return None


def _classify_karplus_type(mol, path: List[int]) -> str:
    """根据 H-C-C-H 路径中间 C-C 键的杂化选择 Karplus 常数."""
    if len(path) < 4:
        return 'aliphatic'
    try:
        from rdkit import Chem as _Chem
    except ImportError:
        return 'aliphatic'
    c1, c2 = path[1], path[2]
    a1 = mol.GetAtomWithIdx(c1)
    a2 = mol.GetAtomWithIdx(c2)
    # 芳香: c 原子在芳香环中
    if a1.GetIsAromatic() and a2.GetIsAromatic():
        return 'aromatic'
    # 烯烃: 任一 c 是 sp2 (= 或 C=C)
    if a1.GetHybridization() == _Chem.HybridizationType.SP2 or \
            a2.GetHybridization() == _Chem.HybridizationType.SP2:
        return 'olefinic'
    return 'aliphatic'


def _karplus_predict(theta_rad: float, A: float, B: float, C: float) -> float:
    """Karplus 方程预测 ³J."""
    cos_theta = float(np.cos(theta_rad))
    return float(A * cos_theta * cos_theta + B * cos_theta + C)


def verify_candidate_3d_j(
    candidate_smiles: str,
    j_constraints: List[dict],
    h_assignments: dict,
    karplus_constants: Optional[dict] = None,
    j_tolerance_hz: float = 1.5,
    max_conflict_hz: float = 3.0,
    num_confs: int = 10,
) -> JVerificationResult:
    """对候选结构用 RDKit 多构象 + Karplus 预测 J, 与实测 J 对比 (Part 2.C).

    算法:
    1. Chem.AddHs + EmbedMultipleConfs + MMFFOptimizeMolecule
    2. 对每个 j_constraint:
       a. 从 h_assignments 找两个 H 原子 idx
       b. BFS 找 H-C-C-H 4 原子路径
       c. 计算 H-C-C-H 二面角
       d. 由路径中间 C-C 杂化选 Karplus 常数
       e. 预测 J, 与实测 j_obs_hz 对比
    3. 取 J 命中率最高的 conformer 为 best

    Args:
        candidate_smiles: 候选 SMILES
        j_constraints: [{h1_label, h2_label, h1_ppm, h2_ppm, j_obs_hz, multiplicity}, ...]
        h_assignments: {dH_ppm: atom_index} 从 FullAssignmentTable 构建
        karplus_constants: {type: (A, B, C)} 覆盖默认
        j_tolerance_hz: matched 判定容差 (Hz, 默认 1.5)
        max_conflict_hz: conflict 判定阈值 (Hz, 默认 3.0)
        num_confs: 嵌入构象数

    Returns:
        JVerificationResult
    """
    result = JVerificationResult()
    if karplus_constants:
        result.karplus_constants = dict(karplus_constants)
    else:
        result.karplus_constants = dict(_DEFAULT_KARPLUS)

    if not j_constraints or not h_assignments:
        result.verification_status = 'skipped'
        return result

    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        from rdkit.Chem import AllChem
    except ImportError:
        result.verification_status = 'skipped'
        result.conflicts.append({'error': 'RDKit 未安装'})
        return result

    mol = Chem.MolFromSmiles(candidate_smiles)
    if mol is None:
        result.verification_status = 'skipped'
        return result

    mol = Chem.AddHs(mol)

    try:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol, numConfs=num_confs, randomSeed=42,
            useBasicKnowledge=True, pruneRmsThresh=0.5)
        if len(conf_ids) == 0:
            result.verification_status = 'fallback_2d'
            return result
    except Exception:
        result.verification_status = 'fallback_2d'
        return result

    for conf_id in conf_ids:
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=int(conf_id))
        except Exception:
            pass

    # 过滤 j_constraints: 两个 H 都能解析到 RDKit atom index
    verifiable: List[tuple] = []  # (idx1, idx2, j_obs, h1_ppm, h2_ppm)
    for jc in j_constraints:
        h1_ppm = round(jc.get('h1_ppm', 0.0), 3)
        h2_ppm = round(jc.get('h2_ppm', 0.0), 3)
        j_obs = jc.get('j_obs_hz', 0.0)
        if j_obs <= 0:
            continue
        idx1 = h_assignments.get(h1_ppm)
        idx2 = h_assignments.get(h2_ppm)
        if idx1 is None or idx2 is None or idx1 == idx2:
            continue
        verifiable.append((idx1, idx2, j_obs, h1_ppm, h2_ppm))

    result.n_total = len(verifiable)
    if result.n_total == 0:
        result.verification_status = 'skipped'
        return result

    # 对每个 conformer 计算命中率
    best_rate = -1.0
    best_conf_idx = 0
    best_matched = 0
    best_conflicts_list: List[dict] = []
    best_rms = 0.0

    for conf_id in conf_ids:
        conf = mol.GetConformer(int(conf_id))
        n_matched = 0
        n_conflict = 0
        conflicts_list: List[dict] = []
        j_errors: List[float] = []

        for idx1, idx2, j_obs, h1_ppm, h2_ppm in verifiable:
            # BFS 找 4 原子路径
            path = _find_h_c_c_h_path(mol, int(idx1), int(idx2))
            if path is None:
                n_conflict += 1
                conflicts_list.append({
                    'h1_ppm': h1_ppm, 'h2_ppm': h2_ppm,
                    'j_obs': j_obs, 'j_pred': None,
                    'theta_deg': None,
                    'reason': 'no_4bond_path',
                })
                continue

            theta_rad = _calculate_dihedral_rad(conf, path)
            if theta_rad is None:
                n_conflict += 1
                conflicts_list.append({
                    'h1_ppm': h1_ppm, 'h2_ppm': h2_ppm,
                    'j_obs': j_obs, 'j_pred': None,
                    'theta_deg': None,
                    'reason': 'dihedral_calc_fail',
                })
                continue

            karplus_type = _classify_karplus_type(mol, path)
            A, B, C = result.karplus_constants.get(
                karplus_type, _DEFAULT_KARPLUS['aliphatic'])
            j_pred = _karplus_predict(theta_rad, A, B, C)
            err = abs(j_pred - j_obs)

            if err <= j_tolerance_hz:
                n_matched += 1
                j_errors.append(err)
            elif err >= max_conflict_hz:
                n_conflict += 1
                conflicts_list.append({
                    'h1_ppm': h1_ppm, 'h2_ppm': h2_ppm,
                    'j_obs': round(j_obs, 2),
                    'j_pred': round(j_pred, 2),
                    'theta_deg': round(np.degrees(theta_rad), 1),
                    'karplus_type': karplus_type,
                    'reason': 'large_j_error',
                })
            else:
                # 中间偏差: 不计 matched 也不计 conflict, 仅记录
                j_errors.append(err)

        rate = n_matched / result.n_total if result.n_total > 0 else 0.0
        if rate > best_rate:
            best_rate = rate
            best_conf_idx = int(conf_id)
            best_matched = n_matched
            best_conflicts_list = conflicts_list
            best_rms = (sum(e**2 for e in j_errors) / len(j_errors))**0.5 if j_errors else 0.0

    result.n_matched = best_matched
    result.n_conflict = len(best_conflicts_list)
    result.rms_j_error_hz = round(best_rms, 3)
    result.conflicts = best_conflicts_list
    result.best_conformer_idx = best_conf_idx

    if best_rate >= 0.8:
        result.verification_status = 'good'
    elif best_rate >= 0.5:
        result.verification_status = 'partial'
    else:
        result.verification_status = 'conflict'

    return result
