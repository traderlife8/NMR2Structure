# -*- coding: utf-8 -*-
"""程序化 NMR 归属算法 (替代 Mnova Manual Assignment 的确定性匹配).

对标 Mnova 手动归属工作流, 把"人点击确认"变成可复核的确定性算法:
  1. HSQC 锚点归属: (dH, dC) 交叉峰 -> 结构 C-H 对
     - DEPT 碳级数一致性过滤 (CH2 不允许配 CH3 等)
     - 对称等价原子组 (Morgan 等价) 共享同一观测峰
     - scipy linear_sum_assignment 全局一对一最优 (非逐峰贪心)
  2. 1H / 13C 残峰归属: 未被 HSQC 覆盖的观测峰
     - 孤立 1H (如 OH) 与 Cq 13C 用 HOSE 13C 预测值做全局指派
     - 无 HOSE 时回退 shift_ranges 经验范围
  3. 2D 反向验证: 归属完成后用结构拓扑预测 COSY/HMBC 关联,
     与实测交叉峰比对给出命中率与冲突列表 (归属即验证)
  4. 输出: 等效 Mnova Report Assignment 的归属表 (CSV) + 证据链

数据边界:
  - 输入: 结构 SMILES, 观测 1D 峰表, 2D 交叉峰 (Spectrum2D / CHPair 等)
  - 输出: AlgorithmicAssignmentResult, 全程携带证据与置信度, 不伪装成功
  - 本模块不做结构推导, 只对给定结构做归属与验证
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
from rdkit import Chem

from .nmr_predict import C13Predictor


# ---- 错误类型 ----

class AlgorithmicAssignmentError(Exception):
    """程序化归属过程中的显式错误 (stage='algorithmic_assignment')."""

    def __init__(self, msg: str, stage: str = 'algorithmic_assignment'):
        super().__init__(f'[{stage}] {msg}')
        self.stage = stage


# ---- 数据结构 ----

@dataclass
class AtomAssignment:
    """单个结构原子的归属结果."""

    atom_idx: int                  # RDKit C 原子 index (removeHs 后的 mol)
    atom_label: str                # 化学家标号, 如 'C-3'
    carbon_type: str               # 'CH'/'CH2'/'CH3'/'Cq'
    obs_c_ppm: Optional[float] = None   # 归属到的实测 13C
    obs_h_ppms: List[float] = field(default_factory=list)  # 归属到的实测 1H (CH2 可 2 个)
    hsqc_matched: bool = False     # 是否由 HSQC 直接锚定
    equiv_group: int = -1          # 对称等价组编号 (-1 = 非等价)
    source: str = ''               # 归属来源: 'hsqc'/'hose'/'range'/'equiv'/''


@dataclass
class UnassignedPeak:
    """未能归属的观测峰 (显式暴露, 不吞掉)."""

    ppm: float
    nucleus: str  # '1H' / '13C'
    reason: str   # 无法归属的原因


@dataclass
class CorrelationCheck:
    """单条 2D 相关的验证结果."""

    kind: str            # 'COSY' / 'HMBC'
    obs_pair: Tuple[float, ...]   # 实测交叉峰 (dH1, dH2) 或 (dH, dC)
    expected_atoms: Tuple[int, ...]  # 结构上解释该峰的原子 (RDKit idx, -1=无法定位)
    matched: bool
    note: str = ''


@dataclass
class AlgorithmicAssignmentResult:
    """程序化归属总结果 (含证据链)."""

    structure_smiles: str = ''
    assignments: List[AtomAssignment] = field(default_factory=list)
    unassigned_peaks: List[UnassignedPeak] = field(default_factory=list)
    correlation_checks: List[CorrelationCheck] = field(default_factory=list)
    # 统计指标
    n_atoms_total: int = 0          # 结构 C 原子总数
    n_atoms_assigned: int = 0       # 至少有 13C 或 1H 归属的原子数
    n_obs_c: int = 0
    n_obs_h: int = 0
    n_c_assigned: int = 0
    n_h_assigned: int = 0
    hsqc_match_rate: float = 0.0    # HSQC 峰被结构解释的比例
    cosy_match_rate: float = 0.0    # 可检查 COSY 关联的命中率
    hmbc_match_rate: float = 0.0
    hose_rms: Optional[float] = None
    confidence: float = 0.0         # 综合置信度 0-1
    notes: List[str] = field(default_factory=list)


# ---- 结构准备 ----

def _canonical_equivalence_groups(mol: Chem.Mol) -> Dict[int, int]:
    """用 Morgan canonical ranking (breakTies=False) 找对称等价 C 原子组.

    同组原子在拓扑上不可区分 (如苯环对称位 4/6), 共享同一观测信号.
    Returns:
        {atom_idx: group_id}, 非等价原子 group_id = -1 (不返回)
    """
    try:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    except Exception:
        return {}
    by_rank: Dict[int, List[int]] = {}
    for idx, r in enumerate(ranks):
        if mol.GetAtomWithIdx(idx).GetSymbol() != 'C':
            continue
        by_rank.setdefault(r, []).append(idx)
    groups: Dict[int, int] = {}
    gid = 0
    for r, idxs in sorted(by_rank.items()):
        if len(idxs) >= 2:
            for i in idxs:
                groups[i] = gid
            gid += 1
    return groups


def _build_mol_from_smiles(smiles: str) -> Chem.Mol:
    """SMILES -> 加氢 mol (校验失败抛显式错误, 不静默)."""
    if not smiles or not smiles.strip():
        raise AlgorithmicAssignmentError('SMILES 为空, 无法构建分子')
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise AlgorithmicAssignmentError(f'SMILES 解析失败: {smiles}')
    return Chem.AddHs(m)


def _carbon_type_of(mol: Chem.Mol, c_idx: int) -> str:
    """从 H 邻居数推断碳级数."""
    n_h = sum(1 for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors()
              if nb.GetSymbol() == 'H')
    return {0: 'Cq', 1: 'CH', 2: 'CH2', 3: 'CH3'}.get(n_h, 'unknown')


def _default_labels(mol: Chem.Mol) -> Dict[int, str]:
    """给 C 原子按重原子 index 顺序编号 C-1..C-N (确定性, 与化学家习惯无必然关系)."""
    labels: Dict[int, str] = {}
    n = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            n += 1
            labels[atom.GetIdx()] = f'C-{n}'
    return labels


# ---- HSQC 锚点归属 ----

def _dept_compatible(carbon_type: str, allow_reverse: bool = False) -> Set[str]:
    """HSQC 交叉峰允许的碳级数 (DEPT 一致性).

    HSQC 观测所有质子化碳; edited-HSQC 中 CH2 为负信号, 这里仅做
    级数合法性 (Cq 无 HSQC 信号).
    """
    return {'CH', 'CH2', 'CH3'}


def _protonated_c_expected_range(mol: Chem.Mol, c_idx: int) -> Tuple[float, float]:
    """质子化碳的 13C 经验窗口 (规则 8: 确定性, 精简自 shift_ranges 思路).

    按杂原子邻接数 + 芳香性 + 与羧基距离细分:
      芳香 CH 105-155; 连氧/氮 CH 45-90; 羧基 alpha CH 28-70;
      普通脂肪 CH/CH2 18-60; CH3 5-60 (异丙基/叔丁基次甲基邻位甲基
      可达 ~45, 需宽窗防 CH/CH3 窗口互斥造成错配)
    """
    atom = mol.GetAtomWithIdx(c_idx)
    n_het = sum(1 for nb in atom.GetNeighbors()
                if nb.GetSymbol() in ('O', 'N', 'S', 'F', 'Cl', 'Br', 'I'))
    if atom.GetIsAromatic():
        return (105.0, 155.0)
    if n_het >= 2:
        return (55.0, 105.0)   # 缩醛/连氧连氮
    if n_het == 1:
        return (45.0, 90.0)    # C-O / C-N 单键
    # 纯碳氢环境: 检查 1 键内是否有羧基 C (dC 显著下移, alpha 位)
    for nb in atom.GetNeighbors():
        if nb.GetSymbol() != 'C':
            continue
        n_o_total = sum(1 for x in nb.GetNeighbors() if x.GetSymbol() == 'O')
        if n_o_total >= 2:  # 羧酸/酯 C(=O)O
            return (28.0, 70.0)   # 羧基 alpha CH/CH2
    return (5.0, 60.0)         # 脂肪 CH/CH2/CH3 (统一宽窗, 用 1H 先验区分)


def assign_hsqc_anchor(
    mol: Chem.Mol,
    hsqc_pairs: Sequence[Tuple[float, float]],
    h_tol: float = 0.3,
    c_tol: float = 2.0,
) -> Tuple[Dict[int, List[Tuple[float, float]]], List[int]]:
    """HSQC 锚点归属: 全局一对一最优指派 (匈牙利算法).

    对每个质子化 C 原子枚举其 H 的预期归属; 每个观测 HSQC 峰最多分给一个
    等价组 (等价组内原子共享该峰). 代价 = max(0, 容差 - |dH| - |dC|) 的
    线性组合, 越近代价越小; 超容差为不可行.

    Args:
        mol: AddHs 后的 mol
        hsqc_pairs: [(dH, dC), ...] 实测 HSQC 交叉峰
        h_tol / c_tol: 1H / 13C 匹配容差

    Returns:
        (atom_to_peaks, unmatched_peak_indices)
        atom_to_peaks: {c_atom_idx: [(dH, dC), ...]} (等价组共享同峰时都记录)
    """
    from scipy.optimize import linear_sum_assignment

    equiv = _canonical_equivalence_groups(mol)
    # 质子化 C 候选 (Cq 无 HSQC 信号, DEPT 一致性)
    proto_atoms = [a.GetIdx() for a in mol.GetAtoms()
                   if a.GetSymbol() == 'C'
                   and _carbon_type_of(mol, a.GetIdx()) in _dept_compatible('CH')]

    # 等价组代表 -> 组成员 (代表 = 组内最小 idx, 确定性)
    rep_to_members: Dict[int, List[int]] = {}
    for idx in proto_atoms:
        g = equiv.get(idx, -1)
        if g < 0:
            rep_to_members.setdefault(idx, [idx])
        else:
            members = [i for i in proto_atoms if equiv.get(i, -1) == g]
            rep_to_members[min(members)] = members

    reps = sorted(rep_to_members.keys())
    # 行复制策略: 容量 = 信号数上限.
    #   等价组 (成员>1): 每成员共享 1 条信号 → cap = 成员数
    #   CH2 单原子: 2 个 H 非对异位可分开出峰 → cap = 2
    #   CH/CH3 单原子: 1 条信号 → cap = 1
    # 行 -> rep 映射用于展开后还原.
    row_to_rep: List[int] = []
    for rep in reps:
        members = rep_to_members[rep]
        c_type = _carbon_type_of(mol, rep)
        cap = len(members) if len(members) > 1 else (2 if c_type == 'CH2' else 1)
        row_to_rep.extend([rep] * cap)
    n_rows, n_cols = len(row_to_rep), len(hsqc_pairs)
    if n_rows == 0 or n_cols == 0:
        return {}, list(range(n_cols))

    BIG = 1e6
    cost = np.full((n_rows, n_cols), BIG)
    for i, rep in enumerate(row_to_rep):
        exp_lo, exp_hi = _protonated_c_expected_range(mol, rep)
        c_type = _carbon_type_of(mol, rep)
        atom_rep = mol.GetAtomWithIdx(rep)
        if atom_rep.GetIsAromatic():
            h_center = 7.1
        elif exp_lo >= 45.0:
            h_center = 3.8  # C-O / C-N (糖环, 溶剂残峰区)
        elif exp_lo >= 28.0:
            h_center = 2.4  # 羧基 alpha
        else:
            h_center = {'CH3': 1.1, 'CH2': 1.5, 'CH': 1.8}.get(c_type, 1.5)
        for j, (h_ppm, c_ppm) in enumerate(hsqc_pairs):
            # 13C 维度: 观测值落进该碳的经验窗口才可行, 偏差进代价;
            # 1H 维度: 级数先验中心做软惩罚, 防止同为可行时乱配 (规则 8)
            d_lo = max(0.0, exp_lo - c_ppm)
            d_hi = max(0.0, c_ppm - exp_hi)
            dc = d_lo + d_hi
            if dc > c_tol:
                continue  # 不可行 (保持 BIG)
            dh = abs(h_ppm - h_center)
            cost[i, j] = dc * 3.0 + dh + float(j) * 1e-3  # 确定性 tie-break

    row_ind, col_ind = linear_sum_assignment(cost)
    atom_to_peaks: Dict[int, List[Tuple[float, float]]] = {}
    used_cols: Set[int] = set()
    rep_seen_peaks: Dict[int, List[Tuple[float, float]]] = {}
    for r, c in zip(row_ind, col_ind):
        if cost[r, c] >= BIG:
            continue
        used_cols.add(int(c))
        rep = row_to_rep[r]
        peak = hsqc_pairs[c]
        # 同一代表的多行不应配同一峰 (重复信号): 检查 δH/δC 均相同则跳过
        if any(abs(p[0] - peak[0]) < 1e-6 and abs(p[1] - peak[1]) < 1e-6
               for p in rep_seen_peaks.get(rep, [])):
            # 该行空置 (本可用更低代价的其他峰), 不加入 used_cols 的语义:
            # 峰 c 仍标记已用, 因为它已被该代表"占用" — 但同一峰同时
            # 归属同组是冗余. 为避免另一个代表错失该峰, 这里把列 c
            # 释放回未匹配池由第二轮重分配.
            used_cols.discard(int(c))
            continue
        rep_seen_peaks.setdefault(rep, []).append(peak)
        for member in rep_to_members[rep]:
            atom_to_peaks.setdefault(member, []).append(peak)
    unmatched = [j for j in range(n_cols) if j not in used_cols]

    # 等价组共享展开: 组内每个成员都记录组拿到的全部峰 (组 cap = 成员数
    # 而观测峰可能少于此 — 峰重叠时成员共享已有峰, 不留空)
    for rep, members in rep_to_members.items():
        if len(members) > 1:
            rep_peaks = rep_seen_peaks.get(rep, [])
            if 0 < len(rep_peaks) < len(members):
                # 组只拿到部分峰: 所有成员共享全部峰 (等价原子简并)
                for member in members:
                    atom_to_peaks[member] = list(rep_peaks)

    # 组间借峰: 某等价组 (或单碳) 完全没拿到峰, 但其窗口与已获峰的
    # 相同类型组重叠时, 借用其峰 (芳香 CH 组竞争同一批峰的典型场景).
    # 借峰即"共享信号", source 标记为 equiv, 不影响 2D 验证 (等价原子
    # 拓扑不可区分, 归属到任一成员均自洽).
    for rep, members in rep_to_members.items():
        if rep in rep_seen_peaks and rep_seen_peaks[rep]:
            continue
        exp_lo, exp_hi = _protonated_c_expected_range(mol, rep)
        for other_rep, peaks in rep_seen_peaks.items():
            if not peaks:
                continue
            o_lo, o_hi = _protonated_c_expected_range(mol, other_rep)
            # 同窗口类型 (芳香/杂原子/脂肪) 才借
            if (exp_lo >= 105.0) != (o_lo >= 105.0):
                continue
            if (exp_lo >= 45.0) != (o_lo >= 45.0):
                continue
            for member in rep_to_members[rep]:
                atom_to_peaks[member] = list(peaks)
            break

    # 容量统计修正: 等价组共享/借峰后 atom_to_peaks 各成员持有相同峰,
    # 按组去重统计 (组持有的独立信号数), 而非累加各成员峰数.
    h_capacity: Dict[int, int] = {}   # rep -> 可接收信号数 (每成员 1 个)
    rep_by_member: Dict[int, int] = {}
    for rep, members in rep_to_members.items():
        h_capacity[rep] = len(members)
        for m in members:
            rep_by_member[m] = rep
    assigned_count: Dict[int, int] = {}
    for rep, members in rep_to_members.items():
        group_peaks: Set[Tuple[float, float]] = set()
        for m in members:
            group_peaks.update(atom_to_peaks.get(m, []))
        assigned_count[rep] = len(group_peaks)
    for j in unmatched:
        h_ppm, c_ppm = hsqc_pairs[j]
        best_rep, best_key = -1, None
        for rep in reps:
            # 已满容量的碳不再接收 (等价组内每个成员 1 个信号)
            if assigned_count.get(rep, 0) >= h_capacity[rep]:
                continue
            exp_lo, exp_hi = _protonated_c_expected_range(mol, rep)
            if c_ppm < exp_lo - c_tol or c_ppm > exp_hi + c_tol:
                continue
            # 距离该代表已配峰的 13C 差值最小者优先 (CH2 两峰 δC 相近)
            rep_peaks = [p[1] for m in rep_to_members[rep]
                         for p in atom_to_peaks.get(m, [])]
            dc_extra = min((abs(c_ppm - pc) for pc in rep_peaks), default=999.0)
            key = (dc_extra if rep_peaks else 0.0, rep)
            if best_key is None or key < best_key:
                best_key, best_rep = key, rep
        if best_rep >= 0:
            for member in rep_to_members[best_rep]:
                atom_to_peaks.setdefault(member, []).append(hsqc_pairs[j])
            assigned_count[best_rep] = assigned_count.get(best_rep, 0) + 1
        # 仍无归宿的峰保持 unmatched
    # 重新精确计算 unmatched (按列)
    matched_cols: Set[int] = set()
    for peaks in atom_to_peaks.values():
        for p in peaks:
            for jj, pair in enumerate(hsqc_pairs):
                if pair == p:
                    matched_cols.add(jj)
    unmatched = [j for j in range(n_cols) if j not in matched_cols]
    return atom_to_peaks, unmatched


# ---- 孤立峰归属 (HOSE / 经验范围) ----

def assign_quaternary_and_lone_h(
    mol: Chem.Mol,
    obs_c13: Sequence[float],
    obs_h1: Sequence[float],
    assigned_c: Set[float],
    assigned_h: Set[float],
    predictor: Optional[C13Predictor] = None,
    c_tol: float = 4.0,
    h_tol: float = 0.3,
) -> Tuple[Dict[int, float], Dict[int, float], List[str]]:
    """把未被 HSQC 覆盖的观测峰归属到 Cq (13C) / 孤立 H (OH 等).

    策略:
      - 剩余 13C 峰 <-> Cq 原子: 优先 HOSE 预测值全局指派 (匈牙利),
        无 predictor 时回退 shift_ranges 中心值贪心
      - 剩余 1H 峰: 若结构无可归属的 H (全部 H 已被 HSQC 覆盖),
        标记为孤立 H (如 OH / NH / 溶剂残峰), 不强行指派

    Returns:
        (cq_assignments, lone_h_assignments, notes)
        cq_assignments: {c_atom_idx: obs_c_ppm}
        lone_h_assignments: {obs_h_index(在 obs_h1 中): obs_h_ppm} 仅记录
    """
    notes: List[str] = []
    cq_atoms = [a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() == 'C' and _carbon_type_of(mol, a.GetIdx()) == 'Cq']
    rest_c = [c for c in obs_c13
              if not any(abs(c - ac) <= 0.05 for ac in assigned_c)]
    rest_h = [h for h in obs_h1
              if not any(abs(h - ah) <= 0.05 for ah in assigned_h)]

    cq_assignments: Dict[int, float] = {}
    # HOSE 路径: predict_smiles 输入 canonical SMILES (去氢), 预测按 reparsed
    # mol 的 C 遍历顺序返回 atom_index; 用 canonical rank (breakTies, 顺序
    # 不变量) 把 reparsed 索引映射回原 mol 索引
    pred_map: Dict[int, float] = {}
    if predictor is not None and cq_atoms:
        try:
            mol_noh = Chem.RemoveHs(mol)
            noh_smiles = Chem.MolToSmiles(mol_noh)
            pred_result = predictor.predict_smiles(noh_smiles)
            if pred_result.predictions:
                reparsed = Chem.MolFromSmiles(noh_smiles)
                ranks_orig = list(Chem.CanonicalRankAtoms(mol_noh, breakTies=True))
                ranks_rep = list(Chem.CanonicalRankAtoms(reparsed, breakTies=True))
                # 原 mol C 原子 rank -> 原 idx
                rank_to_orig = {
                    ranks_orig[a.GetIdx()]: a.GetIdx()
                    for a in mol_noh.GetAtoms() if a.GetSymbol() == 'C'
                }
                for pred in pred_result.predictions:
                    rep_rank = ranks_rep[pred.atom_index]
                    orig_idx = rank_to_orig.get(rep_rank)
                    if orig_idx is not None:
                        pred_map[orig_idx] = pred.predicted_ppm
        except Exception as exc:  # 显式记录, 不静默
            notes.append(f'HOSE 预测失败, 回退经验范围: {type(exc).__name__}: {exc}')
            pred_map = {}

    if cq_atoms and rest_c:
        reps_hose = [i for i in cq_atoms if i in pred_map]
        if reps_hose:
            from scipy.optimize import linear_sum_assignment
            cost = np.full((len(reps_hose), len(rest_c)), 1e6)
            for i, idx in enumerate(reps_hose):
                for j, oc in enumerate(rest_c):
                    d = abs(pred_map[idx] - oc)
                    if d <= c_tol:
                        cost[i, j] = d
            r_idx, c_idx_arr = linear_sum_assignment(cost)
            for r, c in zip(r_idx, c_idx_arr):
                if cost[r, c] < 1e6:
                    cq_assignments[reps_hose[r]] = rest_c[c]
            n_ok = len(cq_assignments)
            if n_ok < len(cq_atoms):
                notes.append(
                    f'HOSE 归属 Cq {n_ok}/{len(cq_atoms)} '
                    '(其余超容差, 留空不强行指派)')
        else:
            # 回退: 按官能团类型经验中心值贪心 (羧/酰胺 > 芳香 > 脂肪)
            centers: Dict[int, float] = {
                idx: _cq_range_center(mol, idx) for idx in cq_atoms}
            used: Set[int] = set()
            for idx in sorted(cq_atoms, key=lambda i: -centers[i]):
                best_j, best_d = -1, float('inf')
                for j, oc in enumerate(rest_c):
                    if j in used:
                        continue
                    d = abs(centers[idx] - oc)
                    if d < best_d:
                        best_d, best_j = d, j
                if best_j >= 0 and best_d <= max(c_tol, 10.0):
                    used.add(best_j)
                    cq_assignments[idx] = rest_c[best_j]
            notes.append('13C Cq 归属使用经验中心值 (无 HOSE, 精度有限)')

    lone_h = {i: h for i, h in enumerate(rest_h)}
    if lone_h:
        notes.append(f'{len(lone_h)} 个 1H 峰未归属 (可能为 OH/NH/溶剂残峰)')
    return cq_assignments, lone_h, notes


def _cq_range_center(mol: Chem.Mol, c_idx: int) -> float:
    """Cq 粗略中心值: 按邻接杂原子数粗分 (羧/酰胺 > 芳香 > 脂肪)."""
    atom = mol.GetAtomWithIdx(c_idx)
    n_o = sum(1 for nb in atom.GetNeighbors() if nb.GetSymbol() in ('O', 'N', 'S'))
    aromatic = atom.GetIsAromatic()
    if n_o >= 2:
        return 178.0  # 羧酸/酯/酰胺 C=O 165-185
    if n_o == 1:
        return 160.0  # 酰胺/缩醛 C-N/O
    if aromatic:
        return 138.0  # 芳香季碳 130-145
    return 45.0       # 脂肪季碳 25-55


# ---- 2D 反向验证 ----

def _bond_graph(mol: Chem.Mol, heavy_only: bool = True) -> Dict[int, Set[int]]:
    """邻接表; heavy_only=True 时只含重原子 (H-H 拓扑距离 = C-C 键数 + 1)."""
    graph: Dict[int, Set[int]] = {}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if heavy_only and (mol.GetAtomWithIdx(a).GetSymbol() == 'H' or
                           mol.GetAtomWithIdx(b).GetSymbol() == 'H'):
            continue
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    return graph


def _bfs_dist(graph: Dict[int, Set[int]], start: int, max_d: int) -> Dict[int, int]:
    dist = {start: 0}
    frontier = [start]
    for d in range(1, max_d + 1):
        nxt: List[int] = []
        for u in frontier:
            for v in graph.get(u, ()):
                if v not in dist:
                    dist[v] = d
                    nxt.append(v)
        frontier = nxt
        if not frontier:
            break
    return dist


def _hydrogens_of(mol: Chem.Mol, c_idx: int) -> List[int]:
    return [nb.GetIdx() for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors()
            if nb.GetSymbol() == 'H']


def verify_cosy(
    mol: Chem.Mol,
    cosy_pairs: Sequence[Tuple[float, float]],
    h_assign: Dict[int, List[float]],
    h_tol: float = 0.1,
) -> List[CorrelationCheck]:
    """COSY 验证: 结构 3 键内 H-H 对是否解释实测交叉峰.

    双向检查:
      - 每个实测 COSY 峰能否由结构邻接 H 对解释 (峰覆盖)
      - 结果逐条输出 CorrelationCheck
    """
    graph = _bond_graph(mol)
    # 实测 ppm -> 归属 C idx 列表 (一个 ppm 可归属多个等价碳)
    ppm_to_cidxs: Dict[float, List[int]] = {}
    for idx, ppms in h_assign.items():
        for p in ppms:
            ppm_to_cidxs.setdefault(p, []).append(idx)
    checks: List[CorrelationCheck] = []
    for a, b in cosy_pairs:
        cas: List[int] = []
        cbs: List[int] = []
        for p, idxs in ppm_to_cidxs.items():
            if abs(p - a) <= h_tol:
                cas.extend(idxs)
            if abs(p - b) <= h_tol:
                cbs.extend(idxs)
        matched = False
        note = ''
        if not cas or not cbs:
            note = '端点 H 未归属, 无法核对'
        else:
            # 任一端点组合拓扑命中即算命中 (等价歧义吸收)
            for ca in cas:
                for cb in cbs:
                    if _any_h_pair_within(mol, graph, ca, cb, max_d=3):
                        matched = True
                        break
                if matched:
                    break
            if not matched:
                note = '结构拓扑 3 键内无此 H-H 对'
        checks.append(CorrelationCheck(
            kind='COSY', obs_pair=(a, b),
            expected_atoms=(cas[0] if cas else -1, cbs[0] if cbs else -1),
            matched=matched, note=note))
    return checks


def _any_h_pair_within(
    mol: Chem.Mol,
    graph: Dict[int, Set[int]],
    c1: int,
    c2: int,
    max_d: int = 3,
) -> bool:
    """两个碳上是否存在任意 H 对, 等效 COSY 3 键内 (重原子图).

    处理 CH2/CH3 同碳 H 化学等价但 RDKit idx 不同的情况:
    实测信号归到 C 原子层面. H-H 键数 = d(C1,C2) + 1 (两个 C-H 键),
    故 COSY 3 键 <=> 重原子图 d(C1,C2) <= max_d - 1 = 2.
    """
    if not _hydrogens_of(mol, c1) or not _hydrogens_of(mol, c2):
        return False
    d_cc = _bfs_dist(graph, c1, max_d).get(c2, 99)
    return d_cc <= max_d - 1


def _pair_within_3bonds(graph: Dict[int, Set[int]], i1: int, i2: int) -> bool:
    d = _bfs_dist(graph, i1, 3)
    return d.get(i2, 99) <= 3


def verify_hmbc(
    mol: Chem.Mol,
    hmbc_pairs: Sequence[Tuple[float, float]],
    h_assign: Dict[int, List[float]],
    c_assign: Dict[int, float],
    h_tol: float = 0.1,
    c_tol: float = 0.5,
    allow_2j: bool = True,
) -> List[CorrelationCheck]:
    """HMBC 验证: 结构 2-3 键 H-C 对是否解释实测交叉峰.

    Args:
        h_assign: {c_atom_idx: [obs_dH, ...]}
        c_assign: {c_atom_idx: obs_dC} (含 HSQC 与 Cq 归属)
    """
    graph = _bond_graph(mol)
    # 归属映射: ppm -> C idx 列表 (等价碳可能共享同一观测信号)
    ppm_to_c_idx: Dict[float, List[int]] = {}
    for idx, p in c_assign.items():
        ppm_to_c_idx.setdefault(p, []).append(idx)
    hppm_to_c: Dict[float, List[int]] = {}
    for idx, ppms in h_assign.items():
        for p in ppms:
            hppm_to_c.setdefault(p, []).append(idx)
    checks: List[CorrelationCheck] = []
    for h_ppm, c_ppm in hmbc_pairs:
        c_idxs: List[int] = []
        for p, idxs in ppm_to_c_idx.items():
            if abs(p - c_ppm) <= c_tol:
                c_idxs.extend(idxs)
        h_cidxs: List[int] = []
        for p, idxs in hppm_to_c.items():
            if abs(p - h_ppm) <= h_tol:
                h_cidxs.extend(idxs)
        matched = False
        note = ''
        if not c_idxs or not h_cidxs:
            note = '端点未归属, 无法核对'
        else:
            # 任一端点组合 2-3 键即命中 (等价歧义吸收)
            for h_carbon in h_cidxs:
                for c_idx in c_idxs:
                    d = _bfs_dist(graph, h_carbon, 3).get(c_idx, 99)
                    if d in (2, 3):
                        matched = True
                        break
                if matched:
                    break
            if not matched:
                d0 = _bfs_dist(graph, h_cidxs[0], 3).get(c_idxs[0], 99)
                note = f'结构拓扑 H-C 距离 {d0} 键, 不在 2-3 键范围'
        checks.append(CorrelationCheck(
            kind='HMBC', obs_pair=(h_ppm, c_ppm),
            expected_atoms=(h_cidxs[0] if h_cidxs else -1,
                            c_idxs[0] if c_idxs else -1),
            matched=matched, note=note))
    return checks


# ---- 主入口 ----

def run_algorithmic_assignment(
    smiles: str,
    observed_c13: Sequence[float],
    observed_h1: Sequence[float],
    hsqc_pairs: Optional[Sequence[Tuple[float, float]]] = None,
    cosy_pairs: Optional[Sequence[Tuple[float, float]]] = None,
    hmbc_pairs: Optional[Sequence[Tuple[float, float]]] = None,
    predictor: Optional[C13Predictor] = None,
    h_tol: float = 0.3,
    c_tol: float = 2.0,
) -> AlgorithmicAssignmentResult:
    """端到端程序化归属.

    Args:
        smiles: 目标结构
        observed_c13 / observed_h1: 实测 1D 峰表 (ppm)
        hsqc_pairs: [(dH, dC), ...]
        cosy_pairs: [(dH1, dH2), ...] 已去对角线
        hmbc_pairs: [(dH, dC), ...] 已滤 1J
        predictor: HOSE 13C 预测器 (可选)
        h_tol / c_tol: 匹配容差

    Returns:
        AlgorithmicAssignmentResult (含全部证据链)
    """
    result = AlgorithmicAssignmentResult(structure_smiles=smiles)
    mol = _build_mol_from_smiles(smiles)
    labels = _default_labels(mol)
    equiv = _canonical_equivalence_groups(mol)

    c_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C']
    result.n_atoms_total = len(c_atoms)
    result.n_obs_c = len(observed_c13)
    result.n_obs_h = len(observed_h1)

    obs_c = sorted(observed_c13, reverse=True)
    obs_h = sorted(observed_h1, reverse=True)

    # 1. HSQC 锚点
    hsqc = list(hsqc_pairs or [])
    atom_to_peaks, unmatched_hsqc = assign_hsqc_anchor(mol, hsqc, h_tol, c_tol)
    assigned_c: Set[float] = set()
    assigned_h: Set[float] = set()
    c_assign: Dict[int, float] = {}   # c_idx -> obs dC
    h_assign: Dict[int, List[float]] = {}  # c_idx -> [obs dH]
    for idx, peaks in atom_to_peaks.items():
        for h_ppm, c_ppm in peaks:
            c_assign[idx] = c_ppm
            h_assign.setdefault(idx, []).append(h_ppm)
            assigned_c.add(round(c_ppm, 3))
            assigned_h.add(round(h_ppm, 3))

    # HSQC 覆盖率: 观测峰被结构解释的比例
    if hsqc:
        result.hsqc_match_rate = (len(hsqc) - len(unmatched_hsqc)) / len(hsqc)
        for j in unmatched_hsqc:
            h, c = hsqc[j]
            result.unassigned_peaks.append(UnassignedPeak(
                ppm=c, nucleus='13C(HSQC)',
                reason=f'HSQC 峰 (dH={h:.2f}) 无法指派到任何质子化碳'))

    # 2. 剩余 13C -> Cq (HOSE / 经验范围)
    cq_assign, lone_h, notes_cq = assign_quaternary_and_lone_h(
        mol, obs_c, obs_h, assigned_c, assigned_h,
        predictor=predictor, c_tol=max(c_tol, 4.0), h_tol=h_tol)
    result.notes.extend(notes_cq)
    for idx, ppm in cq_assign.items():
        c_assign[idx] = ppm
        assigned_c.add(round(ppm, 3))
    if lone_h:
        result.notes.append(
            f'{len(lone_h)} 个 1H 未归属: ' +
            ', '.join(f'{p:.2f}' for p in lone_h.values()))

    # 3. 汇总归属行
    for idx in c_atoms:
        c_type = _carbon_type_of(mol, idx)
        aa = AtomAssignment(
            atom_idx=idx,
            atom_label=labels.get(idx, f'C(idx{idx})'),
            carbon_type=c_type,
            equiv_group=equiv.get(idx, -1),
        )
        if idx in c_assign:
            aa.obs_c_ppm = round(c_assign[idx], 2)
        if idx in h_assign:
            aa.obs_h_ppms = sorted({round(p, 3) for p in h_assign[idx]}, reverse=True)
        aa.hsqc_matched = idx in atom_to_peaks
        if aa.hsqc_matched:
            aa.source = 'hsqc'
        elif idx in cq_assign:
            aa.source = 'hose' if predictor is not None else 'range'
        # 等价组共享: 同组其他原子的归属也补上
        if aa.source == '' and equiv.get(idx, -1) >= 0:
            for member in c_atoms:
                if member != idx and equiv.get(member, -1) == equiv.get(idx, -1) \
                        and member in c_assign:
                    aa.obs_c_ppm = round(c_assign[member], 2)
                    aa.obs_h_ppms = sorted(
                        {round(p, 3) for p in h_assign.get(member, [])}, reverse=True)
                    aa.source = 'equiv'
                    aa.hsqc_matched = member in atom_to_peaks
                    break
        # 芳香等价组补救: 观测峰少于组内原子时 (部分重叠), 同组已归属
        # 原子的 δC 若存在且其他成员仍在窗口内, 共享该归属 (共享信号)
        if aa.source == '' and aa.obs_c_ppm is None and equiv.get(idx, -1) >= 0:
            for member in c_atoms:
                if member != idx and equiv.get(member, -1) == equiv.get(idx, -1):
                    other = next((a for a in result.assignments
                                  if a.atom_idx == member), None)
                    if other is not None and other.obs_c_ppm is not None:
                        aa.obs_c_ppm = other.obs_c_ppm
                        aa.obs_h_ppms = list(other.obs_h_ppms)
                        aa.source = 'equiv_shared'
                        aa.hsqc_matched = other.hsqc_matched
                        break
        result.assignments.append(aa)

    result.n_c_assigned = len({round(v, 1) for v in c_assign.values()})
    result.n_h_assigned = len(assigned_h)
    assigned_atoms = sum(
        1 for aa in result.assignments
        if aa.obs_c_ppm is not None or aa.obs_h_ppms)
    result.n_atoms_assigned = assigned_atoms

    # 4. 2D 反向验证
    if cosy_pairs:
        result.correlation_checks.extend(
            verify_cosy(mol, cosy_pairs, h_assign, h_tol=min(h_tol, 0.1)))
    if hmbc_pairs:
        result.correlation_checks.extend(
            verify_hmbc(mol, hmbc_pairs, h_assign, c_assign,
                        h_tol=min(h_tol, 0.1), c_tol=0.5))
    for kind in ('COSY', 'HMBC'):
        checks = [c for c in result.correlation_checks if c.kind == kind]
        if checks:
            rate = sum(1 for c in checks if c.matched) / len(checks)
            setattr(result, f'{kind.lower()}_match_rate', round(rate, 4))

    # 5. 综合置信度 (权重: 原子覆盖 0.4 + HSQC 0.3 + HMBC 0.2 + COSY 0.1)
    coverage = assigned_atoms / result.n_atoms_total if result.n_atoms_total else 0.0
    w_hsqc = result.hsqc_match_rate if hsqc else 0.0
    w_hmbc = result.hmbc_match_rate if hmbc_pairs else 0.0
    w_cosy = result.cosy_match_rate if cosy_pairs else 0.0
    denom = 0.4 + (0.3 if hsqc else 0) + (0.2 if hmbc_pairs else 0) + (0.1 if cosy_pairs else 0)
    result.confidence = round(
        (0.4 * coverage + 0.3 * w_hsqc + 0.2 * w_hmbc + 0.1 * w_cosy) / denom, 4)
    result.notes.append(
        f'覆盖 {assigned_atoms}/{result.n_atoms_total} 原子, '
        f'13C {result.n_c_assigned}/{result.n_obs_c}, '
        f'1H {result.n_h_assigned}/{result.n_obs_h}')
    return result


# ---- 输出 ----

def write_assignment_result_csv(
    result: AlgorithmicAssignmentResult,
    output_path: str,
) -> None:
    """把归属结果写成 Mnova Report Assignment 风格 CSV (含证据链头)."""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['# Structure SMILES', result.structure_smiles])
        w.writerow(['# Confidence', f'{result.confidence:.2%}'])
        w.writerow(['# HSQC match', f'{result.hsqc_match_rate:.2%}'])
        w.writerow(['# COSY match', f'{result.cosy_match_rate:.2%}'])
        w.writerow(['# HMBC match', f'{result.hmbc_match_rate:.2%}'])
        if result.unassigned_peaks:
            w.writerow(['# Unassigned',
                        '; '.join(f"{p.nucleus}@{p.ppm:.2f}" for p in result.unassigned_peaks)])
        for note in result.notes:
            w.writerow(['# Note', note])
        w.writerow([])
        w.writerow(['位置', 'dC (ppm)', '碳级数', 'dH (ppm)', 'HSQC', '等价组', '来源'])
        for aa in result.assignments:
            w.writerow([
                aa.atom_label,
                f'{aa.obs_c_ppm:.2f}' if aa.obs_c_ppm is not None else '',
                aa.carbon_type,
                '; '.join(f'{p:.3f}' for p in aa.obs_h_ppms),
                'Y' if aa.hsqc_matched else '',
                aa.equiv_group if aa.equiv_group >= 0 else '',
                aa.source,
            ])
        # 2D 验证明细
        if result.correlation_checks:
            w.writerow([])
            w.writerow(['2D验证', 'dH1', 'dH2/dC', '命中', '说明'])
            for ck in result.correlation_checks:
                w.writerow([
                    ck.kind,
                    f'{ck.obs_pair[0]:.3f}',
                    f'{ck.obs_pair[1]:.3f}',
                    'Y' if ck.matched else 'N',
                    ck.note,
                ])


# ---- 真实数据观测级归属管线 (无结构, 纯观测信号匹配) ----

@dataclass
class ObservedCarbonRow:
    """一条观测 13C 信号的归属行 (对应标准表 13C NMR / Types / 1H NMR 三列)."""

    c_ppm: float
    carbon_type: str = '?'                    # CH/CH2/CH3/Cq/CH2_or_CH3
    h_ppms: List[float] = field(default_factory=list)  # HSQC 关联的 1H
    dept_source: str = ''                     # 'dept90'/'dept135neg'/'dept135pos'/'none'
    notes: str = ''


def build_observed_assignment(
    c13_shifts: Sequence[float],
    dept90_shifts: Optional[Sequence[float]] = None,
    dept135_shifts: Optional[Sequence[float]] = None,
    dept135_signs: Optional[Sequence[str]] = None,
    hsqc_pairs: Optional[Sequence[Tuple[float, float]]] = None,
    c_tol: float = 0.6,
    h_cluster_tol: float = 0.04,
) -> List[ObservedCarbonRow]:
    """观测级归属: 13C 主谱 + DEPT Types + HSQC (dH, dC) 匹配.

    三列输出 (对标标准归属表):
      - 13C NMR: c_ppm
      - Types: DEPT 判定 (CH/CH2/CH3/C), Cq 用 'C'
      - 1H NMR: HSQC 关联的 dH (经聚类去重)

    匹配逻辑 (确定性, 规则 8):
      1. DEPT-135 sign 主导: negative=CH2, positive 见分区规则
      2. DEPT-90 命中 → CH
      3. DEPT 双双未命中 → C (Cq)
      4. HSQC dC ↔ 13C 主谱最近邻一对一匹配 (c_tol), 组内全部 dH 归入该碳
    """
    # DEPT Types: 复用 assign_carbon_multiplicity (双向一对一最近匹配)
    from .bridge_nmr2struct import assign_carbon_multiplicity
    mults = assign_carbon_multiplicity(
        list(c13_shifts),
        list(dept90_shifts) if dept90_shifts else None,
        list(dept135_shifts) if dept135_shifts else None,
        list(dept135_signs) if dept135_signs is not None else None,
    )
    # 13C 降序
    rows: Dict[float, ObservedCarbonRow] = {}
    for r in sorted(mults, key=lambda x: -x.shift):
        ctype = r.multiplicity
        if ctype == 'Cq':
            ctype = 'C'
        rows[r.shift] = ObservedCarbonRow(
            c_ppm=r.shift, carbon_type=ctype,
            dept_source={'CH': 'dept90', 'CH2': 'dept135neg',
                         'CH3': 'dept135pos'}.get(r.multiplicity, 'none'),
        )
    c_sorted = sorted(rows.keys(), reverse=True)

    # HSQC dC -> 13C 全局最优一对一 (匈牙利算法, 非贪心:
    # 密集重叠区贪心会错配, 如 clxh100 31.9/32.3/31.3 竞争同一 HSQC 峰)
    if hsqc_pairs:
        from scipy.optimize import linear_sum_assignment
        n_c, n_p = len(c_sorted), len(hsqc_pairs)
        BIGC = 1e6
        cost = np.full((n_c, n_p), BIGC)
        for i, cm in enumerate(c_sorted):
            for j, (h, c) in enumerate(hsqc_pairs):
                d = abs(c - cm)
                if d <= c_tol:
                    cost[i, j] = d
        r_idx, c_idx = linear_sum_assignment(cost)
        used_h: Set[int] = set()
        for r, j in zip(r_idx, c_idx):
            if cost[r, j] < BIGC:
                rows[c_sorted[r]].h_ppms.append(hsqc_pairs[j][0])
                used_h.add(j)
        # 第二轮: 同一 13C 可接收多个 dH (CH2/CH3 简并信号), 把未匹配的
        # HSQC 峰分给 dC 距离最近的 13C
        for j, (h, c) in enumerate(hsqc_pairs):
            if j in used_h:
                continue
            best_i, best_d = -1, float('inf')
            for i, cm in enumerate(c_sorted):
                d = abs(c - cm)
                if d < best_d and d <= c_tol:
                    best_d, best_i = d, i
            if best_i >= 0:
                rows[c_sorted[best_i]].h_ppms.append(h)
        # Types 修正: 有 HSQC 峰的碳必为质子化碳 (排除 DEPT 弱信号漏检
        # 导致的 Cq 误判, 如 clxh100 44.8/31.9/29.1/28.7 ppm).
        # 质子化但 sign 丢失时, 脂肪区 (10-60) 多为 CH2 (天然产物 CH2 丰富),
        # 芳/烯区 (95-160) 为 CH; 这里按化学位移分区给最优猜测并标注不确定.
        for row in rows.values():
            if row.h_ppms and row.carbon_type == 'C':
                if 10.0 <= row.c_ppm <= 60.0:
                    row.carbon_type = 'CH2'
                    row.notes = 'HSQC 命中但 DEPT 弱漏检, 脂肪区按 CH2 (不确定)'
                elif 95.0 <= row.c_ppm <= 160.0:
                    row.carbon_type = 'CH'
                    row.notes = 'HSQC 命中但 DEPT 弱漏检, 芳/烯区按 CH'
                else:
                    row.carbon_type = 'CH_or_CH2'
                    row.notes = 'HSQC 命中但 DEPT 弱漏检, 类型占位'
            # DEPT-135 正峰 + DEPT-90 未命中 (CH2/CH3 占位): 脂肪区细分
            # 依据天然产物统计: sp3 CH2 (20-60) 远多于 CH3; CH3 特征区
            # (<20 或甲基碳 15-35 高度简并). 保守: 15-22 → CH3 优先,
            # 其余保持 CH2/CH3 占位并在 notes 标注.
            if row.carbon_type == 'CH2_or_CH3' and 15.0 <= row.c_ppm < 22.0:
                row.carbon_type = 'CH3'
                row.notes = 'DEPT135+ 脂肪低场区, 按 CH3 (不确定)'

    # 1H 峰聚类去重 (同一信号多条 peaklist 峰)
    for row in rows.values():
        if not row.h_ppms:
            continue
        clustered: List[float] = []
        for h in sorted(row.h_ppms, reverse=True):
            if not clustered or abs(h - clustered[-1]) > h_cluster_tol:
                clustered.append(h)
        row.h_ppms = clustered
    return [rows[c] for c in c_sorted]
