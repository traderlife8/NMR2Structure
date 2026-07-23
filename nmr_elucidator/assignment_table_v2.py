"""化学家 Excel 风格的归属表 (Part 3 Phase 2/3/4/5).

与 `correlation_2d.build_full_assignment_table` 的区别:
  - 按 H 排列 (每个氢一行), 不是按 C
  - 立体氢拆分 (CH2 → α/β 或 a/b), 借助 RDKit 标 stero 推断
  - 多重性 + J 整合 (1H 列显示 'd (5.0)' / 'dd (9.5, 6.5)')
  - COSY / HMBC / ROESY 关联按 H 索引回填
  - 糖基子结构自动识别 (A+B 方案, 见 Phase 5)

对标化学家论文 Table 1 模板:
    No. | 13C NMR | Types | 1H NMR (mult., J in Hz) | 1H-1H COSY | HMBC | ROESY | 备注
"""
from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem


# ---- 数据结构 ----

@dataclass
class AssignmentRowV2:
    """归属表的一行 (按 H 排列, 含立体氢拆分)."""

    # 位置编号 (含立体): '1α', '1β', '7α', '7β', '23a', '23b', '26-O-Glc-1′'
    no_label: str = ''
    # 13C 化学位移 (None 表示孤立 H)
    c_ppm: Optional[float] = None
    # 碳级数: 'CH'/'CH2'/'CH3'/'Cq'/''
    multiplicity: str = ''
    # 1H 化学位移
    h_ppm: Optional[float] = None
    # 1H 多重性 ('s', 'd', 't', 'q', 'dd', 'm' 等, 来自 step1b multiplet_analyses)
    h_mult: str = ''
    # J 耦合常数列表 (Hz, 与 h_mult 配套)
    j_values: List[float] = field(default_factory=list)
    # COSY 关联的氢 (按 H 标号)
    cosy_partners: List[str] = field(default_factory=list)
    # HMBC 关联的碳 (按 C 标号, 如 'H-4/C-2/3/5/6/10')
    hmbc_to_c: List[str] = field(default_factory=list)
    # ROESY 关联的氢 (按 H 标号)
    roesy_partners: List[str] = field(default_factory=list)
    # 备注: 'glycan' (糖基), 'aglycone' (苷元) 等
    notes: str = ''


@dataclass
class AssignmentTableV2:
    """完整归属表 (化学家 Excel 风格)."""

    rows: List[AssignmentRowV2] = field(default_factory=list)
    structure_smiles: str = ''
    confidence: float = 0.0
    unassigned_obs_c: List[float] = field(default_factory=list)
    unassigned_obs_h: List[float] = field(default_factory=list)
    # 化学家元信息 (用户输入, 可选)
    hresims: str = ''
    alpha: str = ''
    ir: str = ''


# ---- 立体氢拆分 ----

def _is_ch2_with_two_hydrogens(mol: Chem.Mol, c_idx: int) -> bool:
    """判断 RDKit 原子是否是 sp3 CH2 且带两个显式 H."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False
    if atom.GetTotalNumHs() != 2:
        return False
    if atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False
    return True


def _stereo_suffix_for_ch2_h(mol: Chem.Mol, c_idx: int, h_idx: int) -> str:
    """对 sp3 CH2 上的某个 H 返回立体后缀 'α'/'β' 或 'a'/'b'.

    策略: 用 RDKit 3D conformer (若已有) 或 2D wedge 推断立体.
    若不确定则按 H 索引顺序给 'a'/'b'.

    返回 'α' / 'β' (化学家常用希腊字母).
    """
    # 简单 fallback: 按 H 索引字典序
    return None  # 调用方负责根据 split 决定


def _split_ch2_to_stereo(
    mol: Chem.Mol,
    c_idx: int,
    h_indices: List[int],
    use_3d: bool = True,
) -> List[Tuple[int, str]]:
    """把 CH2 的两个 H 拆成 'α'/'β' 或 'a'/'b'.

    返回 [(h_idx, stereo_label), ...], stereo_label 是 'α'/'β' 或 'a'/'b'.

    优先级 (Part 3 Phase 8):
      1. 3D conformer (ETKDGv3 + MMFF): 按 Z 坐标判断 (z<0=α, z>0=β).
      2. 2D wedge bond: 按 wedge 方向判断 (朝前/朝后).
      3. Fallback: 按 H 索引顺序分 'a'/'b'.

    Args:
        use_3d: 是否调用 RDKit 生成 3D conformer (耗时但准确).
    """
    if len(h_indices) != 2:
        return [(h, '') for h in h_indices]

    # 1. 3D conformer 优先
    has_3d = mol.GetNumConformers() > 0 and any(
        not all(c == 0.0 for c in conformer.GetAtomPosition(c_idx))
        for conformer in mol.GetConformers()
    )

    if has_3d:
        conf = mol.GetConformer()
        z_vals = []
        for h_idx in h_indices:
            p = conf.GetAtomPosition(h_idx)
            z_vals.append((h_idx, p.z))
        z_vals.sort(key=lambda x: x[1])
        # z 较小者 → α (在环面下), 较大者 → β (在环面上)
        return [(z_vals[0][0], 'α'), (z_vals[1][0], 'β')]

    if use_3d:
        # 尝试生成 3D conformer
        try:
            from rdkit.Chem import AllChem
            m3d = Chem.Mol(mol)  # copy
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            cid_3d = AllChem.EmbedMolecule(m3d, params)
            if cid_3d >= 0:
                try:
                    AllChem.MMFFOptimizeMolecule(m3d, maxIters=200)
                except Exception:
                    pass
                # 用生成的 3D mol 判断
                conf = m3d.GetConformer()
                z_vals = []
                for h_idx in h_indices:
                    p = conf.GetAtomPosition(h_idx)
                    z_vals.append((h_idx, p.z))
                z_vals.sort(key=lambda x: x[1])
                return [(z_vals[0][0], 'α'), (z_vals[1][0], 'β')]
        except Exception:
            pass

    # 2. 2D wedge bond 检查
    wedge_labels = _split_by_wedge(mol, c_idx, h_indices)
    if wedge_labels:
        return wedge_labels

    # 3. Fallback: 按索引字典序分 'a'/'b'
    sorted_h = sorted(h_indices)
    return [(sorted_h[0], 'a'), (sorted_h[1], 'b')]


def _split_by_wedge(
    mol: Chem.Mol,
    c_idx: int,
    h_indices: List[int],
) -> List[Tuple[int, str]]:
    """用 RDKit wedge bond 信息拆 CH2 的两个 H.

    思路: 2D 里有 wedge 时, 'up' bond (实心框) 表示朝前, 'down' (虚线) 朝后.
    在环面上的氢通常标为 α (朝下 / 同侧), 反面为 β.
    这里用简化规则: 有 wedge 的 H 标 α, 无 wedge 的标 β.
    """
    c_atom = mol.GetAtomWithIdx(c_idx)
    wedge_h = None
    non_wedge_h = None
    for h_idx in h_indices:
        bond = mol.GetBondBetweenAtoms(c_idx, h_idx)
        if bond is None:
            continue
        bt = bond.GetBondType()
        # 2D wedge 通常 BondDir = BEGINWEDGE / ENDWEDGE
        bond_dir = bond.GetBondDir()
        if bond_dir in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.ENDWEDGE):
            wedge_h = h_idx
        else:
            non_wedge_h = h_idx
    if wedge_h is not None and non_wedge_h is not None:
        return [(wedge_h, 'α'), (non_wedge_h, 'β')]
    return []


# ---- No. 标号生成 ----

def _assign_atom_labels(mol: Chem.Mol) -> Dict[int, str]:
    """给每个 C 原子分配一个标号 (1, 2, 3, ...).

    顺序: 按 RDKit canonical rank (按 SMILES 输出顺序).
    CH2 上有两个 H 时拆成 1α/1β 或 1a/1b.
    """
    labels: Dict[int, str] = {}
    # 收集所有 C 原子
    c_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == 'C']
    # 按 atom index 排序
    c_atoms.sort(key=lambda a: a.GetIdx())
    n = 1
    for c_atom in c_atoms:
        c_idx = c_atom.GetIdx()
        # 找直接相连的 H (RemoveHs 之前的原子索引: 在 mol 中已经 AddHs 的 H 索引)
        h_indices = []
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                h_indices.append(neighbor.GetIdx())

        if len(h_indices) == 2 and _is_ch2_with_two_hydrogens(mol, c_idx):
            # 拆分 α/β 或 a/b
            splits = _split_ch2_to_stereo(mol, c_idx, h_indices)
            for h_idx, stereo in splits:
                labels[h_idx] = f'{n}{stereo}'
            n += 1
        else:
            # 不拆, 给 C 一个标号; H 也用同样标号
            labels[c_idx] = f'{n}'
            for h_idx in h_indices:
                labels[h_idx] = f'{n}'
            n += 1
    return labels


# ---- 多重性整合 ----

def _format_h_mult(mult: str, j_values: List[float]) -> str:
    """格式化 1H 多重性 + J 耦合.

    输入: mult='d', j_values=[5.0] → 'd (5.0)'
          mult='dd', j_values=[9.5, 6.5] → 'dd (9.5, 6.5)'
          mult='m', j_values=[] → 'm'
    """
    if not mult:
        return ''
    if not j_values:
        return mult
    j_str = ', '.join(f'{j:.1f}' for j in j_values)
    return f'{mult} ({j_str})'


def _lookup_multiplet_for_h(
    h_ppm: float,
    multiplets: List,  # List[MultipletAnalysis]
    h_tol: float = 0.05,
) -> Tuple[str, List[float]]:
    """从 step1b multiplet_analyses 中找最接近 h_ppm 的多重峰.

    返回 (mult, j_values).
    """
    if not multiplets or h_ppm is None:
        return ('', [])
    best = None
    best_diff = float('inf')
    for mp in multiplets:
        diff = abs(mp.h_ppm - h_ppm)
        if diff <= h_tol and diff < best_diff:
            best = mp
            best_diff = diff
    if best is None:
        return ('', [])
    # MultipletAnalysis 字段名: pattern (不是 multiplet_pattern), j_couplings_Hz
    j_values = list(getattr(best, 'j_couplings_Hz', []) or [])
    return (getattr(best, 'pattern', '') or '', j_values)


# ---- HMBC / COSY / ROESY 回填 ----

def _lookup_cosy_partners_for_h(
    h_ppm: float,
    spin_systems: List,  # List[SpinSystem]
    ch_pairs: List,  # List[CHPair]
    h_tol: float = 0.05,
) -> List[str]:
    """从 step2 自旋系统中找 h_ppm 所在的自旋系统, 列出其他 H (按 C-H 对的 assigned_h_label)."""
    if h_ppm is None or not spin_systems or not ch_pairs:
        return []
    # 找 h_ppm 所在的 spin system
    for ss in spin_systems:
        # 用任意一个 H_ppm 作为入口 (SpinSystem.h_ppm_列表 或 .h_ppms)
        ss_ppm = getattr(ss, 'h_ppm_list', None) or getattr(ss, 'h_ppms', None) or []
        # 找入口 H
        entry_h = None
        for sh in ss_ppm:
            if abs(sh - h_ppm) <= h_tol:
                entry_h = sh
                break
        if entry_h is None:
            continue
        # 找同一自旋系统内所有其他 H
        partners = []
        for sh in ss_ppm:
            if abs(sh - h_ppm) > h_tol and abs(sh - entry_h) <= h_tol * 4:  # 同一自旋系统
                # 找最匹配的 CHPair 的 H 标号
                for cp in ch_pairs:
                    if abs(cp.h_ppm - sh) <= h_tol:
                        lbl = cp.assigned_h_label or ''
                        # 去掉 H- 前缀以匹配 v2 no_label
                        if lbl.startswith('H-'):
                            lbl = lbl[2:]
                        if lbl:
                            partners.append(lbl)
                        break
        if partners:
            return partners[:5]  # 限制最多 5 个
    return []


def _lookup_hmbc_for_h(
    h_ppm: float,
    hmbc_connections: List,  # List[HMBCConnection]
    ch_pairs: List,  # List[CHPair]
    h_tol: float = 0.05,
    c_tol: float = 2.0,
) -> List[str]:
    """从 step3 HMBC 中找 h_ppm → C 的远程连接, 格式化为 'H-x/C-y/z/w'."""
    if h_ppm is None or not hmbc_connections or not ch_pairs:
        return []
    c_labels = []
    for conn in hmbc_connections:
        if abs(conn.h_ppm - h_ppm) <= h_tol:
            # 找 c_ppm 对应的 C 标号
            for cp in ch_pairs:
                if abs(cp.c_ppm - conn.c_ppm) <= c_tol:
                    lbl = cp.assigned_c_label or ''
                    if lbl.startswith('C-'):
                        lbl = lbl[2:]
                    if lbl and lbl not in c_labels:
                        c_labels.append(lbl)
                    break
    return c_labels[:5]


def _lookup_roesy_for_h(
    h_ppm: float,
    noe_contacts: List,  # List[NOEContact]
    ch_pairs: List,  # List[CHPair]
    h_tol: float = 0.05,
) -> List[str]:
    """从 step5 NOE 接触中找 h_ppm 的 ROESY 关联 H."""
    if h_ppm is None or not noe_contacts or not ch_pairs:
        return []
    partners = []
    for contact in noe_contacts:
        h1 = getattr(contact, 'h1_ppm', None)
        h2 = getattr(contact, 'h2_ppm', None)
        if h1 is None or h2 is None:
            continue
        other_h = None
        if abs(h1 - h_ppm) <= h_tol:
            other_h = h2
        elif abs(h2 - h_ppm) <= h_tol:
            other_h = h1
        if other_h is None:
            continue
        # 找最匹配的 CHPair 的 H 标号
        for cp in ch_pairs:
            if abs(cp.h_ppm - other_h) <= h_tol:
                lbl = cp.assigned_h_label or ''
                if lbl.startswith('H-'):
                    lbl = lbl[2:]
                if lbl and lbl not in partners:
                    partners.append(lbl)
                break
    return partners[:5]


# ---- 糖基识别 (Phase 5.A 子结构匹配) ----

_GLYCAN_SUBSTRUCTURES = [
    # (糖名, SMARTS, 标号前缀, 期望的 C 数量范围)
    ('β-D-Glc', '[C@@H]1([OX2H])[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O', 'Glc', (6, 7)),
    ('β-D-Gal', '[C@@H]1([OX2H])[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O', 'Gal', (6, 7)),
    ('α-L-Rha', '[C@@H]1([OX2H])[C@H](O)[C@@H](O)[C@H](O)C[C@@H]1O', 'Rha', (6, 7)),
    ('β-D-Xyl', '[C@@H]1([OX2H])[C@H](O)[C@@H](O)[C@H](O)O[C@@H]1O', 'Xyl', (5, 6)),
]

# 糖基 13C 化学位移指纹 (B 方案: 化学位移模式)
_GLYCAN_C13_FINGERPRINT = {
    # 典型 β-D-Glc: 端基碳 ~104-106, 其它碳 ~70-79
    # 字典: 键名 → (化学位移区间, 容差 ppm)
    '端基 C1': (104.0, 106.0, 1.5),
    'C2/C3/C4/C5': (70.0, 80.0, 2.0),
    'C6 (CH2OH)': (60.0, 64.0, 1.5),
}


def _detect_glycan_substructure(mol: Chem.Mol, c_ppm: float) -> Optional[str]:
    """检测分子是否含有糖基子结构 (A 方案).

    返回糖名 ('β-D-Glc' 等), 或 None.
    当前只在分子级别检测, 暂不细化到具体原子.
    """
    if mol is None:
        return None
    smi = Chem.MolToSmiles(mol)
    for name, smarts, prefix, (c_min, c_max) in _GLYCAN_SUBSTRUCTURES:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            continue
        if mol.HasSubstructMatch(pat):
            # 简化: 如果分子含 Glc 模式, 标注 'glycan'
            return name
    return None


def _detect_glycan_fingerprint(
    c_ppm_list: List[float],
) -> int:
    """用 13C 化学位移指纹检测糖基 (B 方案).

    返回匹配的糖基数量.
    """
    matched = 0
    has_c1 = False
    n_other = 0
    has_c6 = False
    for c in c_ppm_list:
        # C1 端基
        if 104.0 - 1.5 <= c <= 106.0 + 1.5:
            has_c1 = True
        # C2-C5 区
        elif 70.0 - 2.0 <= c <= 80.0 + 2.0:
            n_other += 1
        # C6 区
        if 60.0 - 1.5 <= c <= 64.0 + 1.5:
            has_c6 = True
    if has_c1 and n_other >= 3 and has_c6:
        matched = 1
    return matched


# ---- 主入口 ----

def build_assignment_table_v2(
    smiles: str,
    mol_block: str = '',
    observed_c13: Optional[List[float]] = None,
    observed_h1: Optional[List[float]] = None,
    interpretation_2d: Optional[object] = None,  # Interpretation2D
    multiplets: Optional[list] = None,  # List[MultipletAnalysis]
    h_tol: float = 0.05,
    c_tol: float = 2.0,
) -> AssignmentTableV2:
    """构建按 H 排列的归属表 (化学家 Excel 风格).

    Args:
        smiles: 候选 SMILES
        mol_block: mol block (可选, 用于 RDKit 解析; 若空则只用 SMILES)
        observed_c13: 实测 13C 化学位移列表
        observed_h1: 实测 1H 化学位移列表
        interpretation_2d: 2D 解读中间结果 (ch_pairs / spin_systems / hmbc_connections / noe_contacts)
        multiplets: step1b 多重峰分析列表
        h_tol: 1H 匹配容差 (ppm)
        c_tol: 13C 匹配容差 (ppm)
    """
    observed_c13 = observed_c13 or []
    observed_h1 = observed_h1 or []
    multiplets = multiplets or []

    table = AssignmentTableV2(structure_smiles=smiles)

    # 解析分子
    if mol_block:
        try:
            mol = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=False)
        except Exception:
            mol = None
    else:
        mol = None
    if mol is None and smiles:
        try:
            m = Chem.MolFromSmiles(smiles)
            if m is not None:
                mol = Chem.AddHs(m)
        except Exception:
            mol = None
    if mol is None:
        # 解析失败: 直接返回空表
        table.unassigned_obs_c = list(observed_c13)
        table.unassigned_obs_h = list(observed_h1)
        return table

    # 给每个原子分配标号
    labels = _assign_atom_labels(mol)

    # 糖基检测 (A: 子结构匹配, B: 化学位移指纹)
    glycan_name = _detect_glycan_substructure(mol, 0.0)
    glycan_count_fp = _detect_glycan_fingerprint(observed_c13)
    has_glycan = glycan_name is not None or glycan_count_fp > 0

    # 收集 2D 中间数据
    ch_pairs = getattr(interpretation_2d, 'ch_pairs', []) or []
    spin_systems = getattr(interpretation_2d, 'spin_systems', []) or []
    hmbc_connections = getattr(interpretation_2d, 'hmbc_connections', []) or []
    noe_contacts = getattr(interpretation_2d, 'noe_contacts', []) or []
    multiplets_list = list(multiplets or [])

    # 如果没有 2D 数据 (1D 流程), 用实测 ppm 的就近匹配填 13C 和 1H
    # strategy: 对每个 RDKit C 原子, 找 entry 的 spectrum_13c 中按 MOL block 顺序对应的 ppm.
    # 但因为 RDKit mol 保留 H (removeHs=False), atom_idx 不对应 spectrum_13c.
    # 简化: 从 candidate_spectrum_13c (候选 entry 的实测 13C 列表) 按 C 原子 index 顺序分配.
    candidate_13c_by_atom: Dict[int, float] = {}
    candidate_h_by_atom: Dict[int, float] = {}
    if interpretation_2d is None:
        # 传入 observed_c13 / observed_h1 已是实测 ppm, 不需要 entry 内置谱.
        # 这里仅靠 ch_pairs (空) 的 fallback -- 直接按 RDKit mol 顺序给 C 原子分配实测 13C.
        pass

    # 构建 RDKit C 原子 index (移除 H 后) 的顺序
    c_atom_indices_in_order = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_atom_indices_in_order.append(atom.GetIdx())
    # 同样 H 原子顺序 (含立体 H)
    h_atom_indices_in_order = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            h_atom_indices_in_order.append(atom.GetIdx())

    # 给 C 原子分配实测 13C (按 RDKit C index 顺序对应 sorted(observed_c13))
    if not ch_pairs and observed_c13:
        sorted_obs_c = sorted(observed_c13, reverse=True)
        for i, c_idx in enumerate(c_atom_indices_in_order):
            if i < len(sorted_obs_c):
                candidate_13c_by_atom[c_idx] = sorted_obs_c[i]

    # 给 H 原子分配实测 1H (按 RDKit H index 顺序对应 sorted(observed_h1))
    if not ch_pairs and observed_h1:
        sorted_obs_h = sorted(observed_h1, reverse=True)
        for i, h_idx in enumerate(h_atom_indices_in_order):
            if i < len(sorted_obs_h):
                candidate_h_by_atom[h_idx] = sorted_obs_h[i]

    # 构建 H → C 标号, C → 标号 的字典
    h_to_c_label: Dict[int, str] = {}
    h_to_h_label: Dict[int, str] = {}
    for atom_idx, label in labels.items():
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'C':
            h_to_c_label[atom_idx] = label
        elif atom.GetSymbol() == 'H':
            h_to_h_label[atom_idx] = label

    # 按 H 标号排序生成 rows
    rows: List[AssignmentRowV2] = []
    assigned_c_shifts: List[float] = []
    assigned_h_shifts: List[float] = []

    # 给每个 H 原子 (按标号顺序) 生成一行
    h_atoms_sorted = sorted(
        [a for a in mol.GetAtoms() if a.GetSymbol() == 'H'],
        key=lambda a: (
            # 提取标号数字部分排序
            int(re.match(r'(\d+)', labels.get(a.GetIdx(), '999')).group(1))
            if labels.get(a.GetIdx(), '').rstrip('αβab').isdigit()
            or re.match(r'(\d+)', labels.get(a.GetIdx(), '999'))
            else 999,
            labels.get(a.GetIdx(), ''),
        ),
    )

    for h_atom in h_atoms_sorted:
        h_idx = h_atom.GetIdx()
        h_label = labels.get(h_idx, '')
        if not h_label:
            continue

        # 找对应的 C
        c_atom = None
        for neighbor in h_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                c_atom = neighbor
                break
        if c_atom is None:
            continue
        c_idx = c_atom.GetIdx()
        c_label = labels.get(c_idx, '')

        # 找匹配的实测 1H (用 CHPair 中 assigned_h_label 优先, 否则用最近实测)
        matched_h_ppm = None
        for cp in ch_pairs:
            # cp.assigned_h_label 形如 'H-3', v2 no_label 形如 '3', 去掉前缀匹配
            cp_h_label = (cp.assigned_h_label or '').replace('H-', '')
            if cp_h_label and cp_h_label == h_label:
                matched_h_ppm = cp.h_ppm
                break
        if matched_h_ppm is None:
            # 1D fallback: 仅在2D 数据存在时使用；1D 流程下不填 (避免错误关联)
            # Part 3 设计: 1D 流程只填 13C 列, 1H 列留空表示 "需 2D HSQC 关联后才能确定"
            if h_idx in candidate_h_by_atom and interpretation_2d is not None:
                matched_h_ppm = candidate_h_by_atom[h_idx]

        # 找匹配的实测 13C
        matched_c_ppm = None
        for cp in ch_pairs:
            cp_c_label = (cp.assigned_c_label or '').replace('C-', '')
            if cp_c_label and cp_c_label == c_label:
                matched_c_ppm = cp.c_ppm
                break
        # 1D fallback: 用 candidate_13c_by_atom
        if matched_c_ppm is None and c_idx in candidate_13c_by_atom:
            matched_c_ppm = candidate_13c_by_atom[c_idx]

        # 碳级数 (从 RDKit 邻居 H 数推)
        n_h_on_c = sum(1 for nb in c_atom.GetNeighbors() if nb.GetSymbol() == 'H')
        if n_h_on_c == 0:
            mult = 'Cq'
        elif n_h_on_c == 1:
            mult = 'CH'
        elif n_h_on_c == 2:
            mult = 'CH2'
        elif n_h_on_c == 3:
            mult = 'CH3'
        else:
            mult = ''

        # 多重性 + J
        h_mult_str, j_values = _lookup_multiplet_for_h(matched_h_ppm, multiplets, h_tol=h_tol)
        h_mult_full = _format_h_mult(h_mult_str, j_values)

        # COSY / HMBC / ROESY
        cosy = _lookup_cosy_partners_for_h(matched_h_ppm, spin_systems, ch_pairs, h_tol=h_tol)
        hmbc = _lookup_hmbc_for_h(matched_h_ppm, hmbc_connections, ch_pairs,
                                   h_tol=h_tol, c_tol=c_tol)
        roesy = _lookup_roesy_for_h(matched_h_ppm, noe_contacts, ch_pairs, h_tol=h_tol)

        # 糖基标记
        notes = ''
        if has_glycan and matched_c_ppm is not None:
            if 60.0 - 1.5 <= matched_c_ppm <= 64.0 + 1.5 or \
               70.0 - 2.0 <= matched_c_ppm <= 80.0 + 2.0 or \
               104.0 - 1.5 <= matched_c_ppm <= 106.0 + 1.5:
                notes = 'glycan'

        row = AssignmentRowV2(
            no_label=h_label,
            c_ppm=round(matched_c_ppm, 2) if matched_c_ppm is not None else None,
            multiplicity=mult,
            h_ppm=round(matched_h_ppm, 3) if matched_h_ppm is not None else None,
            h_mult=h_mult_full,
            cosy_partners=cosy,
            hmbc_to_c=hmbc,
            roesy_partners=roesy,
            notes=notes,
        )
        rows.append(row)
        if matched_c_ppm is not None:
            assigned_c_shifts.append(matched_c_ppm)
        if matched_h_ppm is not None:
            assigned_h_shifts.append(matched_h_ppm)

    # 添加孤立 C (没有 H 的 Cq)
    c_atoms_in_mol = [a for a in mol.GetAtoms() if a.GetSymbol() == 'C']
    for c_atom in c_atoms_in_mol:
        c_idx = c_atom.GetIdx()
        n_h = sum(1 for nb in c_atom.GetNeighbors() if nb.GetSymbol() == 'H')
        if n_h > 0:
            continue  # 已经在 H 行覆盖
        c_label = labels.get(c_idx, '')
        if not c_label:
            continue
        matched_c_ppm = None
        for cp in ch_pairs:
            cp_c_label = (cp.assigned_c_label or '').replace('C-', '')
            if cp_c_label and cp_c_label == c_label:
                matched_c_ppm = cp.c_ppm
                break
        # 1D fallback
        if matched_c_ppm is None and c_idx in candidate_13c_by_atom:
            matched_c_ppm = candidate_13c_by_atom[c_idx]
        row = AssignmentRowV2(
            no_label=c_label,
            c_ppm=round(matched_c_ppm, 2) if matched_c_ppm is not None else None,
            multiplicity='Cq',
            h_ppm=None,
            h_mult='',
            cosy_partners=[],
            hmbc_to_c=[],
            roesy_partners=[],
            notes='Cq (no H)',
        )
        rows.append(row)
        if matched_c_ppm is not None:
            assigned_c_shifts.append(matched_c_ppm)

    # 按 No. 标号排序
    def _sort_key(r: AssignmentRowV2) -> Tuple[int, str]:
        m = re.match(r'(\d+)', r.no_label)
        n = int(m.group(1)) if m else 999
        # α 在 β 前, a 在 b 前
        suf = r.no_label[len(m.group(0)):] if m else r.no_label
        return (n, suf)

    rows.sort(key=_sort_key)
    table.rows = rows

    # 未归属的观测峰
    table.unassigned_obs_c = [
        c for c in observed_c13
        if not any(abs(c - ac) <= c_tol for ac in assigned_c_shifts)
    ]
    table.unassigned_obs_h = [
        h for h in observed_h1
        if not any(abs(h - ah) <= h_tol for ah in assigned_h_shifts)
    ]

    # 置信度: 已归属 H 占总观测 H 的比例
    if observed_h1:
        table.confidence = min(1.0, len(assigned_h_shifts) / max(len(observed_h1), 1))

    return table


# ---- 写入 CSV ----

def write_assignment_table_v2_csv(
    table: AssignmentTableV2,
    output_path: str,
    rank: int = 0,
    score: Optional[float] = None,
    extra_meta: Optional[Dict[str, str]] = None,
) -> None:
    """把 v2 归属表写成 CSV (贴近化学家 Excel)."""
    extra_meta = extra_meta or {}
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        # 元信息头
        w.writerow(['# Rank', rank])
        if score is not None:
            w.writerow(['# Score', f'{score:.4f}'])
        w.writerow(['# Structure SMILES', table.structure_smiles or ''])
        w.writerow(['# Confidence', f'{table.confidence:.2%}'])
        if table.hresims:
            w.writerow(['# HRESIMS', table.hresims])
        if table.alpha:
            w.writerow(['# [α]', table.alpha])
        if table.ir:
            w.writerow(['# IR', table.ir])
        for k, v in extra_meta.items():
            w.writerow([f'# {k}', v])
        if table.unassigned_obs_c:
            w.writerow(['# Unassigned 13C (ppm)',
                        '; '.join(f'{c:.3f}' for c in table.unassigned_obs_c)])
        if table.unassigned_obs_h:
            w.writerow(['# Unassigned 1H (ppm)',
                        '; '.join(f'{h:.3f}' for h in table.unassigned_obs_h)])
        w.writerow([])
        # 表头 (8 列, 贴近化学家 Excel)
        w.writerow(['No.', '13C NMR', 'Types', '1H NMR (mult., J in Hz)',
                    '1H-1H COSY', 'HMBC', 'ROESY', 'Notes'])
        # 数据行
        for r in table.rows:
            # 1H 列: 化学位移 + 多重性 + J, 格式 "1.16, m"
            if r.h_ppm is not None:
                if r.h_mult:
                    h_cell = f'{r.h_ppm:.2f}, {r.h_mult}'
                else:
                    h_cell = f'{r.h_ppm:.2f}'
            else:
                h_cell = ''
            w.writerow([
                r.no_label,
                f'{r.c_ppm:.1f}' if r.c_ppm is not None else '',
                r.multiplicity,
                h_cell,
                '/'.join(r.cosy_partners),
                '/'.join(r.hmbc_to_c),
                '/'.join(r.roesy_partners),
                r.notes,
            ])