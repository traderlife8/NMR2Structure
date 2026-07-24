# -*- coding: utf-8 -*-
"""从已知 SMILES 生成化学家归属标号.

化学家 Excel 归属表的标号规则 (clxh100 数据归属表示例):
  - 苷元 (非糖) C 原子: 1, 2, 3, ..., N (按 RDKit canonical 顺序)
  - CH2 上的两个 H: α/β (3D 构象判定) 或 a/b (按 RDKit 索引 fallback)
  - 糖基 (Glc/Gal/Xyl 等): 26-O-Glc-1′, 2′, 3′, ..., 6a′, 6b′
  - 端基碳总是 1′, 6′ 是 CH2OH (拆 a/b)

本模块:
  - parse_smiles_to_mol(smiles) -> Chem.Mol (带 Hs)
  - assign_chemist_labels(mol) -> Dict[atom_idx, str]
    返回 RDKit atom index → 化学家标号 的映射

核心步骤:
  1. 检测糖基 (β-D-Glc 等 SMARTS) → 标记哪些 C 是糖
  2. 给非糖 C 按 RDKit canonical 顺序编号 1..N
  3. 给糖 C 按 RDKit canonical 顺序编号 (端基碳=1′)
  4. CH2 上的 H 拆 α/β (3D conformer 优先)
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem


# 糖基识别: 不用 SMARTS (RDKit substruct 阈值敏感), 用环结构启发式
# 含 1 个 O + 5 个 C 的 6 元环 → 吡喃糖环
_GLYCAN_DEFAULT_NAME = 'Glc'


def _find_glycan_rings(mol: Chem.Mol) -> List[Tuple[List[int], str]]:
    """找出分子中所有吡喃糖环 (含 1 O 的 6 元环).

    Returns:
        [(ring_atom_indices, glycan_name), ...]
        ring_atom_indices: 6 个环原子 indices (含 1 个 O + 5 个 C)
        glycan_name: 糖名 (简化: 统一用 'Glc'; 后续可根据子结构细化)
    """
    rings = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) != 6:
            continue
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        n_o = sum(1 for a in ring_atoms if a.GetSymbol() == 'O')
        n_c = sum(1 for a in ring_atoms if a.GetSymbol() == 'C')
        if n_o != 1 or n_c != 5:
            continue
        # 验证: 环上至少有 1 个 C 邻接环外 O (苷元侧链或 CH2OH)
        has_external = False
        for a in ring_atoms:
            if a.GetSymbol() != 'C':
                continue
            for nb in a.GetNeighbors():
                if nb.GetIdx() not in ring and nb.GetSymbol() == 'O':
                    has_external = True
                    break
            if has_external:
                break
        if not has_external:
            continue
        rings.append((list(ring), _GLYCAN_DEFAULT_NAME))
    return rings


@dataclass
class AtomLabel:
    """RDKit 原子 → 化学家标号 + 元信息."""

    c_idx: int = -1                     # RDKit atom index (C 原子, 移除 H 后)
    h_indices: List[int] = field(default_factory=list)  # RDKit H atom indices (AddHs 后)
    chemist_label: str = ''             # 化学家标号, 如 '1', '1α', '26-O-Glc-1′'
    is_glycan: bool = False             # 是否为糖基 C
    glycan_name: str = ''               # 糖名 (Glc/Gal/Xyl/Rha)
    is_anomeric: bool = False           # 是否为端基碳 (1′)
    carbon_type: str = ''               # CH/CH2/CH3/C (从 RDKit H 邻居数推)
    parent_no: str = ''                 # 所属苷元/糖基编号 (如 '1', 'Glc')
    suffix: str = ''                    # α/β/a/b 后缀 (仅 CH2 用)


@dataclass
class ChemistLabels:
    """完整化学家标号映射."""

    mol: object = None                            # RDKit Mol (带 H)
    smiles: str = ''
    c_labels: Dict[int, AtomLabel] = field(default_factory=dict)  # c_idx → AtomLabel
    h_to_label: Dict[int, str] = field(default_factory=dict)     # h_idx → 完整标号 (如 '1α')
    c_to_label: Dict[int, str] = field(default_factory=dict)     # c_idx → C 标号 (如 '1', 'Glc-1′')
    rows_by_label: Dict[str, AtomLabel] = field(default_factory=dict)  # 标号 → AtomLabel
    ordered_labels: List[str] = field(default_factory=list)        # 按化学家顺序排列的标号


def parse_smiles_to_mol(smiles: str, add_hs: bool = True) -> Optional[Chem.Mol]:
    """解析 SMILES → RDKit Mol.

    失败返回 None. 默认 AddHs (让 H 显式存在, 方便 α/β 拆分).
    """
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if add_hs:
            mol = Chem.AddHs(mol)
        return mol
    except Exception:
        return None


def _is_carbohydrate_ring(mol: Chem.Mol, c_idx: int) -> Tuple[bool, str]:
    """检测 c_idx 是否属于吡喃糖环 (含 O 的 6 元环 或 5 元环).

    Returns:
        (is_ring_oxygen_adjacent, glycan_name)
        is_ring_oxygen_adjacent: 是否为糖环上的碳 (邻接环内 O)
        glycan_name: 推测糖名 (Glc/Gal/Xyl/Rha), 未知则空
    """
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetSymbol() != 'C':
        return False, ''

    # 判断该 C 是否属于含 O 的环
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if c_idx not in ring:
            continue
        # 环内是否含 O
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(a.GetSymbol() == 'O' for a in ring_atoms):
            ring_size = len(ring)
            # 6 元环含 O → 吡喃糖; 5 元环含 O → 呋喃糖
            if ring_size == 6:
                # 区分糖: 用 SMARTS (Glc vs Gal vs Xyl 等)
                # 简化: 含 CH2OH (C 邻接 O-H 但不在环内) → 己糖
                for name, smarts, prefix in _GLYCAN_SUBSTRUCTURES:
                    pat = Chem.MolFromSmarts(smarts)
                    if pat and mol.HasSubstructMatch(pat):
                        return True, name
                return True, 'Gly'  # 通用
            elif ring_size == 5:
                return True, 'Fur'
    return False, ''


def _split_ch2_stereo(
    mol: Chem.Mol,
    c_idx: int,
    h_indices: List[int],
    use_3d: bool = True,
) -> List[Tuple[int, str]]:
    """把 CH2 的两个 H 拆成 'α'/'β' 或 'a'/'b' (与 assignment_table_v2 一致)."""
    if len(h_indices) != 2:
        return [(h, '') for h in h_indices]

    # 1. 3D conformer
    has_3d = mol.GetNumConformers() > 0 and any(
        not all(c == 0.0 for c in conformer.GetAtomPosition(c_idx))
        for conformer in mol.GetConformers()
    )

    if has_3d:
        conf = mol.GetConformer()
        z_vals = [(h, conf.GetAtomPosition(h).z) for h in h_indices]
        z_vals.sort(key=lambda x: x[1])
        return [(z_vals[0][0], 'α'), (z_vals[1][0], 'β')]

    if use_3d:
        try:
            from rdkit.Chem import AllChem
            m3d = Chem.Mol(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            cid_3d = AllChem.EmbedMolecule(m3d, params)
            if cid_3d >= 0:
                try:
                    AllChem.MMFFOptimizeMolecule(m3d, maxIters=200)
                except Exception:
                    pass
                conf = m3d.GetConformer()
                z_vals = [(h, conf.GetAtomPosition(h).z) for h in h_indices]
                z_vals.sort(key=lambda x: x[1])
                return [(z_vals[0][0], 'α'), (z_vals[1][0], 'β')]
        except Exception:
            pass

    # 2. 2D wedge
    wedge_labels = []
    wedge_h = non_wedge_h = None
    for h_idx in h_indices:
        bond = mol.GetBondBetweenAtoms(c_idx, h_idx)
        if bond is None:
            continue
        bond_dir = bond.GetBondDir()
        if bond_dir in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.ENDWEDGE):
            wedge_h = h_idx
        else:
            non_wedge_h = h_idx
    if wedge_h is not None and non_wedge_h is not None:
        return [(wedge_h, 'α'), (non_wedge_h, 'β')]

    # 3. Fallback: 按 RDKit 索引 a/b
    sorted_h = sorted(h_indices)
    return [(sorted_h[0], 'a'), (sorted_h[1], 'b')]


def _carbon_type(mol: Chem.Mol, c_idx: int) -> str:
    """根据 RDKit H 邻居数推碳级数."""
    atom = mol.GetAtomWithIdx(c_idx)
    n_h = sum(1 for nb in atom.GetNeighbors() if nb.GetSymbol() == 'H')
    return {0: 'C', 1: 'CH', 2: 'CH2', 3: 'CH3'}.get(n_h, '')


def _find_attachment_position(
    mol: Chem.Mol, glycan_atom_indices: List[int], full_mol: Chem.Mol,
) -> str:
    """糖基通过哪个苷元位置 O 连接到苷元 (如 '26-O-Glc' 中的 '26').

    策略: 找到糖基内连接外部的非环 O, 该 O 的另一个邻居 (非糖) 是苷元 C,
    取该 C 在 full_mol 中的化学家编号.
    """
    full_c_indices = set(a.GetIdx() for a in full_mol.GetAtoms() if a.GetSymbol() == 'C')
    for idx in glycan_atom_indices:
        atom = full_mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != 'O':
            continue
        # 该 O 是否在糖环内 (邻居都在 glycan_atom_indices 内) 还是连接外部
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        outside = [n for n in neighbors if n not in glycan_atom_indices]
        if len(outside) == 1:
            other_idx = outside[0]
            if other_idx in full_c_indices:
                # other_idx 是苷元 C → 它的化学家编号 (在 c_to_label 字典中查询)
                # 注: 调用方在 _assign_chemist_labels 时回填
                return str(other_idx)
    return ''


def assign_chemist_labels(
    smiles: str,
    use_3d: bool = True,
    glycan_substructures: Optional[List[Tuple[str, str, str]]] = None,
) -> ChemistLabels:
    """从 SMILES 生成完整化学家标号映射.

    Args:
        smiles: 已知结构 SMILES (canonical 或任意形式均可)
        use_3d: 是否用 RDKit 3D conformer 拆 CH2 立体氢
        glycan_substructures: 自定义糖基 SMARTS, None 用内置列表

    Returns:
        ChemistLabels 实例. mol/smiles 为空时返回空实例.
    """
    result = ChemistLabels(smiles=smiles)

    mol = parse_smiles_to_mol(smiles, add_hs=True)
    if mol is None:
        return result
    result.mol = mol

    glycan_pats = glycan_substructures or _GLYCAN_SUBSTRUCTURES

    # 1. 检测糖基 (substructure match 找出所有糖基原子 indices)
    # 注: 糖环 C + 环外 O 都在匹配集合中
    glycan_c_indices: set = set()  # 属于糖基的 C 原子 indices
    glycan_rings: List[Tuple[List[int], str]] = []  # [(ring_atoms, glycan_name), ...]
    for name, smarts, _ in glycan_pats:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        for match in matches:
            ring_atoms = list(match)
            # 该糖基的 C 原子
            ring_cs = [i for i in ring_atoms
                       if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
            # 避免重复 (同一糖被匹配多次)
            if any(set(ring_cs) & set(r) for r, _ in glycan_rings):
                continue
            glycan_rings.append((ring_atoms, name))
            glycan_c_indices.update(ring_cs)

    # 2. 找每个糖基的"附着位" (糖内连接外部的 O 的对侧苷元 C)
    # 用于命名: '26-O-Glc' 中的 '26'
    ring_attachment: Dict[int, str] = {}  # glycan ring 起点 C idx → 附着苷元编号
    for ring_atoms, gname in glycan_rings:
        # 找连接外部的 O
        attach_no = _find_attachment_position(mol, ring_atoms, mol)
        # 取端基 C (端基碳是环内与环外 O 直接相连的 C, 但环外 O 不在 ring_atoms 中)
        # 这里简化: 找糖环内第一个邻接 -O- 的 C 作为端基 C
        # SMARTS [OX2H] 是连接苷元的 O, 不在环里. 所以糖环里与外部 O 直接相邻的 C 是端基.
        anomeric_c = -1
        for idx in ring_atoms:
            a = mol.GetAtomWithIdx(idx)
            if a.GetSymbol() != 'C':
                continue
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == 'O' and nb.GetIdx() not in ring_atoms:
                    anomeric_c = idx
                    break
            if anomeric_c >= 0:
                break
        if anomeric_c >= 0:
            ring_attachment[anomeric_c] = attach_no

    # 3. 给非糖 C 按 RDKit canonical 顺序编号 1..N
    non_glycan_c_indices: List[int] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        if atom.GetIdx() in glycan_c_indices:
            continue
        non_glycan_c_indices.append(atom.GetIdx())

    # 按 canonical rank 排序, 缺失时按 RDKit idx 升序
    try:
        ranks = Chem.CanonicalRankAtoms(mol, includeHs=False, breakTies=True)
        non_glycan_c_indices.sort(key=lambda i: (ranks[i], i))
    except Exception:
        non_glycan_c_indices.sort()

    # 给非糖 C 编号 1, 2, 3, ...
    n = 1
    for c_idx in non_glycan_c_indices:
        n_h = sum(1 for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors() if nb.GetSymbol() == 'H')
        h_indices = sorted([nb.GetIdx() for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors()
                            if nb.GetSymbol() == 'H'])
        carbon_type = _carbon_type(mol, c_idx)
        chemist_no = str(n)

        # CH2 拆 α/β
        h_with_suffix = []
        if n_h == 2 and len(h_indices) == 2:
            splits = _split_ch2_stereo(mol, c_idx, h_indices, use_3d=use_3d)
            h_with_suffix = splits  # [(h_idx, suffix), ...]
        else:
            h_with_suffix = [(h, '') for h in h_indices]

        al = AtomLabel(
            c_idx=c_idx,
            h_indices=h_indices,
            chemist_label=chemist_no,
            is_glycan=False,
            carbon_type=carbon_type,
            parent_no='',
            suffix='',
        )
        result.c_labels[c_idx] = al
        result.c_to_label[c_idx] = chemist_no
        for h_idx, suf in h_with_suffix:
            if suf:
                full_lbl = f'{chemist_no}{suf}'
            else:
                full_lbl = chemist_no
            result.h_to_label[h_idx] = full_lbl
            # 在 AtomLabel 内记录 H 的标号 (用 h_to_label 反查时不区分 α/β 时取首)
        n += 1

    # 4. 给糖 C 编号: '26-O-Glc-1′', '2′', ..., '6a′', '6b′'
    for ring_atoms, gname in glycan_rings:
        ring_cs = [i for i in ring_atoms if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
        try:
            ranks = Chem.CanonicalRankAtoms(mol, includeHs=False, breakTies=True)
            ring_cs.sort(key=lambda i: (ranks[i], i))
        except Exception:
            ring_cs.sort()

        # 端基 C: 与环外 O 直接相连的那个
        anomeric_c = -1
        for idx in ring_cs:
            a = mol.GetAtomWithIdx(idx)
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == 'O' and nb.GetIdx() not in ring_atoms:
                    anomeric_c = idx
                    break
            if anomeric_c >= 0:
                break

        # 把端基 C 排到首位
        if anomeric_c in ring_cs:
            ring_cs.remove(anomeric_c)
            ring_cs.insert(0, anomeric_c)

        # 找附着苷元编号
        attach_no = ring_attachment.get(anomeric_c, '')
        if attach_no:
            parent_prefix = f'{attach_no}-O-{gname}'
        else:
            parent_prefix = gname

        # 端基 C 用 '1′', 其余按顺序 2′, 3′, ...
        # 端基 C 一般是 CH (1 个 H), 不拆
        # CH2OH (C6) 拆 a/b
        # 注意: 如果端基是 CH2OH (某些糖如鼠李糖 Rha), 则拆 a/b
        for k, c_idx in enumerate(ring_cs):
            n_h = sum(1 for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors() if nb.GetSymbol() == 'H')
            h_indices = sorted([nb.GetIdx() for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors()
                                if nb.GetSymbol() == 'H'])
            carbon_type = _carbon_type(mol, c_idx)
            is_anomeric = (k == 0)
            if is_anomeric:
                chemist_no = f'{parent_prefix}-1′'
                suffix_main = ''
            else:
                chemist_no = f'{parent_prefix}-{k + 1}′'
                suffix_main = ''
            h_with_suffix = []
            if n_h == 2 and len(h_indices) == 2:
                splits = _split_ch2_stereo(mol, c_idx, h_indices, use_3d=use_3d)
                h_with_suffix = splits
            else:
                h_with_suffix = [(h, '') for h in h_indices]

            al = AtomLabel(
                c_idx=c_idx,
                h_indices=h_indices,
                chemist_label=chemist_no,
                is_glycan=True,
                glycan_name=gname,
                is_anomeric=is_anomeric,
                carbon_type=carbon_type,
                parent_no=parent_prefix,
                suffix='',
            )
            result.c_labels[c_idx] = al
            result.c_to_label[c_idx] = chemist_no
            for h_idx, suf in h_with_suffix:
                if suf:
                    full_lbl = f'{chemist_no}{suf}'
                else:
                    full_lbl = chemist_no
                result.h_to_label[h_idx] = full_lbl

    # 5. 收集所有按"化学家顺序"的标号
    # 顺序: 非糖 C (按 1, 2, 3, ...) → 每个糖基 (按附着苷元编号, 然后端基 1′, 2′, ...)
    ordered = []
    # 非糖 C 的标号
    for c_idx in non_glycan_c_indices:
        al = result.c_labels[c_idx]
        n_h = sum(1 for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors() if nb.GetSymbol() == 'H')
        if n_h >= 2 and len(al.h_indices) >= 2:
            # CH2 → α/β 各自一行
            # 取 H 标号列表
            h_labels = [result.h_to_label.get(h, '') for h in al.h_indices]
            # 按 α 在前 (z 小), β 在后 (z 大). 已经有 h_to_label 的键名相同但顺序未定
            # 这里简化: 按 h_to_label 字符串排序 (α < β < a < b)
            ordered.extend(sorted(set(h_labels)))
        else:
            ordered.append(al.chemist_label)
    # 糖基
    for ring_atoms, gname in glycan_rings:
        ring_cs = [i for i in ring_atoms if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
        try:
            ranks = Chem.CanonicalRankAtoms(mol, includeHs=False, breakTies=True)
            ring_cs.sort(key=lambda i: (ranks[i], i))
        except Exception:
            ring_cs.sort()
        anomeric_c = -1
        for idx in ring_cs:
            a = mol.GetAtomWithIdx(idx)
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == 'O' and nb.GetIdx() not in ring_atoms:
                    anomeric_c = idx
                    break
            if anomeric_c >= 0:
                break
        if anomeric_c in ring_cs:
            ring_cs.remove(anomeric_c)
            ring_cs.insert(0, anomeric_c)
        for c_idx in ring_cs:
            al = result.c_labels[c_idx]
            n_h = sum(1 for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors() if nb.GetSymbol() == 'H')
            if n_h >= 2 and len(al.h_indices) >= 2:
                h_labels = [result.h_to_label.get(h, '') for h in al.h_indices]
                ordered.extend(sorted(set(h_labels)))
            else:
                ordered.append(al.chemist_label)
    result.ordered_labels = ordered
    for lbl in ordered:
        if lbl:
            # 取对应的 AtomLabel (CH2 情况下取首字母的 C label)
            for al in result.c_labels.values():
                if lbl.startswith(al.chemist_label):
                    result.rows_by_label[lbl] = al
                    break

    return result


def format_h_cell(h_ppm: Optional[float], h_mult: str = '', j_values: Optional[List[float]] = None) -> str:
    """格式化 1H 单元格: '1.16, m' / '3.92, dd (9.5, 6.5)'

    Args:
        h_ppm: 1H 化学位移
        h_mult: 多重性字符串 (s/d/t/q/m/dd/...)
        j_values: J 耦合常数列表 (Hz)

    Returns:
        格式化字符串. h_ppm 为 None 时返回 ''.
    """
    if h_ppm is None:
        return ''
    j_values = j_values or []
    s = f'{h_ppm:.2f}'
    if not h_mult and not j_values:
        return s
    if h_mult and j_values:
        j_str = ', '.join(f'{j:.1f}' for j in j_values[:3])
        return f'{s}, {h_mult} ({j_str})'
    if h_mult:
        return f'{s}, {h_mult}'
    # 只有 J
    j_str = ', '.join(f'{j:.1f}' for j in j_values[:3])
    return f'{s}, m ({j_str})'


def format_c_cell(c_ppm: Optional[float]) -> str:
    """格式化 13C 单元格: '38.1' (无 1H 列表示季碳)."""
    if c_ppm is None:
        return ''
    return f'{c_ppm:.1f}'


__all__ = [
    'AtomLabel',
    'ChemistLabels',
    'parse_smiles_to_mol',
    'assign_chemist_labels',
    'format_h_cell',
    'format_c_cell',
]