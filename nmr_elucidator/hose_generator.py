"""HOSE-like 原子环境编码生成器.

基于 RDKit 的原子环境描述, 生成用于 13C 化学位移预测的层次化环境编码.
算法思路与 HOSE (Hierarchically Ordered Spherical Environment) 一致:
对目标碳原子, 按键距 (radius 1~6) 逐层编码周围化学环境.

本模块不依赖外部 hosegen 包, 完全基于 RDKit 实现.
编码格式: 每层用 '/' 分隔, 层内原子按 (原子序数, 键级, 氢数) 排序后拼接.
"""
from __future__ import annotations

from collections import deque
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Mol


# 键级编码
_BOND_ORDER_MAP = {
    Chem.BondType.SINGLE: '1',
    Chem.BondType.DOUBLE: '2',
    Chem.BondType.TRIPLE: '3',
    Chem.BondType.AROMATIC: 'a',
}

# 原子符号缩写 (保持编码紧凑)
_ATOM_SYMBOL_MAP = {
    6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S',
    17: 'Cl', 35: 'Br', 53: 'I', 5: 'B', 14: 'Si',
}


def _atom_token(atom) -> str:
    """生成原子 token: 元素符号 + 氢数 + 芳香性标记."""
    z = atom.GetAtomicNum()
    sym = _ATOM_SYMBOL_MAP.get(z, f'Z{z}')
    nh = atom.GetTotalNumHs()
    arom = '*' if atom.GetIsAromatic() else ''
    return f'{sym}{arom}H{nh}'


def _bond_token(bond) -> str:
    """生成键 token."""
    bt = bond.GetBondType()
    return _BOND_ORDER_MAP.get(bt, '?')


def generate_hose_code(mol: Mol, atom_idx: int, max_radius: int = 6) -> str:
    """为指定原子生成 HOSE-like 环境编码.

    算法:
      1. 从目标原子出发做 BFS, 记录每个原子到目标的键距 (sphere).
      2. 对每一层 (1..max_radius), 收集该层所有原子的 token 和连接键 token.
      3. 层内排序后拼接, 层间用 '/' 分隔.
      4. 若某层为空 (已到分子边界), 后续层省略.

    Args:
        mol: RDKit Mol 对象 (不含显式氢)
        atom_idx: 目标碳原子索引
        max_radius: 最大球层数 (默认 6)

    Returns:
        环境编码字符串, 如 'C*H0=1/C*H1-2,C*H0=1/...'
    """
    if atom_idx < 0 or atom_idx >= mol.GetNumAtoms():
        return ''

    # BFS 计算每个原子到目标的键距
    n_atoms = mol.GetNumAtoms()
    dist = [-1] * n_atoms
    dist[atom_idx] = 0
    queue = deque([atom_idx])
    # 记录到达每个原子的键 (用于编码连接关系)
    parent_bond: Dict[int, int] = {}  # atom_idx -> bond_idx

    while queue:
        cur = queue.popleft()
        if dist[cur] >= max_radius:
            continue
        atom = mol.GetAtomWithIdx(cur)
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtomIdx(cur)
            if dist[nbr] < 0:
                dist[nbr] = dist[cur] + 1
                parent_bond[nbr] = bond.GetIdx()
                queue.append(nbr)

    # 按层收集
    spheres: List[List[Tuple[str, str]]] = []  # [(atom_token, bond_token), ...]
    for r in range(1, max_radius + 1):
        layer_atoms = [i for i in range(n_atoms) if dist[i] == r]
        if not layer_atoms:
            break
        layer_tokens = []
        for ai in layer_atoms:
            a = mol.GetAtomWithIdx(ai)
            at = _atom_token(a)
            # 连接键
            bi = parent_bond.get(ai)
            bt = _bond_token(mol.GetBondWithIdx(bi)) if bi is not None else '?'
            layer_tokens.append((at, bt))
        # 排序保证确定性
        layer_tokens.sort()
        spheres.append(layer_tokens)

    # 拼接编码
    parts = []
    for layer in spheres:
        layer_str = ','.join(f'{at}-{bt}' for at, bt in layer)
        parts.append(layer_str)

    return '/'.join(parts)


def generate_hose_codes_for_carbons(
    mol: Mol, max_radius: int = 6
) -> Dict[int, str]:
    """为分子中所有碳原子生成 HOSE 编码.

    Args:
        mol: RDKit Mol 对象
        max_radius: 最大球层数

    Returns:
        {atom_idx: hose_code} 字典
    """
    codes: Dict[int, str] = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            codes[atom.GetIdx()] = generate_hose_code(mol, atom.GetIdx(), max_radius)
    return codes


def prepare_mol(smiles: str) -> Optional[Mol]:
    """从 SMILES 准备分子对象 (去显式氢, 感知芳香性).

    Args:
        smiles: SMILES 字符串

    Returns:
        RDKit Mol 或 None (解析失败)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.RemoveHs(mol)
    Chem.SetAromaticity(mol)
    return mol
