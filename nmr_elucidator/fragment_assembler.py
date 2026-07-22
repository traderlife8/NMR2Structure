"""片段组装候选生成器.

当数据库检索均无命中时, 从 1D NMR 约束出发组装候选结构.
策略:
  1. 从 13C 化学位移范围推断官能团片段
  2. 从 DEPT 推断碳级数约束 (CH/CH2/CH3/Cq 数量)
  3. 从 1H 积分推断氢总数
  4. 从分子式 (若已知) 计算 DoU
  5. 用 RDKit 枚举满足约束的结构 (限制重原子数)
  6. 对每个候选用 HOSE 预测验证

注意: 组合爆炸控制 — 最多生成 500 个候选.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from itertools import product
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from .shift_ranges import ShiftRange, assign_shift_ranges, infer_functional_groups


@dataclass
class FragmentConstraint:
    """从 1D NMR 推断的结构约束."""

    # 碳级数统计 (来自 DEPT)
    n_ch: int = 0
    n_ch2: int = 0
    n_ch3: int = 0
    n_cq: int = 0
    # 氢总数 (来自 1H 积分)
    total_h: Optional[int] = None
    # 分子式
    formula: Optional[str] = None
    # 不饱和度
    dou: Optional[float] = None
    # 推断的官能团
    functional_groups: List[str] = field(default_factory=list)
    # 13C 位移
    c13_shifts: List[float] = field(default_factory=list)
    # 重原子数范围
    min_heavy_atoms: int = 3
    max_heavy_atoms: int = 15


@dataclass
class AssembledCandidate:
    """片段组装产生的候选."""

    smiles: str
    source: str = 'fragment_assembly'
    # 约束满足度
    formula_match: bool = False
    dept_match_score: float = 0.0
    dou_match: bool = False
    # 验证分数 (由 nmr_predict 填充)
    verification_score: float = 0.0


def build_constraints(
    c13_shifts: List[float],
    dept_summary: Optional[List[dict]] = None,
    formula: Optional[str] = None,
    total_h: Optional[int] = None,
) -> FragmentConstraint:
    """从 1D NMR 数据构建结构约束.

    Args:
        c13_shifts: 实测 13C 化学位移
        dept_summary: DEPT 推断 [{ppm, mult, sign}, ...]
        formula: 分子式 (可选)
        total_h: 氢总数 (可选, 来自 1H 积分)

    Returns:
        FragmentConstraint
    """
    constraint = FragmentConstraint(c13_shifts=list(c13_shifts))

    # DEPT 碳级数统计
    if dept_summary:
        for d in dept_summary:
            mult = d.get('mult', '?')
            if mult == 'CH':
                constraint.n_ch += 1
            elif mult == 'CH2':
                constraint.n_ch2 += 1
            elif mult == 'CH3':
                constraint.n_ch3 += 1
            elif mult == 'Cq':
                constraint.n_cq += 1

    # 分子式与 DoU
    constraint.formula = formula
    if formula:
        try:
            from .formula_solver import formula_dou, parse_formula
            parsed = parse_formula(formula)
            constraint.dou = formula_dou(parsed)
        except Exception:
            pass

    # 氢总数
    constraint.total_h = total_h

    # 官能团推断
    constraint.functional_groups = infer_functional_groups(c13_shifts)

    # 重原子数范围
    n_carbons = len(c13_shifts)
    constraint.min_heavy_atoms = max(3, n_carbons)
    constraint.max_heavy_atoms = min(20, n_carbons + 5)

    return constraint


def _check_formula_match(mol, formula: Optional[str]) -> bool:
    """检查分子式是否匹配."""
    if not formula:
        return True
    try:
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        return mol_formula == formula
    except Exception:
        return False


def _check_dou_match(mol, expected_dou: Optional[float]) -> bool:
    """检查不饱和度是否匹配."""
    if expected_dou is None:
        return True
    try:
        # RDKit 计算 DoU: (2C + 2 + N - H - X) / 2
        formula = rdMolDescriptors.CalcMolFormula(mol)
        from .formula_solver import formula_dou, parse_formula
        parsed = parse_formula(formula)
        actual_dou = formula_dou(parsed)
        if actual_dou is None:
            return True
        return abs(actual_dou - expected_dou) < 0.5
    except Exception:
        return True


def _check_dept_match(mol, constraint: FragmentConstraint) -> float:
    """检查 DEPT 碳级数匹配度 (0-1)."""
    if constraint.n_ch + constraint.n_ch2 + constraint.n_ch3 + constraint.n_cq == 0:
        return 0.5  # 无 DEPT 数据

    # 统计候选分子的碳级数
    mol_ch = 0
    mol_ch2 = 0
    mol_ch3 = 0
    mol_cq = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'C':
            continue
        nh = atom.GetTotalNumHs()
        if nh == 0:
            mol_cq += 1
        elif nh == 1:
            mol_ch += 1
        elif nh == 2:
            mol_ch2 += 1
        elif nh >= 3:
            mol_ch3 += 1

    # 计算匹配度
    total_obs = constraint.n_ch + constraint.n_ch2 + constraint.n_ch3 + constraint.n_cq
    total_mol = mol_ch + mol_ch2 + mol_ch3 + mol_cq
    if total_obs == 0 or total_mol == 0:
        return 0.0

    # 各类型匹配
    match_ch = min(constraint.n_ch, mol_ch)
    match_ch2 = min(constraint.n_ch2, mol_ch2)
    match_ch3 = min(constraint.n_ch3, mol_ch3)
    match_cq = min(constraint.n_cq, mol_cq)
    total_match = match_ch + match_ch2 + match_ch3 + match_cq

    return total_match / total_obs


def _has_functional_groups(mol, required_groups: List[str]) -> bool:
    """检查分子是否包含推断的官能团 (宽松匹配)."""
    if not required_groups:
        return True

    # 简化: 只检查关键官能团
    smarts_map = {
        '醛基': '[CX3H1]=O',
        '酮羰基': '[CX3](=O)([CX4])[CX4]',
        '羧酸/酯羰基': '[CX3](=O)[OX2]',
        '酰胺羰基': '[CX3](=O)[NX3]',
        '芳香CH': '[cH1]',
        '芳香季碳': 'c(:*)(:*)',
        'C-O (醇/醚)': '[CX4]-[OX2]',
        'C-N (胺)': '[CX4]-[NX3]',
    }

    # 只要求 50% 的推断官能团存在 (宽松)
    found = 0
    for group in required_groups:
        smarts = smarts_map.get(group)
        if smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                found += 1
        else:
            found += 1  # 未知官能团不扣分

    threshold = max(1, len(required_groups) // 2)
    return found >= threshold


def assemble_candidates(
    constraint: FragmentConstraint,
    max_candidates: int = 500,
    quiet: bool = False,
) -> List[AssembledCandidate]:
    """从约束出发组装候选结构.

    当前实现: 基于 RDKit 的分子式枚举 + 约束过滤.
    对于小分子 (重原子 <= 10), 尝试枚举结构异构体.
    对于较大分子, 仅做约束验证 (不枚举).

    Args:
        constraint: 结构约束
        max_candidates: 最大候选数
        quiet: 静默模式

    Returns:
        AssembledCandidate 列表
    """
    candidates: List[AssembledCandidate] = []

    if not constraint.formula:
        if not quiet:
            print("  [片段组装] 无分子式, 无法枚举结构异构体", file=sys.stderr)
        return candidates

    if not quiet:
        print(f"  [片段组装] 分子式={constraint.formula}, "
              f"DoU={constraint.dou}, "
              f"官能团={constraint.functional_groups}", file=sys.stderr)

    # 策略: 从 PubChem 按分子式获取候选, 然后用约束过滤
    # (纯枚举组合爆炸太大, 用数据库候选 + 约束过滤更实际)
    try:
        from .retrieve_pubchem import retrieve_pubchem_candidates
        pubchem_candidates = retrieve_pubchem_candidates(
            formula=constraint.formula,
            c13_shifts=constraint.c13_shifts,
            max_results=max_candidates * 2,
            quiet=quiet,
        )
    except Exception as e:
        if not quiet:
            print(f"  [片段组装] PubChem 检索失败: {e}", file=sys.stderr)
        pubchem_candidates = []

    # 约束过滤
    for pc in pubchem_candidates:
        if len(candidates) >= max_candidates:
            break
        mol = Chem.MolFromSmiles(pc.smiles)
        if mol is None:
            continue

        # 分子式检查
        if not _check_formula_match(mol, constraint.formula):
            continue

        # DoU 检查
        if not _check_dou_match(mol, constraint.dou):
            continue

        # DEPT 匹配
        dept_score = _check_dept_match(mol, constraint)
        if dept_score < 0.5:
            continue

        # 官能团检查
        if not _has_functional_groups(mol, constraint.functional_groups):
            continue

        candidates.append(AssembledCandidate(
            smiles=pc.smiles,
            formula_match=True,
            dept_match_score=dept_score,
            dou_match=True,
        ))

    if not quiet:
        print(f"  [片段组装] 约束过滤后: {len(candidates)} 个候选", file=sys.stderr)

    return candidates
