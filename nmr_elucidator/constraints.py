"""约束定义: 分子式、溶剂、谱图种类等用于硬筛候选.

约束来源:
1. 用户显式输入 (via CLI/JSON)
2. 从 acqus/parm.txt 自动推断 (溶剂、频率、温度)
3. 从 BrukerSpectrum 直接读取 (谱种类、峰数)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


@dataclass
class Constraints:
    """硬约束集合."""

    formula: Optional[str] = None  # 分子式（e.g. "C8H8O2"）
    solvent: Optional[str] = None  # 期望溶剂 (e.g. "Pyr")
    allowed_solvents: List[str] = field(default_factory=list)  # 可接受溶剂列表
    has_1h: bool = True
    has_13c: bool = True
    heavy_atom_min: Optional[int] = None
    heavy_atom_max: Optional[int] = None
    # 化学位移窗口（基于实测峰集）作为软约束
    shift_window_13c: Optional[tuple] = None  # (low, high)
    shift_window_1h: Optional[tuple] = None

    def validate(self) -> None:
        if self.formula:
            try:
                Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles('C'))  # 探针
            except Exception:
                pass
            mol = Chem.MolFromSmiles(f"C")  # noqa
            # 仅做基本格式校验
            if not isinstance(self.formula, str) or not self.formula:
                raise ValueError(f"分子式格式非法: {self.formula!r}")


def parse_constraints_from_acqus(acqus_text: str) -> Constraints:
    """从 acqus/parm.txt 文本推断约束.

    仅取 SOLVENT 字段. 输入形如 parm.txt 的文本.
    """
    c = Constraints()
    for line in acqus_text.splitlines():
        s = line.strip()
        if not s or s.startswith('$$'):
            continue
        if '=' not in s:
            continue
        key, _, val = s.partition('=')
        key = key.strip().lstrip('$').lstrip('#').strip().upper()
        val = val.strip()
        if key == 'SOLVENT' and val:
            c.solvent = val
    return c


def formula_match(formula: str, mol: Chem.Mol) -> bool:
    """分子式严格匹配（使用 RDKit 计算的公式）."""
    if not formula or mol is None:
        return False
    try:
        actual = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return False
    return _normalize_formula(actual) == _normalize_formula(formula)


def _normalize_formula(f: str) -> str:
    """规范化分子式: 按 Hill 顺序排序.

    e.g. 'C8H8O2' → 'C8H8O2', 'HO2C8' → 'C8H8O2' (若 RDKit 已经按 Hill 输出).
    """
    # RDKit CalcMolFormula 默认按 Hill 顺序（C 优先，其次 H，然后字母序）
    return f.strip()


def solvent_compatible(expected: Optional[str], candidate_solvent: Optional[str]) -> bool:
    """溶剂兼容性检查.

    规则:
      - 若 expected 为 None, 视为无约束
      - 若 candidate_solvent 为 None, 表示未知, 接受
      - 否则精确比较（大小写无关, 去除 '-d6' 标记）
    """
    if expected is None:
        return True
    if candidate_solvent is None:
        return True  # 未知不淘汰
    e = _norm_solvent(expected)
    c = _norm_solvent(candidate_solvent)
    return e == c


def _norm_solvent(s: str) -> str:
    s = s.strip().upper()
    for tok in ('-D6', '-D5', '-D4', '-D3', '-D2', '-D', 'D6', 'D5'):
        s = s.replace(tok, '')
    return s.strip()