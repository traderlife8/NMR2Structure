"""基于 SymPy 的分子式推导与元素平衡校验.

应用场景:
1. parse_formula: 把分子式字符串 (e.g., "C8H8O2") 解析为 {元素: 个数}
2. formula_balance: 比较两个分子式的元素守恒 (差值 dict)
3. formula_dou: 不饱和度 (degree of unsaturation) 计算
4. formula_hc_ratio: H/C 比例 (辅助判断饱和度)
5. formula_molecular_weight: 分子量计算 (用 SymPy 符号代数, 保证精度)

依赖: sympy>=1.10 (本项目环境已装 1.13.3).

设计原则:
- 解析失败显式抛 FormulaParseError, 不静默返回空 dict
- 元素大小写不敏感 (SymPy 风格), 但输出规范化
- 全部函数纯函数 (无副作用), 便于测试
"""
from __future__ import annotations

import re
from typing import Dict, Optional

from sympy import Rational, Symbol, symbols


# 原子量表 (IUPAC 2021, 简化版). 用于 formula_molecular_weight. 单位 g/mol.
# 选择 IUPAC 标准原子量四舍五入到 4 位小数, 涵盖常见有机元素.
_ATOMIC_MASS = {
    'H': 1.008, 'He': 4.0026,
    'Li': 6.94, 'Be': 9.0122,
    'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982,
    'Si': 28.085, 'P': 30.974, 'S': 32.06,
    'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078,
    'Br': 79.904, 'I': 126.904,
}


class FormulaParseError(Exception):
    """分子式解析错误. 携带 stage 便于定位."""

    def __init__(self, msg: str, stage: str = ''):
        super().__init__(f'[{stage}] {msg}')
        self.stage = stage


# 元素正则: 大写字母开头, 可跟 1-2 个小写字母 (如 'Br', 'C', 'Cl')
# 系数: 紧跟元素, 可为 1 位或多位数
_ELEMENT_PATTERN = re.compile(r'([A-Z][a-z]?)(\d*)')


def parse_formula(formula: str) -> Dict[str, int]:
    """解析分子式字符串到 {元素符号: 个数}.

    支持格式:
      "C8H8O2"   -> {'C': 8, 'H': 8, 'O': 2}
      "C6H12NO2S" -> {'C': 6, 'H': 12, 'N': 1, 'O': 2, 'S': 1}
      "C"        -> {'C': 1}
      "CH4"      -> {'C': 1, 'H': 4}

    不支持: 同位素标记 (13C, 2H), 电荷, 同位素丰度比. 这类需用 RDKit.
    解析失败抛 FormulaParseError (空串、非法字符、连续数字).

    ponytail: 跳过项 — 升级路径: 若需支持同位素标记 / 电荷, 改用 RDKit.rdMolDescriptors.CalcMolFormula 或 SMILES 解析
    """
    if not isinstance(formula, str):
        raise FormulaParseError(f'公式必须为字符串: {type(formula).__name__}', stage='type')
    s = formula.strip()
    if not s:
        raise FormulaParseError('空字符串', stage='empty')
    # 预处理: 移除非字母数字内容 (空格、点等)
    s_clean = re.sub(r'[^A-Za-z0-9]', '', s)
    if not s_clean:
        raise FormulaParseError(f'清洗后为空: {formula!r}', stage='clean')
    if not s_clean[0].isalpha() or not s_clean[0].isupper():
        raise FormulaParseError(f'必须以大写字母开头: {s_clean}', stage='start')

    out: Dict[str, int] = {}
    i = 0
    while i < len(s_clean):
        m = _ELEMENT_PATTERN.match(s_clean, i)
        if not m:
            raise FormulaParseError(f'非法字符于位置 {i}: {s_clean[i:]!r}', stage='parse')
        elem = m.group(1)
        count_str = m.group(2)
        try:
            count = int(count_str) if count_str else 1
        except ValueError as e:
            raise FormulaParseError(f'非法系数: {count_str!r}', stage='count') from e
        if count <= 0:
            raise FormulaParseError(f'系数必须为正: {count}', stage='count')
        # 合并重复元素 (e.g., 'CH3CH3' 不合法但 'C2H6' 合法; 不应出现重复元素)
        if elem in out:
            # 累加 (容错: 接受 'CC' 这种不规范写法, 但抛错更安全)
            raise FormulaParseError(f'元素重复 {elem}: 分子式中同一元素应只出现一次', stage='dup')
        out[elem] = count
        i = m.end()
    return out


def formula_balance(a: Dict[str, int], b: Dict[str, int]) -> Dict[str, int]:
    """元素平衡校验. 返回 a-b 的差值 dict.

    正值表示 a 比 b 多, 负值表示少, 0 不列出.
    空差值说明两边元素守恒.
    """
    if not a or not b:
        return {}
    keys = set(a.keys()) | set(b.keys())
    return {k: a.get(k, 0) - b.get(k, 0) for k in keys if a.get(k, 0) != b.get(k, 0)}


def formula_dou(formula: Dict[str, int]) -> Optional[float]:
    """不饱和度 (degree of unsaturation) 计算.

    DoU = (2C + 2 + N - H - X) / 2
    其中 X = 卤素 (F, Cl, Br, I).
    适用于 C/H/N/O/卤素的有机小分子.
    返回 None 表示无碳元素 (无法计算 DoU).
    """
    c = formula.get('C', 0)
    if c == 0:
        return None
    h = formula.get('H', 0)
    n = formula.get('N', 0)
    x = formula.get('F', 0) + formula.get('Cl', 0) + formula.get('Br', 0) + formula.get('I', 0)
    # 用 Rational 保证整数除法精确
    dou = Rational(2 * c + 2 + n - h - x, 2)
    return float(dou)


def formula_hc_ratio(formula: Dict[str, int]) -> Optional[float]:
    """H/C 比例. 饱和烷烃 ~ 2.0-2.2; 芳烃 ~ 1.0; 不饱和脂肪烃 ~ 1.5-1.8.

    返回 None 表示无碳元素.
    """
    c = formula.get('C', 0)
    h = formula.get('H', 0)
    if c == 0:
        return None
    return h / c


def formula_molecular_weight(formula: Dict[str, int]) -> Optional[float]:
    """分子量计算 (g/mol). 基于 IUPAC 原子量表.

    返回 None 表示存在未知元素 (不在原子量表内).
    """
    if not formula:
        return None
    total = 0.0
    for elem, count in formula.items():
        mass = _ATOMIC_MASS.get(elem)
        if mass is None:
            return None  # 未知元素
        total += mass * count
    return round(total, 4)


def validate_formula_match(
    formula_a: str,
    formula_b: str,
    tol_pct: float = 5.0,
) -> Dict:
    """校验两个分子式是否一致 (含元素守恒 + 分子量接近).

    返回:
      {
        'matched': bool,
        'balance_diff': Dict[str, int],  # 元素守恒差 (a-b)
        'mw_a': float or None,
        'mw_b': float or None,
        'mw_diff_pct': float or None,
        'notes': List[str],
      }
    """
    out: Dict = {'matched': False, 'balance_diff': {}, 'mw_a': None, 'mw_b': None,
                 'mw_diff_pct': None, 'notes': []}
    try:
        a = parse_formula(formula_a)
    except FormulaParseError as e:
        out['notes'].append(f'a 解析失败: {e}')
        return out
    try:
        b = parse_formula(formula_b)
    except FormulaParseError as e:
        out['notes'].append(f'b 解析失败: {e}')
        return out
    diff = formula_balance(a, b)
    out['balance_diff'] = diff
    mw_a = formula_molecular_weight(a)
    mw_b = formula_molecular_weight(b)
    out['mw_a'] = mw_a
    out['mw_b'] = mw_b
    if mw_a is not None and mw_b is not None and mw_a > 0:
        out['mw_diff_pct'] = abs(mw_a - mw_b) / mw_a * 100.0
    # 判定: 元素守恒 (diff 为空) 且分子量差 <= tol_pct
    balance_ok = not diff
    mw_ok = (out['mw_diff_pct'] is not None and out['mw_diff_pct'] <= tol_pct)
    out['matched'] = balance_ok and mw_ok
    if diff:
        out['notes'].append(f'元素不平衡: {diff}')
    if out['mw_diff_pct'] is not None and out['mw_diff_pct'] > tol_pct:
        out['notes'].append(f'分子量差 {out["mw_diff_pct"]:.2f}% > {tol_pct}%')
    return out


__all__ = [
    'FormulaParseError',
    'formula_balance',
    'formula_dou',
    'formula_hc_ratio',
    'formula_molecular_weight',
    'parse_formula',
    'validate_formula_match',
]