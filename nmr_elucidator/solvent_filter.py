"""NMR 溶剂残余峰过滤 (Part 2.B).

常用氘代溶剂在 1H / 13C 维的残余峰位置表. 用于在 1D / 2D 数据准备阶段
排除非样品信号, 避免污染 retrieve 候选库与 step 中间产物.

数据来源: 各大仪器商出厂标准 + 论文常用参考值 (Gottlieb et al. 1997).

例: CDCl3 在 1H 谱中的残余 CHCl3 峰 ≈ 7.26 ppm, 在 13C 谱中 ≈ 77.16 ppm.
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Tuple


# 1H 维溶剂残余峰 (ppm)
# 注意: 每个溶剂列出所有常见残余峰 (如 CDCl3 单一, Pyr 三组)
SOLVENT_RESIDUAL_PEAKS_1H: dict = {
    'CDCl3': [7.26],
    'CHLOROFORM-D': [7.26],
    'CDCl3 + TMS': [7.26, 0.00],
    'DMSO-D6': [2.50],
    'CD3OD': [3.31],
    'METHANOL-D4': [3.31],
    'D2O': [4.79],
    'HOD (D2O)': [4.79],
    'Acetone-D6': [2.05],
    'C6D6': [7.16],
    'BENZENE-D6': [7.16],
    'CD2Cl2': [5.32],
    'DICHLOROMETHANE-D2': [5.32],
    'PYRIDINE-D5': [8.74, 7.58, 7.22],
    'Pyr': [8.74, 7.58, 7.22],
    'PYR': [8.74, 7.58, 7.22],
    'TFA-D': [11.50],
    'THF-D8': [3.58, 1.73],
    'ACETONITRILE-D3': [1.94],
    'CD3CN': [1.94],
    'DMSO + H2O': [2.50, 3.33],
    'H2O (HOD)': [1.56, 2.92],  # D2O 中 HOD 峰
}

# 13C 维溶剂残余峰 (ppm)
SOLVENT_RESIDUAL_PEAKS_13C: dict = {
    'CDCl3': [77.16],
    'CHLOROFORM-D': [77.16],
    'DMSO-D6': [39.52],
    'CD3OD': [49.00],
    'METHANOL-D4': [49.00],
    'Acetone-D6': [29.84, 206.26],
    'C6D6': [128.06],
    'BENZENE-D6': [128.06],
    'CD2Cl2': [53.84],
    'DICHLOROMETHANE-D2': [53.84],
    'PYRIDINE-D5': [149.90, 135.50, 123.50],
    'Pyr': [149.90, 135.50, 123.50],
    'PYR': [149.90, 135.50, 123.50],
    'THF-D8': [67.57, 25.37],
    'ACETONITRILE-D3': [1.32, 118.26],
    'CD3CN': [1.32, 118.26],
}

# 1H 键不区分大小写 + 别名映射
_NUCLEUS_ALIAS = {
    '1H': '1H',
    'PROTON': '1H',
    'PROTONS': '1H',
    '13C': '13C',
    'CARBON': '13C',
    '13CNMR': '13C',
}


def _normalize_solvent_key(solvent: str) -> Optional[str]:
    """归一化溶剂名到表键 (大小写/全角半角无关). 找不到返回 None."""
    if not solvent:
        return None
    s = solvent.strip().upper().replace(' ', '')
    # 直接查完全匹配
    for key in SOLVENT_RESIDUAL_PEAKS_1H.keys():
        if key.upper().replace(' ', '') == s:
            return key
    # 简化搜索: 移除常见后缀 (-D6, -D4, -D5, -D8, -D, -D3)
    aliases = ['-D6', '-D4', '-D5', '-D8', '-D3', '-D', '-D9']
    cleaned = s
    for a in aliases:
        cleaned = cleaned.replace(a, '')
    for key in SOLVENT_RESIDUAL_PEAKS_1H.keys():
        k_upper = key.upper().replace(' ', '')
        for a in aliases:
            k_upper_clean = k_upper.replace(a, '')
        if k_upper_clean == cleaned:
            return key
    return None


def get_solvent_peaks(solvent: str, nucleus: str) -> List[float]:
    """取指定溶剂在指定核 (1H / 13C) 的残余峰位置列表. 无返回 []."""
    nuc_key = _NUCLEUS_ALIAS.get(nucleus.upper(), nucleus)
    solvent_key = _normalize_solvent_key(solvent)
    if not solvent_key:
        return []
    table = (
        SOLVENT_RESIDUAL_PEAKS_1H if nuc_key == '1H'
        else SOLVENT_RESIDUAL_PEAKS_13C
    )
    return list(table.get(solvent_key, []))


def is_solvent_peak(
    ppm: float,
    nucleus: str,
    solvent: str = '',
    tolerance_ppm: float = 0.05,
) -> bool:
    """判断 ppm 是否在已知溶剂残余峰 ±tolerance_ppm 内.

    适用: 1D / 2D 谱中判断一个具体峰是否为溶剂残余.
    即使溶剂名未知, 也会查全表 (用空白溶剂扫描所有常用溶剂).
    """
    nuc_key = _NUCLEUS_ALIAS.get(nucleus.upper(), nucleus)
    table = (
        SOLVENT_RESIDUAL_PEAKS_1H if nuc_key == '1H'
        else SOLVENT_RESIDUAL_PEAKS_13C
    )
    if solvent:
        peaks = get_solvent_peaks(solvent, nuc_key)
    else:
        # 无 solvent 信息时, 全表合并
        peaks = []
        for v in table.values():
            peaks.extend(v)
    for sp in peaks:
        if abs(ppm - sp) <= tolerance_ppm:
            return True
    return False


def filter_peaks_by_solvent(
    peaks: List[Tuple[float, ...]],
    nucleus: str,
    solvent: str = '',
    tolerance_ppm: float = 0.05,
) -> List[Tuple[float, ...]]:
    """过滤掉溶剂峰, 返回非溶剂峰子集.

    输入: peaks 是 tuple 列表. 第 0 列视为 ppm. 其余列透明.
    输出: 移除溶剂峰后的子集.
    """
    out: List[Tuple[float, ...]] = []
    for p in peaks:
        if not p:
            out.append(p)
            continue
        ppm = float(p[0])
        if is_solvent_peak(ppm, nucleus, solvent, tolerance_ppm):
            continue
        out.append(p)
    return out


def annotate_solvent_residual(
    peaks: Sequence[Tuple[float, ...]],
    nucleus: str,
    solvent: str = '',
    tolerance_ppm: float = 0.05,
) -> List[bool]:
    """对每个峰返回一个 bool mask, True 表示溶剂残余峰.

    用途: 给 io_bruker 等模块在解析峰列表后追加 is_solvent_residual 字段.
    """
    out: List[bool] = []
    for p in peaks:
        if not p:
            out.append(False)
            continue
        ppm = float(p[0])
        out.append(is_solvent_peak(ppm, nucleus, solvent, tolerance_ppm))
    return out


def count_solvent_peaks(
    peaks: Iterable[float],
    nucleus: str,
    solvent: str = '',
    tolerance_ppm: float = 0.05,
) -> int:
    """统计 peaks 中的溶剂残余峰数量."""
    n = 0
    for ppm in peaks:
        if is_solvent_peak(float(ppm), nucleus, solvent, tolerance_ppm):
            n += 1
    return n
