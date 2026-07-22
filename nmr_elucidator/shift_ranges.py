"""13C 化学位移经验范围表.

基于文献 (Pretsch et al., Silverstein et al.) 的 13C NMR 化学位移范围,
用于从实测 13C 位移推断官能团片段.

每个条目: (名称, SMARTS, 13C 范围 low, 13C 范围 high, 典型碳类型)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class ShiftRange:
    """化学位移范围条目."""

    name: str  # 官能团名称
    smarts: str  # SMARTS 模式
    ppm_low: float  # 13C 范围下限
    ppm_high: float  # 13C 范围上限
    carbon_type: str  # 典型碳类型: Cq/CH/CH2/CH3
    priority: int = 0  # 优先级 (越高越先匹配)


# 13C 化学位移范围表 (按 ppm 降序排列)
SHIFT_RANGES: List[ShiftRange] = [
    # 羰基区 (160-220 ppm)
    ShiftRange('醛基', '[CX3H1](=O)', 190.0, 220.0, 'CH', 10),
    ShiftRange('酮羰基', '[CX3](=O)([CX4])[CX4]', 195.0, 220.0, 'Cq', 10),
    ShiftRange('羧酸/酯羰基', '[CX3](=O)[OX2]', 160.0, 185.0, 'Cq', 9),
    ShiftRange('酰胺羰基', '[CX3](=O)[NX3]', 160.0, 180.0, 'Cq', 9),
    ShiftRange('酸酐羰基', '[CX3](=O)O[CX3](=O)', 160.0, 175.0, 'Cq', 8),
    ShiftRange('氨基甲酸酯羰基', '[NX3][CX3](=O)[OX2]', 150.0, 165.0, 'Cq', 8),

    # 芳香/烯烃区 (100-160 ppm)
    ShiftRange('芳香碳(连杂原子)', 'c-[OX2, NX2]', 145.0, 165.0, 'Cq', 7),
    ShiftRange('芳香季碳', 'c(:*)(:*)', 125.0, 155.0, 'Cq', 5),
    ShiftRange('芳香CH', '[cH1]', 110.0, 140.0, 'CH', 6),
    ShiftRange('烯烃碳(连O)', '[CX3]=[CX3]-[OX2]', 135.0, 160.0, 'Cq', 6),
    ShiftRange('烯烃CH', '[CX3H1]=[CX3]', 100.0, 150.0, 'CH', 5),
    ShiftRange('烯烃季碳', '[CX3](=[CX3])(*)(*)', 110.0, 155.0, 'Cq', 5),
    ShiftRange('腈基', '[NX1]#[CX2]', 110.0, 125.0, 'Cq', 7),

    # 连杂原子脂肪区 (40-90 ppm)
    ShiftRange('缩醛/缩酮碳', '[CX4]([OX2])[OX2]', 90.0, 110.0, 'CH', 7),
    ShiftRange('糖苷碳(C1)', '[CX4H1]([OX2])-[OX2]', 90.0, 105.0, 'CH', 7),
    ShiftRange('C-O (醇/醚)', '[CX4]-[OX2]', 50.0, 90.0, 'CH', 4),
    ShiftRange('C-N (胺)', '[CX4]-[NX3]', 35.0, 65.0, 'CH', 4),
    ShiftRange('C-X (卤代)', '[CX4]-[F,Cl,Br,I]', 10.0, 80.0, 'CH2', 4),
    ShiftRange('环氧碳', '[CX4]1[OX2][CX4]1', 40.0, 65.0, 'CH', 5),

    # 饱和脂肪区 (0-50 ppm)
    ShiftRange('炔烃碳', '[CX2]#[CX2]', 65.0, 90.0, 'CH', 6),
    ShiftRange('环丙烷碳', '[CX4]1[CX4][CX4]1', -5.0, 15.0, 'CH2', 3),
    ShiftRange('N-甲基', '[NX3]-[CX4H3]', 25.0, 45.0, 'CH3', 3),
    ShiftRange('O-甲基', '[OX2]-[CX4H3]', 45.0, 60.0, 'CH3', 4),
    ShiftRange('S-甲基', '[SX2]-[CX4H3]', 10.0, 25.0, 'CH3', 3),
    ShiftRange('脂肪CH3', '[CX4H3]', 5.0, 30.0, 'CH3', 1),
    ShiftRange('脂肪CH2', '[CX4H2]', 15.0, 45.0, 'CH2', 1),
    ShiftRange('脂肪CH', '[CX4H1]', 20.0, 50.0, 'CH', 1),
    ShiftRange('脂肪季碳', '[CX4H0]', 25.0, 50.0, 'Cq', 2),
]


def assign_shift_ranges(
    c13_shifts: List[float],
    dept_types: Optional[List[str]] = None,
) -> List[Tuple[float, List[ShiftRange]]]:
    """将实测 13C 位移分配到可能的官能团范围.

    Args:
        c13_shifts: 实测 13C 化学位移列表
        dept_types: 对应的 DEPT 碳类型 (可选, 与 c13_shifts 等长)

    Returns:
        [(shift, [匹配的 ShiftRange 列表]), ...] 按 shift 降序
    """
    results: List[Tuple[float, List[ShiftRange]]] = []

    for i, shift in enumerate(c13_shifts):
        dept = dept_types[i] if (dept_types and i < len(dept_types)) else None
        matches: List[ShiftRange] = []

        for sr in SHIFT_RANGES:
            if sr.ppm_low <= shift <= sr.ppm_high:
                # 如果有 DEPT 信息, 进一步过滤
                if dept and dept != '?':
                    if dept == 'Cq' and sr.carbon_type != 'Cq':
                        continue
                    if dept in ('CH', 'CH2', 'CH3') and sr.carbon_type == 'Cq':
                        continue
                    if dept == 'CH' and sr.carbon_type in ('CH2', 'CH3'):
                        continue
                    if dept == 'CH2' and sr.carbon_type in ('CH', 'CH3'):
                        continue
                    if dept == 'CH3' and sr.carbon_type in ('CH', 'CH2'):
                        continue
                matches.append(sr)

        # 按优先级排序
        matches.sort(key=lambda x: x.priority, reverse=True)
        results.append((shift, matches))

    # 按 shift 降序
    results.sort(key=lambda x: x[0], reverse=True)
    return results


def infer_functional_groups(c13_shifts: List[float]) -> List[str]:
    """从 13C 位移推断存在的官能团 (名称列表).

    简化版: 只返回最可能的官能团名称.
    """
    groups: List[str] = []
    seen = set()

    for shift, matches in assign_shift_ranges(c13_shifts):
        for sr in matches:
            if sr.name not in seen:
                seen.add(sr.name)
                groups.append(sr.name)
                break  # 每个位移只取最高优先级的

    return groups
