"""把 nmr_elucidator 读出的 Bruker 实测谱转成 NMR2Struct 单谱图推理输入格式.

NMR2Struct SingleNMRDataset 要求的格式:
  hnmr_file: 文本文件, 每行 "ppm intensity" (空格分隔)
  cnmr_file: 文本文件, 逗号分隔的 ppm 列表 (单行)
  hnmr_shifts / cnmr_shifts: 由 NMR2Struct 提供的 grid pickle

数据源: nmr_elucidator.io_bruker.read_bruker_folder() 读出的 BrukerSpectrum
      同一个项目, 不修改原始谱图

DEPT 整合 (化学位移无畸变极化转移):
  exp 22 pulse=zgdc       13C 主谱 (去耦, 含全部 C)
  exp 23 pulse=deptsp90   DEPT-90 (仅 CH 显示)
  exp 24 pulse=deptsp135  DEPT-135 (CH/CH3 向上, CH2 向下, 季碳消失)
  三者交叉可推断每个 13C 信号的碳类型 (CH / CH2 / CH3 / Cq).
  本模块仅产生 DEPT 推断结果, 实际喂入 NMR2Struct 模型时仅用 H+C 主谱;
  推断的碳类型作为"证据"附加到候选报告中.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .io_bruker import BrukerSpectrum, read_bruker_folder


# DEPT-135 峰方向枚举
DEPT_UNKNOWN = '?'
DEPT_CH = 'CH'        # 三级碳 (一氢)
DEPT_CH2 = 'CH2'      # 二级碳
DEPT_CH3 = 'CH3'      # 一级碳
DEPT_CQ = 'Cq'        # 季碳 (无氢)
DEPT_CH2_or_CH3 = 'CH2/CH3'  # DEPT-135 正峰且 sign 未识别时的占位 (排除 CH 但无法在 CH3 与芳/烯 CH 间区分)
ALL_DEPT = {DEPT_UNKNOWN, DEPT_CH, DEPT_CH2, DEPT_CH3, DEPT_CQ, DEPT_CH2_or_CH3}


@dataclass
class C13WithMultiplicity:
    """单条 13C 化学位移 + 推断的碳类型 (来自 DEPT)."""
    shift: float
    multiplicity: str  # 'CH' / 'CH2' / 'CH3' / 'Cq' / '?'
    source: str        # 'main' / 'dept90' / 'dept135'
    sign: str = 'unknown'  # 'positive' / 'negative' / 'unknown' (DEPT-135 方向)


def _filter_peaks_by_intensity(shifts, intensities, percentile: float = 0.20):
    """保留强度 >= 第 N 百分位的峰 (去噪声). 与 nmr_elucidator.__main__.py 同步."""
    if not shifts or not intensities or len(shifts) <= 5:
        return list(shifts)
    paired = sorted(zip(intensities, shifts), reverse=True)
    n = len(paired)
    cut = max(1, int(n * percentile))
    thresh = paired[cut - 1][0]
    return sorted(p for i, p in paired if i >= thresh)


def _write_hnmr_txt(out_path: str, hnmr_peaks: List[float], h_integrals: Optional[List[tuple]] = None) -> int:
    """写出 HNMR 输入文件: 每行 'ppm intensity' (空格分隔).

    hnmr_peaks: 实测化学位移列表 (ppm, 保留原始精度)
    h_integrals: 可选, [(low_ppm, high_ppm, integral_value), ...] 来自 Bruker
        若 peak 落在 integral 区间内, 用对应积分当强度; 否则 = 1.0
    """
    # 排序 integrals, 二分加速
    # NMR 积分区间端点顺序按显示习惯可能是反的 (从高场到低场, 即大→小),
    # 所以统一取 min/max, 避免漏匹配.
    iv_regions: list[tuple[float, float, float]] = []
    if h_integrals:
        for tup in h_integrals:
            if len(tup) == 3:
                lo, hi, val = tup
            else:
                center, val = tup
                # 没给区间时假定 ±0.05 ppm
                lo, hi = center - 0.05, center + 0.05
            lo, hi = min(lo, hi), max(lo, hi)
            iv_regions.append((float(lo), float(hi), float(val)))
        iv_regions.sort(key=lambda x: x[0])
    lines: list[str] = []
    for ppm in hnmr_peaks:
        intens = 1.0
        if iv_regions:
            for lo, hi, val in iv_regions:
                if lo - 0.01 <= ppm <= hi + 0.01:
                    intens = val
                    break
        lines.append(f'{ppm:.4f} {intens:.6f}')
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or '.', exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    return len(lines)


def _write_cnmr_txt(out_path: str, cnmr_peaks: List[float]) -> int:
    """写出 CNMR 输入文件: 逗号分隔 ppm 列表 (单行)."""
    if not cnmr_peaks:
        cnmr_peaks = []
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or '.', exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(','.join(f'{p:.4f}' for p in sorted(cnmr_peaks)) + '\n')
    return len(cnmr_peaks)


# 容忍窗口: 主谱峰位置 vs DEPT 峰位置最大允许差 (ppm)
DEPT_TOL_PPM = 5.0  # Part 3: 放宽到 5.0 ppm (clxh100 实测 DEPT vs 主谱偏差可达 ±4.5 ppm)


def _match_to_main(main_shifts: List[float], aux_shifts: List[float], tol: float = DEPT_TOL_PPM) -> set:
    """返回 aux_shifts 中与 main_shifts 任意峰距离 <= tol 的子集 (索引)."""
    if not aux_shifts or not main_shifts:
        return set()
    main_sorted = sorted(main_shifts)
    out: set = set()
    for i, p in enumerate(aux_shifts):
        lo, hi = p - tol, p + tol
        # 二分找 main_sorted 中落进 [lo, hi] 的
        import bisect
        left = bisect.bisect_left(main_sorted, lo)
        right = bisect.bisect_right(main_sorted, hi)
        if right > left:
            out.add(i)
    return out


def assign_carbon_multiplicity(
    main_shifts: List[float],
    dept90_shifts: Optional[List[float]] = None,
    dept135_shifts: Optional[List[float]] = None,
    dept135_signs: Optional[List[str]] = None,
    tol: float = DEPT_TOL_PPM,
    dept90_tol: Optional[float] = None,
) -> List[C13WithMultiplicity]:
    """根据 13C 主谱 + DEPT-90 + DEPT-135 推断每条 13C 信号的碳类型.

    优先级:
      1) DEPT-135 peaklist 自带正负号 (推荐): positive -> CH3, negative -> CH2
      2) 化学位移兜底: < 35 ppm -> CH3, 35-95 -> CH2, 95-170 -> CH2/CH3 占位

    规则:
      - 主谱有 & DEPT-90 有:        CH   (DEPT-90 仅显示 CH)
      - 主谱有 & DEPT-90 无 &
              DEPT-135 有 (正峰):   CH3  (DEPT-135 中正峰, 排除 CH)
      - 主谱有 & DEPT-90 无 &
              DEPT-135 有 (负峰):   CH2  (DEPT-135 中负峰, DEPT-90 不显示)
      - 主谱有 & DEPT-90 无 &
              DEPT-135 也无:         Cq  (季碳, DEPT 全不显示)
      - DEPT 数据缺失:              ?

    DEPT-90 容忍度默认比 DEPT-135 更严 (0.15 vs 0.30 ppm), 避免 DEPT-90
    噪声峰误命中, 漏掉 CH/CH2 正确分类.

    返回按化学位移降序 (与主谱同序) 的 C13WithMultiplicity 列表.
    """
    out: list[C13WithMultiplicity] = []
    if not main_shifts:
        return out

    main_sorted = sorted(main_shifts, reverse=True)
    if dept90_tol is None:
        dept90_tol = min(tol, 0.50)  # Part 3: 放宽到 0.5 ppm (不同实验间 ppm 校准差异)

    # 准备各 DEPT 在主谱中的命中索引
    idx90_in_main: dict[int, None] = {}
    if dept90_shifts:
        for p in dept90_shifts:
            for i, m in enumerate(main_sorted):
                if abs(p - m) <= dept90_tol:
                    idx90_in_main[i] = None
                    break
    # DEPT-135 命中: 记录 主谱索引 -> sign
    idx135_in_main: dict[int, str] = {}
    if dept135_shifts:
        signs = dept135_signs if dept135_signs is not None else ['unknown'] * len(dept135_shifts)
        for p, s in zip(dept135_shifts, signs):
            for i, m in enumerate(main_sorted):
                if abs(p - m) <= tol:
                    idx135_in_main[i] = s
                    break

    for i, p in enumerate(main_sorted):
        # DEPT-90 命中即确认为 CH (三级碳/芳 CH/烯 CH).
        # 仅当 DEPT-90 未命中时, 才用 DEPT-135 sign 区分 CH2 vs CH3.
        # 理由: 芳/烯 CH 在 DEPT-135 也是正峰, 与 CH3 不可区分;
        #       但 DEPT-90 仅显示 CH, 所以 DEPT-90 命中即排除 CH2/CH3.
        if i in idx90_in_main:
            mult = DEPT_CH
            sign = 'positive'
        else:
            d135_sign = idx135_in_main.get(i, 'unknown')
            if d135_sign == 'negative':
                mult = DEPT_CH2
                sign = 'negative'
            elif d135_sign == 'positive':
                # DEPT-135 正峰, DEPT-90 无 = 排除 CH
                # 化学位移兜底: < 35 -> CH3, 其它 -> CH2/CH3 占位
                if p < 35.0:
                    mult = DEPT_CH3
                else:
                    mult = DEPT_CH2_or_CH3
                sign = 'positive'
            else:
                # sign 'unknown' 实际混合了三种本质不同的状态, 必须逐一区分:
                #   A) DEPT-135 谱有数据 (dept135_shifts 不为 None) 且该碳未命中
                #      → 物理上 DEPT-90 无 + DEPT-135 无信号 → 排除 CH/CH2/CH3 → Cq 季碳
                #   B) DEPT-135 谱有数据, 该碳命中但 sign 缺失
                #      → 有信号 ≠ 季碳, 但 sign 未知无法区分 CH2/CH3 → CH2_or_CH3 占位
                #   C) DEPT-135 谱数据完全缺失 (调用方未提供)
                #      → 无 DEPT 信息 → 化学位移兜底 (向后兼容旧 API)
                # 之前实现把这三种状态混为一谈, 导致 80 ppm 等被错判 CH2.
                if dept135_shifts is not None and i not in idx135_in_main:
                    mult = DEPT_CQ
                    sign = 'unknown'
                elif i in idx135_in_main:
                    mult = DEPT_CH2_or_CH3
                    sign = 'unknown'
                else:
                    if p < 35.0:
                        mult = DEPT_CH3
                        sign = 'positive'
                    elif p < 95.0:
                        mult = DEPT_CH2
                        sign = 'negative'
                    else:
                        mult = DEPT_CH2_or_CH3
                        sign = 'unknown'
        out.append(C13WithMultiplicity(shift=p, multiplicity=mult, source='main', sign=sign))
    return out


def convert_bruker_folder_to_nmr2struct(
    bruker_folder: str,
    hnmr_out: str,
    cnmr_out: str,
    intensity_percentile: float = 0.20,
    use_dept: bool = True,
) -> dict:
    """把 Bruker 原始谱图转成 NMR2Struct 期望的两个文本文件.

    bruker_folder: 含 21/22/23/24 等子目录的 Bruker 实验目录
    hnmr_out: HNMR 输出 .txt 路径
    cnmr_out: CNMR 输出 .txt 路径
    intensity_percentile: 强度过滤百分位 (默认 0.20, 保留 top 20% 强度的峰)
    use_dept: 是否使用 DEPT 数据做碳类型推断 (写到返回 dict, 不影响 hnmr/cnmr)

    返回元数据 dict: {n_1h, n_13c, exp_h, exp_c, solvent, field_mhz, dept_summary, ...}
    dept_summary: [{ppm, mult}, ...] 按化学位移降序, 供候选报告引用.
    """
    specs = read_bruker_folder(bruker_folder)
    h_spec = next((s for s in specs if s.nucleus == '1H'), None)
    # 碳谱分类
    c_spec = next((s for s in specs if s.nucleus == '13C' and 'DEPT' not in s.nucleus), None)
    d90_spec = next((s for s in specs if s.nucleus == '13C-DEPT' and '90' in (s.pulse_program or '')), None)
    d135_spec = next((s for s in specs if s.nucleus == '13C-DEPT' and '135' in (s.pulse_program or '')), None)
    # 兜底: 如果没有按 pulse program 区分的 DEPT, 就按出现的顺序取前两个 DEPT 实验
    if use_dept and (d90_spec is None or d135_spec is None):
        depts = [s for s in specs if s.nucleus == '13C-DEPT']
        if d90_spec is None and len(depts) >= 1:
            d90_spec = depts[0]
        if d135_spec is None and len(depts) >= 2:
            d135_spec = depts[1]

    h_peaks_raw = [p.shift for p in h_spec.peaks] if h_spec else []
    h_intens_raw = [p.intensity for p in h_spec.peaks] if h_spec else []
    c_peaks_raw = [p.shift for p in c_spec.peaks] if c_spec else []

    # 强度过滤去噪声
    h_shifts = _filter_peaks_by_intensity(h_peaks_raw, h_intens_raw, intensity_percentile)
    c_shifts = _filter_peaks_by_intensity(c_peaks_raw, [p.intensity for p in c_spec.peaks] if c_spec else [], intensity_percentile)

    # 1H 积分 (传给 hnmr.txt 当强度)
    h_integrals: list[tuple] = []
    if h_spec:
        for region in h_spec.integrals:
            h_integrals.append((region.low, region.high, region.value))

    n_h = _write_hnmr_txt(hnmr_out, h_shifts, h_integrals)
    n_c = _write_cnmr_txt(cnmr_out, c_shifts)

    # DEPT 推断 (可选) - 使用全部 13C 峰 (不过滤), 避免丢失弱峰
    dept_summary: list[dict] = []
    if use_dept and c_spec and c_peaks_raw:
        d90_shifts = [p.shift for p in d90_spec.peaks] if d90_spec else None
        d135_shifts = [p.shift for p in d135_spec.peaks] if d135_spec else None
        # 直接从 d135 peaklist 读 sign (intensity > 0 -> positive, < 0 -> negative)
        d135_signs = [p.sign for p in d135_spec.peaks] if d135_spec else None
        c_with_mult = assign_carbon_multiplicity(
            c_peaks_raw, d90_shifts, d135_shifts, d135_signs,
        )
        for item in c_with_mult:
            dept_summary.append({
                'ppm': round(item.shift, 3),
                'mult': item.multiplicity,
                'sign': item.sign,
            })

    return {
        'n_1h_written': n_h,
        'n_13c_written': n_c,
        'exp_1h': h_spec.exp_no if h_spec else None,
        'exp_13c': c_spec.exp_no if c_spec else None,
        'exp_dept90': d90_spec.exp_no if d90_spec else None,
        'exp_dept135': d135_spec.exp_no if d135_spec else None,
        'solvent': h_spec.solvent if h_spec else (c_spec.solvent if c_spec else None),
        'field_mhz': h_spec.field_mhz if h_spec else (c_spec.field_mhz if c_spec else 0.0),
        'hnmr_file': os.path.abspath(hnmr_out),
        'cnmr_file': os.path.abspath(cnmr_out),
        'dept_summary': dept_summary,
        'dept_counts': {
            'CH': sum(1 for d in dept_summary if d['mult'] == DEPT_CH),
            'CH2': sum(1 for d in dept_summary if d['mult'] == DEPT_CH2),
            'CH3': sum(1 for d in dept_summary if d['mult'] == DEPT_CH3),
            'CH2/CH3': sum(1 for d in dept_summary if d['mult'] == DEPT_CH2_or_CH3),
            'Cq': sum(1 for d in dept_summary if d['mult'] == DEPT_CQ),
            '?': sum(1 for d in dept_summary if d['mult'] == DEPT_UNKNOWN),
        },
    }
