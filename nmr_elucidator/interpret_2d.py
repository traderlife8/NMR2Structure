"""2D NMR 逐步解谱引擎.

按化学家真实解谱工作流 (对标 MestReNova_解谱指南.md 5.4 节) 逐步解读 2D NMR:
  步骤 1: HSQC → 建立 C-H 直接连接对 ("拼图块")
  步骤 2: COSY → 从已知 H 出发, 沿耦合链建立自旋系统
  步骤 3: HMBC → 从已归属 H 出发, 找 2-3 键外 C, 建立片段间连接
  步骤 4: 推断季碳 (13C 中无 HSQC 对应的峰)

每步产出结构化中间结果并写入独立 CSV 文件, 供化学家复核.

依赖: numpy, correlation_2d (峰提取), io_bruker_2d (Spectrum2D)
"""
from __future__ import annotations

import csv
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from .io_bruker_2d import Spectrum2D


# ---- 数据结构 ----

@dataclass
class CHPair:
    """HSQC 建立的一个 C-H 直接连接对 (一个'拼图块')."""

    h_ppm: float
    c_ppm: float
    multiplicity: str = 'unknown'  # 'CH'/'CH2'/'CH3'/'unknown'
    h_integral: Optional[float] = None  # 1H 积分值 (对应 H 个数)
    assigned_c_label: Optional[str] = None  # 归属标签 (如 'C-3')
    assigned_h_label: Optional[str] = None  # 归属标签 (如 'H-3')
    j_couplings: List[float] = field(default_factory=list)  # 该 H 的 J 耦合 (Hz)
    h_neighbors: List[float] = field(default_factory=list)  # 邻位 H 的 δH (ppm)
    # Part 2.B: 溶剂残余峰标记. 该 C-H 对若来自溶剂残余峰, 不进入归属表.
    is_solvent_residual: bool = False


@dataclass
class SpinSystem:
    """COSY 建立的一个氢自旋系统 (一组耦合相连的 H)."""

    h_ppms: List[float] = field(default_factory=list)
    cosy_links: List[Tuple[float, float]] = field(default_factory=list)
    label: Optional[str] = None  # 如 'spin-A'


@dataclass
class HMBCConnection:
    """HMBC 建立的一个远程 H-C 连接 (片段间桥梁)."""

    h_ppm: float
    c_ppm: float
    bond_type: str = 'unknown'  # '2J'/'3J'/'2J/3J'/'long_range'/'unknown'
    connects_fragments: Optional[str] = None  # 如 'spin-A -> Cq-170.5'
    inferred_distance: int = 0  # 推断键数 (2/3)
    inference_basis: str = ''   # 'topology'/'unresolved'
    # Part 2.A: 是否为 ¹J 直接连接残留. 默认为 False (远程相关) 在
    # 过滤后. True 表示可能是 ¹J 残留但未被滤掉 (阈值边缘).
    is_oneJ_residue: bool = False


@dataclass
class OverlapRegion:
    """一个重叠区的分析结果 (对标 5.9 流程)."""

    h_ppm_center: float          # 重叠区中心 δH
    integral_total: float        # 该区域总积分
    n_peaks_in_13c_slice: int    # ¹³C 切片中检出的峰数
    c_peaks: List[float] = field(default_factory=list)  # 切片中各峰的 δC
    c_peak_intensities: List[float] = field(default_factory=list)  # 各峰相对强度
    resolved_pairs: List[CHPair] = field(default_factory=list)  # 分解后的 C-H 对
    status: str = 'pending'      # 'resolved'/'partial'/'unresolvable'
    note: str = ''               # 分析说明


@dataclass
class MultipletAnalysis:
    """对单个 CHPair 的 1H 多重峰分析结果 (对标 TopSpin Manual).

    包含多重性、子峰数、J 耦合、半高宽、强度均匀性、质量等级.
    用于逐 C-H 对精细化分析, 区别于 _step1b_resolve_overlapping 的
    区域分解.
    """

    h_ppm: float                # CHPair 的 1H 化学位移
    c_ppm: float = 0.0          # 对应 13C 化学位移
    pattern: str = ''           # 's'/'d'/'t'/'q'/'m'/''
    num_subpeaks: int = 0
    j_couplings_Hz: List[float] = field(default_factory=list)
    halfwidth_Hz: float = 0.0   # 估算子峰半高宽
    intensity_uniformity: float = 0.0  # 0-1: 子峰强度剖面均匀性
    quality_grade: str = ''     # 'Excellent'/'Good'/'Poor'/'Very Poor'/'N/A'
    multiplet_center_ppm: float = 0.0  # 加权中心
    subpeak_ppms: List[float] = field(default_factory=list)
    subpeak_intensities: List[float] = field(default_factory=list)
    note: str = ''              # 拟合说明
    # 关联到 CHPair (仅内部使用, 不输出)
    c_label: str = ''
    h_label: str = ''


@dataclass
class NOEContact:
    """NOESY/ROESY 检测到的一个空间近邻 H-H 对."""

    h1_ppm: float
    h2_ppm: float
    intensity: float = 0.0
    confidence: str = ''         # 'strong'/'medium'/'weak'
    expected_distance: float = 0.0  # 期望距离 (A), 基于强度-距离映射
    is_symmetric: bool = False   # 是否为对角线对称点


@dataclass
class AssignmentValidation:
    """归属完整性核对结果 (对标 5.4 节第六步核对)."""

    total_observed_c: int = 0
    assigned_c: int = 0
    unassigned_c: List[float] = field(default_factory=list)
    total_observed_h: int = 0
    assigned_h: int = 0
    unassigned_h: List[float] = field(default_factory=list)
    unassigned_qc: List[float] = field(default_factory=list)  # 无 HSQC 对应的 13C
    warnings: List[str] = field(default_factory=list)
    is_complete: bool = False


@dataclass
class Interpretation2D:
    """2D 解谱的完整中间结果 (化学家可审查)."""

    ch_pairs: List[CHPair] = field(default_factory=list)
    spin_systems: List[SpinSystem] = field(default_factory=list)
    tocsy_spin_systems: List[SpinSystem] = field(default_factory=list)
    hmbc_connections: List[HMBCConnection] = field(default_factory=list)
    overlap_regions: List[OverlapRegion] = field(default_factory=list)
    multiplet_analyses: List[MultipletAnalysis] = field(default_factory=list)
    noe_contacts: List[NOEContact] = field(default_factory=list)
    assignment_validation: Optional[AssignmentValidation] = None
    n_cq_inferred: int = 0
    cq_shifts: List[float] = field(default_factory=list)
    cq_evidence: List[str] = field(default_factory=list)  # 每个季碳的推断依据
    notes: List[str] = field(default_factory=list)


# ---- 步骤 1: HSQC → C-H 对 ----

def _step1_hsqc_to_ch_pairs(
    spectra_2d: Dict[str, Spectrum2D],
    dept_summary: List[dict],
    proton_integrals: List[tuple] = None,
    multiplets: List = None,
    solvent: str = '',
    h_tol: float = 0.3,
    c_tol: float = 2.5,
    min_snr: float = 20.0,
) -> List[CHPair]:
    """从 HSQC 谱提取 C-H 直接连接对, 结合 DEPT 和 1H 积分确定碳级数.

    HSQC F1=13C, F2=1H. 每个交叉峰 = 一个直接 C-H 连接.
    碳级数来源优先级: DEPT 推断 > 1H 积分推断 > 化学位移经验 > 'unknown'.
    若提供 multiplets (1H 多重峰分析结果), 会按 h_ppm 匹配将 J 耦合填入 CHPair.

    Part 2.B: 若提供 solvent, HSQC 中位于溶剂残余 (h, c) 位置的交叉峰会被标记
    为 is_solvent_residual, 仍在列表中返回 (供人工复核), 但不阻塞下游处理.

    Args:
        proton_integrals: [(center_ppm, integral_value), ...] 1H 积分数据
        multiplets: [Multiplet, ...] 1H 多重峰分析结果 (J 耦合)
        solvent: 溶剂名 (e.g. 'CDCl3', 'DMSO-D6') 用于标记溶剂残余峰
    """
    from .correlation_2d import extract_hsqc_pairs

    if 'HSQC' not in spectra_2d:
        return []

    hsqc_pairs_raw = extract_hsqc_pairs(spectra_2d['HSQC'], min_snr=min_snr)
    if not hsqc_pairs_raw:
        return []

    # 构建 DEPT 查找: 13C_ppm -> multiplicity
    dept_by_c: Dict[float, str] = {}
    for d in dept_summary:
        ppm = d.get('ppm')
        mult = d.get('mult', 'unknown')
        if ppm is not None:
            dept_by_c[round(ppm, 2)] = mult

    # 构建积分查找: 用于匹配 h_ppm 到积分区间
    # proton_integrals: [(center_ppm, value), ...]
    integrals_sorted = sorted(proton_integrals or [], key=lambda x: x[0])

    def _find_integral(h_ppm: float) -> Optional[float]:
        """找到包含 h_ppm 的积分区间的值."""
        if not integrals_sorted:
            return None
        # 找最近的积分中心 (容差 0.5 ppm)
        best_val = None
        best_dist = 0.5
        for center, val in integrals_sorted:
            dist = abs(h_ppm - center)
            if dist < best_dist:
                best_dist = dist
                best_val = val
        return best_val

    def _integral_to_mult(integral: Optional[float]) -> str:
        """根据积分值推断碳级数."""
        if integral is None:
            return 'unknown'
        # 积分值四舍五入到整数 H 个数
        n_h = round(integral)
        if n_h <= 0:
            return 'unknown'
        if n_h == 1:
            return 'CH'
        if n_h == 2:
            return 'CH2'
        if n_h == 3:
            return 'CH3'
        # 4+ H 可能是多个重叠峰, 不推断
        return 'unknown'

    # 构建 multiplet 查找: center_ppm -> j_couplings (容差 h_tol)
    mult_by_ppm: Dict[float, List[float]] = {}
    if multiplets:
        for m in multiplets:
            mult_by_ppm[round(m.center_ppm, 3)] = list(m.j_couplings or [])

    def _find_j_couplings(h_ppm: float) -> List[float]:
        """查找 h_ppm 所在 multiplet 的 J 耦合. 找不到返回 []."""
        if not mult_by_ppm:
            return []
        for mp_ppm, j_vals in mult_by_ppm.items():
            if abs(h_ppm - mp_ppm) <= h_tol:
                return list(j_vals)
        return []

    ch_pairs: List[CHPair] = []
    # Part 2.B: 溶剂残余交叉峰检测
    from .solvent_filter import is_solvent_peak as _is_sol
    n_solvent_filtered = 0
    for h_ppm, c_ppm in hsqc_pairs_raw:
        # 用 DEPT 确定碳级数 (容差内匹配)
        mult = 'unknown'
        for dept_c, dept_m in dept_by_c.items():
            if abs(c_ppm - dept_c) <= c_tol:
                mult = dept_m
                break

        # 查找 1H 积分
        h_integral = _find_integral(h_ppm)

        # 如果 DEPT 未确定, 用积分推断
        if mult == 'unknown' and h_integral is not None:
            mult = _integral_to_mult(h_integral)

        # 查找该 H 的 J 耦合 (按 h_ppm 匹配 multiplet)
        j_vals = _find_j_couplings(round(h_ppm, 3))

        # Part 2.B: 检查是否为溶剂残余 HSQC 交叉峰
        is_sol_residual = False
        if solvent:
            if (_is_sol(h_ppm, '1H', solvent, 0.05)
                    and _is_sol(c_ppm, '13C', solvent, 1.0)):
                is_sol_residual = True
                n_solvent_filtered += 1

        ch_pairs.append(CHPair(
            h_ppm=round(h_ppm, 3),
            c_ppm=round(c_ppm, 3),
            multiplicity=mult,
            h_integral=round(h_integral, 3) if h_integral is not None else None,
            j_couplings=j_vals,
            is_solvent_residual=is_sol_residual,
        ))

    # 按 13C 位移降序排列 (高场到低场, 与化学家习惯一致)
    ch_pairs.sort(key=lambda p: p.c_ppm, reverse=True)

    # 编号
    for i, pair in enumerate(ch_pairs, 1):
        pair.assigned_c_label = f'C-{i}'
        pair.assigned_h_label = f'H-{i}'

    return ch_pairs


# ---- 步骤 1b: 重叠峰分解 (对标 5.9 流程) ----

def _step1b_resolve_overlapping(
    hsqc_spectrum: Spectrum2D,
    ch_pairs: List[CHPair],
    integral_threshold: float = 4.0,
    h_slice_width: float = 0.15,
    c_min_distance: float = 2.0,
    h_cluster_tol: float = 0.3,
    multiplets: List = None,
) -> List[OverlapRegion]:
    """对标 5.9 流程: 对重叠区做 ¹³C 切片 + 分峰检测 + 分解.

    算法:
    1. 筛选 h_integral >= threshold 的 C-H 对
    2. 按 δH 聚类 (间距 < h_cluster_tol 的归为同一重叠区)
    3. 对每个重叠区:
       a. 在 HSQC 数据中提取 δH ± h_slice_width 的窄带
       b. 沿 ¹H 维求和 → 得到 ¹³C 切片 (trace)
       c. 用 scipy.signal.find_peaks 检测切片中的峰
       d. 峰数 = 该区域实际 C-H 对数
       e. 每个 C-H 对的 H 个数 = total_integral / n_peaks
    f. 分解出的 CHPair 同时从 multiplets 继承 J 耦合 (按 h_ppm 匹配)

    Args:
        hsqc_spectrum: HSQC 2D 谱数据
        ch_pairs: 已提取的 C-H 对列表
        integral_threshold: 积分阈值 (默认 4.0)
        h_slice_width: ¹H 切片宽度 (ppm, 默认 ±0.15)
        c_min_distance: ¹³C 峰最小间距 (ppm, 默认 2.0)
        h_cluster_tol: ¹H 聚类容差 (ppm, 默认 0.3)
        multiplets: [Multiplet, ...] 1H 多重峰 (J 耦合)

    Returns:
        重叠区分析结果列表
    """
    import numpy as np

    # 筛选积分 >= threshold 的 C-H 对
    high_pairs = [
        p for p in ch_pairs
        if p.h_integral is not None and p.h_integral >= integral_threshold
    ]
    if not high_pairs:
        return []

    # 按 δH 聚类 (间距 < h_cluster_tol 的归为同一重叠区)
    high_pairs_sorted = sorted(high_pairs, key=lambda p: p.h_ppm)
    clusters: List[List[CHPair]] = []
    current_cluster: List[CHPair] = [high_pairs_sorted[0]]
    for p in high_pairs_sorted[1:]:
        if abs(p.h_ppm - current_cluster[-1].h_ppm) <= h_cluster_tol:
            current_cluster.append(p)
        else:
            clusters.append(current_cluster)
            current_cluster = [p]
    clusters.append(current_cluster)

    # 对每个重叠区做 ¹³C 切片分析
    f1_ppm = hsqc_spectrum.f1_ppm_scale  # ¹³C
    f2_ppm = hsqc_spectrum.f2_ppm_scale  # ¹H
    data = hsqc_spectrum.data

    # 计算 ¹³C 轴分辨率 (ppm/点)
    if len(f1_ppm) > 1:
        c_ppm_per_pt = abs(f1_ppm[1] - f1_ppm[0])
    else:
        c_ppm_per_pt = 1.0
    min_dist_pts = max(1, int(c_min_distance / c_ppm_per_pt))

    regions: List[OverlapRegion] = []
    for cluster in clusters:
        # 计算重叠区中心和总积分
        h_center = sum(p.h_ppm for p in cluster) / len(cluster)
        integral_total = max(p.h_integral for p in cluster)  # 取最大积分

        # 提取 ¹H 窄带: δH_center ± h_slice_width
        h_lo = h_center - h_slice_width
        h_hi = h_center + h_slice_width
        f2_mask = (f2_ppm >= min(h_lo, h_hi)) & (f2_ppm <= max(h_lo, h_hi))

        if not f2_mask.any():
            regions.append(OverlapRegion(
                h_ppm_center=round(h_center, 3),
                integral_total=round(integral_total, 2),
                n_peaks_in_13c_slice=0,
                status='unresolvable',
                note='¹H 切片窗口内无数据',
            ))
            continue

        # 沿 ¹H 维 (axis=1) 求和 → ¹³C 切片
        trace_13c = np.sum(data[:, f2_mask], axis=1)

        # 检峰: 只检测正峰 (取绝对值处理负峰)
        trace_abs = np.abs(trace_13c)
        noise_level = np.median(trace_abs) * 3
        peak_height = max(noise_level, np.max(trace_abs) * 0.05)

        try:
            from scipy.signal import find_peaks
            peaks_idx, props = find_peaks(
                trace_abs, distance=min_dist_pts, height=peak_height)
        except ImportError:
            # scipy 未安装时回退: 用简单峰值检测
            peaks_idx = _simple_find_peaks(trace_abs, min_dist_pts, peak_height)
            props = {'peak_heights': trace_abs[peaks_idx]} if len(peaks_idx) > 0 else {'peak_heights': np.array([])}

        # 提取峰位置和强度
        c_peaks = [round(float(f1_ppm[i]), 2) for i in peaks_idx]
        c_intensities = [float(trace_abs[i]) for i in peaks_idx]
        n_peaks = len(c_peaks)

        # 归一化强度
        max_int = max(c_intensities) if c_intensities else 1.0
        c_intensities_norm = [round(v / max_int, 3) for v in c_intensities]

        # 判定状态
        expected_h = round(integral_total)
        if n_peaks >= expected_h:
            status = 'resolved'
            note = f'¹³C 切片检出 {n_peaks} 个峰, 可解释积分 {integral_total:.1f}H'
        elif n_peaks > 1:
            status = 'partial'
            note = f'¹³C 切片检出 {n_peaks} 个峰, 少于预期 {expected_h}H, 部分重叠'
        else:
            status = 'unresolvable'
            note = f'¹³C 切片仅检出 {n_peaks} 个峰, 可能是同碳等价 H 或真正简并'

        # 构建分解后的 C-H 对
        resolved: List[CHPair] = []
        # 准备 multiplet J 查找 (同 _step1)
        mult_by_ppm_res: Dict[float, List[float]] = {}
        if multiplets:
            for m in multiplets:
                mult_by_ppm_res[round(m.center_ppm, 3)] = list(m.j_couplings or [])

        def _find_j_res(h_ppm: float) -> List[float]:
            for mp_ppm, j_vals in mult_by_ppm_res.items():
                if abs(h_ppm - mp_ppm) <= h_cluster_tol:
                    return list(j_vals)
            return []

        if n_peaks > 0:
            h_per_peak = integral_total / n_peaks
            h_round = round(h_center, 3)
            j_for_h = _find_j_res(h_round)
            for c_pk in c_peaks:
                mult = 'unknown'
                n_h = round(h_per_peak)
                if n_h == 1:
                    mult = 'CH'
                elif n_h == 2:
                    mult = 'CH2'
                elif n_h == 3:
                    mult = 'CH3'
                resolved.append(CHPair(
                    h_ppm=h_round,
                    c_ppm=c_pk,
                    multiplicity=mult,
                    h_integral=round(h_per_peak, 2),
                    j_couplings=j_for_h,
                ))

        regions.append(OverlapRegion(
            h_ppm_center=round(h_center, 3),
            integral_total=round(integral_total, 2),
            n_peaks_in_13c_slice=n_peaks,
            c_peaks=c_peaks,
            c_peak_intensities=c_intensities_norm,
            resolved_pairs=resolved,
            status=status,
            note=note,
        ))

    return regions


def _simple_find_peaks(arr, min_dist: int, height: float):
    """scipy 不可用时的简单峰值检测回退."""
    import numpy as np
    peaks = []
    for i in range(1, len(arr) - 1):
        if arr[i] > arr[i - 1] and arr[i] >= arr[i + 1] and arr[i] >= height:
            if not peaks or (i - peaks[-1]) >= min_dist:
                peaks.append(i)
            elif arr[i] > arr[peaks[-1]]:
                peaks[-1] = i
    return np.array(peaks, dtype=int)


# ---- 步骤 1b (新): 对每个 CHPair 做 1H 多重峰精细化分析 ----

def _step1b_multiplet_analysis_for_ch_pairs(
    ch_pairs: List['CHPair'],
    peaklist_path: Optional[str] = None,
    field_mhz: float = 0.0,
    window_ppm: float = 0.15,
) -> List[MultipletAnalysis]:
    """对每个 CHPair 做精细化的 1H 多重峰分析 (对标 TopSpin Manual).

    与 _step1b_resolve_overlapping (重叠峰分解, 按 H 簇) 不同:
    本函数对每个 CHPair 的 h_ppm 单独取窗口内的 Peak1D 做精细化拟合,
    输出多重性、子峰数、J 耦合、半高宽、强度均匀性、质量等级.

    Args:
        ch_pairs: HSQC C-H 对列表
        peaklist_path: 1H peaklist.xml 路径; None 或不存在则跳过
        field_mhz: 谱仪 1H 频率 (MHz)
        window_ppm: 取峰窗口半宽 (ppm, 默认 ±0.15)

    Returns:
        MultipletAnalysis 列表 (每个非季碳 CHPair 一个)
    """
    from .io_bruker import analyze_multiplet_at_ppm

    if not peaklist_path or field_mhz <= 0:
        return []

    analyses: List[MultipletAnalysis] = []
    for p in ch_pairs:
        # 季碳 / 无 HSQC 直接 H → 跳过 (无 1H 多重峰可分析)
        if p.multiplicity in ('Cq', 'unknown') and (p.h_integral is None or p.h_integral <= 0):
            continue

        try:
            info = analyze_multiplet_at_ppm(
                peaklist_path=peaklist_path,
                center_ppm=p.h_ppm,
                field_mhz=field_mhz,
                window_ppm=window_ppm,
            )
        except (OSError, ValueError, KeyError):
            continue

        analyses.append(MultipletAnalysis(
            h_ppm=p.h_ppm,
            c_ppm=p.c_ppm,
            pattern=info.get('pattern', ''),
            num_subpeaks=info.get('num_subpeaks', 0),
            j_couplings_Hz=info.get('j_couplings_Hz', []),
            halfwidth_Hz=info.get('halfwidth_Hz', 0.0),
            intensity_uniformity=info.get('intensity_uniformity', 0.0),
            quality_grade=info.get('quality_grade', ''),
            multiplet_center_ppm=info.get('center_ppm', 0.0),
            subpeak_ppms=info.get('subpeak_ppms', []),
            subpeak_intensities=info.get('subpeak_intensities', []),
            note=info.get('note', ''),
            c_label=p.assigned_c_label or '',
            h_label=p.assigned_h_label or '',
        ))
    return analyses


# ---- 步骤 2: COSY → 自旋系统 ----

def _step2_cosy_spin_systems(
    spectra_2d: Dict[str, Spectrum2D],
    ch_pairs: List[CHPair],
    multiplets: List = None,
    h_tol: float = 0.3,
    min_snr: float = 15.0,
    j_threshold: float = 4.0,
) -> List[SpinSystem]:
    """从 COSY 谱建立氢自旋系统.

    策略: 用 Union-Find 将 COSY 相关的 H 归入同一系统.
    仅保留含 >= 2 个 H 的系统 (单峰无耦合信息).

    若有 multiplets (J 耦合数据), 对每个 COSY 交叉峰 (dHa, dHb) 验证:
    仅当 Ha 的 multiplet 中含 |J| > j_threshold 时才视为真实邻位耦合.
    """
    from .correlation_2d import extract_cosy_pairs

    if 'COSY' not in spectra_2d:
        return []

    cosy_pairs = extract_cosy_pairs(spectra_2d['COSY'], min_snr=min_snr)
    if not cosy_pairs:
        return []

    # J 耦合验证: 若有 multiplets 数据, 过滤掉无 J 支持的 COSY 交叉峰
    if multiplets:
        # 构建 multiplet 查找: center_ppm -> j_couplings
        mult_by_ppm: Dict[float, List[float]] = {}
        for m in multiplets:
            mult_by_ppm[round(m.center_ppm, 3)] = m.j_couplings

        def _has_j_support(h_ppm: float) -> bool:
            """Check if h_ppm has any J coupling > threshold."""
            for mp_ppm, j_vals in mult_by_ppm.items():
                if abs(h_ppm - mp_ppm) <= h_tol:
                    return any(abs(j) > j_threshold for j in j_vals)
            return True  # 无 multiplet 数据时不过滤

        # 保留至少一侧有 J 支持的 COSY 对
        cosy_pairs = [(a, b) for a, b in cosy_pairs
                      if _has_j_support(round(a, 3)) or _has_j_support(round(b, 3))]
        if not cosy_pairs:
            return []

    # 收集所有出现在 COSY 中的 H 位移
    all_h: Set[float] = set()
    for a, b in cosy_pairs:
        all_h.add(round(a, 3))
        all_h.add(round(b, 3))

    # Union-Find
    parent: Dict[float, float] = {h: h for h in all_h}

    def find(x: float) -> float:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: float, b: float):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in cosy_pairs:
        a_r, b_r = round(a, 3), round(b, 3)
        if a_r in parent and b_r in parent:
            union(a_r, b_r)

    # 分组
    groups: Dict[float, List[float]] = {}
    for h in all_h:
        root = find(h)
        if root not in groups:
            groups[root] = []
        groups[root].append(h)

    # 构建 SpinSystem (仅 >= 2 个 H 的组)
    systems: List[SpinSystem] = []
    label_idx = 0
    for root, members in sorted(groups.items(), key=lambda x: -len(x[1])):
        if len(members) < 2:
            continue
        label_idx += 1
        label = f'spin-{chr(64 + label_idx)}'  # spin-A, spin-B, ...
        members_sorted = sorted(members, reverse=True)
        # 收集系统内的 COSY links
        links = [(a, b) for a, b in cosy_pairs
                 if round(a, 3) in members and round(b, 3) in members]
        systems.append(SpinSystem(
            h_ppms=members_sorted,
            cosy_links=links,
            label=label,
        ))

    return systems


# ---- 步骤 3: HMBC → 远程连接 ----

def _step3_hmbc_connections(
    spectra_2d: Dict[str, Spectrum2D],
    ch_pairs: List[CHPair],
    cq_shifts: List[float],
    h_tol: float = 0.4,
    c_tol: float = 3.0,
    min_snr: float = 10.0,
    max_connections: int = 100,
) -> List[HMBCConnection]:
    """从 HMBC 谱建立远程 H-C 连接.

    策略:
    - 对每个 HMBC 峰 (h_ppm, c_ppm):
      - 确认 h_ppm 对应一个已知 C-H 对 (即 H 已归属)
      - 确认 c_ppm 不是该 H 自身的直接连接 C (排除 1J)
      - 推测键类型: c_ppm 在已知 C-H 对的 C 列表中 → 3J; 在 Cq 中 → 2J/3J
    """
    from .correlation_2d import extract_hmbc_pairs

    if 'HMBC' not in spectra_2d:
        return []

    # Part 2.A: 传入 HSQC 配对过滤 ¹J 残留. ch_pairs (h_ppm, c_ppm)
    # 来自 HSQC 视为 ¹J 真值. HMBC 中相同 (h, c) 视为 ¹J 残留删除.
    hsqc_known = [(p.h_ppm, p.c_ppm) for p in ch_pairs]
    hmbc_pairs, n_oneJ_filtered = extract_hmbc_pairs(
        spectra_2d['HMBC'], min_snr=min_snr,
        hsqc_pairs=hsqc_known if hsqc_known else None,
    )
    if not hmbc_pairs:
        return []

    # 已知 C-H 对中的 H 和 C 位移集合
    known_h = {round(p.h_ppm, 2) for p in ch_pairs}
    known_c_from_hsqc = {round(p.c_ppm, 2) for p in ch_pairs}
    # 直接连接映射: h_ppm -> 其直接连接的 c_ppm
    direct_c_for_h: Dict[float, float] = {
        round(p.h_ppm, 2): round(p.c_ppm, 2) for p in ch_pairs
    }

    connections: List[HMBCConnection] = []
    seen: Set[Tuple[float, float]] = set()

    for h_ppm, c_ppm in hmbc_pairs:
        if len(connections) >= max_connections:
            break
        h_r, c_r = round(h_ppm, 2), round(c_ppm, 2)

        # H 必须是已归属的 (出现在 HSQC 中)
        h_matched = any(abs(h_r - kh) <= h_tol for kh in known_h)
        if not h_matched:
            continue

        # 排除 1J: c_ppm 不能是该 H 的直接连接 C
        for kh, dc in direct_c_for_h.items():
            if abs(h_r - kh) <= h_tol and abs(c_r - dc) <= c_tol:
                break
        else:
            # 未命中排除条件, 保留
            key = (h_r, c_r)
            if key in seen:
                continue
            seen.add(key)

            # 推测键类型
            bond_type = 'unknown'
            c_is_cq = any(abs(c_r - cq) <= c_tol for cq in cq_shifts)
            c_is_known = any(abs(c_r - kc) <= c_tol for kc in known_c_from_hsqc)
            if c_is_cq:
                bond_type = '2J/3J'  # 到季碳通常是 2J 或 3J
            elif c_is_known:
                bond_type = '3J'  # 到另一个 CH 通常是 3J

            connections.append(HMBCConnection(
                h_ppm=round(h_ppm, 3),
                c_ppm=round(c_ppm, 3),
                bond_type=bond_type,
            ))

    return connections


# ---- 步骤 4: 推断季碳 ----

def _step4_infer_quaternary_carbons(
    carbon_shifts_1d: List[float],
    ch_pairs: List[CHPair],
    dept_summary: List[dict] = None,
    c_tol: float = 2.5,
) -> List[float]:
    """推断季碳: 综合 DEPT 和 HSQC 两种证据.

    优先级:
    1. DEPT 明确判定为 Cq 的 → 直接信任 (DEPT 中 Cq 无信号是物理事实)
    2. 无 DEPT 信息时, 13C 有峰但无 HSQC 对应 → 推断为 Cq

    Args:
        carbon_shifts_1d: 1D 13C 观测位移列表
        ch_pairs: HSQC C-H 对列表
        dept_summary: DEPT 推断结果 [{ppm, mult, sign}, ...]
        c_tol: HSQC 匹配容差 (ppm)
    """
    cq: List[float] = []

    # 证据 1: DEPT 明确判定为 Cq 的碳 (最可靠)
    dept_cq_shifts: List[float] = []
    if dept_summary:
        dept_cq_shifts = [d['ppm'] for d in dept_summary if d.get('mult') == 'Cq']

    # 证据 2: 13C 有峰但无 HSQC 对应
    hsqc_c_shifts = [p.c_ppm for p in ch_pairs] if ch_pairs else []
    no_hsqc_shifts: List[float] = []
    if hsqc_c_shifts:
        for c_obs in carbon_shifts_1d:
            matched = any(abs(c_obs - hc) <= c_tol for hc in hsqc_c_shifts)
            if not matched:
                no_hsqc_shifts.append(c_obs)

    # 合并: DEPT Cq 优先, 再补充无 HSQC 匹配的 (去重)
    # 返回 (ppm, 依据) 元组列表, 便于 CSV 输出区分证据来源
    seen: set = set()
    results: List[tuple] = []  # [(ppm, evidence_str), ...]
    for c in dept_cq_shifts:
        cr = round(c, 3)
        if cr not in seen:
            seen.add(cr)
            results.append((cr, 'DEPT判定为Cq'))
    for c in no_hsqc_shifts:
        cr = round(c, 3)
        # 避免与 DEPT Cq 重复 (容差内视为同一个峰)
        if not any(abs(cr - s) <= c_tol for s in seen):
            seen.add(cr)
            results.append((cr, '13C有峰但无HSQC对应'))

    results.sort(key=lambda x: x[0], reverse=True)
    return results


# ---- 步骤 2b: TOCSY 长程自旋系统 ----

def _step2_tocsy_spin_systems(
    spectra_2d: Dict[str, Spectrum2D],
    ch_pairs: List[CHPair],
    h_tol: float = 0.6,
    min_snr: float = 15.0,
) -> List[SpinSystem]:
    """从 TOCSY 谱建立长程氢自旋系统.

    与 COSY 同算法 (Union-Find), 但容差放大 (0.6 ppm),
    覆盖整个自旋系统 (通常 4-7 个 H).
    """
    from .correlation_2d import extract_cosy_pairs

    if 'TOCSY' not in spectra_2d:
        return []

    # TOCSY 谱的峰提取复用 COSY 逻辑 (同为 1H x 1H)
    tocsy_pairs = extract_cosy_pairs(spectra_2d['TOCSY'], min_snr=min_snr)
    if not tocsy_pairs:
        return []

    # Union-Find
    all_h: Set[float] = set()
    for a, b in tocsy_pairs:
        all_h.add(round(a, 3))
        all_h.add(round(b, 3))

    parent: Dict[float, float] = {h: h for h in all_h}

    def find(x: float) -> float:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: float, b: float):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in tocsy_pairs:
        a_r, b_r = round(a, 3), round(b, 3)
        if a_r in parent and b_r in parent:
            union(a_r, b_r)

    groups: Dict[float, List[float]] = {}
    for h in all_h:
        root = find(h)
        if root not in groups:
            groups[root] = []
        groups[root].append(h)

    systems: List[SpinSystem] = []
    label_idx = 0
    for root, members in sorted(groups.items(), key=lambda x: -len(x[1])):
        if len(members) < 2:
            continue
        label_idx += 1
        label = f'tocsy-{chr(64 + label_idx)}'
        members_sorted = sorted(members, reverse=True)
        links = [(a, b) for a, b in tocsy_pairs
                 if round(a, 3) in members and round(b, 3) in members]
        systems.append(SpinSystem(
            h_ppms=members_sorted,
            cosy_links=links,
            label=label,
        ))

    return systems


def _fill_h_neighbors_from_spin_systems(
    ch_pairs: List['CHPair'],
    spin_systems: List[SpinSystem],
    tocsy_systems: Optional[List[SpinSystem]] = None,
    h_match_tol: float = 0.05,
) -> None:
    """回填 CHPair.h_neighbors: 根据自旋系统的 cosy_links 找出每个 H 的邻位 H ppm.

    对每个 CHPair, 在 spin_systems 中查找其 h_ppm 所在的系统, 然后用
    cosy_links (键对) 找出所有通过 COSY 耦合相连的 H ppm. 去重 + 排序.

    TOCSY 系统仅用于补充远程相邻 (TOCSY 覆盖 COSY 看不到的耦合).

    Args:
        ch_pairs: 已生成的 C-H 对列表 (原地修改)
        spin_systems: COSY 自旋系统
        tocsy_systems: TOCSY 自旋系统 (可选)
        h_match_tol: h_ppm 匹配容差 (ppm, 默认 0.05)
    """
    if not spin_systems and not tocsy_systems:
        return

    # 构建 H ppm -> 邻位 H ppm 集合 的查找表
    neighbors_map: Dict[float, Set[float]] = {}

    def _add_neighbor(a: float, b: float):
        ra = round(a, 3)
        rb = round(b, 3)
        if ra == rb:
            return
        neighbors_map.setdefault(ra, set()).add(rb)
        neighbors_map.setdefault(rb, set()).add(ra)

    for ss in spin_systems or []:
        for a, b in ss.cosy_links:
            _add_neighbor(a, b)
    for ss in tocsy_systems or []:
        for a, b in ss.cosy_links:
            _add_neighbor(a, b)

    if not neighbors_map:
        return

    def _find_neighbors(h_ppm: float) -> List[float]:
        rh = round(h_ppm, 3)
        # 直接匹配
        if rh in neighbors_map:
            return sorted(neighbors_map[rh], reverse=True)
        # 模糊匹配 (容差内最近的 H)
        best_key = None
        best_dist = h_match_tol
        for k in neighbors_map:
            d = abs(rh - k)
            if d <= best_dist:
                best_dist = d
                best_key = k
        if best_key is not None:
            return sorted(neighbors_map[best_key], reverse=True)
        return []

    for p in ch_pairs:
        if p.h_neighbors:
            continue
        p.h_neighbors = _find_neighbors(p.h_ppm)


# ---- 步骤 3b: HMBC 2J/3J 键数推断 ----

def _annotate_hmbc_bond_type(
    hmbc_connections: List[HMBCConnection],
    ch_pairs: List[CHPair],
    spin_systems: List[SpinSystem],
) -> List[HMBCConnection]:
    """基于自旋系统拓扑推断 HMBC 是 2J 还是 3J.

    规则:
    - 若 H 的 HSQC C 在同一自旋系统 -> 距离 1-2 -> likely 2J
    - 若 H 与 C 在不同自旋系统 -> 跨系统 -> likely 3J
    - 自旋系统不确定 -> long_range
    """
    # 构建 H ppm -> 所属自旋系统 label 的映射
    h_to_spin: Dict[float, str] = {}
    for ss in spin_systems:
        for h in ss.h_ppms:
            h_to_spin[round(h, 3)] = ss.label or ''

    # 构建 C ppm -> 其直接连接 H 的自旋系统
    c_to_spin: Dict[float, str] = {}
    for p in ch_pairs:
        h_r = round(p.h_ppm, 3)
        if h_r in h_to_spin:
            c_to_spin[round(p.c_ppm, 2)] = h_to_spin[h_r]

    for conn in hmbc_connections:
        h_r = round(conn.h_ppm, 3)
        c_r = round(conn.c_ppm, 2)

        h_spin = h_to_spin.get(h_r, '')
        c_spin = c_to_spin.get(c_r, '')

        if not h_spin or not c_spin:
            conn.bond_type = 'long_range'
            conn.inference_basis = 'unresolved'
            continue

        if h_spin == c_spin:
            # 同一自旋系统内: 距离近, likely 2J
            conn.bond_type = '2J'
            conn.inferred_distance = 2
            conn.inference_basis = 'topology'
        else:
            # 跨自旋系统: 距离远, likely 3J
            conn.bond_type = '3J'
            conn.inferred_distance = 3
            conn.inference_basis = 'topology'

    return hmbc_connections


# ---- 步骤 5: NOESY 空间近邻提取 ----

def _step5_extract_noe_contacts(
    spectra_2d: Dict[str, Spectrum2D],
    diagonal_exclusion: float = 0.1,
    min_snr: float = 10.0,
) -> List[NOEContact]:
    """从 NOESY/ROESY 谱提取空间近邻 H-H 对.

    算法:
    1. 排除主对角线附近 (+/- diagonal_exclusion ppm)
    2. 检测交叉峰
    3. 对称化: (h1, h2) 和 (h2, h1) 合并
    4. 按强度分级: strong (top 30%) / medium / weak (bottom 30%)
    5. 强度-距离映射: strong < 2.5A, medium 2.5-3.5A, weak 3.5-5.0A
    """
    # 查找 NOESY 或 ROESY 谱
    noe_spectrum = spectra_2d.get('NOESY') or spectra_2d.get('ROESY')
    if noe_spectrum is None:
        return []

    if not noe_spectrum.peaks:
        return []

    # 过滤对角线附近的峰
    filtered: List[Tuple[float, float, float]] = []  # (h1, h2, intensity)
    for peak in noe_spectrum.peaks:
        h1, h2 = peak.f2_ppm, peak.f1_ppm  # NOESY: F1=1H, F2=1H
        if abs(h1 - h2) < diagonal_exclusion:
            continue  # 排除对角线
        filtered.append((round(h1, 3), round(h2, 3), peak.intensity))

    if not filtered:
        return []

    # 对称化: 合并 (h1, h2) 和 (h2, h1)
    merged: Dict[Tuple[float, float], float] = {}
    for h1, h2, inten in filtered:
        key = (min(h1, h2), max(h1, h2))
        if key in merged:
            merged[key] = max(merged[key], inten)
        else:
            merged[key] = inten

    # 按强度分级
    intensities = sorted(merged.values(), reverse=True)
    n = len(intensities)
    if n == 0:
        return []
    strong_thresh = intensities[max(0, int(n * 0.3) - 1)] if n >= 3 else intensities[0]
    weak_thresh = intensities[min(n - 1, int(n * 0.7))] if n >= 3 else intensities[-1]

    contacts: List[NOEContact] = []
    for (h1, h2), inten in sorted(merged.items(), key=lambda x: -x[1]):
        if inten >= strong_thresh:
            conf = 'strong'
            exp_dist = 2.0  # < 2.5 A
        elif inten <= weak_thresh:
            conf = 'weak'
            exp_dist = 4.5  # 3.5-5.0 A
        else:
            conf = 'medium'
            exp_dist = 3.0  # 2.5-3.5 A

        contacts.append(NOEContact(
            h1_ppm=h1,
            h2_ppm=h2,
            intensity=inten,
            confidence=conf,
            expected_distance=exp_dist,
            is_symmetric=True,
        ))

    return contacts


# ---- 步骤 7: 归属完整性核对 ----

def validate_assignment(
    interpretation: Interpretation2D,
    observed_c_shifts: List[float],
    observed_h_shifts: List[float],
    dept_summary: List[dict] = None,
    h_tol: float = 0.3,
    c_tol: float = 2.5,
) -> AssignmentValidation:
    """核对归属完整性 (对标 5.4 节第六步核对规则).

    检查项:
    1. 每个非季碳 (CH/CH2/CH3) 都应有 HSQC 对应
    2. 每个观察到的 H 都应能找到归属
    3. 总积分与观察 H 数一致性
    """
    result = AssignmentValidation()
    result.total_observed_c = len(observed_c_shifts)
    result.total_observed_h = len(observed_h_shifts)

    ch_pairs = interpretation.ch_pairs
    hsqc_c_shifts = [p.c_ppm for p in ch_pairs]
    hsqc_h_shifts = [p.h_ppm for p in ch_pairs]

    # 检查 1: 13C 中非季碳是否有 HSQC 对应
    result.assigned_c = len(hsqc_c_shifts)
    for c_obs in observed_c_shifts:
        matched = any(abs(c_obs - hc) <= c_tol for hc in hsqc_c_shifts)
        if not matched:
            result.unassigned_c.append(round(c_obs, 3))
            # 检查是否为 DEPT 判定的季碳
            is_dept_cq = False
            if dept_summary:
                is_dept_cq = any(
                    abs(c_obs - d.get('ppm', -999)) <= c_tol
                    and d.get('mult') == 'Cq'
                    for d in dept_summary
                )
            if not is_dept_cq:
                result.warnings.append(
                    f'13C {c_obs:.2f} ppm 无 HSQC 对应且非 DEPT 季碳, 可能漏归属')

    # 检查 2: 1H 中是否有峰无 HSQC 对应
    result.assigned_h = len(hsqc_h_shifts)
    for h_obs in observed_h_shifts:
        matched = any(abs(h_obs - hh) <= h_tol for hh in hsqc_h_shifts)
        if not matched:
            result.unassigned_h.append(round(h_obs, 3))
            # 活泼氢通常在 1-5 ppm 宽峰, 或 > 10 ppm
            if h_obs > 10.0 or (1.0 < h_obs < 5.0):
                result.warnings.append(
                    f'1H {h_obs:.2f} ppm 无 HSQC 对应, 可能为活泼氢 (NH/OH)')

    # 判断完整性
    result.is_complete = (len(result.unassigned_c) == 0 and len(result.unassigned_h) == 0)

    return result


# ---- 顶层接口 ----

def interpret_2d_sequential(
    spectra_2d: Dict[str, Spectrum2D],
    carbon_shifts_1d: List[float],
    proton_shifts_1d: List[float],
    dept_summary: List[dict],
    proton_integrals: List[tuple] = None,
    multiplets: List = None,
    h1_peaklist_path: Optional[str] = None,
    h1_field_mhz: float = 0.0,
    solvent: str = '',
    quiet: bool = False,
) -> Interpretation2D:
    """按化学家顺序逐步解读 2D NMR, 返回结构化中间结果.

    步骤:
    1. HSQC → C-H pairs (用 DEPT + 1H 积分确定 CH/CH2/CH3)
    2. COSY/TOCSY → 从 C-H pairs 中的 H 出发, 沿耦合链建立自旋系统
    3. HMBC → 从已归属的 H 出发, 找 2-3 键外的 C, 建立片段间连接 + 2J/3J 推断
    4. 推断季碳: 13C 中无 HSQC 对应的 = Cq
    5. NOESY → 提取空间近邻 H-H 对
    7. 归属完整性核对

    Args:
        spectra_2d: {experiment_type: Spectrum2D} 字典
        carbon_shifts_1d: 1D 13C 观测位移列表
        proton_shifts_1d: 1D 1H 观测位移列表
        dept_summary: DEPT 推断 [{ppm, mult, sign}, ...]
        proton_integrals: [(center_ppm, integral_value), ...] 1H 积分数据
        multiplets: [Multiplet, ...] 1H 多重峰 (全局峰聚类结果)
        h1_peaklist_path: 1H peaklist.xml 路径 (用于精细化多重峰分析)
        h1_field_mhz: 谱仪 1H 频率 (MHz)
        solvent: 溶剂名 (Part 2.B: 标记 HSQC 溶剂残余 HSQC 交叉峰)
        quiet: 不打印进度

    Returns:
        Interpretation2D 完整中间结果
    """
    result = Interpretation2D()

    if not spectra_2d:
        result.notes.append('无 2D 数据, 跳过逐步解读')
        return result

    # 步骤 1: HSQC → C-H 对 (结合 DEPT + 1H 积分推断碳级数, 填充 J 耦合)
    ch_pairs = _step1_hsqc_to_ch_pairs(
        spectra_2d, dept_summary, proton_integrals, multiplets, solvent=solvent)
    result.ch_pairs = ch_pairs

    # 步骤 1b: 重叠峰分解 (对标 5.9 流程: ¹³C 切片 + 分峰检测)
    if 'HSQC' in spectra_2d:
        overlap_regions = _step1b_resolve_overlapping(
            spectra_2d['HSQC'], ch_pairs, multiplets=multiplets)
        result.overlap_regions = overlap_regions
        # 用分解结果修正 C-H 对的 multiplicity
        for region in overlap_regions:
            if region.status == 'resolved' and region.resolved_pairs:
                h_per_peak = region.integral_total / region.n_peaks_in_13c_slice
                for p in ch_pairs:
                    if (p.h_integral is not None
                            and p.h_integral >= 4.0
                            and abs(p.h_ppm - region.h_ppm_center) <= 0.3):
                        # 修正: 每个 C-H 对的实际 H 个数
                        p.h_integral = round(h_per_peak, 2)
                        n_h = round(h_per_peak)
                        if n_h == 1:
                            p.multiplicity = 'CH'
                        elif n_h == 2:
                            p.multiplicity = 'CH2'
                        elif n_h == 3:
                            p.multiplicity = 'CH3'

    # 步骤 1b (新): 对每个 CHPair 做精细化 1H 多重峰分析 (对标 TopSpin Manual)
    if h1_peaklist_path and h1_field_mhz > 0:
        multiplet_analyses = _step1b_multiplet_analysis_for_ch_pairs(
            ch_pairs, peaklist_path=h1_peaklist_path, field_mhz=h1_field_mhz)
        result.multiplet_analyses = multiplet_analyses
        # 用精细化结果更新 CHPair 的 j_couplings (覆盖全局聚类的旧值)
        for ma in multiplet_analyses:
            if not ma.j_couplings_Hz:
                continue
            for p in ch_pairs:
                if (abs(p.h_ppm - ma.h_ppm) <= 0.05
                        and ma.pattern in ('d', 't', 'q', 'm')):
                    p.j_couplings = list(ma.j_couplings_Hz)
                    break

    # 步骤 4 先做 (步骤 3 需要季碳列表)
    cq_results = _step4_infer_quaternary_carbons(carbon_shifts_1d, ch_pairs, dept_summary)
    cq_shifts = [ppm for ppm, _ in cq_results]
    result.cq_shifts = cq_shifts
    result.cq_evidence = [ev for _, ev in cq_results]
    result.n_cq_inferred = len(cq_shifts)

    # 步骤 2: COSY → 自旋系统
    spin_systems = _step2_cosy_spin_systems(spectra_2d, ch_pairs, multiplets)
    result.spin_systems = spin_systems

    # 步骤 2b: TOCSY → 长程自旋系统 (若有 TOCSY 谱)
    tocsy_systems = _step2_tocsy_spin_systems(spectra_2d, ch_pairs)
    result.tocsy_spin_systems = tocsy_systems

    # 回填 CHPair.h_neighbors: 根据自旋系统的 cosy_links 找出每个 H 的邻位 H ppm
    _fill_h_neighbors_from_spin_systems(ch_pairs, spin_systems, tocsy_systems)

    # 步骤 3: HMBC → 远程连接
    hmbc_connections = _step3_hmbc_connections(
        spectra_2d, ch_pairs, cq_shifts)
    # 步骤 3b: HMBC 2J/3J 键数推断
    all_spin_systems = spin_systems + tocsy_systems
    _annotate_hmbc_bond_type(hmbc_connections, ch_pairs, all_spin_systems)
    result.hmbc_connections = hmbc_connections

    # 步骤 5: NOESY → 空间近邻 H-H 对
    noe_contacts = _step5_extract_noe_contacts(spectra_2d)
    result.noe_contacts = noe_contacts

    # 步骤 7: 归属完整性核对
    validation = validate_assignment(
        result, carbon_shifts_1d, proton_shifts_1d, dept_summary)
    result.assignment_validation = validation

    # 控制台仅输出一行摘要 (详细结果见 CSV 文件)
    if not quiet:
        n_overlap = len(result.overlap_regions)
        n_resolved = sum(1 for r in result.overlap_regions if r.status == 'resolved')
        overlap_str = f' / 重叠峰 {n_overlap} 区(分解 {n_resolved})' if n_overlap else ''
        noe_str = f' / NOE {len(noe_contacts)} 对' if noe_contacts else ''
        tocsy_str = f' / TOCSY {len(tocsy_systems)} 系统' if tocsy_systems else ''
        print(f"  2D 解读完成: HSQC {len(ch_pairs)} 对 / "
              f"COSY {len(spin_systems)} 系统{tocsy_str} / "
              f"HMBC {len(hmbc_connections)} 连接 / "
              f"季碳 {len(cq_shifts)} 个{overlap_str}{noe_str}",
              file=sys.stderr)

    # 汇总注释
    result.notes.append(
        f'HSQC: {len(ch_pairs)} 个 C-H 对; '
        f'COSY: {len(spin_systems)} 个自旋系统; '
        f'TOCSY: {len(tocsy_systems)} 个; '
        f'HMBC: {len(hmbc_connections)} 个远程连接; '
        f'季碳: {len(cq_shifts)} 个; '
        f'NOE: {len(noe_contacts)} 个空间近邻; '
        f'归属完整性: {"OK" if validation.is_complete else "有警告"}'
    )

    return result


# ---- 中间结果文件输出 (化学家复核用) ----

def _draw_hsqc_contour_all_pairs(
    hsqc_spectrum: Spectrum2D,
    ch_pairs: List[CHPair],
    output_dir: str,
    integral_threshold: float = 4.0,
    c_window: float = 5.0,
    h_window: float = 0.5,
) -> List[str]:
    """为所有 C-H 对绘制 Mnova 风格 HSQC 等高线图, 按重叠/非重叠分文件夹.

    输出目录结构:
      step1_hsqc_ch_pairs/overlap/   ← 积分 >= threshold 的重叠峰
      step1_hsqc_ch_pairs/normal/    ← 积分 < threshold 或无积分的非重叠峰

    图表包含:
    - 中央: 2D 等高线图 (蓝色正峰, 红色负峰)
    - 上方: 1H 投影 (F2 维求和)
    - 左侧: 13C 投影 (F1 维求和)
    - 右侧: 13C 切片 (带检峰标注)
    - 标记: 峰位置 + 十字线

    Args:
        hsqc_spectrum: HSQC 2D 谱数据
        ch_pairs: C-H 对列表 (全部)
        output_dir: 输出目录
        integral_threshold: 积分阈值, >= 此值归入 overlap 文件夹 (默认 4.0)
        c_window: 13C 轴窗口 (ppm, 默认 ±5)
        h_window: 1H 轴窗口 (ppm, 默认 ±0.5)

    Returns:
        生成的图片文件路径列表
    """
    import matplotlib
    matplotlib.use('Agg')  # 非交互后端
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import numpy as np

    if not ch_pairs:
        return []

    # 创建两个子目录: 重叠峰 / 非重叠峰
    contour_base = os.path.join(output_dir, 'step1_hsqc_ch_pairs')
    overlap_dir = os.path.join(contour_base, 'overlap')
    normal_dir = os.path.join(contour_base, 'normal')
    os.makedirs(overlap_dir, exist_ok=True)
    os.makedirs(normal_dir, exist_ok=True)

    f1_ppm = hsqc_spectrum.f1_ppm_scale  # 13C
    f2_ppm = hsqc_spectrum.f2_ppm_scale  # 1H
    data = hsqc_spectrum.data

    written: List[str] = []
    for p in ch_pairs:
        # 判断是否为重叠峰, 选择输出目录
        is_overlap = (p.h_integral is not None and p.h_integral >= integral_threshold)
        contour_dir = overlap_dir if is_overlap else normal_dir
        # 提取局部区域
        c_min, c_max = p.c_ppm - c_window, p.c_ppm + c_window
        h_min, h_max = p.h_ppm - h_window, p.h_ppm + h_window

        # 找到对应的索引范围
        f1_mask = (f1_ppm >= min(c_min, c_max)) & (f1_ppm <= max(c_min, c_max))
        f2_mask = (f2_ppm >= min(h_min, h_max)) & (f2_ppm <= max(h_min, h_max))

        if not f1_mask.any() or not f2_mask.any():
            continue

        f1_idx = np.where(f1_mask)[0]
        f2_idx = np.where(f2_mask)[0]
        local_data = data[f1_idx[0]:f1_idx[-1]+1, f2_idx[0]:f2_idx[-1]+1]
        local_f1 = f1_ppm[f1_idx[0]:f1_idx[-1]+1]  # 13C
        local_f2 = f2_ppm[f2_idx[0]:f2_idx[-1]+1]  # 1H

        if local_data.size == 0:
            continue

        # 计算 1D 投影
        proj_1h = np.sum(local_data, axis=0)   # 沿 13C 求和 -> 1H 投影
        proj_13c = np.sum(local_data, axis=1)  # 沿 1H 求和 -> 13C 投影

        # ¹³C 切片: 在峰位置 δH ± 0.15 ppm 的窄带求和
        h_slice_w = 0.15
        slice_mask = (local_f2 >= p.h_ppm - h_slice_w) & (local_f2 <= p.h_ppm + h_slice_w)
        if slice_mask.any():
            trace_13c = np.sum(local_data[:, slice_mask], axis=1)
        else:
            trace_13c = proj_13c  # 回退到全投影

        # 在切片中检峰
        trace_abs = np.abs(trace_13c)
        noise_lev = np.median(trace_abs) * 3
        pk_height = max(noise_lev, np.max(trace_abs) * 0.05)
        if len(local_f1) > 1:
            c_ppm_pt = abs(local_f1[1] - local_f1[0])
        else:
            c_ppm_pt = 1.0
        min_d = max(1, int(2.0 / c_ppm_pt))
        try:
            from scipy.signal import find_peaks as _fp
            pk_idx, _ = _fp(trace_abs, distance=min_d, height=pk_height)
        except ImportError:
            pk_idx = _simple_find_peaks(trace_abs, min_d, pk_height)

        # 创建 Mnova 风格布局: 上方 1H + 左侧 13C + 中央 2D + 右侧 ¹³C 切片
        fig = plt.figure(figsize=(10, 7))
        gs = GridSpec(2, 3, width_ratios=[1, 5, 1.2], height_ratios=[1, 5],
                     hspace=0.05, wspace=0.05)

        ax_2d = fig.add_subplot(gs[1, 1])    # 中央 2D 等高线
        ax_1h = fig.add_subplot(gs[0, 1])    # 上方 1H 投影
        ax_13c = fig.add_subplot(gs[1, 0])   # 左侧 13C 投影
        ax_slice = fig.add_subplot(gs[1, 2]) # 右侧 ¹³C 切片

        # --- 中央 2D 等高线图 (蓝色正峰 + 红色虚线负峰, Mnova 风格) ---
        max_val = np.max(np.abs(local_data))
        min_val = np.min(local_data)
        if max_val > 0:
            # 正峰等高线 (蓝色, Mnova 风格)
            pos_levels = []
            level = max_val * 0.08
            while level < max_val * 0.95:
                pos_levels.append(level)
                level *= 1.4
            if pos_levels:
                ax_2d.contour(local_f2, local_f1, local_data,
                             levels=pos_levels, colors='#0066CC', linewidths=0.7)
            # 负峰等高线 (红色虚线, CH2 在相敏 HSQC 中为负)
            if min_val < 0:
                neg_levels = []
                level = min_val * 0.08
                while level > min_val * 0.95:
                    neg_levels.append(level)
                    level *= 1.4
                neg_levels.sort()
                if neg_levels:
                    ax_2d.contour(local_f2, local_f1, local_data,
                                 levels=neg_levels, colors='#CC0000',
                                 linewidths=0.7, linestyles='--')

        # 标记峰位置 (十字线)
        ax_2d.axvline(p.h_ppm, color='gray', linewidth=0.5, linestyle=':', alpha=0.6)
        ax_2d.axhline(p.c_ppm, color='gray', linewidth=0.5, linestyle=':', alpha=0.6)
        ax_2d.plot(p.h_ppm, p.c_ppm, 'r+', markersize=14, markeredgewidth=2)

        ax_2d.set_xlabel('δH (ppm)', fontsize=10)
        ax_2d.set_ylabel('δC (ppm)', fontsize=10)
        ax_2d.invert_xaxis()  # NMR 惯例: 高场在右
        ax_2d.invert_yaxis()

        # --- 上方 1H 投影 ---
        ax_1h.fill_between(local_f2, proj_1h, alpha=0.3, color='#0066CC')
        ax_1h.plot(local_f2, proj_1h, color='#0066CC', linewidth=0.8)
        ax_1h.axvline(p.h_ppm, color='red', linewidth=1, linestyle='--', alpha=0.7)
        ax_1h.set_xlim(ax_2d.get_xlim())
        ax_1h.invert_xaxis()
        ax_1h.set_yticks([])
        # 标题: 含积分信息 (无积分时省略)
        if p.h_integral is not None:
            title_str = (f"{p.assigned_c_label}/{p.assigned_h_label}: "
                        f"δC {p.c_ppm:.2f}, δH {p.h_ppm:.3f} | "
                        f"1H积分={p.h_integral:.1f} (~{round(p.h_integral)}H) | "
                        f"碳级数={p.multiplicity}")
        else:
            title_str = (f"{p.assigned_c_label}/{p.assigned_h_label}: "
                        f"δC {p.c_ppm:.2f}, δH {p.h_ppm:.3f} | "
                        f"碳级数={p.multiplicity}")
        ax_1h.set_title(title_str, fontsize=10, fontweight='bold')

        # --- 左侧 13C 投影 ---
        ax_13c.fill_betweenx(local_f1, proj_13c, alpha=0.3, color='#0066CC')
        ax_13c.plot(proj_13c, local_f1, color='#0066CC', linewidth=0.8)
        ax_13c.axhline(p.c_ppm, color='red', linewidth=1, linestyle='--', alpha=0.7)
        ax_13c.set_ylim(ax_2d.get_ylim())
        ax_13c.invert_yaxis()
        ax_13c.set_xticks([])

        # --- 右侧 ¹³C 切片 (带检峰标注) ---
        ax_slice.fill_betweenx(local_f1, trace_abs, alpha=0.2, color='#CC6600')
        ax_slice.plot(trace_abs, local_f1, color='#CC6600', linewidth=1.0)
        # 标注检出的峰 (红色三角)
        if len(pk_idx) > 0:
            ax_slice.plot(trace_abs[pk_idx], local_f1[pk_idx],
                         'rv', markersize=8, label=f'{len(pk_idx)} 峰')
            for idx in pk_idx:
                ax_slice.annotate(f'{local_f1[idx]:.1f}',
                                 (trace_abs[idx], local_f1[idx]),
                                 fontsize=7, color='red',
                                 xytext=(3, 0), textcoords='offset points')
        ax_slice.set_ylim(ax_2d.get_ylim())
        ax_slice.invert_yaxis()
        ax_slice.set_yticks([])
        ax_slice.set_xlabel('强度', fontsize=8)
        ax_slice.set_title(f'13C切片\n({len(pk_idx)}峰)', fontsize=8)

        # 保存 (文件名含积分信息)
        if p.h_integral is not None:
            fname = f"{p.assigned_c_label}_{p.assigned_h_label}_integral{int(p.h_integral)}.png"
        else:
            fname = f"{p.assigned_c_label}_{p.assigned_h_label}.png"
        fpath = os.path.join(contour_dir, fname)
        fig.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        written.append(fpath)

    return written


def write_step_files(
    result: Interpretation2D,
    output_dir: str,
    spectra_2d: Dict[str, Spectrum2D] = None,
) -> List[str]:
    """将四步中间结果写入独立 CSV 文件, 供化学家在 Excel 中复核.

    文件命名带步骤序号:
      step1_hsqc_ch_pairs.csv
      step1_hsqc_ch_pairs/overlap/  (积分>=4 的重叠峰等高线图)
      step1_hsqc_ch_pairs/normal/   (非重叠峰等高线图)
      step1b_multiplet_analysis.csv  (每个 CHPair 的 1H 多重峰精细化分析)
      step2_cosy_spin_systems.csv
      step3_hmbc_connections.csv
      step4_quaternary_carbons.csv

    Args:
        result: interpret_2d_sequential 的返回值
        output_dir: 输出目录 (已存在)
        spectra_2d: 2D 谱数据字典 (用于绘制等高线图)

    Returns:
        实际写入的文件路径列表
    """
    written: List[str] = []

    # 步骤 1: HSQC C-H 对 (含 1H 积分值 + J 耦合)
    path1 = os.path.join(output_dir, 'step1_hsqc_ch_pairs.csv')
    with open(path1, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['编号', 'C_label', 'H_label', 'dH_ppm', 'dC_ppm', '1H积分', 'H个数',
                    '碳级数', 'J_couplings_Hz', 'h_neighbors_ppm',
                    'is_solvent_residual'])  # Part 2.B
        for i, p in enumerate(result.ch_pairs, 1):
            h_int = f'{p.h_integral:.3f}' if p.h_integral is not None else ''
            n_h = str(round(p.h_integral)) if p.h_integral is not None else ''
            j_str = '; '.join(f'{j:.1f}' for j in p.j_couplings) if p.j_couplings else ''
            hn_str = '; '.join(f'{h:.3f}' for h in p.h_neighbors) if p.h_neighbors else ''
            w.writerow([i, p.assigned_c_label, p.assigned_h_label,
                        f'{p.h_ppm:.4f}', f'{p.c_ppm:.4f}', h_int, n_h,
                        p.multiplicity, j_str, hn_str,
                        'Y' if p.is_solvent_residual else ''])
    written.append(path1)

    # 步骤 1 附加: 所有 C-H 对等高线图 (分 overlap/normal 两个子文件夹)
    if spectra_2d and 'HSQC' in spectra_2d:
        contour_files = _draw_hsqc_contour_all_pairs(
            spectra_2d['HSQC'], result.ch_pairs, output_dir)
        written.extend(contour_files)

    # 步骤 1b (新): 对每个 CHPair 做精细化 1H 多重峰分析 (对标 TopSpin Manual)
    if result.multiplet_analyses:
        path1b = os.path.join(output_dir, 'step1b_multiplet_analysis.csv')
        with open(path1b, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['编号', 'C_label', 'H_label', 'dH_ppm', 'dC_ppm',
                        '碳级数', 'multiplet_pattern', '子峰数',
                        'J1_Hz', 'J2_Hz', 'J3_Hz',
                        '多重峰中心_ppm', '半高宽_Hz',
                        '强度均匀性(0-1)', '质量等级', '说明'])
            for i, ma in enumerate(result.multiplet_analyses, 1):
                j = ma.j_couplings_Hz + [''] * 3  # 填到 3 列
                w.writerow([
                    i, ma.c_label, ma.h_label,
                    f'{ma.h_ppm:.4f}', f'{ma.c_ppm:.4f}',
                    ma.pattern, ma.num_subpeaks,
                    f'{j[0]:.2f}' if j[0] != '' else '',
                    f'{j[1]:.2f}' if j[1] != '' else '',
                    f'{j[2]:.2f}' if j[2] != '' else '',
                    f'{ma.multiplet_center_ppm:.4f}',
                    f'{ma.halfwidth_Hz:.2f}',
                    f'{ma.intensity_uniformity:.3f}',
                    ma.quality_grade,
                    ma.note,
                ])
        written.append(path1b)
    elif result.overlap_regions:
        # fallback: 旧的重叠峰分解分析 (无 multiplet 分析结果时)
        path1b = os.path.join(output_dir, 'step1b_overlap_regions.csv')
        with open(path1b, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['区域序号', 'dH_center', '总积分', '13C切片峰数',
                        '13C峰位置', '峰强度(归一化)', '分解后H/峰', '状态', '说明'])
            for i, r in enumerate(result.overlap_regions, 1):
                c_peaks_str = '; '.join(f'{c:.2f}' for c in r.c_peaks)
                c_int_str = '; '.join(f'{v:.3f}' for v in r.c_peak_intensities)
                h_per = f'{r.integral_total / r.n_peaks_in_13c_slice:.2f}' if r.n_peaks_in_13c_slice > 0 else ''
                w.writerow([i, f'{r.h_ppm_center:.4f}', f'{r.integral_total:.2f}',
                            r.n_peaks_in_13c_slice, c_peaks_str, c_int_str,
                            h_per, r.status, r.note])
        written.append(path1b)

    # 步骤 2: COSY 自旋系统
    path2 = os.path.join(output_dir, 'step2_cosy_spin_systems.csv')
    with open(path2, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['系统编号', 'label', 'H_ppm_列表', 'H_个数', 'COSY_links'])
        for ss in result.spin_systems:
            h_str = '; '.join(f'{h:.4f}' for h in ss.h_ppms)
            links_str = '; '.join(f'({a:.3f},{b:.3f})' for a, b in ss.cosy_links)
            w.writerow([ss.label, ss.label, h_str, len(ss.h_ppms), links_str])
    written.append(path2)

    # 步骤 3: HMBC 远程连接 (完整, 不截断, 含 bond_type)
    path3 = os.path.join(output_dir, 'step3_hmbc_connections.csv')
    with open(path3, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        # Part 2.A: 新增 is_oneJ_residue 列 (应全 False, 仅含远程相关)
        w.writerow(['序号', 'dH_ppm', 'dC_ppm', '键类型', '推断键数',
                    '推断依据', '连接描述', 'is_oneJ_residue'])
        for i, conn in enumerate(result.hmbc_connections, 1):
            w.writerow([i, f'{conn.h_ppm:.4f}', f'{conn.c_ppm:.4f}',
                        conn.bond_type, conn.inferred_distance or '',
                        conn.inference_basis, conn.connects_fragments or '',
                        'Y' if conn.is_oneJ_residue else ''])
    written.append(path3)

    # 步骤 4: 季碳
    path4 = os.path.join(output_dir, 'step4_quaternary_carbons.csv')
    with open(path4, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['序号', 'dC_ppm', '推断依据'])
        for i, c in enumerate(result.cq_shifts, 1):
            ev = result.cq_evidence[i - 1] if i - 1 < len(result.cq_evidence) else ''
            w.writerow([i, f'{c:.4f}', ev])
    written.append(path4)

    # 步骤 2b: TOCSY 自旋系统 (若有)
    if result.tocsy_spin_systems:
        path2t = os.path.join(output_dir, 'step2_tocsy_spin_systems.csv')
        with open(path2t, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['系统编号', 'label', 'H_ppm_列表', 'H_个数', 'TOCSY_links'])
            for ss in result.tocsy_spin_systems:
                h_str = '; '.join(f'{h:.4f}' for h in ss.h_ppms)
                links_str = '; '.join(f'({a:.3f},{b:.3f})' for a, b in ss.cosy_links)
                w.writerow([ss.label, ss.label, h_str, len(ss.h_ppms), links_str])
        written.append(path2t)

    # 步骤 5: NOE 空间近邻 (若有 NOESY/ROESY)
    if result.noe_contacts:
        path5 = os.path.join(output_dir, 'step5_noe_spatial_contacts.csv')
        with open(path5, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['序号', 'H1_ppm', 'H2_ppm', '强度', '置信度', '期望距离_A'])
            for i, nc in enumerate(result.noe_contacts, 1):
                w.writerow([i, f'{nc.h1_ppm:.4f}', f'{nc.h2_ppm:.4f}',
                            f'{nc.intensity:.3f}', nc.confidence,
                            f'{nc.expected_distance:.1f}'])
        written.append(path5)

    # 步骤 7: 归属完整性核对
    if result.assignment_validation:
        av = result.assignment_validation
        path7 = os.path.join(output_dir, 'step7_assignment_validation.csv')
        with open(path7, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['类型', 'd_ppm', '警告'])
            for c in av.unassigned_c:
                w.writerow(['未归属C', f'{c:.4f}', '无HSQC对应'])
            for h in av.unassigned_h:
                w.writerow(['未归属H', f'{h:.4f}', '可能为活泼氢'])
            for warn in av.warnings:
                w.writerow(['警告', '', warn])
            # 汇总行
            w.writerow(['汇总', '', f'C: {av.assigned_c}/{av.total_observed_c} 已归属; '
                        f'H: {av.assigned_h}/{av.total_observed_h} 已归属; '
                        f'完整性: {"OK" if av.is_complete else "有警告"}'])
        written.append(path7)

    return written
