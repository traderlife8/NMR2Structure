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
    bond_type: str = 'unknown'  # '2J'/'3J'/'unknown'
    connects_fragments: Optional[str] = None  # 如 'spin-A -> Cq-170.5'


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
class Interpretation2D:
    """2D 解谱的完整中间结果 (化学家可审查)."""

    ch_pairs: List[CHPair] = field(default_factory=list)
    spin_systems: List[SpinSystem] = field(default_factory=list)
    hmbc_connections: List[HMBCConnection] = field(default_factory=list)
    overlap_regions: List[OverlapRegion] = field(default_factory=list)
    n_cq_inferred: int = 0
    cq_shifts: List[float] = field(default_factory=list)
    cq_evidence: List[str] = field(default_factory=list)  # 每个季碳的推断依据
    notes: List[str] = field(default_factory=list)


# ---- 步骤 1: HSQC → C-H 对 ----

def _step1_hsqc_to_ch_pairs(
    spectra_2d: Dict[str, Spectrum2D],
    dept_summary: List[dict],
    proton_integrals: List[tuple] = None,
    h_tol: float = 0.3,
    c_tol: float = 2.5,
    min_snr: float = 20.0,
) -> List[CHPair]:
    """从 HSQC 谱提取 C-H 直接连接对, 结合 DEPT 和 1H 积分确定碳级数.

    HSQC F1=13C, F2=1H. 每个交叉峰 = 一个直接 C-H 连接.
    碳级数来源优先级: DEPT 推断 > 1H 积分推断 > 化学位移经验 > 'unknown'.

    Args:
        proton_integrals: [(center_ppm, integral_value), ...] 1H 积分数据
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

    ch_pairs: List[CHPair] = []
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

        ch_pairs.append(CHPair(
            h_ppm=round(h_ppm, 3),
            c_ppm=round(c_ppm, 3),
            multiplicity=mult,
            h_integral=round(h_integral, 3) if h_integral is not None else None,
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

    Args:
        hsqc_spectrum: HSQC 2D 谱数据
        ch_pairs: 已提取的 C-H 对列表
        integral_threshold: 积分阈值 (默认 4.0)
        h_slice_width: ¹H 切片宽度 (ppm, 默认 ±0.15)
        c_min_distance: ¹³C 峰最小间距 (ppm, 默认 2.0)
        h_cluster_tol: ¹H 聚类容差 (ppm, 默认 0.3)

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
        if n_peaks > 0:
            h_per_peak = integral_total / n_peaks
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
                    h_ppm=round(h_center, 3),
                    c_ppm=c_pk,
                    multiplicity=mult,
                    h_integral=round(h_per_peak, 2),
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


# ---- 步骤 2: COSY → 自旋系统 ----

def _step2_cosy_spin_systems(
    spectra_2d: Dict[str, Spectrum2D],
    ch_pairs: List[CHPair],
    h_tol: float = 0.3,
    min_snr: float = 15.0,
) -> List[SpinSystem]:
    """从 COSY 谱建立氢自旋系统.

    策略: 用 Union-Find 将 COSY 相关的 H 归入同一系统.
    仅保留含 >= 2 个 H 的系统 (单峰无耦合信息).
    """
    from .correlation_2d import extract_cosy_pairs

    if 'COSY' not in spectra_2d:
        return []

    cosy_pairs = extract_cosy_pairs(spectra_2d['COSY'], min_snr=min_snr)
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

    hmbc_pairs = extract_hmbc_pairs(spectra_2d['HMBC'], min_snr=min_snr)
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


# ---- 顶层接口 ----

def interpret_2d_sequential(
    spectra_2d: Dict[str, Spectrum2D],
    carbon_shifts_1d: List[float],
    proton_shifts_1d: List[float],
    dept_summary: List[dict],
    proton_integrals: List[tuple] = None,
    quiet: bool = False,
) -> Interpretation2D:
    """按化学家顺序逐步解读 2D NMR, 返回结构化中间结果.

    步骤:
    1. HSQC → C-H pairs (用 DEPT + 1H 积分确定 CH/CH2/CH3)
    2. COSY → 从 C-H pairs 中的 H 出发, 沿 COSY 相关建立自旋系统
    3. HMBC → 从已归属的 H 出发, 找 2-3 键外的 C, 建立片段间连接
    4. 推断季碳: 13C 中无 HSQC 对应的 = Cq

    Args:
        spectra_2d: {experiment_type: Spectrum2D} 字典
        carbon_shifts_1d: 1D 13C 观测位移列表
        proton_shifts_1d: 1D 1H 观测位移列表
        dept_summary: DEPT 推断 [{ppm, mult, sign}, ...]
        proton_integrals: [(center_ppm, integral_value), ...] 1H 积分数据
        quiet: 不打印进度

    Returns:
        Interpretation2D 完整中间结果
    """
    result = Interpretation2D()

    if not spectra_2d:
        result.notes.append('无 2D 数据, 跳过逐步解读')
        return result

    # 步骤 1: HSQC → C-H 对 (结合 DEPT + 1H 积分推断碳级数)
    ch_pairs = _step1_hsqc_to_ch_pairs(spectra_2d, dept_summary, proton_integrals)
    result.ch_pairs = ch_pairs

    # 步骤 1b: 重叠峰分解 (对标 5.9 流程: ¹³C 切片 + 分峰检测)
    if 'HSQC' in spectra_2d:
        overlap_regions = _step1b_resolve_overlapping(
            spectra_2d['HSQC'], ch_pairs)
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

    # 步骤 4 先做 (步骤 3 需要季碳列表)
    cq_results = _step4_infer_quaternary_carbons(carbon_shifts_1d, ch_pairs, dept_summary)
    cq_shifts = [ppm for ppm, _ in cq_results]
    result.cq_shifts = cq_shifts
    result.cq_evidence = [ev for _, ev in cq_results]
    result.n_cq_inferred = len(cq_shifts)

    # 步骤 2: COSY → 自旋系统
    spin_systems = _step2_cosy_spin_systems(spectra_2d, ch_pairs)
    result.spin_systems = spin_systems

    # 步骤 3: HMBC → 远程连接
    hmbc_connections = _step3_hmbc_connections(
        spectra_2d, ch_pairs, cq_shifts)
    result.hmbc_connections = hmbc_connections

    # 控制台仅输出一行摘要 (详细结果见 CSV 文件)
    if not quiet:
        n_overlap = len(result.overlap_regions)
        n_resolved = sum(1 for r in result.overlap_regions if r.status == 'resolved')
        overlap_str = f' / 重叠峰 {n_overlap} 区(分解 {n_resolved})' if n_overlap else ''
        print(f"  2D 解读完成: HSQC {len(ch_pairs)} 对 / "
              f"COSY {len(spin_systems)} 系统 / "
              f"HMBC {len(hmbc_connections)} 连接 / "
              f"季碳 {len(cq_shifts)} 个{overlap_str}",
              file=sys.stderr)

    # 汇总注释
    result.notes.append(
        f'HSQC: {len(ch_pairs)} 个 C-H 对; '
        f'COSY: {len(spin_systems)} 个自旋系统; '
        f'HMBC: {len(hmbc_connections)} 个远程连接; '
        f'季碳: {len(cq_shifts)} 个'
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

        # --- 中央 2D 等高线图 ---
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

    # 步骤 1: HSQC C-H 对 (含 1H 积分值)
    path1 = os.path.join(output_dir, 'step1_hsqc_ch_pairs.csv')
    with open(path1, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['编号', 'C_label', 'H_label', 'dH_ppm', 'dC_ppm', '1H积分', 'H个数', '碳级数'])
        for i, p in enumerate(result.ch_pairs, 1):
            h_int = f'{p.h_integral:.3f}' if p.h_integral is not None else ''
            n_h = str(round(p.h_integral)) if p.h_integral is not None else ''
            w.writerow([i, p.assigned_c_label, p.assigned_h_label,
                        f'{p.h_ppm:.4f}', f'{p.c_ppm:.4f}', h_int, n_h, p.multiplicity])
    written.append(path1)

    # 步骤 1 附加: 所有 C-H 对等高线图 (分 overlap/normal 两个子文件夹)
    if spectra_2d and 'HSQC' in spectra_2d:
        contour_files = _draw_hsqc_contour_all_pairs(
            spectra_2d['HSQC'], result.ch_pairs, output_dir)
        written.extend(contour_files)

    # 步骤 1b: 重叠峰分解分析结果
    if result.overlap_regions:
        path1b = os.path.join(output_dir, 'step1b_overlap_analysis.csv')
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

    # 步骤 3: HMBC 远程连接 (完整, 不截断)
    path3 = os.path.join(output_dir, 'step3_hmbc_connections.csv')
    with open(path3, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['序号', 'dH_ppm', 'dC_ppm', '键类型', '连接描述'])
        for i, conn in enumerate(result.hmbc_connections, 1):
            w.writerow([i, f'{conn.h_ppm:.4f}', f'{conn.c_ppm:.4f}',
                        conn.bond_type, conn.connects_fragments or ''])
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

    return written
