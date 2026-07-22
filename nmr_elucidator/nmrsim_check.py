"""nmrsim 反向校验层.

角色：
- nmrsim 本身不能从 SMILES 预测化学位移，它是从「化学位移 + 耦合常数」生成谱的工具
- 反向校验的价值在于：给定一组化学位移 + 耦合表，生成模拟谱图，与实测谱做形状比对

本模块功能：
1. simulate_peaks_to_spectrum: 用 nmrsim 把 (shift, J_Hz, intensity) 列表生成谱图数据
2. compare_spectra: 把候选预测峰列表（来自 retrieve）转成模拟谱，与实测谱图数据对比

注意：
- 谱图数据通常从 Bruker 1r 文件读取；这里只接受 (shift, intensity) 列表
- 1r 二进制文件解析较为复杂，本版本仅基于峰表做对比，不读取 1r
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class SimulatedPeak:
    """用于反向校验的 1H 峰."""
    shift: float  # ppm
    intensity: float = 1.0
    j_hz: float = 0.0  # 多重峰的 J 耦合常数（0 表示单峰）
    multiplicity: str = 's'


def simulate_peaks_to_spectrum(
    peaks: List[SimulatedPeak],
    sweep_width_ppm: float = 12.0,
    n_points: int = 4096,
    field_mhz: float = 500.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """用 nmrsim 生成 1H 模拟谱.

    返回 (ppm_axis, intensity_array). ppm_axis 降序（与 Bruker 一致）.
    """
    # nmrsim 0.6 API: Multiplet + Spectrum
    from nmrsim import Multiplet
    try:
        from nmrsim import Spectrum
    except ImportError:
        # 兼容老版本
        from nmrsim.spectrum import Spectrum

    if not peaks:
        return np.array([]), np.array([])

    # 生成 ppm 轴（降序）
    ppm = np.linspace(sweep_width_ppm, 0.0, n_points)
    # Hz 轴
    hz = ppm * field_mhz

    multiplets = []
    for p in peaks:
        # nmrsim Multiplet 接受 J 列表（Hz）
        # shift 单位 ppm
        if p.j_hz > 0:
            # 用一阶近似：n+1 重峰, intensity 总和 = p.intensity
            mult = int(round(p.intensity))  # 假设整数积分
            j_list = [p.j_hz] * max(1, mult - 1)
            m = Multiplet(p.shift, j_list, peak_multiplicity=mult)
        else:
            m = Multiplet(p.shift, [], peak_multiplicity=1)
        multiplets.append(m)
    try:
        spec = Spectrum(multiplets)
    except Exception:
        return ppm, np.zeros_like(ppm)
    # nmrsim 返回强度，按 Hz 取样
    # Spectrum.peaklist() 返回 list of (ppm, intensity)
    # 我们需要连续谱 → 把每个 multiplet 在 hz 上做高斯加宽
    arr = np.zeros_like(hz)
    for m in multiplets:
        c_hz = m.v * field_mhz  # 中心 Hz
        for split_ppm, split_int in m.peaklist():
            s_hz = split_ppm * field_mhz
            # 高斯加宽（半高宽 1 Hz）
            sigma = 0.5
            arr += split_int * np.exp(-((hz - s_hz) ** 2) / (2 * sigma ** 2))
    return ppm, arr


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """两个等长谱图的余弦相似度."""
    if a.size == 0 or b.size == 0:
        return 0.0
    if a.shape != b.shape:
        n = min(a.size, b.size)
        a, b = a[:n], b[:n]
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def build_predicted_peaks_from_entry(
    entry,  # NmrShiftDb2Entry
    nucleus: str = '1H',
    min_intensity: float = 0.1,
) -> List[SimulatedPeak]:
    """从 nmrshiftdb2 条目里抽 (shift, intensity) 转 SimulatedPeak 列表.

    注意 nmrshiftdb2 里 <Spectrum> 字段的 integral 是归一化值，
    与实测的相对积分可能不同。这里只取峰值位置.
    """
    if nucleus == '13C':
        src = entry.spectrum_13c
    else:
        src = entry.spectrum_1h
    out: List[SimulatedPeak] = []
    for sig in src:
        if abs(sig.shift) < min_intensity:
            continue
        out.append(SimulatedPeak(shift=sig.shift, intensity=1.0, j_hz=0.0,
                                 multiplicity=sig.multiplicity or 's'))
    return out


def check_candidate_with_nmrsim(
    observed_peaks_ppm: List[float],
    predicted_peaks_ppm: List[float],
    sweep_width_ppm: float = 12.0,
    n_points: int = 4096,
    field_mhz: float = 500.0,
) -> float:
    """用 nmrsim 模拟两组峰，转成谱图后计算余弦相似度.

    返回 0..1 之间, 越大越好. 用于辅助判断候选谱与实测谱的整体形状兼容性.

    注：当前用 1H 简化模拟. 真实情况下应考虑耦合常数 (J).
    """
    obs_sim = [SimulatedPeak(shift=p, intensity=1.0, j_hz=0.0) for p in observed_peaks_ppm]
    pred_sim = [SimulatedPeak(shift=p, intensity=1.0, j_hz=0.0) for p in predicted_peaks_ppm]
    _, a = simulate_peaks_to_spectrum(obs_sim, sweep_width_ppm, n_points, field_mhz)
    _, b = simulate_peaks_to_spectrum(pred_sim, sweep_width_ppm, n_points, field_mhz)
    return cosine_similarity(a, b)