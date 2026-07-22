"""读取 Bruker 已处理谱图 (pdata/1/1r) 与重采样到 ML 网格.

设计要点:
- 1r 是 TopSpin 处理后的最终实谱 (int32/float, 已做 FT + 相位 + 基线)
  相比 peaklist.xml 的竹签, 1r 保留完整线形, 是画"真实谱"和喂 ML 模型的最佳来源
- 用 nmrglue.bruker.read_pdata 读, 用 procs 的 SF/SW_p/OFFSET/SI 计算 ppm 轴
- 只读原始数据, 不改写

用途:
1. read_processed_spectrum / read_processed_folder: 读真实谱曲线, 供报告绘图
2. resample_to_grid: 把真实谱重采样到均匀 ppm 网格
   - NMR2Struct: HNMR 28000 维 [-2,12], CNMR 80 维 [3.4,231.3]
   - ChefNMR:    h_10k 10000 维 [-2,10],  c_80 80 维 [3.42,231.3]
   两个 ML 库都吃这种网格向量, 重采样真实谱比竹签装箱保真度更高
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np


class ProcessedReadError(Exception):
    """读取已处理谱图失败. 携带路径与所处阶段, 便于定位."""

    def __init__(self, message: str, path: str = '', stage: str = ''):
        super().__init__(message)
        self.path = path
        self.stage = stage


@dataclass
class ProcessedSpectrum:
    """单条已处理实谱 (来自 pdata/1/1r)."""

    exp_no: int
    nucleus: str                 # '1H' / '13C' / '13C-DEPT' / '?'
    ppm: np.ndarray              # ppm 轴 (与 intensity 等长, 通常高场在右)
    intensity: np.ndarray        # 真实处理后强度 (已按 NC_proc 缩放)
    sf_mhz: float = 0.0          # 观测频率 MHz
    sw_hz: float = 0.0           # 谱宽 Hz
    pdata_dir: str = ''

    def ppm_range(self) -> tuple:
        """返回 (min_ppm, max_ppm)."""
        return float(self.ppm.min()), float(self.ppm.max())


def read_processed_spectrum(pdata_dir: str, exp_no: int = -1, nucleus: str = '?') -> ProcessedSpectrum:
    """读单个 pdata/1 目录的 1r 实谱.

    pdata_dir: 形如 .../<exp>/pdata/1
    exp_no / nucleus: 由调用方 (通常 io_bruker 已解析) 传入, 仅作标注
    失败时抛 ProcessedReadError, 明确阶段.
    """
    try:
        import nmrglue as ng
    except ImportError as e:
        raise ProcessedReadError(
            'nmrglue 未安装. 安装: python -m pip install nmrglue -i https://pypi.tuna.tsinghua.edu.cn/simple',
            path=pdata_dir, stage='import-nmrglue') from e

    r1 = os.path.join(pdata_dir, '1r')
    if not os.path.isfile(r1):
        raise ProcessedReadError(f'未找到 1r 文件: {r1}', path=pdata_dir, stage='check-1r')

    try:
        dic, data = ng.bruker.read_pdata(pdata_dir)
    except Exception as e:  # nmrglue 内部多种异常, 统一包装并保留上下文
        raise ProcessedReadError(f'nmrglue 读取 pdata 失败: {e}', path=pdata_dir, stage='read-pdata') from e

    if data is None or getattr(data, 'ndim', 0) != 1:
        raise ProcessedReadError(
            f'期望 1D 谱, 实得 shape={getattr(data, "shape", None)}',
            path=pdata_dir, stage='validate-dim')

    # ppm 轴: 用 nmrglue 的单位换算 (基于 procs 的 SF/SW/OFFSET/SI)
    try:
        udic = ng.bruker.guess_udic(dic, data)
        u0 = udic[0]
        uc = ng.fileiobase.unit_conversion(
            u0['size'], u0['complex'], u0['sw'], u0['obs'], u0['car'])
        ppm = np.asarray(uc.ppm_scale(), dtype=np.float64)
    except Exception as e:
        raise ProcessedReadError(f'ppm 轴计算失败: {e}', path=pdata_dir, stage='ppm-axis') from e

    intensity = np.asarray(data, dtype=np.float64)
    if ppm.shape[0] != intensity.shape[0]:
        raise ProcessedReadError(
            f'ppm 轴与强度长度不一致: {ppm.shape[0]} vs {intensity.shape[0]}',
            path=pdata_dir, stage='align')

    procs = dic.get('procs', {}) if isinstance(dic, dict) else {}
    sf = float(procs.get('SF', 0.0) or 0.0)
    sw = float(procs.get('SW_p', 0.0) or 0.0)

    return ProcessedSpectrum(
        exp_no=exp_no, nucleus=nucleus, ppm=ppm, intensity=intensity,
        sf_mhz=sf, sw_hz=sw, pdata_dir=pdata_dir,
    )


def read_processed_folder(folder: str, specs: Optional[List] = None) -> Dict[int, ProcessedSpectrum]:
    """读一个化合物目录下所有实验的 1r 实谱.

    folder: 含 21/22/... 数字子目录的 Bruker 目录
    specs:  可选, io_bruker.read_bruker_folder 的结果, 用于取每个 exp 的核种标注
            (避免在此重复核种识别逻辑)
    返回 {exp_no: ProcessedSpectrum}. 单个实验读失败会跳过并记录到 stderr, 不中断整体.
    """
    import sys

    if not os.path.isdir(folder):
        raise ProcessedReadError(f'原始数据目录不存在: {folder}', path=folder, stage='check-folder')

    # exp_no -> nucleus 标注 (来自已解析的 specs)
    # 注意: io_bruker 的 exp_no 是 str (如 "21"), 这里统一转 int 以匹配 folder 数字子目录
    nucleus_by_exp: Dict[int, str] = {
        int(s.exp_no): s.nucleus for s in (specs or []) if str(s.exp_no).isdigit()
    }

    out: Dict[int, ProcessedSpectrum] = {}
    for name in sorted(os.listdir(folder)):
        if not name.isdigit():
            continue
        pdata_dir = os.path.join(folder, name, 'pdata', '1')
        if not os.path.isfile(os.path.join(pdata_dir, '1r')):
            continue
        exp_no = int(name)
        try:
            ps = read_processed_spectrum(
                pdata_dir, exp_no=exp_no, nucleus=nucleus_by_exp.get(exp_no, '?'))
            out[exp_no] = ps
        except ProcessedReadError as e:
            print(f'WARN: exp {exp_no} 1r 读取失败 (stage={e.stage}): {e}', file=sys.stderr)
    return out


def resample_to_grid(
    ppm: np.ndarray,
    intensity: np.ndarray,
    grid_min: float,
    grid_max: float,
    n_points: int,
    normalize: bool = True,
) -> np.ndarray:
    """把真实谱重采样到均匀 ppm 网格 (供 ML 模型输入).

    ppm / intensity: 真实谱 (长度一致)
    grid_min/grid_max/n_points: 目标网格 (如 ChefNMR c_80 = 3.42, 231.3, 80)
    normalize: True 时把网格向量归一化到 max=1 (多数 ML 模型期望)
    网格范围外的区域填 0. 返回长度 n_points 的 float32 向量.
    """
    if ppm.shape[0] != intensity.shape[0]:
        raise ValueError(f'ppm 与 intensity 长度不一致: {ppm.shape[0]} vs {intensity.shape[0]}')
    if n_points <= 0:
        raise ValueError(f'n_points 必须为正: {n_points}')

    # np.interp 要求 xp 单调递增, 故按 ppm 升序排序
    order = np.argsort(ppm)
    xp = ppm[order]
    fp = intensity[order]

    grid = np.linspace(grid_min, grid_max, n_points)
    # 网格范围外用 left/right=0 填充, 避免边界外插造成假信号
    resampled = np.interp(grid, xp, fp, left=0.0, right=0.0)

    # 负强度 (基线波动/相位残留) 截断为 0, ML 网格通常只表达正峰
    resampled = np.clip(resampled, 0.0, None)

    if normalize:
        peak = resampled.max()
        if peak > 0:
            resampled = resampled / peak
    return resampled.astype(np.float32)
