"""Bruker 2D NMR 谱图读取与峰检测.

只读取不改写。原始数据放在 NMR/<folder>/<expNo>/ 下：
  - pdata/1/2rr          2D 处理后二进制数据 (实部)
  - pdata/1/procs        F2 处理参数
  - pdata/1/proc2s       F1 处理参数
  - acqus                F2 采集参数
  - acqu2s               F1 采集参数

输出 Spectrum2D：统一表示 COSY/HSQC/HMBC/ROESY/NOESY/TOCSY 2D 谱，
含元数据 + ppm 轴 + 自动检测的交叉峰列表。

依赖: nmrglue (可选, 未装时静默返回空结果, 不中断主流程)

对标: MestReNova 解谱指南 第五章 (二维谱关联 COSY/HSQC/HMBC/NOESY)
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---- 数据结构 ----

@dataclass
class Peak2D:
    """2D 谱单条交叉峰."""

    f1_ppm: float   # F1 (间接维) 化学位移 ppm
    f2_ppm: float   # F2 (直接维) 化学位移 ppm
    intensity: float = 0.0
    snr: Optional[float] = None  # 信噪比 (MAD-based)


@dataclass
class Spectrum2D:
    """单个 2D NMR 实验."""

    exp_no: str                    # e.g. "25"
    experiment_type: str           # 'COSY'/'HSQC'/'HMBC'/'ROESY'/'NOESY'/'TOCSY'
    f1_nucleus: str                # e.g. '13C' or '1H' (间接维核种)
    f2_nucleus: str                # e.g. '1H' (直接维核种)
    f1_ppm_scale: np.ndarray       # F1 ppm 轴
    f2_ppm_scale: np.ndarray       # F2 ppm 轴
    data: np.ndarray               # 2D 强度矩阵 (F1 x F2)
    peaks: List[Peak2D] = field(default_factory=list)
    solvent: str = ''
    field_mhz: float = 0.0
    pulse_program: str = ''
    source_path: str = ''
    name: str = ''


# ---- 错误类型 ----

class Bruker2DReadError(Exception):
    """2D Bruker 读取过程中的显式错误."""

    def __init__(self, msg: str, path: str = '', stage: str = ''):
        super().__init__(f'[{stage}] {msg} (path={path})')
        self.path = path
        self.stage = stage


# ---- 2D 实验类型检测 ----

# 2D 脉冲程序关键词白名单 (小写匹配)
_2D_PULSE_KEYWORDS = (
    'cosy', 'hmbc', 'hsqc', 'hmqc', 'noesy', 'roesy',
    'tocsy', 'mlev', 'dipsi', 'inv4',
)


def is_2d_pulse_program(pulse_program: str) -> bool:
    """判断脉冲程序是否为 2D 实验 (关键词匹配)."""
    pp = pulse_program.lower().strip()
    return any(kw in pp for kw in _2D_PULSE_KEYWORDS)


def detect_2d_experiment_type(
    pulse_program: str,
    f1_nucleus: str = '',
    f2_nucleus: str = '',
) -> str:
    """从脉冲程序和核种推断 2D 实验类型.

    参考 lucy-ng 的 _detect_experiment_type 逻辑.
    返回: 'HSQC'/'HMBC'/'COSY'/'TOCSY'/'NOESY'/'ROESY'
    失败时返回空串.
    """
    pp = pulse_program.lower().strip()

    # 顺序重要: HMBC 在 HSQC 之前 (部分脉冲程序名含两者)
    if 'hmbc' in pp:
        return 'HMBC'
    # inv4 系列: 含 lr/lplr/lrnd 为 HMBC, 否则为 HSQC
    if 'inv4' in pp:
        if 'lr' in pp or 'lplr' in pp or 'lrnd' in pp:
            return 'HMBC'
        return 'HSQC'
    if 'hsqc' in pp:
        return 'HSQC'
    # HMQC 等同 HSQC (一键异核关联)
    if 'hmqc' in pp:
        return 'HSQC'
    if 'cosy' in pp:
        return 'COSY'
    if 'tocsy' in pp or 'mlev' in pp or 'dipsi' in pp:
        return 'TOCSY'
    if 'noesy' in pp:
        return 'NOESY'
    if 'roesy' in pp:
        return 'ROESY'

    # 兜底: 根据核种判断同核/异核
    if f1_nucleus and f2_nucleus:
        if f1_nucleus != f2_nucleus:
            # 异核: 无法确定是 HSQC 还是 HMBC, 默认 HSQC
            return 'HSQC'
    return ''


# ---- 2D 数据读取 ----

def _strip_brackets(value: str) -> str:
    """去掉 Bruker 参数值中的尖括号."""
    if isinstance(value, str) and value.startswith('<') and value.endswith('>'):
        return value[1:-1]
    return value


def _get_acqus_param(dic: dict, key: str, default=None):
    """从 acqus 字典安全取参数."""
    try:
        value = dic['acqus'][key]
        return _strip_brackets(value) if isinstance(value, str) else value
    except (KeyError, TypeError):
        return default


def _get_acqu2s_param(dic: dict, key: str, default=None):
    """从 acqu2s 字典安全取参数 (F1 间接维)."""
    try:
        value = dic['acqu2s'][key]
        return _strip_brackets(value) if isinstance(value, str) else value
    except (KeyError, TypeError):
        return default


def read_bruker_2d_experiment(exp_dir: str, snr_floor: float = 5.0) -> Spectrum2D:
    """读取单个 Bruker 2D 实验目录.

    exp_dir 形如: '.../NMR/clxh100/25'
    snr_floor: 峰检测 SNR 阈值 (MAD-based, 默认 5.0)

    失败时抛 Bruker2DReadError.
    """
    try:
        import nmrglue as ng
    except ImportError as e:
        raise Bruker2DReadError(
            'nmrglue 未安装, 无法读取 2D 数据. '
            '安装: python -m pip install nmrglue -i https://pypi.tuna.tsinghua.edu.cn/simple',
            path=exp_dir, stage='import-nmrglue',
        ) from e

    exp_dir = os.path.abspath(exp_dir)
    if not os.path.isdir(exp_dir):
        raise Bruker2DReadError(f'实验目录不存在: {exp_dir}', exp_dir, stage='exp-dir')

    pdata_dir = os.path.join(exp_dir, 'pdata', '1')
    if not os.path.isdir(pdata_dir):
        raise Bruker2DReadError(f'pdata/1 目录不存在', pdata_dir, stage='pdata-dir')

    # 读取 2D 处理后数据
    try:
        dic, data = ng.bruker.read_pdata(pdata_dir)
    except Exception as e:
        raise Bruker2DReadError(f'nmrglue 读取 2D pdata 失败: {e}', pdata_dir, stage='read-pdata')

    if data is None or getattr(data, 'ndim', 0) != 2:
        raise Bruker2DReadError(
            f'期望 2D 谱, 实得 ndim={getattr(data, "ndim", None)}',
            pdata_dir, stage='validate-dim',
        )

    # 读取采集参数 (acqus + acqu2s)
    try:
        acqus_dic, _ = ng.bruker.read(exp_dir)
        dic.update(acqus_dic)
    except Exception:
        # acqu2s 可能不存在 (某些同核 2D), 不致命
        pass

    # 提取核种
    f2_nucleus = _get_acqus_param(dic, 'NUC1', '') or ''
    f1_nucleus = _get_acqu2s_param(dic, 'NUC1', '') or ''
    # 兜底: 同核实验 (COSY/ROESY) F1 核种 = F2 核种
    if not f1_nucleus:
        f1_nucleus = f2_nucleus

    # 脉冲程序
    pulse_program = _get_acqus_param(dic, 'PULPROG', '') or ''

    # 实验类型
    experiment_type = detect_2d_experiment_type(pulse_program, f1_nucleus, f2_nucleus)
    if not experiment_type:
        raise Bruker2DReadError(
            f'无法识别 2D 实验类型 (PULPROG={pulse_program})',
            exp_dir, stage='detect-type',
        )

    # 溶剂
    solvent = _get_acqus_param(dic, 'SOLVENT', '') or ''

    # 频率 (F2 主核)
    field_mhz = 0.0
    sfo1 = _get_acqus_param(dic, 'SFO1')
    if sfo1 is not None:
        try:
            field_mhz = float(sfo1)
        except (TypeError, ValueError):
            pass

    # 计算 ppm 轴
    try:
        udic = ng.bruker.guess_udic(dic, data)
        # F1 (dim=0, 间接维)
        uc_f1 = ng.fileiobase.uc_from_udic(udic, dim=0)
        f1_ppm_scale = np.asarray(uc_f1.ppm_scale(), dtype=np.float64)
        # F2 (dim=1, 直接维)
        uc_f2 = ng.fileiobase.uc_from_udic(udic, dim=1)
        f2_ppm_scale = np.asarray(uc_f2.ppm_scale(), dtype=np.float64)
    except Exception as e:
        raise Bruker2DReadError(f'ppm 轴计算失败: {e}', pdata_dir, stage='ppm-axis')

    # 转实数 (部分 2D 数据是复数)
    data_real = np.real(np.asarray(data, dtype=np.complex128)).astype(np.float64)

    # 自动峰检测
    peaks = pick_peaks_2d(data_real, f1_ppm_scale, f2_ppm_scale, snr_floor=snr_floor)

    exp_no = os.path.basename(exp_dir)
    name = os.path.basename(os.path.dirname(exp_dir))

    return Spectrum2D(
        exp_no=exp_no,
        experiment_type=experiment_type,
        f1_nucleus=f1_nucleus,
        f2_nucleus=f2_nucleus,
        f1_ppm_scale=f1_ppm_scale,
        f2_ppm_scale=f2_ppm_scale,
        data=data_real,
        peaks=peaks,
        solvent=solvent,
        field_mhz=field_mhz,
        pulse_program=pulse_program,
        source_path=exp_dir,
        name=name,
    )


# ---- 2D 峰检测 ----

def _compute_2d_noise_sigma(data: np.ndarray) -> float:
    """计算 2D 谱的 MAD-based 噪声标准差.

    MAD = median(|data - median(data)|), sigma = 1.4826 * MAD.
    相比角落 std 估计, MAD 对强交叉峰更鲁棒.
    """
    if data.size == 0:
        return 1.0
    med = float(np.median(data))
    mad = float(np.median(np.abs(data - med)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma == 0.0:
        max_abs = float(np.max(np.abs(data)))
        return (0.05 * max_abs / 5.0) if max_abs > 0 else 1.0
    return sigma


def pick_peaks_2d(
    data: np.ndarray,
    f1_ppm_scale: np.ndarray,
    f2_ppm_scale: np.ndarray,
    snr_floor: float = 5.0,
    min_separation: Tuple[int, int] = (3, 3),
) -> List[Peak2D]:
    """2D 谱自动峰检测 (MAD-SNR 阈值 + nmrglue connected 算法).

    Args:
        data: 2D 强度矩阵 (F1 x F2)
        f1_ppm_scale: F1 ppm 轴
        f2_ppm_scale: F2 ppm 轴
        snr_floor: SNR 阈值 (默认 5.0, 即峰强度 >= 5 * noise_sigma)
        min_separation: 最小峰间距 (点数, F1/F2)

    Returns:
        检测到的交叉峰列表 (按 F1 ppm 降序, F2 ppm 降序)
    """
    try:
        import nmrglue as ng
    except ImportError:
        return []

    sigma = _compute_2d_noise_sigma(data)
    abs_threshold = snr_floor * sigma

    peaks: List[Peak2D] = []

    try:
        # 正峰检测
        result = ng.peakpick.pick(
            data,
            pthres=abs_threshold,
            nthres=None,
            msep=min_separation,
            algorithm='connected',
            est_params=True,
            cluster=True,
            table=True,
        )
        for row in result:
            y_idx = int(row['Y_AXIS'])
            x_idx = int(row['X_AXIS'])
            if y_idx >= len(f1_ppm_scale) or x_idx >= len(f2_ppm_scale):
                continue
            f1_ppm = float(f1_ppm_scale[y_idx])
            f2_ppm = float(f2_ppm_scale[x_idx])
            intensity = float(data[y_idx, x_idx])
            peaks.append(Peak2D(
                f1_ppm=f1_ppm,
                f2_ppm=f2_ppm,
                intensity=intensity,
                snr=abs(intensity) / sigma if sigma > 0 else None,
            ))
    except Exception as e:
        # 检峰失败不致命, 返回空列表
        print(f'WARN: 2D 峰检测失败: {e}', file=sys.stderr)
        return []

    # 按 F1 降序, F2 降序排列
    peaks.sort(key=lambda p: (-p.f1_ppm, -p.f2_ppm))
    return peaks


# ---- 顶层接口 ----

def read_bruker_2d_folder(
    folder: str,
    snr_floor: float = 5.0,
    quiet: bool = False,
) -> Dict[str, Spectrum2D]:
    """读取样本目录下所有 2D NMR 实验.

    识别逻辑:
      1. 子目录 pdata/1/ 下有 '2rr' 文件 (2D 处理数据)
      2. 或 acqus 中 PULPROG 匹配 2D 脉冲程序白名单

    Args:
        folder: 样本目录 (如 .../NMR/clxh100)
        snr_floor: 峰检测 SNR 阈值
        quiet: 静默模式

    Returns:
        {experiment_type: Spectrum2D} 字典.
        若无 2D 数据则返回空 dict (触发 1D 回退).
        同类型多个实验取第一个.
    """
    folder = os.path.abspath(folder)
    if not os.path.isdir(folder):
        if not quiet:
            print(f'WARN: 2D 读取 - 目录不存在: {folder}', file=sys.stderr)
        return {}

    out: Dict[str, Spectrum2D] = {}

    for entry in sorted(os.listdir(folder)):
        exp_dir = os.path.join(folder, entry)
        if not os.path.isdir(exp_dir):
            continue
        if not entry.isdigit():
            continue

        pdata_dir = os.path.join(exp_dir, 'pdata', '1')
        # 判据 1: 有 2rr 文件 (2D 处理后数据)
        has_2rr = os.path.isfile(os.path.join(pdata_dir, '2rr'))
        if not has_2rr:
            continue

        # 判据 2: 排除已有 peaklist.xml 的 1D 实验 (不太可能有 2rr, 但防御)
        # 2D 实验通常没有 peaklist.xml (TopSpin 不自动为 2D 生成)

        try:
            spec = read_bruker_2d_experiment(exp_dir, snr_floor=snr_floor)
        except Bruker2DReadError as e:
            if not quiet:
                print(f'WARN: 2D exp {entry} 读取失败 (stage={e.stage}): {e}',
                      file=sys.stderr)
            continue
        except Exception as e:
            if not quiet:
                print(f'WARN: 2D exp {entry} 未知错误: {e}', file=sys.stderr)
            continue

        # 同类型取第一个
        if spec.experiment_type not in out:
            out[spec.experiment_type] = spec
            if not quiet:
                print(f'  2D exp {entry}: {spec.experiment_type} '
                      f'({spec.f1_nucleus}-{spec.f2_nucleus}) '
                      f'peaks={len(spec.peaks)} '
                      f'pulprog={spec.pulse_program}',
                      file=sys.stderr)

    return out
