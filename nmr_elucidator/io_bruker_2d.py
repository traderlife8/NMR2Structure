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


# ---- 2D 投影与切片读取（产物 B/C） ----

@dataclass
class ProjectionPeak:
    """从 1D 投影上检测到的峰 (辅助 13C 谱表增强)."""

    shift: float  # ppm (主轴)
    intensity: float = 0.0
    snr: Optional[float] = None
    source: str = ''  # 'hmbc_col_proj' / 'hmbc_row_proj' / 'hsqc_row_slice'


@dataclass
class RowSliceResult:
    """HSQC 单行 (13C ppm) 切片的读取结果."""

    hsqc_exp: str                # 原始 HSQC exp
    row_pprocno2: int            # used_from 中 PROCNO2 (行号)
    f1_ppm: float                # 该行的 13C ppm (从 ppm 轴推算)
    ppm_scale: np.ndarray         # F2 ppm 轴
    intensity: np.ndarray         # F2 强度轴
    detected_peaks: List[ProjectionPeak] = field(default_factory=list)
    note: str = ''


def _read_1d_projection_via_nmrglue(pdata_dir: str) -> Optional[Tuple[np.ndarray, np.ndarray, dict]]:
    """用 nmrglue 读一个 1D 投影数据 (如 HMBC 列投影).

    Args:
        pdata_dir: 含 1r + proc/procs 的目录

    Returns:
        (ppm_scale, intensity, dic) 或 None
    """
    if not os.path.isdir(pdata_dir):
        return None
    one_r = os.path.join(pdata_dir, '1r')
    if not os.path.isfile(one_r):
        return None
    try:
        import nmrglue as ng
        dic, data = ng.bruker.read_pdata(pdata_dir)
    except Exception:
        return None
    if data is None:
        return None
    data_arr = np.asarray(data)
    if data_arr.ndim == 2:
        # 复数或双实数，取绝对值或第一列
        data_arr = np.real(data_arr[:, 0]) if data_arr.shape[1] > 0 else np.real(data_arr[0])
    else:
        data_arr = np.real(data_arr)
    # 重建 ppm 轴
    try:
        udic = ng.bruker.guess_udic(dic, data)
        uc = ng.fileiobase.uc_from_udic(udic, dim=0)
        ppm_scale = np.asarray(uc.ppm_scale(), dtype=np.float64)
    except Exception:
        return None
    # 截断对齐
    n = min(len(ppm_scale), len(data_arr))
    return ppm_scale[:n], data_arr[:n], dic


def _build_f1_ppm_from_acqu2s(
    acqu2s_path: str,
    n_points: int,
) -> Optional[np.ndarray]:
    """从 acqu2s / proc2s 重建 F1 ppm 轴.

    用于投影 / 切片文件没有正确 udic 的场景 (如 HMBC 的 pdata/99 列投影).

    公式: ppm[i] = OFFSET + SIGF1/2 - SIGF1 * i / (n-1) (传统 Bruker ppm 公式)
           注: Bruker 的 ppm 轴是 高场 → 低场 (从大到小), 所以 ppm[0] = OFFSET + SIGF1/2.

    支持字段 (按优先级):
      - OFFSET (proc2s 或 acqu2s)
      - SIGF1 (proc2s) 或 SW (acqu2s)

    Args:
        acqu2s_path: acqu2s 或 proc2s 文件路径
        n_points: F1 维点数

    Returns:
        np.ndarray (F1 ppm) 或 None
    """
    if not os.path.isfile(acqu2s_path):
        return None
    try:
        info = {}
        with open(acqu2s_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                # 同时匹配 OFFSET / SW / SIGF1
                for key in ('##$OFFSET=', '##$SW=', '##$SW_h=',
                             '##$SIGF1=', '##$SIGF2=', '##$SF=',
                             '##$BF1=', '##$BF2=', '##$SI=', '##$STSI='):
                    if line.startswith(key):
                        try:
                            info[key.replace('##$', '').replace('=', '')] = float(line[len(key):])
                        except ValueError:
                            pass
        offset = info.get('OFFSET')
        sw = info.get('SIGF1') or info.get('SW')
        if offset is None or sw is None:
            return None
        if n_points <= 1:
            return np.array([offset], dtype=np.float64)
        idx = np.arange(n_points)
        ppm = offset + sw / 2.0 - sw * idx / (n_points - 1)
        return ppm.astype(np.float64)
    except Exception:
        return None


def read_hmbc_f1_projection(
    hmbc_exp_dir: str,
    snr_floor: float = 5.0,
    prefer: str = 'col',  # 'col' = pdata/99 (column), 'row' = pdata/999 (row subset)
) -> List[ProjectionPeak]:
    """从 HMBC 投影提取 13C 峰 (产物 B).

    HMBC 的 F1 (13C) 维投影位于:
      - pdata/99  : positive projection of all 2048 columns (全列投影 = 13C 维投影)
      - pdata/999 : positive projection of rows 298-334 (某区域行投影)

    算法: nmrglue 读投影数据 + 重建 ppm 轴 + find_peaks 检峰.
    返回 ProjectionPeak 列表 (按强度递减).

    Args:
        hmbc_exp_dir: HMBC 实验目录 (如 .../clxh100/27)
        snr_floor: SNR 阈值 (默认 5.0)
        prefer: 优先选择 'col' (pdata/99) 还是 'row' (pdata/999)

    Returns:
        [(dC_ppm, intensity, snr), ...] 按强度降序
    """
    candidates = []
    if prefer == 'col':
        order = ['99', '999']
    else:
        order = ['999', '99']
    for pprocno in order:
        pdata_dir = os.path.join(hmbc_exp_dir, 'pdata', pprocno)
        if os.path.isdir(pdata_dir):
            candidates.append((pprocno, pdata_dir))

    if not candidates:
        return []

    pprocno, pdata_dir = candidates[0]
    src_label = f'hmbc_{pprocno}_proj'
    result = _read_1d_projection_via_nmrglue(pdata_dir)
    if result is None:
        return []
    ppm_scale, intensity, _ = result

    # 检查 ppm_scale 数值是否合理 (13C 谱区 0-250 ppm).
    # 若 > 1000 ppm 则表示 nmrglue 用错了维度, fallback 到从 acqu2s/proc2s 重算.
    bruker_offset = 0.0  # Bruker OFFSET (proc 设定值)
    o1_ppm_tms = 0.0     # O1 ppm (TMS 参考) = (SFO1 - SF) / SF * 1e6
    if ppm_scale.size > 0 and (np.max(np.abs(ppm_scale)) > 1000.0 or np.min(ppm_scale) < -500.0):
        acqu2s = os.path.join(hmbc_exp_dir, 'acqu2s')
        proc2s = os.path.join(hmbc_exp_dir, 'pdata', '1', 'proc2s')
        rebuilt = None
        bruker_offset = 0.0
        o1_ppm_tms = 0.0
        sigf1_for_correction = 0.0
        for cand_path in (acqu2s, proc2s):
            # 同时读取 OFFSET 和 SFO1/SF 以便后续校正
            offset_val = 0.0
            sfo1_val = 0.0
            sf_val = 0.0
            sigf1_val = 0.0
            sw_val = 0.0
            try:
                with open(cand_path, 'r', encoding='utf-8', errors='replace') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('##$OFFSET='):
                            try:
                                offset_val = float(line[len('##$OFFSET='):])
                            except ValueError:
                                pass
                        elif line.startswith('##$SFO1='):
                            try:
                                sfo1_val = float(line[len('##$SFO1='):])
                            except ValueError:
                                pass
                        elif line.startswith('##$SF='):
                            try:
                                sf_val = float(line[len('##$SF='):])
                            except ValueError:
                                pass
                        elif line.startswith('##$SIGF1='):
                            try:
                                sigf1_val = float(line[len('##$SIGF1='):])
                            except ValueError:
                                pass
                        elif line.startswith('##$SW=') and not line.startswith('##$SW_'):
                            try:
                                sw_val = float(line[len('##$SW='):])
                            except ValueError:
                                pass
            except OSError:
                offset_val = sfo1_val = sf_val = sigf1_val = sw_val = 0.0
            rebuilt = _build_f1_ppm_from_acqu2s(cand_path, len(intensity))
            if rebuilt is not None:
                bruker_offset = offset_val
                if sfo1_val > 0 and sf_val > 0:
                    o1_ppm_tms = (sfo1_val - sf_val) / sf_val * 1e6
                sigf1_for_correction = sigf1_val if sigf1_val > 0 else sw_val
                break
        if rebuilt is not None and len(rebuilt) == len(ppm_scale):
            # 校正公式 (经验证): Bruker ppm = TMS ppm + SIGF1/2
            #   => TMS ppm = Bruker ppm - SIGF1/2
            # 经验证 (clxh100 HMBC, OFFSET=225.98, SIGF1=228.63):
            #   Bruker 287.553 -> TMS 173.23 (匹配 main 173.32) ✓
            #   Bruker 264.087 -> TMS 149.77 (匹配 main 150.11) ✓
            #   Bruker 184.748 -> TMS 70.43  (匹配 main 70.62)  ✓
            sw_half_ppm = sigf1_for_correction / 2.0
            ppm_scale = rebuilt - sw_half_ppm

    # 计算 MAD-based sigma
    med = float(np.median(intensity))
    mad = float(np.median(np.abs(intensity - med)))
    sigma = 1.4826 * mad
    if sigma <= 0:
        max_abs = float(np.max(np.abs(intensity))) if intensity.size > 0 else 0
        sigma = (0.05 * max_abs / 5.0) if max_abs > 0 else 1.0

    peaks: List[ProjectionPeak] = []
    try:
        from scipy.signal import find_peaks
        # 投影数据已是求和增强, 信号可能较强, 阈值可以略低
        abs_thresh = max(snr_floor * sigma * 0.5, 0.0)
        # 估计最小峰间距 (按 ppm 轴分辨率)
        if len(ppm_scale) > 1:
            ppm_per_pt = abs(float(ppm_scale[1] - ppm_scale[0]))
        else:
            ppm_per_pt = 1.0
        # 最小峰间距 1 ppm (13C 上是合理的)
        min_dist_pts = max(1, int(1.0 / max(ppm_per_pt, 1e-6)))
        pk_idx, props = find_peaks(
            np.abs(intensity), height=abs_thresh, distance=min_dist_pts,
        )
        for idx in pk_idx:
            try:
                idx_i = int(idx)
            except (TypeError, ValueError):
                continue
            if idx_i < 0 or idx_i >= len(ppm_scale):
                continue
            shift = float(ppm_scale[idx_i])
            inten = float(intensity[idx_i])
            peaks.append(ProjectionPeak(
                shift=round(shift, 3),
                intensity=round(inten, 3),
                snr=round(abs(inten) / sigma, 2) if sigma > 0 else None,
                source=src_label,
            ))
    except ImportError:
        # scipy 不可用时退化为简单峰值检测
        for i in range(1, len(intensity) - 1):
            if intensity[i] > intensity[i - 1] and intensity[i] >= intensity[i + 1]:
                if abs(intensity[i]) > snr_floor * sigma:
                    peaks.append(ProjectionPeak(
                        shift=round(float(ppm_scale[i]), 3),
                        intensity=round(float(intensity[i]), 3),
                        snr=round(float(intensity[i]) / sigma, 2) if sigma > 0 else None,
                        source=src_label,
                    ))
    except Exception:
        return []

    # 按强度绝对值降序
    peaks.sort(key=lambda p: -abs(p.intensity))
    return peaks


def read_hsqc_row_slice(
    hsqc_exp_dir: str,
    pprocno: str,
    row_idx: Optional[int] = None,
    snr_floor: float = 5.0,
) -> Optional[RowSliceResult]:
    """读 HSQC 单行切片 (产物 C).

    Args:
        hsqc_exp_dir: HSQC 实验目录 (如 .../clxh100/26)
        pprocno: pdata 子目录名 (如 '1000' 或 '1001')
        row_idx: 原始 2D 矩阵的行号 (从 used_from 的 PROCNO2 读取); 若 None
                则从 pdata/<pprocno>/used_from 自动读取
        snr_floor: SNR 阈值

    Returns:
        RowSliceResult 或 None (若切片不存在 / 数据读取失败)
    """
    pdata_dir = os.path.join(hsqc_exp_dir, 'pdata', pprocno)
    if row_idx is None:
        used_from_path = os.path.join(pdata_dir, 'used_from')
        row_idx = read_used_from_procno2(used_from_path)
    result = _read_1d_projection_via_nmrglue(pdata_dir)
    if result is None:
        return None
    ppm_scale, intensity, dic = result

    # 检查 ppm_scale 数值是否合理 (HSQC F2 是 1H 谱区, 0-15 ppm).
    # 若 > 1000 则 fallback 从 acqu2s 重算 (F2 维) 或重读 2D 取 F2 ppm 轴.
    if ppm_scale.size > 0 and (np.max(np.abs(ppm_scale)) > 1000.0 or np.min(ppm_scale) < -500.0):
        # 切片是 F2 维投影, 重读 2D 矩阵取 F2 ppm 轴
        try:
            import nmrglue as ng
            full_pdata = os.path.join(hsqc_exp_dir, 'pdata', '1')
            full_dic, full_data = ng.bruker.read_pdata(full_pdata)
            udic = ng.bruker.guess_udic(full_dic, full_data)
            uc_f2 = ng.fileiobase.uc_from_udic(udic, dim=1)
            f2_ppm_scale = np.asarray(uc_f2.ppm_scale(), dtype=np.float64)
            if len(f2_ppm_scale) == len(ppm_scale):
                ppm_scale = f2_ppm_scale
        except Exception:
            pass

    # 从 used_from 推断该行的 13C ppm (通过原始 2D 矩阵的行号)
    # 用 dic + nmrglue 的 guess_udic 计算 ppm_scale, 用行号索引
    f1_ppm = float('nan')
    if row_idx is not None and row_idx > 0:
        try:
            import nmrglue as ng
            # 重新读原始 2D 矩阵以拿到 F1 ppm 轴
            full_pdata = os.path.join(hsqc_exp_dir, 'pdata', '1')
            full_dic, full_data = ng.bruker.read_pdata(full_pdata)
            udic = ng.bruker.guess_udic(full_dic, full_data)
            uc_f1 = ng.fileiobase.uc_from_udic(udic, dim=0)
            f1_ppm_scale = np.asarray(uc_f1.ppm_scale(), dtype=np.float64)
            if 0 < row_idx < len(f1_ppm_scale):
                f1_ppm = float(f1_ppm_scale[row_idx])
        except Exception:
            pass

    # 找峰
    med = float(np.median(intensity))
    mad = float(np.median(np.abs(intensity - med)))
    sigma = 1.4826 * mad
    if sigma <= 0:
        max_abs = float(np.max(np.abs(intensity))) if intensity.size > 0 else 0
        sigma = (0.05 * max_abs / 5.0) if max_abs > 0 else 1.0

    detected: List[ProjectionPeak] = []
    try:
        from scipy.signal import find_peaks
        abs_thresh = max(snr_floor * sigma * 0.5, 0.0)
        pk_idx, _ = find_peaks(
            np.abs(intensity), height=abs_thresh,
        )
        for idx in pk_idx:
            try:
                idx_i = int(idx)
            except (TypeError, ValueError):
                continue
            if idx_i < 0 or idx_i >= len(ppm_scale):
                continue
            detected.append(ProjectionPeak(
                shift=round(float(ppm_scale[idx_i]), 3),
                intensity=round(float(intensity[idx_i]), 3),
                snr=round(float(intensity[idx_i]) / sigma, 2) if sigma > 0 else None,
                source=f'hsqc_row_{row_pprocno2}_slice',
            ))
    except Exception:
        pass

    detected.sort(key=lambda p: -abs(p.intensity))
    return RowSliceResult(
        hsqc_exp=os.path.basename(hsqc_exp_dir),
        row_pprocno2=row_idx if row_idx is not None else -1,
        f1_ppm=round(f1_ppm, 3),
        ppm_scale=ppm_scale,
        intensity=intensity,
        detected_peaks=detected,
    )


def discover_hsqc_row_slices(hsqc_exp_dir: str) -> List[Tuple[str, str]]:
    """探测 HSQC 实验下所有 pdata 子目录 (用于发现行切片).

    Returns:
        [(pprocno, used_from 路径), ...] 例如 [('1000', '...pdata/1000/used_from'), ...]
    """
    pdata_root = os.path.join(hsqc_exp_dir, 'pdata')
    if not os.path.isdir(pdata_root):
        return []
    out = []
    for entry in sorted(os.listdir(pdata_root)):
        if entry == '1':  # 主谱目录跳过
            continue
        sub = os.path.join(pdata_root, entry)
        if not os.path.isdir(sub):
            continue
        used_from = os.path.join(sub, 'used_from')
        if os.path.isfile(used_from):
            out.append((entry, used_from))
    return out


def read_used_from_procno2(used_from_path: str) -> Optional[int]:
    """从 used_from 文件提取 PROCNO2 (HSQC 行切片的关键索引)."""
    if not os.path.isfile(used_from_path):
        return None
    try:
        with open(used_from_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##$PROCNO2='):
                    val = line[len('##$PROCNO2='):].strip()
                    return int(val)
    except (OSError, ValueError):
        return None
    return None


# ---- 产物 F: 2D 谱溶剂一致性校验 ----

def validate_2d_solvent_consistency(
    spectra_2d: Dict[str, Spectrum2D],
) -> dict:
    """检查所有 2D 实验的 SOLVENT 一致性.

    Args:
        spectra_2d: {experiment_type: Spectrum2D} 字典

    Returns:
        dict with keys:
          - n_total: int (检查的 2D 实验数)
          - unique_solvents: List[str] (唯一溶剂列表)
          - consistent: bool (所有 SOLVENT 是否一致)
          - details: [{exp_type, exp_no, solvent, nucleus}, ...]
          - warning: str (若 inconsistent, 给出告警文本)
    """
    out = {
        'n_total': 0,
        'unique_solvents': [],
        'consistent': True,
        'details': [],
        'warning': '',
    }
    if not spectra_2d:
        return out
    solvents = []
    for etype, spec in spectra_2d.items():
        sol = (spec.solvent or '').strip()
        out['n_total'] += 1
        out['details'].append({
            'exp_type': etype,
            'exp_no': spec.exp_no,
            'solvent': sol,
            'nucleus_pair': f'{spec.f1_nucleus}-{spec.f2_nucleus}',
        })
        if sol:
            solvents.append(sol)
    out['unique_solvents'] = sorted(set(solvents))
    out['consistent'] = len(out['unique_solvents']) <= 1
    if not out['consistent'] and len(out['unique_solvents']) > 1:
        out['warning'] = (
            f'2D 实验溶剂不一致: {out["unique_solvents"]}. '
            f'可能影响 HSQC/HMBC 化学位移标度对齐. '
            f'建议: 优先使用溶剂一致的 2D 实验, 或在 derive 流程中按溶剂分组.'
        )
    return out


def merge_13c_with_hmbc_projection(
    main_13c_shifts: List[float],
    hmbc_proj_peaks: List[ProjectionPeak],
    merge_tol: float = 0.1,
) -> dict:
    """合并主谱 13C 峰表 + HMBC 投影峰 (产物 B 落盘准备).

    Args:
        main_13c_shifts: 主谱 13C 化学位移 (来自 exp 22 peaklist.xml)
        hmbc_proj_peaks: HMBC 投影检测的 13C 峰 (来自 read_hmbc_f1_projection)
        merge_tol: 合并容差 (ppm, 默认 0.1)

    Returns:
        dict with keys:
          - merged: [(dC_ppm, source, intensity), ...] 按 dC 降序
            source = 'main' / 'hmbc_proj_only'
          - main_count: int
          - proj_only_count: int
          - intersection_count: int
    """
    main_sorted = sorted(set(round(c, 3) for c in main_13c_shifts), reverse=True)
    used_main = set()
    merged = []
    intersection = 0
    proj_only = 0
    for proj_p in hmbc_proj_peaks:
        d_proj = round(proj_p.shift, 3)
        # 找最接近的 main 峰
        best_main = None
        best_diff = float('inf')
        for m in main_sorted:
            diff = abs(d_proj - m)
            if diff < best_diff:
                best_diff = diff
                best_main = m
        if best_main is not None and best_diff <= merge_tol:
            used_main.add(best_main)
            merged.append((best_main, 'main+hmbc_proj', proj_p.intensity))
            intersection += 1
        else:
            merged.append((d_proj, 'hmbc_proj_only', proj_p.intensity))
            proj_only += 1
    # 补充 main 中未出现在 proj 的
    for m in main_sorted:
        if m not in used_main:
            merged.append((m, 'main_only', None))
    # 按 dC 降序
    merged.sort(key=lambda x: -x[0])
    return {
        'merged': merged,
        'main_count': len(main_sorted),
        'proj_only_count': proj_only,
        'intersection_count': intersection,
    }


def write_13c_merged_csv(
    merge_result: dict,
    output_path: str,
) -> None:
    """把合并结果落盘 (产物 B 的主要 CSV 产物)."""
    import csv as _csv

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or '.', exist_ok=True)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        w = _csv.writer(f)
        w.writerow(['# 13C 峰表 (主谱 ∪ HMBC F1 投影)'])
        w.writerow(['# 生成时间', str(__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S'))])
        w.writerow([])
        w.writerow(['主谱 (13C) 峰数', merge_result.get('main_count', 0)])
        w.writerow(['HMBC 投影独有峰数', merge_result.get('proj_only_count', 0)])
        w.writerow(['交集 (主谱 + 投影) 峰数', merge_result.get('intersection_count', 0)])
        w.writerow([])
        w.writerow(['dC_ppm', 'source', 'intensity'])
        for dC, src, inten in merge_result.get('merged', []):
            w.writerow([f'{dC:.3f}', src, '' if inten is None else f'{inten:.3f}'])


def write_hsqc_row_slices_csv(
    hsqc_exp_dir: str,
    output_dir: str,
) -> List[str]:
    """扫描 HSQC 实验下所有 pdata 子目录的 used_from, 自动读行切片并落盘 (产物 C).

    Returns:
        实际写入的 CSV 文件路径列表
    """
    written: List[str] = []
    slices = discover_hsqc_row_slices(hsqc_exp_dir)
    if not slices:
        return written

    os.makedirs(output_dir, exist_ok=True)
    import csv as _csv
    for pprocno, used_from_path in slices:
        row_idx = read_used_from_procno2(used_from_path)
        if row_idx is None:
            continue
        result = read_hsqc_row_slice(hsqc_exp_dir, pprocno, row_idx=row_idx)
        if result is None:
            continue
        out_csv = os.path.join(output_dir, 'hsqc_row_{}.csv'.format(pprocno))
        try:
            with open(out_csv, 'w', newline='', encoding='utf-8') as f:
                w = _csv.writer(f)
                w.writerow(['# HSQC 单行切片 (产物 C)'])
                w.writerow(['# HSQC exp', result.hsqc_exp])
                w.writerow(['# PROCNO2 (行号)', result.row_pprocno2])
                w.writerow(['# 推算 13C ppm', f'{result.f1_ppm:.3f}'])
                w.writerow(['# 检测峰数', len(result.detected_peaks)])
                w.writerow([])
                w.writerow(['H_ppm', 'intensity', 'snr'])
                for p in result.detected_peaks:
                    w.writerow([
                        f'{p.shift:.4f}',
                        f'{p.intensity:.3f}',
                        '' if p.snr is None else f'{p.snr:.2f}',
                    ])
            written.append(out_csv)
        except OSError:
            continue
    return written
