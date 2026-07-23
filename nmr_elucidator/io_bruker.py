"""Bruker 原始谱图读取.

只读取不改写。原始数据放在 NMR/<folder>/<expNo>/ 下：
  - pdata/1/peaklist.xml  峰表
  - pdata/1/parm.txt      采集/处理参数
  - pdata/1/integrals.txt 1H 积分区间（仅氢谱）
  - pdata/1/1r            二进制谱图（nmrglue 可读）
  - acqu / acqus          采集参数
  - title                 实验标题

输出 BrukerSpectrum：统一表示氢谱或碳谱，含元数据 + 峰表 + 积分。
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
from xml.etree import ElementTree as ET


# ---- 数据结构 ----

@dataclass
class Peak:
    """单条峰记录."""

    shift: float  # ppm
    intensity: float = 0.0  # 原始强度（peaklist.xml 的 intensity）
    region_low: Optional[float] = None  # peaklist region F1 lower
    region_high: Optional[float] = None  # peaklist region F1 upper
    # DEPT-135 方向. peaklist 默认无方向信息, 由相对强度推断:
    #   'positive'  = 正峰 (CH, CH3)
    #   'negative'  = 负峰 (CH2 in DEPT-135)
    #   'unknown'   = 1H 谱, DEPT-90, 或无法判断
    sign: str = 'unknown'


@dataclass
class IntegralRegion:
    """1H 积分区间."""

    low: float  # ppm
    high: float  # ppm
    value: float  # 积分值


@dataclass
class BrukerSpectrum:
    """单个 Bruker 实验（1D 谱）."""

    exp_no: str  # e.g. "21"
    name: str  # e.g. "clxh3"
    solvent: str  # e.g. "Pyr", "CDCl3"
    nucleus: str  # "1H" or "13C"
    field_mhz: float  # 频率 (MHz)
    temperature_k: Optional[float]
    pulse_program: str  # e.g. "zg", "zgdc", "dept"
    peaks: List[Peak] = field(default_factory=list)
    integrals: List[IntegralRegion] = field(default_factory=list)
    source_path: str = ""
    # 一些诊断信息
    parm_path: Optional[str] = None
    peaklist_path: Optional[str] = None

    def shift_values(self) -> List[float]:
        return [p.shift for p in self.peaks]


# ---- 错误类型 ----

class BrukerReadError(Exception):
    """Bruker 读取过程中任何显式错误."""

    def __init__(self, msg: str, path: str = "", stage: str = ""):
        super().__init__(f"[{stage}] {msg} (path={path})")
        self.path = path
        self.stage = stage


# ---- parm.txt 解析 ----

def _parse_parm(path: str) -> dict:
    """解析 Bruker parm.txt 或 acqus 文件.

    支持两种格式:
      A) 空格分隔（TopSpin parm.txt）:
           NAME              clxh3
           SOLVENT             Pyr
      B) 等号分隔（TopSpin ##header = value 风格）:
           ##SOLVENT= Pyr
           ##$SFO1= 125.78 MHz

    返回 {key: value}; 多个相同 key 取最后一个（与 TopSpin 行为一致）.
    """
    out: dict = {}
    if not os.path.exists(path):
        return out
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            text = f.read()
    except OSError as e:
        raise BrukerReadError(f"无法读取文件: {e}", path, stage="parm-open")
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith('$$'):
            continue
        # 跳过标题行
        if s.startswith('=' * 3) or s.startswith('-' * 3):
            continue
        # 跳过 "Current Data Parameters" 之类无值的标题
        if ' - ' in s and '=' not in s:
            continue
        # 跳过 "Current ..." 之类（首词是 Current 且无 '='）
        if s.startswith('Current') and '=' not in s:
            continue
        key = ''
        val = ''
        if '=' in s:
            # 等号格式
            head, _, tail = s.partition('=')
            # lstrip 接收字符集合, 不是字符串前缀; 'lstrip("#$")' 一次性去掉所有 #/$
            # 原 'lstrip("$").lstrip("#")' 顺序错: 去 $ 时不动 (# 不在字符集), 去 # 时残留 $
            # 导致 '##$SFO1' 变成 '$SFO1', 首字符非字母被 isalpha() 过滤掉 → 键丢失
            key = head.strip().lstrip('#$').strip().upper()
            val = tail.strip()
        else:
            # 空格格式：第一段为 key, 其余拼接为 val
            parts = s.split(None, 1)
            if len(parts) < 2:
                continue
            key = parts[0].strip().upper()
            val = parts[1].strip()
        # 键的合法性: 首字符必须为字母, 其余字母/数字/_/.<>
        if not key or not key[0].isalpha():
            continue
        if key:
            out[key] = val  # 后写覆盖前写 = 取最后一个
    return out


def _nucleus_from_parm(p: dict) -> str:
    """从 parm 推断核种. 优先 NUC1, 缺则取 NUCLEUS."""
    n = p.get('NUC1') or p.get('NUCLEUS') or p.get('NUC2') or ''
    n = n.strip().upper()
    # 13C, 1H, 15N, 31P 等
    if n.startswith('1H'):
        return '1H'
    if n.startswith('13C') or n.startswith('13 C'):
        return '13C'
    if n:
        return n
    return ''


def _strip_unit(s: str) -> str:
    """去掉末尾单位 Hz/MHz/sec/usec/K/W. 不影响数字部分."""
    import re
    return re.sub(r'\s*(MHz|Hz|usec|sec|K|W)\s*$', '', s).strip()


def _field_mhz_from_parm(p: dict) -> float:
    """从 parm 推断共振频率 MHz. 优先 SFO2 (1H 频率), 缺则 SFO1.

    SFO 值形如 '125.7848246 MHz' 或 '500.1827508 MHz'.
    """
    # 优先 1H 频率 (SFO2), 缺则用 SFO1 (主核)
    for k in ('SFO2', 'SFO1', 'BF1', 'SF'):
        v = p.get(k, '')
        if not v:
            continue
        try:
            f = float(_strip_unit(v))
            # SFO 单位通常是 MHz；若 < 50 则视为 Hz → 转 MHz
            if f < 50:
                f = f * 1e-6 if f < 50 else f / 1e6  # 防御
                # 简化: < 50 视为 Hz
                f = f * 1e-6
            return f
        except (TypeError, ValueError):
            continue
    return 0.0


def _temperature_from_parm(p: dict) -> Optional[float]:
    v = p.get('TE') or p.get('TEMPERATURE')
    if not v:
        return None
    try:
        # 形如 "298.2 K" 或 "298.2"
        return float(_strip_unit(v))
    except (TypeError, ValueError):
        return None


# 1H 脉冲程序白名单（不穷举，作为兜底）
_1H_PULSES = {'zg', 'zg30', 'zgpr', 'zgesgp', 'zgesgppe', 'noesygppr1d',
              'watergate', 'excitation', 'cpmg', 'led', 'cpmgpr1d'}

# 13C 脉冲程序白名单
_13C_PULSES = {'zgdc', 'zgpg', 'zgpg30', 'cpmg', 'zgig', 'zgig30',
               'dept', 'dept135', 'dept45', 'deptsp90', 'deptsp135',
               'deptsp45', 'apt', 'jmod', 'zgdcpr', 'invig', 'zgdcwfg',
               'mlev'}


def _detect_pulse_role(p: dict) -> str:
    """区分碳谱/氢谱/DEPT. 失败返回空串, 调用方应回退到 NUC1 推断."""
    pp = (p.get('PULPROG') or '').lower().strip()
    if not pp:
        return ''
    if pp in _13C_PULSES:
        if pp.startswith('dept'):
            return '13C-DEPT'
        return '13C'
    if pp in _1H_PULSES:
        return '1H'
    # 部分匹配
    if pp.startswith('dept'):
        return '13C-DEPT'
    if pp.startswith(('zgdc', 'zgpg', 'zgig')):
        return '13C'
    if pp.startswith(('zg', 'noesy', 'cosy')):
        return '1H'
    return ''


# ---- peaklist.xml 解析 ----

def _parse_peaklist_xml(path: str) -> List[Peak]:
    """解析 Bruker peaklist.xml.

    元素结构：
      <PeakList1DHeader>
        <PeakPickDetails>F1=8.946ppm, F2=6.581ppm, ...</PeakPickDetails>
      </PeakList1DHeader>
      <Peak1D F1="8.7093" intensity="4.61" type="0"/>
    """
    if not os.path.exists(path):
        return []
    try:
        tree = ET.parse(path)
    except ET.ParseError as e:
        raise BrukerReadError(f"peaklist.xml 解析失败: {e}", path, stage="peaklist-xml")
    root = tree.getroot()
    out: List[Peak] = []
    # PeakPickDetails 提供 region F1/F2
    for plist in root.iter('PeakList1D'):
        header = plist.find('PeakList1DHeader')
        region_low = region_high = None
        if header is not None:
            det = header.findtext('PeakPickDetails', default='')
            # 形如: F1=8.946ppm, F2=6.581ppm, MI=0.54cm, MAXI=20.59cm, PC=1.000
            for part in det.split(','):
                part = part.strip()
                if part.startswith('F1='):
                    region_low = _to_ppm(part[3:])
                elif part.startswith('F2='):
                    region_high = _to_ppm(part[3:])
        for peak in plist.iter('Peak1D'):
            f1_text = peak.attrib.get('F1', '')
            if not f1_text:
                continue
            try:
                shift = float(f1_text)
            except ValueError:
                continue
            try:
                intensity = float(peak.attrib.get('intensity', '0'))
            except ValueError:
                intensity = 0.0
            out.append(Peak(
                shift=shift, intensity=intensity,
                region_low=region_low, region_high=region_high,
                sign=_infer_peak_sign(intensity, nucleus_hint='', pulse_hint=''),
            ))
    return out


def _infer_peak_sign(intensity: float, nucleus_hint: str = '', pulse_hint: str = '') -> str:
    """根据 peak intensity 推断 DEPT-135 峰方向.

    Bruker DEPT-135 peaklist 保留正负号:
      intensity > 0 -> positive (CH, CH3)
      intensity < 0 -> negative (CH2)
      intensity == 0 -> unknown
    1H 谱与 DEPT-90 的 intensity 总是 >= 0, sign 标 'unknown'.
    """
    if intensity > 1e-6:
        return 'positive'
    if intensity < -1e-6:
        return 'negative'
    return 'unknown'


def _to_ppm(s: str) -> Optional[float]:
    """'8.946ppm' -> 8.946. 也兼容无单位."""
    s = s.strip()
    for unit in ('ppm', 'PPM'):
        if s.endswith(unit):
            s = s[:-len(unit)].strip()
            break
    try:
        return float(s)
    except ValueError:
        return None


# ---- integrals.txt 解析 ----

def _parse_integrals(path: str) -> List[IntegralRegion]:
    if not os.path.exists(path):
        return []
    out: List[IntegralRegion] = []
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
    except OSError as e:
        raise BrukerReadError(f"无法读取积分文件: {e}", path, stage="integrals-open")
    for line in lines:
        s = line.strip()
        # 形如: '    1      7.165      7.025         3.80669'
        parts = s.split()
        if len(parts) < 4:
            continue
        # 第一列是 Number; 第二/三列是 ppm 范围; 第四列是积分值
        try:
            low = float(parts[1])
            high = float(parts[2])
            value = float(parts[3])
        except ValueError:
            continue
        out.append(IntegralRegion(low=low, high=high, value=value))
    return out


# ---- 顶层接口 ----

def read_bruker_experiment(exp_dir: str) -> BrukerSpectrum:
    """读取单个 Bruker 实验目录.

    exp_dir 形如: '.../NMR/clxh3/21'
    """
    exp_dir = os.path.abspath(exp_dir)
    if not os.path.isdir(exp_dir):
        raise BrukerReadError(f"实验目录不存在: {exp_dir}", exp_dir, stage="exp-dir")

    exp_no = os.path.basename(exp_dir)
    name = os.path.basename(os.path.dirname(exp_dir))

    # 找 parm.txt 或 acqus
    parm_candidates = [
        os.path.join(exp_dir, 'pdata', '1', 'parm.txt'),
        os.path.join(exp_dir, 'acqus'),
    ]
    parm = {}
    parm_path = ''
    for c in parm_candidates:
        if os.path.exists(c):
            parm_path = c
            parm = _parse_parm(c)
            if parm:
                break

    # 核种推断优先级：NUC1 > PULSE_PROG 推断
    nucleus = _nucleus_from_parm(parm)
    pp = parm.get('PULPROG', '').strip()
    role = _detect_pulse_role(parm)
    if nucleus == '13C' and role.startswith('13C-DEPT'):
        nucleus = '13C-DEPT'
    elif role == '13C-DEPT':
        nucleus = '13C-DEPT'
    elif role == '13C' and not nucleus:
        nucleus = '13C'
    elif role == '1H' and not nucleus:
        nucleus = '1H'

    # 溶剂：可能没有 SOLVENT 字段时从 title 推断
    solvent = parm.get('SOLVENT', '').strip()
    title_path = os.path.join(exp_dir, 'pdata', '1', 'title')
    if not solvent and os.path.exists(title_path):
        try:
            with open(title_path, 'r', encoding='utf-8', errors='replace') as f:
                t = f.read().strip().upper()
            for cand in ('CDCL3', 'PYR', 'CD3OD', 'D2O', 'DMSO-D6', 'ACETONE-D6', 'C6D6'):
                if cand in t.replace(' ', '').replace('-', ''):
                    solvent = cand
                    break
        except OSError:
            pass

    # 峰表
    peaklist_path = os.path.join(exp_dir, 'pdata', '1', 'peaklist.xml')
    peaks = _parse_peaklist_xml(peaklist_path)

    # 积分（仅氢谱有意义）
    integrals_path = os.path.join(exp_dir, 'pdata', '1', 'integrals.txt')
    integrals = _parse_integrals(integrals_path) if nucleus == '1H' else []

    if not peaks:
        raise BrukerReadError(
            f"未在 peaklist.xml 中找到任何峰（{peaklist_path}）",
            peaklist_path, stage="peaklist-empty",
        )

    return BrukerSpectrum(
        exp_no=exp_no,
        name=name,
        solvent=solvent,
        nucleus=nucleus,
        field_mhz=_field_mhz_from_parm(parm),
        temperature_k=_temperature_from_parm(parm),
        pulse_program=pp,
        peaks=peaks,
        integrals=integrals,
        source_path=exp_dir,
        parm_path=parm_path,
        peaklist_path=peaklist_path,
    )


def read_bruker_folder(folder: str) -> List[BrukerSpectrum]:
    """读取一个样本目录（如 .../clxh3）下的所有 1D 实验.

    子目录名（数字）会被当作 expNo. 任何子目录里只要存在 pdata/1/peaklist.xml 就读取.
    """
    folder = os.path.abspath(folder)
    if not os.path.isdir(folder):
        raise BrukerReadError(f"样本目录不存在: {folder}", folder, stage="folder")
    out: List[BrukerSpectrum] = []
    for entry in sorted(os.listdir(folder)):
        exp_dir = os.path.join(folder, entry)
        if not os.path.isdir(exp_dir):
            continue
        # 必须有 peaklist.xml 才是有效 1D 实验
        peaklist = os.path.join(exp_dir, 'pdata', '1', 'peaklist.xml')
        if not os.path.exists(peaklist):
            continue
        out.append(read_bruker_experiment(exp_dir))
    return out


# ---- scipy.signal.find_peaks 二次找峰 (可选增强) ----

def enhance_peaks_with_scipy(
    pdata_dir: str,
    existing_peaks: List[Peak],
    tol_ppm: float = 0.1,
    min_distance_pts: int = 5,
    prominence: float = 0.05,
    height_min: float = 0.0,
) -> List[Peak]:
    """从 1r 实谱用 scipy.signal.find_peaks 二次找峰, 与 peaklist.xml 合并去重.

    设计目标: peaklist.xml 是 Bruker TopSpin 自动峰检测, 可能漏检弱峰/重叠峰;
    本函数从 1r 实谱 (连续谱线) 补充, 与现有峰位相差 <tol_ppm 视为同峰, 丢弃.

    依赖 (可选, 未装时静默返回原峰列表):
      - scipy.signal.find_peaks
      - nmrglue (io_processed.read_processed_spectrum)

    失败路径 (scipy 未装 / nmrglue 未装 / 1r 缺失) 静默返回原峰列表, 不抛错,
    保证 read_bruker_experiment 主路径不受影响.

    ponytail: 跳过项 — 升级路径: 若需支持重叠峰分峰, 改用 scipy.signal.peak_widths + 拟合 (curve_fit).
    """
    try:
        from scipy.signal import find_peaks
    except ImportError:
        return list(existing_peaks)
    try:
        from .io_processed import ProcessedReadError, read_processed_spectrum
    except ImportError:
        return list(existing_peaks)
    try:
        ps = read_processed_spectrum(pdata_dir, nucleus='?')
    except ProcessedReadError:
        return list(existing_peaks)
    except Exception:
        return list(existing_peaks)
    intensity = ps.intensity
    ppm = ps.ppm
    if intensity is None or len(intensity) == 0:
        return list(existing_peaks)
    peak_max = float(intensity.max())
    if peak_max <= 0:
        return list(existing_peaks)
    norm = intensity / peak_max
    # 默认参数: 最小峰间隔 + 突出度阈值. 可由调用方覆盖.
    kwargs: dict = {'distance': max(1, int(min_distance_pts))}
    if prominence is not None:
        kwargs['prominence'] = float(prominence)
    if height_min is not None:
        kwargs['height'] = float(height_min)
    try:
        indices, _ = find_peaks(norm, **kwargs)
    except Exception:
        return list(existing_peaks)
    if len(indices) == 0:
        return list(existing_peaks)
    # 去重合并: 与 existing_peaks 比, |ppm delta| < tol_ppm 视为同峰, 丢弃新峰.
    existing_shifts = [p.shift for p in existing_peaks]
    merged = list(existing_peaks)
    for idx in indices:
        try:
            idx_i = int(idx)
        except (TypeError, ValueError):
            continue
        if idx_i < 0 or idx_i >= len(ppm):
            continue
        s = float(ppm[idx_i])
        if any(abs(s - es) < tol_ppm for es in existing_shifts):
            continue  # 与原 peaklist 重复
        merged.append(Peak(
            shift=s,
            intensity=float(norm[idx_i]),
            sign='unknown',
        ))
    # 按化学位移降序 (与原 _parse_peaklist_xml 一致)
    merged.sort(key=lambda p: p.shift, reverse=True)
    return merged


# ---- 1H 多重峰分析 (J 耦合) ----

@dataclass
class Multiplet:
    """化学家分析后的一个 1H 多重峰.

    来源: Bruker peaklist.xml 中的 <Multiplet> 标签 (TopSpin 手动多重峰分析),
    或从峰间距估算 J 值.
    """

    center_ppm: float       # 化学位移中心 (ppm)
    integral: float = 0.0   # 积分 (H 个数)
    j_couplings: List[float] = field(default_factory=list)  # J 耦合常数 (Hz)
    multiplicity: str = ''  # 's'/'d'/'t'/'q'/'m' 等
    pattern_type: str = ''  # 'singlet'/'doublet'/'triplet'/'quartet'/'multiplet'


def parse_multiplets_from_peaklist(
    peaklist_path: str,
    field_mhz: float = 0.0,
) -> List[Multiplet]:
    """从 Bruker peaklist.xml 解析多重峰信息 (J 耦合).

    Bruker TopSpin 的多重峰分析结果可能以三种形式存储:
    1. <Multiplet> 标签 (TopSpin 3.x/4.x 手动分析)
    2. <Peak1D> 带 multiplicity/J 属性 (部分版本)
    3. (自动推断) 从 <Peak1D> 峰间距反推 J 耦合, 作为 fallback

    若文件中无任何多重峰信息, 静默返回空列表 (不报错, 不阻断主流程).

    Args:
        peaklist_path: peaklist.xml 文件路径
        field_mhz: 谱仪 1H 共振频率 (MHz); 仅在策略 1/2 失败后,
                   用于自动推断 J 耦合. 0 = 未知 (跳过推断)

    Returns:
        Multiplet 列表 (可能为空)
    """
    import re

    if not os.path.isfile(peaklist_path):
        return []

    try:
        with open(peaklist_path, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
    except OSError:
        return []

    multiplets: List[Multiplet] = []

    # 策略 1: 查找 <Multiplet> 标签 (TopSpin 3.x/4.x)
    # 格式: <Multiplet F1="7.26" F2="7.24" intensity="1.0" J="7.2,7.2" multiplicity="t">
    mult_pattern = re.compile(
        r'<Multiplet\b([^>]*)/?>',
        re.IGNORECASE,
    )
    for m in mult_pattern.finditer(content):
        attrs_str = m.group(1)
        try:
            # 解析属性
            f1_m = re.search(r'F1="([^"]+)"', attrs_str)
            f2_m = re.search(r'F2="([^"]+)"', attrs_str)
            int_m = re.search(r'intensity="([^"]+)"', attrs_str)
            j_m = re.search(r'J="([^"]+)"', attrs_str)
            mult_m = re.search(r'multiplicity="([^"]+)"', attrs_str)

            if f1_m:
                f1 = float(f1_m.group(1))
                f2 = float(f2_m.group(1)) if f2_m else f1
                center = (f1 + f2) / 2.0
            else:
                continue

            integral = float(int_m.group(1)) if int_m else 0.0

            # J 值: 逗号或分号分隔
            j_vals: List[float] = []
            if j_m:
                for jv in re.split(r'[,;\s]+', j_m.group(1)):
                    try:
                        j_vals.append(float(jv))
                    except ValueError:
                        pass

            multiplicity = mult_m.group(1) if mult_m else ''

            multiplets.append(Multiplet(
                center_ppm=round(center, 4),
                integral=integral,
                j_couplings=j_vals,
                multiplicity=multiplicity,
                pattern_type=_multiplicity_to_pattern(multiplicity),
            ))
        except (ValueError, AttributeError):
            continue  # 容错: 单个 multiplet 解析失败不影响其他

    # 策略 2: 查找带 multiplicity 属性的 <Peak1D> (部分 TopSpin 版本)
    if not multiplets:
        peak_mult_pattern = re.compile(
            r'<Peak1D\b([^>]*multiplicity[^>]*)/?>',
            re.IGNORECASE,
        )
        for m in peak_mult_pattern.finditer(content):
            attrs_str = m.group(1)
            try:
                f1_m = re.search(r'F1="([^"]+)"', attrs_str)
                int_m = re.search(r'intensity="([^"]+)"', attrs_str)
                j_m = re.search(r'J="([^"]+)"', attrs_str)
                mult_m = re.search(r'multiplicity="([^"]+)"', attrs_str)

                if not f1_m:
                    continue

                center = float(f1_m.group(1))
                integral = float(int_m.group(1)) if int_m else 0.0

                j_vals = []
                if j_m:
                    for jv in re.split(r'[,;\s]+', j_m.group(1)):
                        try:
                            j_vals.append(float(jv))
                        except ValueError:
                            pass

                multiplicity = mult_m.group(1) if mult_m else ''

                multiplets.append(Multiplet(
                    center_ppm=round(center, 4),
                    integral=integral,
                    j_couplings=j_vals,
                    multiplicity=multiplicity,
                    pattern_type=_multiplicity_to_pattern(multiplicity),
                ))
            except (ValueError, AttributeError):
                continue

    # 策略 3 (fallback): 从 <Peak1D> 峰间距自动推断
    if not multiplets and field_mhz > 0:
        multiplets = infer_multiplets_from_peaks(peaklist_path, field_mhz=field_mhz)

    return multiplets


def infer_multiplets_from_peaks(
    peaklist_path: str,
    field_mhz: float = 0.0,
    cluster_gap_hz: float = 50.0,
    min_peak_intensity_ratio: float = 0.05,
) -> List[Multiplet]:
    """从 peaklist.xml 的 <Peak1D> 列表自动检测多重峰并反推 J 耦合.

    适用场景: 原始 Bruker 数据未做手工多重峰分析(peaklist.xml 无 <Multiplet>
    标签, 也无 multiplicity 属性)时, 作为 fallback 从峰间距估算 J 值.

    算法:
      1. 解析所有 Peak1D (F1=ppm, intensity)
      2. 按 ppm 降序排序 (NMR 惯例)
      3. 聚类: 相邻峰间距 < cluster_gap_hz Hz (默认 50 Hz, 覆盖 99% 的
         多重峰跨度) 归为同一簇
      4. 过滤噪声: 簇总强度 < 所有峰中位数 × min_peak_intensity_ratio
         的簇视为噪声, 丢弃
      5. 簇内峰数推断 multiplicity (1=s, 2=d, 3=t, 4=q, 5+=m)
      6. J = 簇内相邻峰平均间距 × field_mhz (Hz)
         若 field_mhz <= 0, 无法计算 J, 返回空 J 列表
      7. center_ppm = 簇内加权平均 (权重 = intensity)
      8. integral = 簇总强度

    Args:
        peaklist_path: peaklist.xml 路径
        field_mhz: 谱仪 1H 共振频率 (MHz); 0 表示未知, J 无法计算
        cluster_gap_hz: 同一多重峰内允许的最大相邻间距 (Hz)
        min_peak_intensity_ratio: 簇总强度占全部峰中位数的最低比例 (噪声门)

    Returns:
        Multiplet 列表 (峰数不足时可能为空, 不报错)
    """
    import re

    if not os.path.isfile(peaklist_path):
        return []

    try:
        with open(peaklist_path, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
    except OSError:
        return []

    # 解析所有 Peak1D
    peaks: List[tuple] = []  # (ppm, intensity)
    for m in re.finditer(r'<Peak1D\b([^/>]*)/?>', content, re.IGNORECASE):
        f1_m = re.search(r'F1="([^"]+)"', m.group(1))
        int_m = re.search(r'intensity="([^"]+)"', m.group(1))
        if not f1_m:
            continue
        try:
            ppm = float(f1_m.group(1))
            intensity = float(int_m.group(1)) if int_m else 0.0
        except (ValueError, TypeError):
            continue
        peaks.append((ppm, intensity))

    if len(peaks) < 2:
        return []

    # 按 ppm 降序排序 (NMR 习惯)
    peaks.sort(key=lambda x: -x[0])

    # 噪声门: 总峰强度的中位数
    intensities_sorted = sorted(p[1] for p in peaks)
    median_intensity = intensities_sorted[len(intensities_sorted) // 2]
    min_total_intensity = median_intensity * min_peak_intensity_ratio

    # Hz/ppm 转换系数
    hz_per_ppm = field_mhz if field_mhz > 0 else 0.0

    # 聚类阈值: cluster_gap_hz 转 ppm (无频率时用 0.3 ppm 作为寤底)
    cluster_gap_ppm = (cluster_gap_hz / hz_per_ppm) if hz_per_ppm > 0 else 0.3

    # 贪心聚类: 相邻峰 ppm 距离 < cluster_gap_ppm 归同一簇
    clusters: List[List[tuple]] = []
    current: List[tuple] = [peaks[0]]
    for p in peaks[1:]:
        prev_ppm = current[-1][0]
        if abs(p[0] - prev_ppm) <= cluster_gap_ppm:
            current.append(p)
        else:
            clusters.append(current)
            current = [p]
    clusters.append(current)

    multiplets: List[Multiplet] = []
    for cluster in clusters:
        n_peaks = len(cluster)
        if n_peaks < 1:
            continue

        total_intensity = sum(p[1] for p in cluster)

        # 噪声过滤: 单峰且低于门限才丢 (不丢多重峰)
        if n_peaks == 1 and total_intensity < min_total_intensity:
            continue

        # center = 加权平均
        if total_intensity > 0:
            center = sum(p[0] * p[1] for p in cluster) / total_intensity
        else:
            center = sum(p[0] for p in cluster) / n_peaks

        # 簇内相邻间距 (ppm) → J 值 (Hz)
        diffs_ppm = [abs(cluster[i][0] - cluster[i + 1][0])
                     for i in range(n_peaks - 1)]
        j_couplings: List[float] = []
        if hz_per_ppm > 0 and diffs_ppm:
            j_couplings = [round(d * hz_per_ppm, 2) for d in diffs_ppm]

        # 多重性推断
        if n_peaks == 1:
            multiplicity, pattern_type = 's', 'singlet'
        elif n_peaks == 2:
            multiplicity, pattern_type = 'd', 'doublet'
        elif n_peaks == 3:
            multiplicity, pattern_type = 't', 'triplet'
        elif n_peaks == 4:
            multiplicity, pattern_type = 'q', 'quartet'
        else:
            multiplicity, pattern_type = 'm', 'multiplet'

        multiplets.append(Multiplet(
            center_ppm=round(center, 4),
            integral=round(total_intensity, 3),
            j_couplings=j_couplings,
            multiplicity=multiplicity,
            pattern_type=pattern_type,
        ))

    # 按 ppm 降序输出
    multiplets.sort(key=lambda m: -m.center_ppm)
    return multiplets


def analyze_multiplet_at_ppm(
    peaklist_path: str,
    center_ppm: float,
    field_mhz: float = 0.0,
    window_ppm: float = 0.15,
    cluster_gap_hz: float = 30.0,
) -> dict:
    """在以 center_ppm 为中心的窗口内做精细化多重峰分析 (对标 TopSpin Manual).

    用途: 对某个 C-H 对 (CHPair) 的 1H 位移, 取附近 peaklist 峰, 拟合该
    H 处的多重峰结构. 输出结果可与 TopSpin "Multiple Peak Analysis"
    报告对应.

    算法:
      1. 取 center_ppm ± window_ppm 范围内的 Peak1D
      2. 在窗口内贪心聚类 (相邻峰间距 < cluster_gap_hz Hz) - 取含峰数最多的簇
      3. 簇内峰数 → multiplicity (s/d/t/q/m)
      4. J = 簇内相邻峰平均间距 × field_mhz
      5. halfwidth_Hz ≈ J_median / 3 (NMR 经验: 子峰半高宽约为 J 间距的 1/3)
      6. intensity_uniformity: 子峰强度剖面与理论 1:1/1:2:1/1:3:3:1 比值的吻合度
      7. quality_grade: Excellent (uniform>0.9), Good (>0.7), Poor (>0.5), Very Poor (<0.5)

    Args:
        peaklist_path: peaklist.xml 路径
        center_ppm: CHPair 的 1H 化学位移 (ppm)
        field_mhz: 谱仪 1H 频率 (MHz); 0 表示未知, J 无法计算
        window_ppm: 取峰窗口半宽 (ppm, 默认 ±0.15)
        cluster_gap_hz: 簇内峰间允许的最大间距 (Hz, 默认 30)

    Returns:
        dict 包含:
          - center_ppm: 实际加权中心
          - pattern: 's'/'d'/'t'/'q'/'m'/''
          - num_subpeaks: int
          - j_couplings_Hz: List[float]
          - halfwidth_Hz: float (估算)
          - intensity_uniformity: 0-1 (子峰强度均匀性)
          - quality_grade: 'Excellent'/'Good'/'Poor'/'Very Poor'/''
          - subpeak_ppms: List[float]
          - subpeak_intensities: List[float]
          - note: str (说明)
    """
    import re

    result = {
        'center_ppm': round(center_ppm, 4),
        'pattern': '',
        'num_subpeaks': 0,
        'j_couplings_Hz': [],
        'halfwidth_Hz': 0.0,
        'intensity_uniformity': 0.0,
        'quality_grade': '',
        'subpeak_ppms': [],
        'subpeak_intensities': [],
        'note': '',
    }

    if not os.path.isfile(peaklist_path):
        result['note'] = 'peaklist 不存在'
        return result

    try:
        with open(peaklist_path, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
    except OSError:
        result['note'] = 'peaklist 读取失败'
        return result

    # 解析窗口内 Peak1D
    peaks: List[tuple] = []  # (ppm, intensity)
    for m in re.finditer(r'<Peak1D\b([^/>]*)/?>', content, re.IGNORECASE):
        f1_m = re.search(r'F1="([^"]+)"', m.group(1))
        int_m = re.search(r'intensity="([^"]+)"', m.group(1))
        if not f1_m:
            continue
        try:
            ppm = float(f1_m.group(1))
            intensity = float(int_m.group(1)) if int_m else 0.0
        except (ValueError, TypeError):
            continue
        if abs(ppm - center_ppm) <= window_ppm:
            peaks.append((ppm, intensity))

    if not peaks:
        result['note'] = f'窗口 ±{window_ppm:.2f} ppm 内无 Peak1D'
        return result

    # 按 ppm 降序
    peaks.sort(key=lambda x: -x[0])

    hz_per_ppm = field_mhz if field_mhz > 0 else 0.0
    cluster_gap_ppm = (cluster_gap_hz / hz_per_ppm) if hz_per_ppm > 0 else 0.1

    # 贪心聚类
    clusters: List[List[tuple]] = []
    current: List[tuple] = [peaks[0]]
    for p in peaks[1:]:
        if abs(p[0] - current[-1][0]) <= cluster_gap_ppm:
            current.append(p)
        else:
            clusters.append(current)
            current = [p]
    clusters.append(current)

    # 取含峰数最多的簇 (主多重峰)
    cluster = max(clusters, key=lambda c: len(c))
    n_peaks = len(cluster)

    # center = 加权平均
    total_intensity = sum(p[1] for p in cluster)
    if total_intensity > 0:
        center = sum(p[0] * p[1] for p in cluster) / total_intensity
    else:
        center = sum(p[0] for p in cluster) / n_peaks

    # 子峰 ppm + intensity (取实际位置)
    subpeak_ppms = [round(p[0], 4) for p in cluster]
    subpeak_intensities = [round(p[1], 3) for p in cluster]

    # J 值
    diffs_ppm = [abs(cluster[i][0] - cluster[i + 1][0])
                 for i in range(n_peaks - 1)]
    j_couplings: List[float] = []
    if hz_per_ppm > 0 and diffs_ppm:
        j_couplings = [round(d * hz_per_ppm, 2) for d in diffs_ppm]

    # 多重性推断
    if n_peaks == 1:
        pattern = 's'
    elif n_peaks == 2:
        pattern = 'd'
    elif n_peaks == 3:
        pattern = 't'
    elif n_peaks == 4:
        pattern = 'q'
    else:
        pattern = 'm'

    # 强度均匀性: 与理论 1:1 / 1:2:1 / 1:3:3:1 比值的吻合度
    intensity_uniformity = 0.0
    if n_peaks >= 2:
        # 理论比值 (二项式系数, n-1 次偶合)
        # n=2 (d): 1:1; n=3 (t): 1:2:1; n=4 (q): 1:3:3:1; n>=5 (m): 不估计
        theoretical = _binomial_ratios(n_peaks)
        # 实际归一化强度剖面 (取正顺序, 使强度从低到高再匹配)
        actual = [p[1] for p in cluster]
        max_a = max(actual) if max(actual) > 0 else 1.0
        actual_norm = [a / max_a for a in actual]
        # 计算与理论比值的吻合度 (1 - 平均相对偏差)
        if theoretical and len(actual_norm) == len(theoretical):
            deviations = []
            for a, t in zip(actual_norm, theoretical):
                if t > 0:
                    deviations.append(abs(a - t) / t)
            if deviations:
                mean_dev = sum(deviations) / len(deviations)
                intensity_uniformity = round(max(0.0, 1.0 - mean_dev), 3)

    # 半高宽估算 (Hz)
    halfwidth_hz = 0.0
    if j_couplings:
        # NMR 经验: 子峰半高宽约为 J 的 1/3 (避免重叠但保持可辨)
        halfwidth_hz = round(sum(j_couplings) / len(j_couplings) / 3.0, 2)

    # 质量等级
    if intensity_uniformity >= 0.9:
        quality = 'Excellent'
    elif intensity_uniformity >= 0.7:
        quality = 'Good'
    elif intensity_uniformity >= 0.5:
        quality = 'Poor'
    elif intensity_uniformity > 0:
        quality = 'Very Poor'
    else:
        quality = 'N/A'  # 仅 s 或 m (无法评估)

    result.update({
        'center_ppm': round(center, 4),
        'pattern': pattern,
        'num_subpeaks': n_peaks,
        'j_couplings_Hz': j_couplings,
        'halfwidth_Hz': halfwidth_hz,
        'intensity_uniformity': intensity_uniformity,
        'quality_grade': quality,
        'subpeak_ppms': subpeak_ppms,
        'subpeak_intensities': subpeak_intensities,
        'note': f'{n_peaks}个子峰, J={j_couplings}Hz',
    })
    return result


def _binomial_ratios(n: int) -> List[float]:
    """理论偶合多重峰的强度比 (二项式系数)."""
    if n <= 0:
        return []
    # n 个子峰 → n-1 个偶合 → 系数 = C(n-1, k), k=0..n-1
    from math import comb
    coefs = [comb(n - 1, k) for k in range(n)]
    max_c = max(coefs) if coefs else 1
    return [c / max_c for c in coefs]


def _multiplicity_to_pattern(mult: str) -> str:
    """将简写多重性转为完整名称."""
    mapping = {
        's': 'singlet',
        'd': 'doublet',
        't': 'triplet',
        'q': 'quartet',
        'm': 'multiplet',
        'dd': 'doublet of doublets',
        'dt': 'doublet of triplets',
        'td': 'triplet of doublets',
        'tt': 'triplet of triplets',
        'br': 'broad',
    }
    return mapping.get(mult.lower().strip(), mult)