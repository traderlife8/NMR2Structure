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