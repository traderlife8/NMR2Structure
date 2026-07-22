"""nmrshiftdb2 SDF 数据读取.

设计要点：
- 用 bytes 直接读，避免 RDKit SDMolSupplier 的中文路径 bug（C++ 层 file open 问题）
- 自解析 MOL block + SDF tag，再用 RDKit MolFromMolBlock 解析单个分子
- 提取关键 tag：
    <nmrshiftdb2 ID>          数据库 ID
    <INChI> / <INChI key>     InChI 字符串
    <Spectrum 13C N>           13C 化学位移表（"shift;int;mult;atom|..."）
    <Spectrum 1H N>            1H 化学位移表
    <Solvent>                  溶剂（形如 "0:CDCl3 "）
    <Field Strength [MHz]>     频率
    <Temperature [K]>          温度

性能说明：
- 单个文件 150 MB（58187 记录）
- 默认提供流式迭代（stream_nmrshiftdb2_sdf），不一次性读全文件
- 提供 read_nmrshiftdb2_sdf 一次性读取，返回 list（用于 < 1 MB 测试场景）
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from rdkit import Chem


# ---- 数据结构 ----

@dataclass
class SpectrumSignal:
    """单条化学位移信号.

    来自 <Spectrum 13C N> / <Spectrum 1H N>, 原始形如:
        shift;integral;multiplicity;atom_idx  之间 '|' 分隔
    shift 字段单位 ppm (13C) 或 ppm (1H)
    multiplicity 单字符: S=单峰, D=双峰, T=三重峰, Q=四重峰, ...
    """

    shift: float
    multiplicity: str  # 'S', 'D', 'T', 'Q', 'm', ...
    atom_index: Optional[int]  # SDF 内的原子索引；None 表示不确定
    raw: str = ""


@dataclass
class NmrShiftDb2Entry:
    """nmrshiftdb2 中的一条记录."""

    rec_id: Optional[int]  # <nmrshiftdb2 ID>
    smiles: str
    canonical_smiles: str
    mol_block: str  # 原始 MOL block 字符串（不含 $$$$）
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    solvent: Optional[str] = None
    field_mhz: Optional[float] = None
    temperature_k: Optional[float] = None
    spectrum_13c: List[SpectrumSignal] = field(default_factory=list)
    spectrum_1h: List[SpectrumSignal] = field(default_factory=list)
    extra_tags: Dict[str, str] = field(default_factory=dict)

    @property
    def has_13c(self) -> bool:
        return bool(self.spectrum_13c)

    @property
    def has_1h(self) -> bool:
        return bool(self.spectrum_1h)

    def atom_count(self) -> int:
        """MOL block 中原子数（不依赖 RDKit 解析）."""
        # MOL block 第 4 行: 'aa bb ...'
        lines = self.mol_block.splitlines()
        for line in lines:
            s = line.strip()
            # 原子/键计数行典型形如 ' 18 19  0  0  0  0  0  0  0  0999 V2000'
            parts = s.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                # V2000 标识
                if 'V2000' in s or 'V3000' in s:
                    return int(parts[0])
        return 0

    def mol(self):
        """用 RDKit 解析分子对象（每次新建，调用方负责 sanitize/removeHs）.

        返回 None 表示解析失败.
        """
        try:
            return Chem.MolFromMolBlock(self.mol_block, sanitize=True, removeHs=False)
        except Exception:
            return None


# ---- 错误 ----

class SdfReadError(Exception):
    def __init__(self, msg: str, path: str = "", stage: str = ""):
        super().__init__(f"[{stage}] {msg} (path={path})")
        self.path = path
        self.stage = stage


# ---- 内部解析函数 ----

_RECORD_SEP = b'$$$$\r\n'


def _split_records(raw: bytes) -> Iterator[bytes]:
    """按 $$$$\\r\\n 切分. 不要求末尾有 $$$$."""
    if not raw:
        return
    # 用 split 会丢失分隔符，我们手动扫描
    start = 0
    while True:
        idx = raw.find(_RECORD_SEP, start)
        if idx == -1:
            # 最后一段
            tail = raw[start:]
            if tail.strip():
                yield tail
            return
        chunk = raw[start:idx]
        if chunk.strip():
            yield chunk
        start = idx + len(_RECORD_SEP)


def _parse_tags(record_text: str) -> Tuple[Optional[str], Dict[str, str]]:
    """从 SDF 记录文本中提取 SDF tag 数据部分.

    返回 (m mol_block 字符串, dict[tag, value]).
    MOL block 以 'M  END' 结束; 之后的所有 '> <TAG>' 块为 SDF tag.
    """
    # 找 'M  END'
    m_end_marker = 'M  END'
    idx = record_text.find(m_end_marker)
    if idx == -1:
        # 没有 M  END 视为无效
        raise SdfReadError(f"记录缺少 'M  END' 终止符")
    # MOL block 包含 M  END 这一行
    mol_block = record_text[: idx + len(m_end_marker)]
    # tag 部分：从 M  END 之后开始
    tag_section = record_text[idx + len(m_end_marker):]

    tags: Dict[str, str] = {}
    # 块结构:
    #   > <TAG1>\n
    #   value1\n
    #   \n
    #   > <TAG2>\n
    #   value2\n
    #   \n
    # 也可以 value 一行或多行；空行结束当前 tag.
    lines = tag_section.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()
        if s.startswith('>') and '<' in s:
            # 提取 <TAG>
            try:
                tag_name = s.split('<', 1)[1].split('>', 1)[0]
            except IndexError:
                i += 1
                continue
            # 收集值（直到空行或下一个 '>')
            i += 1
            value_lines: List[str] = []
            while i < n:
                t = lines[i]
                if t.strip().startswith('>'):
                    break
                if t.strip() == '':
                    # 遇到空行：表示 tag 结束
                    i += 1
                    break
                value_lines.append(t)
                i += 1
            tags[tag_name] = '\n'.join(value_lines).strip()
        else:
            i += 1
    return mol_block, tags


def _parse_spectrum_value(value: str) -> List[SpectrumSignal]:
    """解析 <Spectrum 13C N> / <Spectrum 1H N> 的值.

    形如:
      23.3;0.0Q;8|23.3;0.0Q;9|23.5;0.0Q;5|...
    字段之间 ';' 分隔: shift;integral;multiplicity;atom_idx
    各信号之间 '|' 分隔
    """
    if not value:
        return []
    out: List[SpectrumSignal] = []
    for token in value.split('|'):
        token = token.strip()
        if not token:
            continue
        parts = token.split(';')
        # 至少需要 shift 字段
        try:
            shift = float(parts[0])
        except (ValueError, IndexError):
            continue
        multiplicity = ''
        atom_idx = None
        # integral 在 parts[1], 但格式是 "0.0Q" → integral=0.0, mult=Q
        if len(parts) >= 2:
            s2 = parts[1].strip()
            # 取末尾字符作为 multiplicity
            if s2 and s2[-1].isalpha():
                multiplicity = s2[-1].upper()
            else:
                multiplicity = ''
        if len(parts) >= 3:
            try:
                atom_idx = int(parts[2].strip())
            except (ValueError, IndexError):
                atom_idx = None
        out.append(SpectrumSignal(
            shift=shift,
            multiplicity=multiplicity,
            atom_index=atom_idx,
            raw=token,
        ))
    return out


def _parse_first_spectrum(tags: Dict[str, str], prefix: str) -> List[SpectrumSignal]:
    """从 tags 中挑出第一个 <Spectrum {prefix} 0> 或 <Spectrum {prefix} 1> ..."""
    candidates = sorted(
        (k for k in tags.keys()
         if k.startswith(f'Spectrum {prefix} ') and k[len(f'Spectrum {prefix} '):].isdigit()),
        key=lambda k: int(k[len(f'Spectrum {prefix} '):]),
    )
    if not candidates:
        return []
    # 优先使用 0（默认常规谱），缺则用 1（DEPT/子谱）
    for k in (f'Spectrum {prefix} 0',):
        if k in tags:
            return _parse_spectrum_value(tags[k])
    return _parse_spectrum_value(tags[candidates[0]])


def _solvent_clean(value: str) -> Optional[str]:
    """<Solvent> 形如 '0:CDCl3 '  → 'CDCl3'."""
    if not value:
        return None
    v = value.strip()
    if ':' in v:
        v = v.split(':', 1)[1]
    return v.strip() or None


def _float_first_token(value: str) -> Optional[float]:
    if not value:
        return None
    try:
        return float(value.strip().split()[0])
    except (ValueError, IndexError):
        return None


def _int_first_token(value: str) -> Optional[int]:
    if not value:
        return None
    try:
        return int(value.strip().split()[0])
    except (ValueError, IndexError):
        return None


def _to_smiles(mol_block: str) -> Tuple[str, str]:
    """从 MOL block 转 SMILES. 返回 (canonical, isomeric).

    None 表示解析失败.
    """
    try:
        mol = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=False)
    except Exception:
        return '', ''
    if mol is None:
        return '', ''
    try:
        canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return '', ''
    return canon or '', isomeric or ''


def _record_to_entry(record_bytes: bytes) -> Optional[NmrShiftDb2Entry]:
    """单条记录 bytes → NmrShiftDb2Entry. 解析失败返回 None."""
    try:
        text = record_bytes.decode('utf-8', errors='replace')
    except Exception:
        return None
    try:
        mol_block, tags = _parse_tags(text)
    except SdfReadError:
        return None
    canon, isomeric = _to_smiles(mol_block)
    if not isomeric:
        # 没解析出 SMILES 视为无效记录
        return None
    return NmrShiftDb2Entry(
        rec_id=_int_first_token(tags.get('nmrshiftdb2 ID', '')),
        smiles=isomeric,
        canonical_smiles=canon,
        mol_block=mol_block,
        inchi=tags.get('INChI') or None,
        inchikey=tags.get('INChI key') or None,
        solvent=_solvent_clean(tags.get('Solvent', '')),
        field_mhz=_float_first_token(tags.get('Field Strength [MHz]', '')),
        temperature_k=_float_first_token(tags.get('Temperature [K]', '')),
        spectrum_13c=_parse_first_spectrum(tags, '13C'),
        spectrum_1h=_parse_first_spectrum(tags, '1H'),
        extra_tags=tags,
    )


# ---- 顶层接口 ----

def stream_nmrshiftdb2_sdf(path: str) -> Iterator[NmrShiftDb2Entry]:
    """流式读取 nmrshiftdb2 SDF. 每次 yield 一条 NmrShiftDb2Entry.

    跳过 SMILES 解析失败的记录（不抛错）.
    """
    p = Path(path)
    if not p.exists():
        raise SdfReadError(f"SDF 文件不存在: {path}", path, stage="sdf-open")
    if not p.is_file():
        raise SdfReadError(f"路径不是文件: {path}", path, stage="sdf-open")
    # 用二进制读（避免 RDKit C++ 中文路径 bug）
    with open(p, 'rb') as f:
        raw = f.read()
    total = 0
    ok = 0
    for rec_bytes in _split_records(raw):
        total += 1
        entry = _record_to_entry(rec_bytes)
        if entry is not None:
            ok += 1
            yield entry
    # 在迭代结束后不打印统计，由调用方决定如何报告


def read_nmrshiftdb2_sdf(path: str) -> List[NmrShiftDb2Entry]:
    """一次性读全（适用于测试或小样本，不适用于 58187 记录）."""
    return list(stream_nmrshiftdb2_sdf(path))