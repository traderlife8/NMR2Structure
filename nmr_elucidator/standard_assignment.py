# -*- coding: utf-8 -*-
"""读取化学家 Excel 标准归属表 (已知化合物流程).

提供:
  - StandardAssignmentRow / StandardAssignmentTable: 标准化数据结构
  - read_standard_assignment_table(path): 主读取函数
  - find_standard_assignment(folder): 在目录中自动发现数据归属表
  - find_known_smiles(folder): 在同目录读取 smiles.txt
  - extract_smiles_from_xlsx(path): 从 xlsx 元信息区提取 SMILES

化学家 Excel 标准格式 (clxh100数据归属表.xlsx 范本):
  Row 1:  No. | 13C NMR | Types | 1H NMR (mult., J in Hz) | 1H-1H COSY | HMBC | ROESY
  Row 2+: 数据行, 含 13C 化学位移 + 碳级数 (CH/CH2/CH3/C) + 1H 字符串 (位移 + 多重性 + J)

设计要点:
  - 直接用 openpyxl 而非 pandas, 避免特殊字符 (如 '′') / 糖基标号 (26-O-Glc-1′) 的解析问题
  - 自适应表头位置 (含 'No.' 的行) + 列名变体
  - 类型归一化: Cq/C/C → 'C'; CH2/CH3 → 保留
  - 1H 字符串解析: 提取首个数值作为化学位移, 剩余文本作为 multiplicity (含 J)
  - 不依赖 RDKit / numpy (纯标准库 + openpyxl)
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple


# 标准归属表可能的文件名模式 (优先级从高到低)
# 例: clxh100数据归属表.xlsx, clxh100数据归属.xlsx, clxh100归属表.xlsx, 归属表.xlsx
_STANDARD_TABLE_PATTERNS = [
    '*数据归属表*.xlsx',
    '*数据归属*.xlsx',
    '*归属表*.xlsx',
    '*归属*.xlsx',
    'assignment_table.xlsx',
]

_SMILES_FILE_CANDIDATES = ['smiles.txt', '*.smi', '*.smiles']


@dataclass
class StandardAssignmentRow:
    """标准归属表的一行 (已知化合物 ground truth)."""

    no_label: str = ''                       # 原子标号 (如 "1α", "Glc-1′")
    c_ppm: Optional[float] = None            # 13C 化学位移
    carbon_type: str = ''                    # CH/CH2/CH3/C (归一化)
    c_raw: str = ''                          # 13C 原始字符串
    h_ppm: Optional[float] = None            # 1H 化学位移 (从字符串中提取的首个数值)
    h_full_text: str = ''                    # 1H 整段文本 (含多重性 + J)
    h_mult: str = ''                         # 多重性 (s/d/t/q/m/dd/...)
    j_values: List[float] = field(default_factory=list)  # J 耦合常数 (Hz)
    xlsx_row: int = 0                        # 在原 xlsx 中的实际行号 (1-based)
    is_quaternary: bool = False              # 是否为季碳行 (无 1H 列)


@dataclass
class StandardAssignmentTable:
    """完整的标准归属表 (化学家 ground truth)."""

    rows: List[StandardAssignmentRow] = field(default_factory=list)
    file_path: str = ''                      # xlsx 路径
    sheet_name: str = ''                     # sheet 名
    header_row_idx: int = 0                  # 表头所在行 (1-based, 用于程序输出 xlsx 写回)
    known_smiles: str = ''                   # 同目录的 smiles.txt (可选)
    n_rows: int = 0                          # 数据行数 (过滤后)
    n_carbon: int = 0                        # 含 13C 的行数
    n_proton: int = 0                        # 含 1H 的行数
    n_quaternary: int = 0                    # 季碳行数


# ============ Excel 解析 ============

def _find_header_row(ws, max_rows: int = 30, max_cols: int = 8) -> int:
    """在 worksheet 中定位表头行 (含 'No.'/'No'/'No,'/'编号' 等标识的行).

    返回 1-based 行号. 未找到时返回 -1.
    """
    for r in range(1, min(max_rows + 1, ws.max_row + 1)):
        for c in range(1, min(max_cols + 1, ws.max_column + 1)):
            v = ws.cell(row=r, column=c).value
            if v is None:
                continue
            s = str(v).strip()
            # 表头标识: 'No.', 'No', 'No ', '编号', 'Position' (容忍常见变体)
            if s in ('No.', 'No', 'No ', 'No\t', '编号', 'Position'):
                return r
            if s.lower().startswith('no.') or s.lower().startswith('no ') or s.lower() == 'no':
                return r
    return -1


def _build_col_map(ws, header_row: int) -> dict:
    """根据表头行内容建立列名 → 列号的映射.

    列名变体容忍:
      - 13C: '13C NMR', '13C', 'δC', 'δ13C'
      - Types: 'Types', 'Type', 'C-type', 'Multiplicity', '碳类型'
      - 1H: '1H NMR (mult., J in Hz)', '1H NMR', '1H', 'δH', 'δ1H'
      - COSY: '1H-1H COSY', 'COSY'
      - HMBC: 'HMBC', 'H→C'
      - ROESY: 'ROESY', 'NOESY'
    """
    col_map = {}
    for c in range(1, ws.max_column + 1):
        v = ws.cell(row=header_row, column=c).value
        if v is None:
            continue
        s = str(v).strip()
        if not s:
            continue

        # 顺序敏感: 'COSY' 包含 '1H', 必须先判 COSY 再判 1H
        # 顺序敏感: 'HMBC' 含 'C', 但 '13C' 必须先判
        if s.startswith('No') or s in ('编号', 'Position'):
            if 'No' not in col_map:
                col_map['No'] = c
        elif '13C' in s or 'δC' in s or 'δ13C' in s or s == 'C':
            if '13C' not in col_map:
                col_map['13C'] = c
        elif s in ('Types', 'Type', 'C-type', 'Multiplicity', '碳类型'):
            col_map['Types'] = c
        elif 'COSY' in s:
            col_map['COSY'] = c
        elif 'HMBC' in s or 'H→C' in s:
            col_map['HMBC'] = c
        elif 'ROESY' in s or 'NOESY' in s:
            col_map['ROESY'] = c
        elif '1H' in s and 'COSY' not in s:
            if '1H' not in col_map:
                col_map['1H'] = c
        elif s in ('1H NMR', 'δH', 'δ1H'):
            if '1H' not in col_map:
                col_map['1H'] = c
    return col_map


def _normalize_type(t) -> str:
    """规范化碳类型: Cq → C, 去空格, 大写. 空值保留为空字符串."""
    if t is None:
        return ''
    s = str(t).strip()
    if not s:
        return ''
    # Cq / CQ / q / C → 'C'
    if s.upper() in ('CQ', 'Q', 'C', 'C-Q'):
        return 'C'
    # CH2/CH3 / CH2ORCH3 → 保留原文 (已含分隔)
    return s.upper().replace(' ', '')


# 1H 字符串解析: 提取数值 + 多重性 + J
# 例: '3.82 m' / '1.16, m' / '4.54, m' / '4.54 d (7.7)' / '5.74, m'
#     '3.92, dd (9.5, 6.5)' / '0.94, d (6.7)' / '2.11, s' / '2.63, m (2H)'
_H_SHIFT_RE = re.compile(r'^\s*([\d.]+)\s*,?\s*(.*)$')
_H_J_RE = re.compile(r'\(?\s*([\d.]+)\s*\)?')


def _parse_1h_text(h_text) -> Tuple[Optional[float], str, str, List[float]]:
    """解析 1H NMR 字符串.

    Returns:
        (h_ppm, full_text, h_mult, j_values)
        h_ppm: 首个浮点数 (化学位移)
        full_text: 整段文本 (位移 + 逗号 + 多重性 + J)
        h_mult: 多重性 (s/d/t/q/m/dd/dt/dq/td/...)
        j_values: J 耦合常数列表 (Hz)
    """
    if h_text is None:
        return None, '', '', []
    full = str(h_text).strip()
    if not full:
        return None, '', '', []

    m = _H_SHIFT_RE.match(full)
    if not m:
        return None, full, '', []

    ppm_str = m.group(1)
    try:
        h_ppm = float(ppm_str)
    except ValueError:
        return None, full, '', []

    rest = m.group(2).strip()

    # 提取多重性: 第一个非数字的英文 token
    # 'm', 'd (7.7)', 'dd (9.5, 6.5)', 'm (2H)', 's'
    h_mult = ''
    mult_m = re.match(r'([a-zA-Z]+)', rest)
    if mult_m:
        h_mult = mult_m.group(1).lower()
        # 标准化 'CH3' 等异常值
        if h_mult not in ('s', 'd', 't', 'q', 'm', 'dd', 'dt', 'dq', 'td',
                          'tt', 'ddd', 'ddt', 'dtd', 'dtt', 'qdd'):
            h_mult = ''  # 未知 → 空

    # 提取 J 值: 在括号内找所有数值
    j_values = []
    j_matches = _H_J_RE.findall(rest)
    for jm in j_matches:
        try:
            jv = float(jm)
            # J 值应在 0-30 Hz 之间, 排除 '2H' 等质子数标记
            if 0.1 < jv < 30:
                j_values.append(jv)
        except ValueError:
            pass

    return h_ppm, full, h_mult, j_values


# ============ SMILES 提取 (从 xlsx 元信息区) ============

def extract_smiles_from_xlsx(path: str | Path) -> str:
    """从 xlsx 表头上方元信息区提取 SMILES.

    扫描 Row 1..header_row-1, 寻找含 'SMILES' 关键字的行, 取后续列的非空字符串作为值.
    找不到时返回空字符串.

    注: 当前程序输出 xlsx **不再包含顶部元信息行** (元信息独立到 meta.txt),
    此函数主要用于兼容旧版 xlsx, 新版请优先调用 read_assignment_meta() 获取 SMILES.
    """
    from openpyxl import load_workbook
    try:
        wb = load_workbook(str(path), data_only=True, read_only=True)
    except Exception:
        return ''

    try:
        ws = wb.active
        header_idx = _find_header_row(ws)
        scan_max = header_idx - 1 if header_idx > 0 else min(5, ws.max_row)
        for r in range(1, max(1, scan_max + 1)):
            for c in range(1, min(4, ws.max_column + 1)):
                v = ws.cell(row=r, column=c).value
                if v is None:
                    continue
                s = str(v).strip()
                if 'SMILES' in s.upper():
                    # 取同行后续列的非空值
                    for cc in range(c + 1, ws.max_column + 1):
                        vv = ws.cell(row=r, column=cc).value
                        if vv is None:
                            continue
                        s2 = str(vv).strip()
                        if len(s2) > 5:
                            return s2
        return ''
    finally:
        wb.close()


# ============ meta.txt 元信息读取 ============

def read_assignment_meta(path: str | Path) -> Dict[str, str]:
    """读取归属表元信息文件 meta.txt (UTF-8).

    文件格式 (与 write_assignment_meta 输出对齐):
        # NMR 归属表元信息 (与 assignment_table.xlsx 分离存放)
        # 生成时间: YYYY-MM-DD HH:MM:SS

        Solvent: CDCl3
        Formula: C33H52O11
        Field_MHz: 500.18
        HRESIMS: [M+Na]+ 647.3407

    Args:
        path: meta.txt 路径

    Returns:
        字典 {字段名: 字段值}. 文件不存在时返回空字典.
    """
    path = Path(path)
    result: Dict[str, str] = {}
    if not path.is_file():
        return result
    try:
        text = path.read_text(encoding='utf-8')
    except (OSError, UnicodeDecodeError):
        return result
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        if '(无元信息)' in s:
            continue
        # 格式: 'Key: Value' (Value 可含空格)
        if ':' in s:
            key, _, val = s.partition(':')
            key = key.strip()
            val = val.strip()
            if key:
                result[key] = val
    return result


# ============ 主读取函数 ============

def read_standard_assignment_table(
    path: str | Path,
    known_smiles: Optional[str] = None,
) -> StandardAssignmentTable:
    """读取化学家 Excel 标准归属表.

    Args:
        path: xlsx 文件路径
        known_smiles: 可选, 已知的 SMILES (从同目录 smiles.txt). 也可后续通过 table.known_smiles 设置.

    Returns:
        StandardAssignmentTable 实例

    Raises:
        FileNotFoundError: 文件不存在
        ValueError: 文件格式不合法 (找不到表头 / No 列 / 13C 列)
    """
    from openpyxl import load_workbook

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'标准归属表不存在: {path}')

    wb = load_workbook(str(path), data_only=True, read_only=True)
    try:
        ws = wb.active
        sheet_name = ws.title or ''

        header_idx = _find_header_row(ws)
        if header_idx < 0:
            raise ValueError(
                f'未找到表头行 (含 "No." 的行). 文件: {path}, max_row={ws.max_row}'
            )

        col_map = _build_col_map(ws, header_idx)
        if 'No' not in col_map:
            raise ValueError(
                f'表头缺少必要列 (No.). 文件: {path}, '
                f'实际识别列: {list(col_map.keys())}'
            )
        if '13C' not in col_map:
            raise ValueError(
                f'表头缺少必要列 (13C). 文件: {path}, '
                f'实际识别列: {list(col_map.keys())}'
            )

        table = StandardAssignmentTable(
            file_path=str(path),
            sheet_name=sheet_name,
            header_row_idx=header_idx,
            known_smiles=known_smiles or extract_smiles_from_xlsx(path),
        )

        # 遍历数据行
        for r in range(header_idx + 1, ws.max_row + 1):
            no_val = ws.cell(row=r, column=col_map['No']).value
            if no_val is None:
                continue
            no_str = str(no_val).strip()
            if not no_str:
                continue
            # 过滤表尾注释 / 子标题行
            if no_str.startswith('No') or no_str.startswith('注') \
                    or no_str in ('STRUCTURE', '未归属', 'Notes', '注释', '说明'):
                continue
            # 过滤元信息残留 (出现关键字 NMR / ppm / 归属)
            if 'NMR' in no_str or '归属' in no_str or 'ppm' in no_str:
                continue

            # 13C 数值
            c_raw = ws.cell(row=r, column=col_map['13C']).value
            c_str = str(c_raw).strip() if c_raw is not None else ''
            try:
                c_ppm = float(c_str)
            except (ValueError, TypeError):
                c_ppm = None

            # 碳类型
            c_type_raw = (ws.cell(row=r, column=col_map['Types']).value
                          if 'Types' in col_map else None)
            carbon_type = _normalize_type(c_type_raw)

            # 1H 字符串
            h_raw = (ws.cell(row=r, column=col_map['1H']).value
                     if '1H' in col_map else None)
            h_ppm, h_full, h_mult, j_values = _parse_1h_text(h_raw)

            is_quat = (h_ppm is None) and (carbon_type == 'C')

            row = StandardAssignmentRow(
                no_label=no_str,
                c_ppm=c_ppm,
                carbon_type=carbon_type,
                c_raw=c_str,
                h_ppm=h_ppm,
                h_full_text=h_full,
                h_mult=h_mult,
                j_values=j_values,
                xlsx_row=r,
                is_quaternary=is_quat,
            )
            table.rows.append(row)

        # 统计
        table.n_rows = len(table.rows)
        table.n_carbon = sum(1 for r in table.rows if r.c_ppm is not None)
        table.n_proton = sum(1 for r in table.rows if r.h_ppm is not None)
        table.n_quaternary = sum(1 for r in table.rows if r.is_quaternary)

        return table
    finally:
        wb.close()


# ============ 自动发现 ============

def find_standard_assignment(folder: str | Path) -> Optional[Path]:
    """在 folder 目录中自动发现"标准归属表".

    匹配优先级:
      1. 含 '数据归属表' 的 xlsx (优先级最高, 如 clxh100数据归属表.xlsx)
      2. 含 '数据归属' 的 xlsx
      3. 含 '归属表' 的 xlsx
      4. 含 '归属' 的 xlsx
      5. assignment_table.xlsx (回退, 与程序输出同名时存在歧义)

    Returns:
        找到的 xlsx 路径, 或 None.
    """
    folder = Path(folder)
    if not folder.is_dir():
        return None

    found: list[Path] = []
    for pattern in _STANDARD_TABLE_PATTERNS:
        try:
            matches = sorted(folder.glob(pattern))
        except OSError:
            matches = []
        for m in matches:
            if m.is_file() and m.suffix.lower() == '.xlsx':
                if m not in found:
                    found.append(m)

    if not found:
        return None
    return found[0]


def find_known_smiles(folder: str | Path) -> str:
    """在 folder 目录中查找 SMILES 文件, 返回 SMILES 字符串 (无则为空).

    优先级:
      1. smiles.txt
      2. *.smi / *.smiles (按字典序取第一)
    """
    folder = Path(folder)
    if not folder.is_dir():
        return ''

    for name in _SMILES_FILE_CANDIDATES:
        try:
            if '*' in name:
                matches = sorted(folder.glob(name))
                for m in matches:
                    if m.is_file():
                        txt = m.read_text(encoding='utf-8').strip()
                        if txt:
                            return txt
            else:
                p = folder / name
                if p.is_file():
                    txt = p.read_text(encoding='utf-8').strip()
                    if txt:
                        return txt
        except (OSError, UnicodeDecodeError):
            continue
    return ''


def discover_known_compound(
    compound_name: str,
    real_structure_dir: str | Path,
) -> Optional[Tuple[Path, str]]:
    """在 real_structure_dir/<compound_name>/ 下检测"已知化合物"资源.

    Args:
        compound_name: 化合物名 (从 --folder 末段推断, 如 'clxh100')
        real_structure_dir: 真实结构根目录 (如 E:\\NMR2Structure\\真实结构)

    Returns:
        (标准归属表路径, 已知 SMILES) 或 None (未发现).
    """
    real_dir = Path(real_structure_dir)
    sub_dir = real_dir / compound_name
    if not sub_dir.is_dir():
        return None

    table_path = find_standard_assignment(sub_dir)
    if table_path is None:
        return None

    known_smiles = find_known_smiles(sub_dir)
    return (table_path, known_smiles)


# ============ 便捷方法 ============

def get_obs_peaks(table: StandardAssignmentTable) -> Tuple[List[float], List[float]]:
    """从标准归属表提取实测 13C / 1H 化学位移列表 (用作 ground truth 验证).

    Returns:
        (carbon_shifts, proton_shifts), 已去重 + 排序(降序).
    """
    cs = sorted({r.c_ppm for r in table.rows if r.c_ppm is not None}, reverse=True)
    hs = sorted({r.h_ppm for r in table.rows if r.h_ppm is not None}, reverse=True)
    return cs, hs


def get_row_by_no(table: StandardAssignmentTable, no_label: str) -> Optional[StandardAssignmentRow]:
    """按 no_label (区分大小写) 查找行."""
    for r in table.rows:
        if r.no_label == no_label:
            return r
    return None


__all__ = [
    'StandardAssignmentRow',
    'StandardAssignmentTable',
    'read_standard_assignment_table',
    'find_standard_assignment',
    'find_known_smiles',
    'extract_smiles_from_xlsx',
    'read_assignment_meta',
    'discover_known_compound',
    'get_obs_peaks',
    'get_row_by_no',
    'compare_program_to_standard',
    'ComparisonResult',
    'write_matched_xlsx',
    'write_joint_xlsx',
    'write_comparison_report',
]


# ============ 程序输出归属表对比 ============

@dataclass
class ComparisonResult:
    """程序输出 xlsx 与标准归属表的对比结果."""

    matches: list = field(default_factory=list)             # [(prog_row_dict, std_row, c_diff, h_diff, type_ok)]
    unmatched_program: list = field(default_factory=list)   # 程序行未匹配
    unmatched_standard: list = field(default_factory=list)  # 标准行未匹配
    type_mismatches: list = field(default_factory=list)
    prog_smiles: str = ''                                   # 程序输出 SMILES (从 meta.txt SMILES 字段读)
    std_smiles: str = ''                                    # 标准归属表对应的 SMILES
    smiles_equal: Optional[bool] = None                    # RDKit canonical SMILES 判等 (None=未启用)
    prog_xlsx_row: int = 0                                  # 程序输出 xlsx 表头行 (1-based, 新版 = 1)
    prog_meta: dict = field(default_factory=dict)           # 从 meta.txt 读到的元信息 (Solvent / Formula / Field_MHz / HRESIMS ...)
    match_rate: float = 0.0                                 # 总行匹配率
    verdict: str = ''                                       # '一致' / '部分一致' / '不一致'
    n_prog: int = 0
    n_prog_c: int = 0
    n_prog_h: int = 0
    n_std: int = 0
    n_std_c: int = 0
    n_std_h: int = 0
    c13_matched: int = 0                                    # 含 13C 的匹配对数


def _load_program_table(prog_xlsx_path: str | Path,
                         prog_meta_path: str | Path | None = None):
    """读取程序输出归属表 xlsx (复用 standard_assignment 的 _find_header_row / _build_col_map).

    新版 (Part 3.x): 程序输出 xlsx **不含** 顶部元信息行, 表头在 row 1.
    元信息 (Solvent / Formula / Field_MHz / HRESIMS / SMILES 等) 单独存在 meta.txt 中.

    Args:
        prog_xlsx_path: 程序输出 xlsx 路径
        prog_meta_path: 可选 meta.txt 路径. 提供时 SMILES 等元信息从该文件读取.

    Returns:
        (rows, prog_smiles, prog_meta_dict, header_idx, sheet_name)
    """
    from openpyxl import load_workbook

    path = Path(prog_xlsx_path)
    if not path.exists():
        raise FileNotFoundError(f'程序输出 xlsx 不存在: {path}')

    wb = load_workbook(str(path), data_only=True, read_only=True)
    try:
        ws = wb.active
        sheet_name = ws.title or ''
        header_idx = _find_header_row(ws)
        if header_idx < 0:
            raise ValueError(f'程序输出 xlsx 未找到表头: {path}')
        col_map = _build_col_map(ws, header_idx)
        if 'No' not in col_map or '13C' not in col_map:
            raise ValueError(f'程序输出 xlsx 表头缺少列 (No/13C): {path}')

        rows = []
        for r in range(header_idx + 1, ws.max_row + 1):
            no_val = ws.cell(row=r, column=col_map['No']).value
            if no_val is None:
                continue
            no_str = str(no_val).strip()
            if not no_str or no_str.startswith('No') \
                    or no_str in ('STRUCTURE', '未归属', 'Notes', '注释'):
                continue
            if 'NMR' in no_str or '归属' in no_str or 'ppm' in no_str:
                continue
            c_raw = ws.cell(row=r, column=col_map['13C']).value
            c_str = str(c_raw).strip() if c_raw is not None else ''
            try:
                c_ppm = float(c_str)
            except (ValueError, TypeError):
                c_ppm = None
            c_type_raw = (ws.cell(row=r, column=col_map['Types']).value
                          if 'Types' in col_map else None)
            carbon_type = _normalize_type(c_type_raw)
            h_raw = (ws.cell(row=r, column=col_map['1H']).value
                     if '1H' in col_map else None)
            h_ppm, _, _, _ = _parse_1h_text(h_raw)
            rows.append({
                'no': no_str,
                'c_ppm': c_ppm,
                'type': carbon_type,
                'h_ppm': h_ppm,
                'xlsx_row': r,
            })

        # 元信息从 meta.txt 读取 (新流程); 若 meta.txt 不存在, 尝试从 xlsx 顶部读 (旧版兼容)
        prog_meta: dict = {}
        if prog_meta_path is not None:
            prog_meta = read_assignment_meta(prog_meta_path)
        if not prog_meta:
            # 回退: 尝试从 xlsx 顶部读 (旧版兼容)
            legacy_smiles = extract_smiles_from_xlsx(path)
            if legacy_smiles:
                prog_meta['SMILES'] = legacy_smiles
        prog_smiles = prog_meta.get('SMILES', '')

        return rows, prog_smiles, prog_meta, header_idx, sheet_name
    finally:
        wb.close()


def _joint_match(prog_rows, std_table, c_tol: float, h_tol: float,
                type_placeholder_ok: bool = True):
    """以 (13C, 1H, 碳类型) 联合约束做 Hungarian 一对一最优匹配.

    与 util/compare_tables.py 的 joint_match 逻辑一致, 但此处 std_rows 用
    StandardAssignmentRow 而不是 DataFrame.

    type_placeholder_ok: True 时, 程序输出的 'CH2/CH3' 与标准的 CH2/CH3 兼容
    """
    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError:
        # scipy 不可用 → 退化为贪心匹配 (按 13C 差最小)
        return _greedy_joint_match(prog_rows, std_table.rows, c_tol, h_tol,
                                   type_placeholder_ok)

    n1, n2 = len(prog_rows), len(std_table.rows)
    if n1 == 0 or n2 == 0:
        return [], list(prog_rows), list(std_table.rows)

    BIG = 1e6
    TYPE_PENALTY = 0.5
    cost = [[BIG] * n2 for _ in range(n1)]

    def _type_compatible(t1: str, t2: str) -> bool:
        """两个碳类型是否兼容 (含 CH2/CH3 占位符兼容)."""
        if not t1 or not t2:
            return True
        if t1 == t2:
            return True
        if type_placeholder_ok:
            if t1 == 'CH2/CH3' and t2 in ('CH2', 'CH3'):
                return True
            if t2 == 'CH2/CH3' and t1 in ('CH2', 'CH3'):
                return True
        return False

    for i, pr in enumerate(prog_rows):
        c1, h1, t1 = pr['c_ppm'], pr['h_ppm'], pr['type']
        for j, sr in enumerate(std_table.rows):
            c2, h2, t2 = sr.c_ppm, sr.h_ppm, sr.carbon_type

            # 13C 必须满足容差
            if c1 is not None and c2 is not None:
                if abs(c1 - c2) > c_tol:
                    continue
                c_diff = abs(c1 - c2)
            elif (c1 is None) != (c2 is None):
                continue
            else:
                c_diff = 0.0

            # 1H 容差
            if h1 is not None and h2 is not None:
                if abs(h1 - h2) > h_tol:
                    continue
                h_diff = abs(h1 - h2)
            else:
                h_diff = 0.0

            d = (c_diff / c_tol if c_tol > 0 else 0) \
                + (h_diff / h_tol if h_tol > 0 else 0)
            if not _type_compatible(t1, t2):
                # 改进点 #32: dC > 2.0 且类型不一致时禁止匹配 (避免远距离误配)
                if c_diff > 2.0:
                    continue
                d += TYPE_PENALTY
            cost[i][j] = d

    row_ind, col_ind = linear_sum_assignment(cost)
    matches, m1, m2 = [], set(), set()
    for r, c in zip(row_ind, col_ind):
        if cost[r][c] >= BIG:
            continue
        c1 = prog_rows[r]['c_ppm']
        c2 = std_table.rows[c].c_ppm
        h1 = prog_rows[r]['h_ppm']
        h2 = std_table.rows[c].h_ppm
        t1 = prog_rows[r]['type']
        t2 = std_table.rows[c].carbon_type
        c_diff = abs(c1 - c2) if (c1 is not None and c2 is not None) else None
        h_diff = abs(h1 - h2) if (h1 is not None and h2 is not None) else None
        type_ok = _type_compatible(t1, t2)
        matches.append((prog_rows[r], std_table.rows[c], c_diff, h_diff, type_ok))
        m1.add(r)
        m2.add(c)

    unmatched_prog = [prog_rows[i] for i in range(n1) if i not in m1]
    unmatched_std = [std_table.rows[i] for i in range(n2) if i not in m2]
    return matches, unmatched_prog, unmatched_std


def _greedy_joint_match(prog_rows, std_rows, c_tol, h_tol,
                        type_placeholder_ok: bool = True):
    """无 scipy 时的贪心回退: 按 13C 差最小优先, 一对一匹配."""

    def _type_compatible(t1: str, t2: str) -> bool:
        if not t1 or not t2:
            return True
        if t1 == t2:
            return True
        if type_placeholder_ok:
            if t1 == 'CH2/CH3' and t2 in ('CH2', 'CH3'):
                return True
            if t2 == 'CH2/CH3' and t1 in ('CH2', 'CH3'):
                return True
        return False

    BIG = 1e6
    candidates = []
    for i, pr in enumerate(prog_rows):
        for j, sr in enumerate(std_rows):
            c1, c2 = pr['c_ppm'], sr.c_ppm
            h1, h2 = pr['h_ppm'], sr.h_ppm
            if c1 is not None and c2 is not None and abs(c1 - c2) > c_tol:
                continue
            if h1 is not None and h2 is not None and abs(h1 - h2) > h_tol:
                continue
            c_d = abs(c1 - c2) if (c1 is not None and c2 is not None) else 0
            h_d = abs(h1 - h2) if (h1 is not None and h2 is not None) else 0
            candidates.append((c_d + h_d, i, j))
    candidates.sort()
    used_p, used_s, matches = set(), set(), []
    for score, i, j in candidates:
        if score >= BIG or i in used_p or j in used_s:
            continue
        pr, sr = prog_rows[i], std_rows[j]
        t1, t2 = pr['type'], sr.carbon_type
        type_ok = _type_compatible(t1, t2)
        c1, c2 = pr['c_ppm'], sr.c_ppm
        h1, h2 = pr['h_ppm'], sr.h_ppm
        c_diff = abs(c1 - c2) if (c1 is not None and c2 is not None) else None
        h_diff = abs(h1 - h2) if (h1 is not None and h2 is not None) else None
        matches.append((pr, sr, c_diff, h_diff, type_ok))
        used_p.add(i)
        used_s.add(j)
    unmatched_prog = [pr for i, pr in enumerate(prog_rows) if i not in used_p]
    unmatched_std = [sr for j, sr in enumerate(std_rows) if j not in used_s]
    return matches, unmatched_prog, unmatched_std


def compare_program_to_standard(
    prog_xlsx_path: str | Path,
    std_table: StandardAssignmentTable,
    prog_meta_path: str | Path | None = None,
    c_tol: float = 2.5,
    h_tol: float = 1.0,
    use_smiles: bool = False,
    type_placeholder_ok: bool = True,
) -> ComparisonResult:
    """对比程序输出 xlsx 与标准归属表.

    Args:
        prog_xlsx_path: 程序输出 xlsx 路径 (化学家 Excel 格式)
        std_table: 已加载的标准归属表
        prog_meta_path: 可选 meta.txt 路径. 提供时元信息从该文件读取.
        c_tol: 13C 匹配容差 (默认 0.5 ppm, 覆盖实际推导偏差)
        h_tol: 1H 匹配容差 (默认 0.1 ppm)
        use_smiles: 是否用 RDKit canonical SMILES 做化合物结构判等 (复核级别)
        type_placeholder_ok: 是否接受 'CH2/CH3' 占位符与 CH2/CH3 兼容

    Returns:
        ComparisonResult 实例

    约束说明:
      - SMILES 仅在 use_smiles=True 时用于结构判等 (canonical SMILES 一致性),
        不参与行级 ppm 匹配
      - 推导过程的化合物结构验证不依赖 SMILES, 只依赖谱图推断
      - 标准归属表本身的标号 (化学家标号) 保持不变
      - 程序输出 xlsx 的元信息 (Solvent / Formula / Field_MHz / HRESIMS 等)
        从 meta.txt 读取, 不再嵌入 xlsx 顶部
    """
    prog_rows, prog_smiles, prog_meta, prog_header_idx, prog_sheet \
        = _load_program_table(prog_xlsx_path, prog_meta_path=prog_meta_path)

    matches, unmatched_prog, unmatched_std = _joint_match(
        prog_rows, std_table, c_tol, h_tol,
        type_placeholder_ok=type_placeholder_ok,
    )

    type_mismatches = [m for m in matches if not m[4]]
    c13_matched = sum(1 for m in matches if m[2] is not None)

    max_rows = max(len(prog_rows), len(std_table.rows))
    match_rate = (len(matches) / max_rows * 100) if max_rows > 0 else 0.0
    if match_rate >= 80 and len(type_mismatches) == 0:
        verdict = '一致 (两个归属表代表同一化合物)'
    elif match_rate >= 60:
        verdict = '部分一致 (可能同一化合物但归属有差异)'
    else:
        verdict = '不一致 (两个归属表代表不同化合物)'

    smiles_equal = None
    if use_smiles and prog_smiles and std_table.known_smiles:
        try:
            from rdkit import Chem
            mk = Chem.MolFromSmiles(std_table.known_smiles)
            mp = Chem.MolFromSmiles(prog_smiles)
            if mk and mp:
                smiles_equal = (Chem.MolToSmiles(mk) == Chem.MolToSmiles(mp))
        except ImportError:
            smiles_equal = None

    return ComparisonResult(
        matches=matches,
        unmatched_program=unmatched_prog,
        unmatched_standard=unmatched_std,
        type_mismatches=type_mismatches,
        prog_smiles=prog_smiles,
        std_smiles=std_table.known_smiles or '',
        smiles_equal=smiles_equal,
        prog_xlsx_row=prog_header_idx,
        prog_meta=prog_meta,
        match_rate=match_rate,
        verdict=verdict,
        n_prog=len(prog_rows),
        n_prog_c=sum(1 for r in prog_rows if r['c_ppm'] is not None),
        n_prog_h=sum(1 for r in prog_rows if r['h_ppm'] is not None),
        n_std=std_table.n_rows,
        n_std_c=std_table.n_carbon,
        n_std_h=std_table.n_proton,
        c13_matched=c13_matched,
    )


def write_matched_xlsx(
    prog_xlsx_path: str | Path,
    out_xlsx_path: str | Path,
    comparison: ComparisonResult,
    placeholder: str = '—',
) -> None:
    """在程序输出 xlsx 基础上加一列 "标准归属No.", 另存为 out_xlsx_path.

    Args:
        prog_xlsx_path: 原程序输出 xlsx
        out_xlsx_path: 输出 xlsx (建议带 _matched 后缀)
        comparison: compare_program_to_standard 返回值
        placeholder: 未匹配行的占位符
    """
    from openpyxl import load_workbook

    src = Path(prog_xlsx_path)
    dst = Path(out_xlsx_path)
    if not src.exists():
        raise FileNotFoundError(f'程序输出 xlsx 不存在: {src}')

    wb = load_workbook(str(src))
    ws = wb.active
    new_col = ws.max_column + 1
    header_row = comparison.prog_xlsx_row
    ws.cell(row=header_row, column=new_col, value='标准归属No.')

    # prog_xlsx_row (1-based) → std no_label
    prog_to_std_no: dict[int, str] = {}
    for pr, sr, _, _, _ in comparison.matches:
        std_no = str(sr.no_label).strip()
        if std_no in ('', '.', 'None'):
            continue
        prog_to_std_no[pr['xlsx_row']] = std_no

    for r in range(header_row + 1, ws.max_row + 1):
        if r in prog_to_std_no:
            ws.cell(row=r, column=new_col, value=prog_to_std_no[r])
        else:
            # 判断是否为空行/注释行 (避免破坏表结构)
            no_val = ws.cell(row=r, column=1).value
            c_val = ws.cell(row=r, column=2).value
            if no_val is None or str(no_val).strip().startswith('No') \
                    or (c_val is not None and 'NMR' in str(c_val)):
                continue
            ws.cell(row=r, column=new_col, value=placeholder)

    dst.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(dst))


def write_joint_xlsx(
    prog_xlsx_path: str | Path,
    std_table: StandardAssignmentTable,
    comparison: ComparisonResult,
    out_xlsx_path: str | Path,
    placeholder: str = '—',
) -> None:
    """生成双 sheet 的左右拼接 xlsx (推荐输出格式).

    **Sheet 1 — '按标准表为准'**:
      - 每行 = 标准归属表的一行 (行序与标准表完全一致)
      - 列结构: [标准表全部列] + '拼接' + [程序表全部列] + ΔdC / ΔdH / 类型一致
      - 未匹配的标准行右侧程序列填 placeholder, Δ列留空

    **Sheet 2 — '按程序表为准'**:
      - 每行 = 程序输出的一行 (行序与程序 xlsx 完全一致)
      - 列结构: [程序表全部列] + '拼接' + [标准表全部列] + ΔdC / ΔdH / 类型一致
      - 未匹配的程序行右侧标准列填 placeholder, Δ列留空

    此格式取代旧的 write_matched_xlsx() 单 sheet + 単列产品.

    Args:
        prog_xlsx_path: 原程序输出 xlsx (提供列结构 / 列名)
        std_table: 标准归属表实例 (提供表数据 + xlsx_row + 标准列定义)
        comparison: compare_program_to_standard 返回值
        out_xlsx_path: 输出 xlsx 路径
        placeholder: 未匹配单元格的占位符 (默认 '—')
    """
    from openpyxl import Workbook, load_workbook
    from openpyxl.styles import Font, Border, Side

    out_path = Path(out_xlsx_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    src = Path(prog_xlsx_path)
    if not src.exists():
        raise FileNotFoundError(f'程序输出 xlsx 不存在: {src}')

    # 读程序 xlsx 列定义 (默认表头在 row 1)
    src_wb = load_workbook(str(src), data_only=True)
    src_ws = src_wb.active
    prog_header_row = comparison.prog_xlsx_row  # 1-based
    prog_col_count = src_ws.max_column
    prog_headers = [
        (str(src_ws.cell(row=prog_header_row, column=c).value or '').strip()
         or f'col{c}')
        for c in range(1, prog_col_count + 1)
    ]
    # 先收集程序侧的逐行原始数据 (后用)
    prog_rows_raw: list[list] = []
    prog_row_xlsx_ids: list[int] = []  # 每行在源 xlsx 中的实际行号
    for r in range(prog_header_row + 1, src_ws.max_row + 1):
        # 与原 compare 逻辑一致: 过滤空/注释行 但仍占据一个数据行位
        no_val = src_ws.cell(row=r, column=1).value
        c_val = src_ws.cell(row=r, column=2).value
        if no_val is None:
            continue
        no_str = str(no_val).strip()
        if not no_str or no_str.startswith('No') \
                or no_str in ('STRUCTURE', '未归属', 'Notes', '注释'):
            continue
        if c_val is not None and 'NMR' in str(c_val):
            continue
        row_data = [src_ws.cell(row=r, column=c).value
                    for c in range(1, prog_col_count + 1)]
        prog_rows_raw.append(row_data)
        prog_row_xlsx_ids.append(r)
    src_wb.close()

    # 读标准 xlsx 列定义 + 逐行原始数据
    std_wb = load_workbook(str(std_table.file_path), data_only=True)
    std_ws = std_wb.active
    std_header_row = std_table.header_row_idx
    std_col_count = std_ws.max_column
    std_headers = [
        (str(std_ws.cell(row=std_header_row, column=c).value or '').strip()
         or f'col{c}')
        for c in range(1, std_col_count + 1)
    ]
    # 映射 std.xlsx_row → 原始行数据
    std_row_raw_by_xrow: dict[int, list] = {}
    for r in range(std_header_row + 1, std_ws.max_row + 1):
        no_val = std_ws.cell(row=r, column=1).value
        if no_val is None:
            continue
        no_str = str(no_val).strip()
        if not no_str or no_str.startswith('No') \
                or no_str in ('STRUCTURE', '未归属', 'Notes', '注释'):
            continue
        c_val = std_ws.cell(row=r, column=2).value
        if c_val is not None and ('NMR' in str(c_val) or '归属' in str(c_val)):
            continue
        std_row_raw_by_xrow[r] = [
            std_ws.cell(row=r, column=c).value
            for c in range(1, std_col_count + 1)
        ]
    std_wb.close()

    # 建立匹配映射
    # 1) std.no_label → (prog_row_raw, cd, hd, ok)
    std_label_to_match: dict[str, tuple] = {}
    for pr, sr, cd, hd, ok in comparison.matches:
        sl = str(sr.no_label).strip()
        if sl and sl not in ('', '.', 'None'):
            # 找 prog raw 行
            try:
                pidx = prog_row_xlsx_ids.index(pr['xlsx_row'])
                prog_raw = prog_rows_raw[pidx]
            except ValueError:
                prog_raw = None
            std_label_to_match[sl] = (prog_raw, cd, hd, ok)

    # 2) prog.xlsx_row → (std_row_raw, cd, hd, ok)
    prog_xlsx_to_match: dict[int, tuple] = {}
    for pr, sr, cd, hd, ok in comparison.matches:
        std_raw = std_row_raw_by_xrow.get(sr.xlsx_row)
        prog_xlsx_to_match[pr['xlsx_row']] = (std_raw, cd, hd, ok)

    thin = Side(style='thin', color='CCCCCC')
    border = Border(top=thin, bottom=thin, left=thin, right=thin)
    bold = Font(bold=True)
    sep_marker = '拼接'
    diff_cols = ['ΔdC (ppm)', 'ΔdH (ppm)', '类型一致']

    wb = Workbook()

    # ===== Sheet 1: 按标准表为准 =====
    ws1 = wb.active
    ws1.title = '按标准表为准'
    headers1 = std_headers + [sep_marker] + prog_headers + diff_cols
    for c, h in enumerate(headers1, 1):
        cell = ws1.cell(row=1, column=c, value=h)
        cell.font = bold
        cell.border = border
    # 以标准表为行基
    for i, sr in enumerate(std_table.rows):
        row_num = i + 2
        # 左侧: 标准列 (从原 xlsx 复制)
        std_raw = std_row_raw_by_xrow.get(sr.xlsx_row)
        if std_raw is not None:
            for c, v in enumerate(std_raw, 1):
                ws1.cell(row=row_num, column=c, value=v)
        # 中间分隔
        sep_col = len(std_headers) + 1
        ws1.cell(row=row_num, column=sep_col, value='|')
        # 右侧: 程序侧 (matched) 或 placeholder
        match = std_label_to_match.get(str(sr.no_label).strip())
        right_offset = sep_col
        if match is not None:
            prog_raw, cd, hd, ok = match
            if prog_raw is not None:
                for c_off, v in enumerate(prog_raw):
                    ws1.cell(row=row_num, column=right_offset + 1 + c_off, value=v)
            # 末尾 Δ列
            dc_col = right_offset + len(prog_headers) + 1
            dh_col = dc_col + 1
            ok_col = dh_col + 1
            ws1.cell(row=row_num, column=dc_col,
                     value=f'{cd:.2f}' if cd is not None else '')
            ws1.cell(row=row_num, column=dh_col,
                     value=f'{hd:.3f}' if hd is not None else '')
            ws1.cell(row=row_num, column=ok_col,
                     value=('OK' if ok else 'DIFF'))
        else:
            for c_off in range(len(prog_headers) + len(diff_cols)):
                ws1.cell(row=row_num, column=right_offset + 1 + c_off,
                         value=placeholder)

    # ===== Sheet 2: 以程序表为准 =====
    ws2 = wb.create_sheet('按程序表为准')
    headers2 = prog_headers + [sep_marker] + std_headers + diff_cols
    for c, h in enumerate(headers2, 1):
        cell = ws2.cell(row=1, column=c, value=h)
        cell.font = bold
        cell.border = border
    for i, prog_raw in enumerate(prog_rows_raw):
        row_num = i + 2
        # 左侧: 程序列 (直接复制)
        for c, v in enumerate(prog_raw, 1):
            ws2.cell(row=row_num, column=c, value=v)
        # 中间分隔
        sep_col = len(prog_headers) + 1
        ws2.cell(row=row_num, column=sep_col, value='|')
        # 右侧: 标准侧
        prog_xrow = prog_row_xlsx_ids[i]
        match = prog_xlsx_to_match.get(prog_xrow)
        right_offset = sep_col
        if match is not None:
            std_raw, cd, hd, ok = match
            if std_raw is not None:
                for c_off, v in enumerate(std_raw):
                    ws2.cell(row=row_num, column=right_offset + 1 + c_off, value=v)
            dc_col = right_offset + len(std_headers) + 1
            dh_col = dc_col + 1
            ok_col = dh_col + 1
            ws2.cell(row=row_num, column=dc_col,
                     value=f'{cd:.2f}' if cd is not None else '')
            ws2.cell(row=row_num, column=dh_col,
                     value=f'{hd:.3f}' if hd is not None else '')
            ws2.cell(row=row_num, column=ok_col,
                     value=('OK' if ok else 'DIFF'))
        else:
            for c_off in range(len(std_headers) + len(diff_cols)):
                ws2.cell(row=row_num, column=right_offset + 1 + c_off,
                         value=placeholder)

    # 自动列宽
    for ws in (ws1, ws2):
        for c in range(1, ws.max_column + 1):
            col_letter = ws.cell(row=1, column=c).column_letter
            if c <= ws.max_column:
                ws.column_dimensions[col_letter].width = 14

    wb.save(str(out_path))


def write_comparison_report(
    out_path: str | Path,
    prog_xlsx_path: str | Path,
    std_table: StandardAssignmentTable,
    comparison: ComparisonResult,
    c_tol: float = 2.5,
    h_tol: float = 1.0,
    prog_meta_path: str | Path | None = None,
) -> None:
    """把对比结果写入 UTF-8 文本报告.

    Args:
        out_path: 输出报告路径
        prog_xlsx_path: 程序输出 xlsx 路径 (仅用于报告头打印)
        std_table: 标准归属表实例
        comparison: compare_program_to_standard 返回值
        c_tol: 13C 匹配容差 (仅报告打印用)
        h_tol: 1H 匹配容差 (仅报告打印用)
        prog_meta_path: 可选 meta.txt 路径. 提供时重新读取 meta 信息放入报告顶部.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # 如果 caller 没传 prog_meta_path, 但 comparison.prog_meta 已有数据, 优先用后者
    prog_meta = dict(comparison.prog_meta or {})
    if prog_meta_path is not None and not prog_meta:
        try:
            prog_meta = read_assignment_meta(prog_meta_path)
        except Exception:
            pass

    lines: list[str] = []
    lines.append('=' * 80)
    lines.append('  归属表化合物结构一致性对比报告 (已知化合物流程)')
    lines.append('=' * 80)
    lines.append(f'程序输出: {prog_xlsx_path}')
    lines.append(f'标准归属: {std_table.file_path}')
    if prog_meta_path is not None:
        lines.append(f'程序 meta: {prog_meta_path}')
    lines.append(f'生成时间: {__import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    lines.append('')

    # 程序输出元信息 (从 meta.txt 读, 与化学家 Excel 分离存放)
    if prog_meta:
        lines.append('-' * 80)
        lines.append('【1】程序输出元信息 (meta.txt)')
        lines.append('-' * 80)
        for key, val in prog_meta.items():
            lines.append(f'  {key}: {val}')
        lines.append('')
        section2_idx = 2
    else:
        section2_idx = 1

    # 化合物结构信息
    sec_idx = section2_idx
    lines.append('-' * 80)
    lines.append(f'【{sec_idx}】化合物结构信息')
    lines.append('-' * 80)
    lines.append(f'  真实结构 SMILES (标准归属): {std_table.known_smiles or "(未提供)"}')
    lines.append(f'  程序输出 SMILES: {comparison.prog_smiles or "(未包含)"}')
    if comparison.smiles_equal is True:
        lines.append('  RDKit canonical SMILES 判等: 一致 ✓')
    elif comparison.smiles_equal is False:
        lines.append('  RDKit canonical SMILES 判等: 不一致 ✗')
    else:
        lines.append('  RDKit canonical SMILES 判等: (未启用)')
    lines.append('')

    # 联合匹配
    sec_idx += 1
    lines.append('-' * 80)
    lines.append(f'【{sec_idx}】逐行联合匹配 (13C ±{c_tol} ppm, 1H ±{h_tol} ppm, 碳类型)')
    lines.append('-' * 80)
    lines.append(f'  程序输出: {comparison.n_prog} 行 ({comparison.n_prog_c} 个 13C, '
                 f'{comparison.n_prog_h} 个 1H)')
    lines.append(f'  标准归属: {comparison.n_std} 行 ({comparison.n_std_c} 个 13C, '
                 f'{comparison.n_std_h} 个 1H)')
    lines.append(f'  成功匹配: {len(comparison.matches)} 对')
    lines.append(f'  碳类型不一致: {len(comparison.type_mismatches)}/{len(comparison.matches)}')
    lines.append('')
    if comparison.matches:
        lines.append(f"  {'程序No':>12} {'程序13C':>8} {'类型':>4} {'程序1H':>7}"
                     f" | {'标准No':>12} {'标准13C':>8} {'类型':>4} {'标准1H':>7}"
                     f" | {'dC':>5} {'dH':>5} {'T':>4}")
        lines.append('  ' + '-' * 110)
        # 按 13C 偏差降序
        for pr, sr, cd, hd, ok in sorted(comparison.matches, key=lambda x: -(x[2] or 0)):
            c1 = f'{pr["c_ppm"]:.1f}' if pr['c_ppm'] is not None else '-'
            c2 = f'{sr.c_ppm:.1f}' if sr.c_ppm is not None else '-'
            h1 = f'{pr["h_ppm"]:.2f}' if pr['h_ppm'] is not None else '-'
            h2 = f'{sr.h_ppm:.2f}' if sr.h_ppm is not None else '-'
            cd_s = f'{cd:.2f}' if cd is not None else '-'
            hd_s = f'{hd:.3f}' if hd is not None else '-'
            flag = 'OK' if ok else 'DIFF'
            lines.append(f"  {pr['no']!s:>12} {c1:>8} {pr['type']:>4} {h1:>7}"
                         f" | {sr.no_label!s:>12} {c2:>8} {sr.carbon_type:>4} {h2:>7}"
                         f" | {cd_s:>5} {hd_s:>5} {flag:>4}")
    lines.append('')

    if comparison.unmatched_program:
        lines.append(f'  --- 程序输出中未匹配 ({len(comparison.unmatched_program)} 行) ---')
        for pr in comparison.unmatched_program:
            c = f"13C={pr['c_ppm']:.1f}" if pr['c_ppm'] is not None else '13C=-'
            h = f"1H={pr['h_ppm']:.2f}" if pr['h_ppm'] is not None else '1H=-'
            lines.append(f"    {pr['no']!s:>12}  {c}  {pr['type']:>4}  {h}")
        lines.append('')

    if comparison.unmatched_standard:
        lines.append(f'  --- 标准归属中未匹配 ({len(comparison.unmatched_standard)} 行) ---')
        for sr in comparison.unmatched_standard:
            c = f"13C={sr.c_ppm:.1f}" if sr.c_ppm is not None else '13C=-'
            h = f"1H={sr.h_ppm:.2f}" if sr.h_ppm is not None else '1H=-'
            lines.append(f"    {sr.no_label!s:>12}  {c}  {sr.carbon_type:>4}  {h}")
        lines.append('')

    # 综合结论
    sec_idx += 1
    lines.append('-' * 80)
    lines.append(f'【{sec_idx}】综合结论')
    lines.append('-' * 80)
    max_rows = max(comparison.n_prog, comparison.n_std)
    lines.append(f'  总行匹配率: {len(comparison.matches)}/{max_rows} = {comparison.match_rate:.1f}%')
    lines.append(f'  13C 匹配对数: {comparison.c13_matched} '
                 f'(程序 {comparison.n_prog_c} vs 标准 {comparison.n_std_c})')
    lines.append(f'  碳类型不一致: {len(comparison.type_mismatches)}/{len(comparison.matches)}')
    lines.append('')
    lines.append(f'  判定结果: {comparison.verdict}')
    lines.append('')
    lines.append('=' * 80)

    out.write_text('\n'.join(lines) + '\n', encoding='utf-8')