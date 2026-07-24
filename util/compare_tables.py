# -*- coding: utf-8 -*-
"""归属表对比脚本: 判断两个归属表代表的化合物结构是否一致.

核心逻辑:
  1. 以 (13C, 1H, 碳类型) 联合约束做 Hungarian 最优一对一匹配 (默认)
  2. 可选 --smiles 启用 RDKit SMILES 结构判等

默认仅比较 NMR 归属表, 不进行化合物结构判等.
行序可以不同, 列数/元信息行数/表头位置均可自适应.

默认行为 (变更):
  在程序输出的归属表 xlsx 同目录下另存为两个文件:
    - assignment_table_matched.xlsx  (复制原表 + 新增 "化学家No." 列, 未匹配填 "—")
    - comparison_report.txt          (文本报告, 含综合结论与匹配明细)

用法:
    python compare_tables.py
    python compare_tables.py --prog /path/to/output.xlsx
    python compare_tables.py --ref /path/to/reference.xlsx
    python compare_tables.py --smiles    # 额外启用 SMILES 结构判等
    python compare_tables.py --no-write  # 仅打印到 stdout, 不写文件
    python compare_tables.py --in-place  # 原地修改原 xlsx (不开新文件)
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')

# ============ 配置 ============
# 程序输出 (默认在脚本同目录, 可用 --prog 覆盖)
FILE_PROGRAM = 'assignment_table.xlsx'
# 真实结构归属表 (化学家 MestReNova 解谱结果)
FILE_REFERENCE = r'E:\NMR2Structure\真实结构\clxh100\clxh100数据归属.xlsx'
# 真实结构 SMILES (由 cdxml_to_smiles.py 从 ChemDraw 生成)
FILE_REF_SMILES = r'E:\NMR2Structure\真实结构\clxh100\smiles.txt'

C13_TOLERANCE = 0.2   # 13C 匹配容差 (ppm)
H1_TOLERANCE = 0.05   # 1H 匹配容差 (ppm)


# ============ Excel 解析 ============

def _find_header_row(raw: pd.DataFrame) -> int:
    """在原始 DataFrame 中定位表头行 (含 'No.' 的行)."""
    for i in range(min(30, len(raw))):
        for j in range(min(4, raw.shape[1])):
            val = raw.iat[i, j]
            if pd.notna(val) and str(val).strip().startswith('No'):
                return i
    raise ValueError('未找到表头行 (含 "No." 的行), 请检查文件格式')


def _build_col_map(header_row: pd.Series) -> dict:
    """根据表头内容建立列名映射, 自适应列数和列序.

    注意: '1H-1H COSY' 同时包含 '1H' 和 'COSY', 必须先判 COSY 再判 1H.
    """
    col_map = {}
    for j, val in enumerate(header_row):
        if pd.isna(val):
            continue
        s = str(val).strip()
        if s.startswith('No'):
            col_map['No'] = j
        elif '13C' in s:
            col_map['13C'] = j
        elif s in ('Types', 'Type'):
            col_map['Types'] = j
        elif 'COSY' in s:
            col_map['COSY'] = j
        elif 'HMBC' in s:
            col_map['HMBC'] = j
        elif 'ROESY' in s or 'NOESY' in s:
            col_map['ROESY'] = j
        elif 'Note' in s:
            col_map['Notes'] = j
        elif '1H' in s and '1H' not in col_map:
            col_map['1H'] = j
    return col_map


def _extract_smiles_above_header(raw: pd.DataFrame, header_idx: int) -> str | None:
    """从表头上方的元信息区提取 SMILES."""
    for i in range(header_idx):
        for j in range(min(4, raw.shape[1])):
            val = raw.iat[i, j]
            if pd.notna(val) and 'SMILES' in str(val).upper():
                # SMILES 值可能在同一行的后续列
                for k in range(j + 1, raw.shape[1]):
                    v = raw.iat[i, k]
                    if pd.notna(v) and len(str(v).strip()) > 5:
                        return str(v).strip()
    return None


def _extract_1h_shift(h_str) -> float | None:
    """从 1H NMR 字符串提取首个化学位移值.

    支持: '3.82, m' / '1.16, m (2H)' / '4.54 d (7.7)' / 纯数值
    """
    if pd.isna(h_str):
        return None
    s = str(h_str).strip()
    m = re.match(r'([\d.]+)', s)
    return float(m.group(1)) if m else None


def _normalize_type(t) -> str:
    """规范化碳类型字符串: Cq→C, 去空格, 大写."""
    if pd.isna(t) or str(t).strip() == '':
        return ''
    s = str(t).strip()
    if s in ('Cq', 'CQ', 'q', 'Q', 'C'):
        return 'C'
    return s.upper()


def load_table(filepath: str | Path) -> tuple[pd.DataFrame, str | None, int]:
    """加载归属表 xlsx, 返回 (数据 DataFrame, SMILES 或 None, 表头所在行号).

    自适应: 表头位置 / 列数 / 列序 / 元信息行数.

    DataFrame 中每行附带 '_xlsx_row' 列, 表示该行在原 xlsx 中的实际行号
    (从 1 开始, 用于后续写回原文件时定位).
    """
    raw = pd.read_excel(filepath, header=None)
    header_idx = _find_header_row(raw)
    col_map = _build_col_map(raw.iloc[header_idx])
    smiles = _extract_smiles_above_header(raw, header_idx)

    if 'No' not in col_map or '13C' not in col_map:
        raise ValueError(f'表头缺少必要列 (No/13C), 实际识别到: {list(col_map.keys())}')

    # 取表头下方的数据行 (记录原 xlsx 行号 = header_idx + 1 + 局部行号, header_idx 本身 0-based, xlsx 行号从 1 开始)
    data = raw.iloc[header_idx + 1:].reset_index(drop=True).copy()
    data['_xlsx_row'] = data.index + header_idx + 2  # +2 = +1 表头行偏移, +1 1-based

    # 按映射提取列
    result = pd.DataFrame()
    result['_xlsx_row'] = data['_xlsx_row']
    result['No'] = data.iloc[:, col_map['No']]
    result['13C_raw'] = data.iloc[:, col_map['13C']]
    result['Types'] = (data.iloc[:, col_map['Types']]
                       if 'Types' in col_map else pd.Series([''] * len(data)))
    result['1H_raw'] = (data.iloc[:, col_map['1H']]
                        if '1H' in col_map else pd.Series([''] * len(data)))

    # 过滤: 去掉空行 / 表尾注释行
    result = result.dropna(subset=['No'])
    result = result[~result['No'].astype(str).str.startswith('No')]
    result = result[~result['No'].astype(str).isin(['STRUCTURE', '未归属', ''])]
    # 去掉 13C 列为非数值且非空的行 (如注释文字)
    result = result[~result['13C_raw'].astype(str).str.contains('NMR|ppm|归属', na=False)]

    # 数值化
    result['13C_num'] = pd.to_numeric(result['13C_raw'], errors='coerce')
    result['1H_shift'] = result['1H_raw'].apply(_extract_1h_shift)
    result['Type_norm'] = result['Types'].apply(_normalize_type)

    return result.reset_index(drop=True), smiles, header_idx


# ============ 联合匹配 ============

def joint_match(df1: pd.DataFrame, df2: pd.DataFrame,
                c_tol: float, h_tol: float):
    """以 (13C, 1H, 碳类型) 联合约束做 Hungarian 最优一对一匹配.

    规则:
      - 两行都有 13C: 13C 差 <= c_tol 为必要条件
      - 两行都有 1H: 1H 差 <= h_tol 为必要条件
      - 碳类型不一致: 代价加惩罚 (但不完全禁止, 因 DEPT 推断可能有误)
      - 一行有 13C 另一行没有: 不可匹配

    Returns:
        matches: [(row1, row2, c_diff, h_diff, type_ok)]
        unmatched_1: DataFrame
        unmatched_2: DataFrame
    """
    from scipy.optimize import linear_sum_assignment

    n1, n2 = len(df1), len(df2)
    if n1 == 0 or n2 == 0:
        return [], df1, df2

    BIG = 1e6
    TYPE_PENALTY = 0.5  # 类型不一致的额外代价 (归一化单位)

    cost = np.full((n1, n2), BIG)

    for i in range(n1):
        c1 = df1.at[i, '13C_num']
        h1 = df1.at[i, '1H_shift']
        t1 = df1.at[i, 'Type_norm']
        for j in range(n2):
            c2 = df2.at[j, '13C_num']
            h2 = df2.at[j, '1H_shift']
            t2 = df2.at[j, 'Type_norm']

            # 13C: 双方都有值时必须满足容差
            if pd.notna(c1) and pd.notna(c2):
                c_diff = abs(c1 - c2)
                if c_diff > c_tol:
                    continue
            elif pd.notna(c1) != pd.notna(c2):
                # 一方有 13C 另一方没有 → 不匹配
                continue
            else:
                c_diff = 0.0  # 双方都无 13C (纯 1H 行)

            # 1H: 双方都有值时必须满足容差
            if pd.notna(h1) and pd.notna(h2):
                h_diff = abs(h1 - h2)
                if h_diff > h_tol:
                    continue
            else:
                h_diff = 0.0

            # 归一化代价
            c_norm = c_diff / c_tol if c_tol > 0 else 0
            h_norm = h_diff / h_tol if h_tol > 0 else 0
            d = c_norm + h_norm

            # 碳类型惩罚
            if t1 and t2 and t1 != t2:
                d += TYPE_PENALTY

            cost[i, j] = d

    row_ind, col_ind = linear_sum_assignment(cost)

    matches = []
    matched_1 = set()
    matched_2 = set()
    for r, c in zip(row_ind, col_ind):
        if cost[r, c] >= BIG:
            continue
        c1 = df1.at[r, '13C_num']
        c2 = df2.at[c, '13C_num']
        h1 = df1.at[r, '1H_shift']
        h2 = df2.at[c, '1H_shift']
        c_diff = abs(c1 - c2) if pd.notna(c1) and pd.notna(c2) else None
        h_diff = abs(h1 - h2) if pd.notna(h1) and pd.notna(h2) else None
        t1 = df1.at[r, 'Type_norm']
        t2 = df2.at[c, 'Type_norm']
        type_ok = (t1 == t2) or (not t1) or (not t2)
        matches.append((df1.iloc[r], df2.iloc[c], c_diff, h_diff, type_ok))
        matched_1.add(r)
        matched_2.add(c)

    unmatched_1 = df1[~df1.index.isin(matched_1)]
    unmatched_2 = df2[~df2.index.isin(matched_2)]
    return matches, unmatched_1, unmatched_2


# ============ 输出 ============

def print_section(title: str):
    print()
    print('=' * 80)
    print(f'  {title}')
    print('=' * 80)


# ============ 写回 xlsx 与报告 ============

def write_matched_xlsx(src_path: Path, dst_path: Path,
                       header_row_1based: int,
                       prog_to_ref_no: dict[int, str],
                       placeholder: str = '—') -> None:
    """在原程序输出 xlsx 基础上加一列 "化学家No.", 另存为 dst_path.

    参数:
        src_path: 原 xlsx 路径
        dst_path: 输出 xlsx 路径
        header_row_1based: 原 xlsx 中表头所在行号 (1-based)
        prog_to_ref_no: 程序输出行号 (1-based) → 化学家归属表 No 值的映射.
                        未匹配的行不在 dict 中, 在写入时填 placeholder.
    """
    from openpyxl import load_workbook

    wb = load_workbook(src_path)
    ws = wb.active

    new_col = ws.max_column + 1
    # 表头行 (header_row_1based) 写列名
    ws.cell(row=header_row_1based, column=new_col, value='化学家No.')

    # 数据行: 从 header_row_1based + 1 开始, 写到 ws.max_row
    for r in range(header_row_1based + 1, ws.max_row + 1):
        # 只填充被 load_table 识别为数据行的行 (即 prog_to_ref_no 键或
        # 实际是 df_prog['_xlsx_row'] 的取值集合). 未匹配行不写入映射
        # 时也写入 placeholder, 保持行序对齐.
        if r in prog_to_ref_no:
            ws.cell(row=r, column=new_col, value=prog_to_ref_no[r])
        else:
            # 判断是否为空行/注释行 (与 load_table 过滤逻辑一致):
            # No. 列为空 或 以 No 开头 或 包含 NMR/ppm/归属 等则跳过
            no_val = ws.cell(row=r, column=1).value
            c13_val = ws.cell(row=r, column=2).value
            if no_val is None or str(no_val).strip().startswith('No') \
                    or (c13_val is not None and 'NMR' in str(c13_val)):
                # 元信息行 / 注释行: 不填任何东西
                continue
            ws.cell(row=r, column=new_col, value=placeholder)

    wb.save(dst_path)


def write_joint_xlsx(prog_path: Path, ref_path: Path,
                     df_prog: pd.DataFrame, df_ref: pd.DataFrame,
                     matches: list,
                     out_path: Path,
                     placeholder: str = '—') -> None:
    """生成双 sheet 的左右拼接 xlsx. 取代旧的 write_matched_xlsx(单 sheet).

    Sheet 1 '按标准表为准': 每行 = df_ref 的一行, 匹配上的 df_prog 行接右侧.
    Sheet 2 '按程序表为准': 每行 = df_prog 的一行, 匹配上的 df_ref 行接右侧.

    参数:
        prog_path: 程序输出 xlsx (提供列定义 / 原始数据)
        ref_path:  标准/参考 xlsx (提供列定义 / 原始数据)
        df_prog:   load_table(prog_path) 的结果 DataFrame (供快速定位原 xlsx_row)
        df_ref:    load_table(ref_path) 的结果 DataFrame
        matches:   joint_match() 返回的 matches 列表
        out_path:  输出 xlsx 路径
        placeholder: 未匹配单元格的占位符
    """
    from openpyxl import Workbook, load_workbook
    from openpyxl.styles import Font, Border, Side

    thin = Side(style='thin', color='CCCCCC')
    border = Border(top=thin, bottom=thin, left=thin, right=thin)
    bold = Font(bold=True)
    sep_marker = '拼接'
    diff_cols = ['ΔdC (ppm)', 'ΔdH (ppm)', '类型一致']

    # 读取程序 xlsx (源)
    src_wb = load_workbook(str(prog_path), data_only=True)
    src_ws = src_wb.active
    prog_header_row = _find_header_row(pd.read_excel(prog_path, header=None)) + 1
    prog_col_count = src_ws.max_column
    prog_headers = [
        (str(src_ws.cell(row=prog_header_row, column=c).value or '').strip()
         or f'col{c}')
        for c in range(1, prog_col_count + 1)
    ]
    # 程序逐行原始数据 (仅保留 load_table 识别的数据行)
    prog_row_xlsx_ids: list[int] = []
    prog_rows_raw: list[list] = []
    for xrow in df_prog['_xlsx_row'].tolist():
        row = [src_ws.cell(row=int(xrow), column=c).value
               for c in range(1, prog_col_count + 1)]
        prog_row_xlsx_ids.append(int(xrow))
        prog_rows_raw.append(row)
    src_wb.close()

    # 读取标准 xlsx (源)
    ref_wb = load_workbook(str(ref_path), data_only=True)
    ref_ws = ref_wb.active
    ref_header_row = _find_header_row(pd.read_excel(ref_path, header=None)) + 1
    ref_col_count = ref_ws.max_column
    ref_headers = [
        (str(ref_ws.cell(row=ref_header_row, column=c).value or '').strip()
         or f'col{c}')
        for c in range(1, ref_col_count + 1)
    ]
    ref_row_raw_by_xrow: dict[int, list] = {}
    for xrow in df_ref['_xlsx_row'].tolist():
        ref_row_raw_by_xrow[int(xrow)] = [
            ref_ws.cell(row=int(xrow), column=c).value
            for c in range(1, ref_col_count + 1)
        ]
    ref_wb.close()

    # 映射: (1) ref No. → (prog row raw, cd, hd, ok) (2) prog xrow → (ref row raw, cd, hd, ok)
    # matches 元素: (r1, r2, c_diff, h_diff, type_ok) 其中 r1=prog, r2=ref
    ref_no_to_match: dict[str, tuple] = {}
    for r1, r2, cd, hd, ok in matches:
        ref_no = str(r2['No']).strip()
        if ref_no in ('', '.', 'None'):
            continue
        try:
            pidx = prog_row_xlsx_ids.index(int(r1['_xlsx_row']))
            prog_raw = prog_rows_raw[pidx]
        except ValueError:
            prog_raw = None
        ref_no_to_match[ref_no] = (prog_raw, cd, hd, ok)

    prog_xrow_to_match: dict[int, tuple] = {}
    for r1, r2, cd, hd, ok in matches:
        ref_raw = ref_row_raw_by_xrow.get(int(r2['_xlsx_row']))
        prog_xrow_to_match[int(r1['_xlsx_row'])] = (ref_raw, cd, hd, ok)

    wb = Workbook()

    # ===== Sheet 1: 按标准表为准 =====
    ws1 = wb.active
    ws1.title = '按标准表为准'
    headers1 = ref_headers + [sep_marker] + prog_headers + diff_cols
    for c, h in enumerate(headers1, 1):
        cell = ws1.cell(row=1, column=c, value=h)
        cell.font = bold
        cell.border = border
    # 逐行: 标准表为基
    for i, (_, row) in enumerate(df_ref.iterrows()):
        row_num = i + 2
        ref_xrow = int(row['_xlsx_row'])
        ref_raw = ref_row_raw_by_xrow.get(ref_xrow)
        if ref_raw is not None:
            for c, v in enumerate(ref_raw, 1):
                ws1.cell(row=row_num, column=c, value=v)
        sep_col = len(ref_headers) + 1
        ws1.cell(row=row_num, column=sep_col, value='|')
        match = ref_no_to_match.get(str(row['No']).strip())
        right_offset = sep_col
        if match is not None:
            prog_raw, cd, hd, ok = match
            if prog_raw is not None:
                for c_off, v in enumerate(prog_raw):
                    ws1.cell(row=row_num, column=right_offset + 1 + c_off, value=v)
            dc_col = right_offset + len(prog_headers) + 1
            ws1.cell(row=row_num, column=dc_col,
                     value=f'{cd:.2f}' if cd is not None else '')
            ws1.cell(row=row_num, column=dc_col + 1,
                     value=f'{hd:.3f}' if hd is not None else '')
            ws1.cell(row=row_num, column=dc_col + 2,
                     value='OK' if ok else 'DIFF')
        else:
            for c_off in range(len(prog_headers) + len(diff_cols)):
                ws1.cell(row=row_num, column=right_offset + 1 + c_off,
                         value=placeholder)

    # ===== Sheet 2: 按程序表为准 =====
    ws2 = wb.create_sheet('按程序表为准')
    headers2 = prog_headers + [sep_marker] + ref_headers + diff_cols
    for c, h in enumerate(headers2, 1):
        cell = ws2.cell(row=1, column=c, value=h)
        cell.font = bold
        cell.border = border
    for i, (_, row) in enumerate(df_prog.iterrows()):
        row_num = i + 2
        prog_xrow = int(row['_xlsx_row'])
        idx = prog_row_xlsx_ids.index(prog_xrow) if prog_xrow in prog_row_xlsx_ids else i
        prog_raw = prog_rows_raw[idx]
        for c, v in enumerate(prog_raw, 1):
            ws2.cell(row=row_num, column=c, value=v)
        sep_col = len(prog_headers) + 1
        ws2.cell(row=row_num, column=sep_col, value='|')
        match = prog_xrow_to_match.get(prog_xrow)
        right_offset = sep_col
        if match is not None:
            ref_raw, cd, hd, ok = match
            if ref_raw is not None:
                for c_off, v in enumerate(ref_raw):
                    ws2.cell(row=row_num, column=right_offset + 1 + c_off, value=v)
            dc_col = right_offset + len(ref_headers) + 1
            ws2.cell(row=row_num, column=dc_col,
                     value=f'{cd:.2f}' if cd is not None else '')
            ws2.cell(row=row_num, column=dc_col + 1,
                     value=f'{hd:.3f}' if hd is not None else '')
            ws2.cell(row=row_num, column=dc_col + 2,
                     value='OK' if ok else 'DIFF')
        else:
            for c_off in range(len(ref_headers) + len(diff_cols)):
                ws2.cell(row=row_num, column=right_offset + 1 + c_off,
                         value=placeholder)

    for ws in (ws1, ws2):
        for c in range(1, ws.max_column + 1):
            col_letter = ws.cell(row=1, column=c).column_letter
            ws.column_dimensions[col_letter].width = 14

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(out_path))


def write_report(dst_path: Path, lines: list[str]) -> None:
    """把报告行写入文件 (UTF-8)."""
    dst_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def build_report_lines(file_prog: Path, file_ref: Path,
                       df_prog: pd.DataFrame, df_ref: pd.DataFrame,
                       smiles_prog: str | None, smiles_ref: str | None,
                       known_smiles: str | None,
                       use_smiles: bool,
                       matches: list, unmatched_prog: pd.DataFrame,
                       unmatched_ref: pd.DataFrame,
                       type_mismatches: list) -> tuple[list[str], str, str, dict]:
    """构造报告文本行, 同时返回判定与匹配率信息.

    Returns: (lines, verdict, match_rate_str, prog_to_ref_no)
    """
    n_prog_c = df_prog['13C_num'].notna().sum()
    n_ref_c = df_ref['13C_num'].notna().sum()
    n_prog_h = df_prog['1H_shift'].notna().sum()
    n_ref_h = df_ref['1H_shift'].notna().sum()

    max_rows = max(len(df_prog), len(df_ref))
    match_rate = len(matches) / max_rows * 100 if max_rows > 0 else 0
    c13_matched_with_c = sum(1 for _, _, cd, _, _ in matches if cd is not None)

    # 构造 prog_xlsx_row -> ref No. 映射
    # 过滤 No='.' 等异常空行 (这些不是真实归属编号)
    prog_to_ref_no: dict[int, str] = {}
    for r1, r2, cd, hd, ok in matches:
        prog_xlsx_row = int(r1['_xlsx_row'])
        ref_no_str = str(r2['No']).strip()
        if ref_no_str in ('', '.', 'None'):
            # 异常空行 → 当作未匹配
            continue
        prog_to_ref_no[prog_xlsx_row] = ref_no_str

    lines: list[str] = []
    lines.append('=' * 80)
    lines.append('  归属表化合物结构一致性对比报告')
    lines.append('=' * 80)
    lines.append(f'程序输出: {file_prog}')
    lines.append(f'真实结构: {file_ref}')
    lines.append(f'生成时间: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}')
    lines.append('')

    # ---- 【1】化合物结构信息 ----
    lines.append('-' * 80)
    lines.append('【1】化合物结构信息')
    lines.append('-' * 80)
    if known_smiles:
        lines.append(f'  真实结构 SMILES: {known_smiles}')
    else:
        lines.append('  真实结构 SMILES: (未提供)')
    if smiles_prog:
        lines.append(f'  程序输出 SMILES: {smiles_prog}')
    else:
        lines.append('  程序输出 SMILES: (未包含)')
    if use_smiles:
        from rdkit import Chem
        can_eq = None
        try:
            if known_smiles and smiles_prog:
                mk = Chem.MolFromSmiles(known_smiles)
                mp = Chem.MolFromSmiles(smiles_prog)
                if mk and mp:
                    can_eq = (Chem.MolToSmiles(mk) == Chem.MolToSmiles(mp))
                    lines.append(f'  RDKit canonical SMILES: 一致={can_eq}')
                    if not can_eq:
                        lines.append(f'  InChIKey 真实: {Chem.inchi.MolToInchiKey(mk)}')
                        lines.append(f'  InChIKey 程序: {Chem.inchi.MolToInchiKey(mp)}')
        except ImportError:
            lines.append('  (RDKit 未安装, 跳过结构判等)')
    else:
        lines.append('  SMILES 结构判等: 未启用 (传 --smiles 启用)')
    lines.append('')

    # ---- 【2】逐行联合匹配 ----
    lines.append('-' * 80)
    lines.append(f'【2】逐行联合匹配 (13C ±{C13_TOLERANCE} ppm, 1H ±{H1_TOLERANCE} ppm, 碳类型)')
    lines.append('-' * 80)
    lines.append(f'  程序输出: {len(df_prog)} 行 ({n_prog_c} 个 13C, {n_prog_h} 个 1H)')
    lines.append(f'  真实结构: {len(df_ref)} 行 ({n_ref_c} 个 13C, {n_ref_h} 个 1H)')
    lines.append(f'  成功匹配: {len(matches)} 对')
    lines.append(f'  碳类型不一致: {len(type_mismatches)}/{len(matches)}')
    lines.append('')
    if matches:
        lines.append(f"  {'程序No':>10} {'程序13C':>8} {'程序Type':>8} {'程序1H':>7}"
                     f" | {'真实No':>10} {'真实13C':>8} {'真实Type':>8} {'真实1H':>7}"
                     f" | {'dC':>5} {'dH':>5} {'类型':>4}")
        lines.append('  ' + '-' * 100)
        for r1, r2, cd, hd, ok in sorted(matches, key=lambda x: -(x[2] or 0)):
            c1 = f"{r1['13C_num']:.1f}" if pd.notna(r1['13C_num']) else '-'
            c2 = f"{r2['13C_num']:.1f}" if pd.notna(r2['13C_num']) else '-'
            h1 = f"{r1['1H_shift']:.2f}" if pd.notna(r1['1H_shift']) else '-'
            h2 = f"{r2['1H_shift']:.2f}" if pd.notna(r2['1H_shift']) else '-'
            cd_s = f'{cd:.2f}' if cd is not None else '-'
            hd_s = f'{hd:.3f}' if hd is not None else '-'
            flag = 'OK' if ok else 'DIFF'
            lines.append(f"  {str(r1['No']):>10} {c1:>8} {r1['Type_norm']:>8} {h1:>7}"
                         f" | {str(r2['No']):>10} {c2:>8} {r2['Type_norm']:>8} {h2:>7}"
                         f" | {cd_s:>5} {hd_s:>5} {flag:>4}")
    lines.append('')

    if len(unmatched_prog) > 0:
        lines.append(f'  --- 程序输出中未匹配 ({len(unmatched_prog)} 行) ---')
        for _, row in unmatched_prog.iterrows():
            c = f"13C={row['13C_num']:.1f}" if pd.notna(row['13C_num']) else '13C=-'
            h = f"1H={row['1H_shift']:.2f}" if pd.notna(row['1H_shift']) else '1H=-'
            lines.append(f"    {str(row['No']):>10}  {c}  {row['Type_norm']:>4}  {h}")
        lines.append('')
    if len(unmatched_ref) > 0:
        lines.append(f'  --- 真实结构中未匹配 ({len(unmatched_ref)} 行) ---')
        for _, row in unmatched_ref.iterrows():
            c = f"13C={row['13C_num']:.1f}" if pd.notna(row['13C_num']) else '13C=-'
            h = f"1H={row['1H_shift']:.2f}" if pd.notna(row['1H_shift']) else '1H=-'
            lines.append(f"    {str(row['No']):>10}  {c}  {row['Type_norm']:>4}  {h}")
        lines.append('')

    # ---- 【3】综合结论 ----
    lines.append('-' * 80)
    lines.append('【3】综合结论')
    lines.append('-' * 80)
    lines.append(f'  总行匹配率: {len(matches)}/{max_rows} = {match_rate:.1f}%')
    lines.append(f'  13C 匹配对数: {c13_matched_with_c} (程序 {n_prog_c} vs 真实 {n_ref_c})')
    lines.append(f'  碳类型不一致: {len(type_mismatches)}/{len(matches)}')
    lines.append('')
    if match_rate >= 80 and len(type_mismatches) == 0:
        verdict = '一致 ✓ — 两个归属表代表同一化合物'
    elif match_rate >= 60:
        verdict = '部分一致 △ — 可能同一化合物但归属有差异'
    else:
        verdict = '不一致 ✗ — 两个归属表代表不同化合物'
    lines.append(f'  判定结果: {verdict}')
    lines.append('')
    lines.append('=' * 80)

    match_rate_str = f'{len(matches)}/{max_rows} = {match_rate:.1f}%'
    return lines, verdict, match_rate_str, prog_to_ref_no


def main():
    import argparse
    parser = argparse.ArgumentParser(description='归属表化合物结构一致性对比')
    parser.add_argument('--prog', default=None, help='程序输出 xlsx 路径')
    parser.add_argument('--ref', default=None, help='真实结构 xlsx 路径')
    parser.add_argument('--smiles', action='store_true',
                        help='启用 SMILES 结构判等 (默认仅做 NMR 归属表比较)')
    parser.add_argument('--no-write', action='store_true',
                        help='仅打印到 stdout, 不写文件')
    parser.add_argument('--in-place', action='store_true',
                        help='原地修改原 xlsx (不开新文件, 会覆盖原 assignment_table.xlsx)')
    args = parser.parse_args()

    base_dir = Path(__file__).parent
    file_prog = Path(args.prog) if args.prog else base_dir / FILE_PROGRAM
    file_ref = Path(args.ref) if args.ref else Path(FILE_REFERENCE)

    if not file_prog.exists():
        print(f'错误: 找不到程序输出文件 {file_prog}')
        return
    if not file_ref.exists():
        print(f'错误: 找不到真实结构文件 {file_ref}')
        return

    # 加载
    df_prog, smiles_prog, header_idx_prog = load_table(file_prog)
    df_ref, smiles_ref, _ = load_table(file_ref)

    print(f'程序输出: {file_prog.name} ({len(df_prog)} 行)')
    print(f'真实结构: {file_ref.name} ({len(df_ref)} 行)')

    # ---- 化合物结构信息 (可选) ----
    print_section('【1】化合物结构信息')

    # 参考 SMILES: 始终报告文件位置 (不默认做结构判等)
    ref_smiles_path = Path(FILE_REF_SMILES)
    known_smiles = None
    if ref_smiles_path.exists():
        known_smiles = ref_smiles_path.read_text(encoding='utf-8').strip()
        print(f'  真实结构 SMILES 文件: {ref_smiles_path}')
    else:
        print(f'  真实结构 SMILES 文件不存在: {ref_smiles_path} (需先运行 cdxml_to_smiles.py)')

    if smiles_prog:
        print(f'  程序输出 SMILES: {smiles_prog}')
    else:
        print('  程序输出未包含 SMILES 元信息')

    if not args.smiles:
        print('\n  [默认跳过 SMILES 结构判等]')
        print('  启用方法: python compare_tables.py --smiles')
    else:
        print('\n  [已启用 SMILES 结构判等]')
        # RDKit 结构判等: 程序输出 vs 真实结构
        try:
            from rdkit import Chem
            if known_smiles and smiles_prog:
                mol_known = Chem.MolFromSmiles(known_smiles)
                mol_prog = Chem.MolFromSmiles(smiles_prog)
                if mol_known and mol_prog:
                    can_known = Chem.MolToSmiles(mol_known)
                    can_prog = Chem.MolToSmiles(mol_prog)
                    print(f'  规范化 (真实): {can_known}')
                    print(f'  规范化 (程序): {can_prog}')
                    if can_known == can_prog:
                        print('  >>> 结构判等: 一致 ✓ (同一化合物)')
                    else:
                        # 进一步用 InChIKey 比较 (忽略立体)
                        ik_known = Chem.inchi.MolToInchiKey(mol_known)
                        ik_prog = Chem.inchi.MolToInchiKey(mol_prog)
                        if ik_known == ik_prog:
                            print('  >>> 结构判等: 一致 ✓ (InChIKey 相同, 仅立体差异)')
                        else:
                            print(f'  >>> 结构判等: 不一致 ✗ (不同化合物)')
                            print(f'      InChIKey 真实: {ik_known}')
                            print(f'      InChIKey 程序: {ik_prog}')
                else:
                    print('  (RDKit 无法解析其中一个 SMILES, 跳过结构判等)')
            elif known_smiles and not smiles_prog:
                print('  (程序输出无 SMILES, 无法判等)')
            elif not known_smiles:
                print('  (真实结构 SMILES 文件不存在, 无法判等)')
        except ImportError:
            print('  (未安装 RDKit, 跳过结构判等)')

    # ---- 联合匹配 ----
    print_section(f'【2】逐行联合匹配 (13C ±{C13_TOLERANCE} ppm, 1H ±{H1_TOLERANCE} ppm, 碳类型)')

    matches, unmatched_prog, unmatched_ref = joint_match(
        df_prog, df_ref, C13_TOLERANCE, H1_TOLERANCE)

    n_prog_c = df_prog['13C_num'].notna().sum()
    n_ref_c = df_ref['13C_num'].notna().sum()
    n_prog_h = df_prog['1H_shift'].notna().sum()
    n_ref_h = df_ref['1H_shift'].notna().sum()

    print(f'\n  程序输出: {len(df_prog)} 行 ({n_prog_c} 个 13C, {n_prog_h} 个 1H)')
    print(f'  真实结构: {len(df_ref)} 行 ({n_ref_c} 个 13C, {n_ref_h} 个 1H)')
    print(f'  成功匹配: {len(matches)} 对')

    type_mismatches = [(r1, r2, cd, hd) for r1, r2, cd, hd, ok in matches if not ok]
    print(f'  碳类型不一致: {len(type_mismatches)}/{len(matches)}')

    # 匹配明细 (按 13C 偏差降序)
    if matches:
        print(f"\n  {'程序No':>10} {'程序13C':>8} {'程序Type':>8} {'程序1H':>7}"
              f" | {'真实No':>10} {'真实13C':>8} {'真实Type':>8} {'真实1H':>7}"
              f" | {'dC':>5} {'dH':>5} {'类型':>4}")
        print('  ' + '-' * 100)
        for r1, r2, cd, hd, ok in sorted(
                matches, key=lambda x: -(x[2] or 0)):
            c1 = f"{r1['13C_num']:.1f}" if pd.notna(r1['13C_num']) else '-'
            c2 = f"{r2['13C_num']:.1f}" if pd.notna(r2['13C_num']) else '-'
            h1 = f"{r1['1H_shift']:.2f}" if pd.notna(r1['1H_shift']) else '-'
            h2 = f"{r2['1H_shift']:.2f}" if pd.notna(r2['1H_shift']) else '-'
            cd_s = f'{cd:.2f}' if cd is not None else '-'
            hd_s = f'{hd:.3f}' if hd is not None else '-'
            flag = 'OK' if ok else 'DIFF'
            print(f"  {str(r1['No']):>10} {c1:>8} {r1['Type_norm']:>8} {h1:>7}"
                  f" | {str(r2['No']):>10} {c2:>8} {r2['Type_norm']:>8} {h2:>7}"
                  f" | {cd_s:>5} {hd_s:>5} {flag:>4}")

    # 未匹配
    if len(unmatched_prog) > 0:
        print(f"\n  --- 程序输出中未匹配 ({len(unmatched_prog)} 行) ---")
        for _, row in unmatched_prog.iterrows():
            c = f"13C={row['13C_num']:.1f}" if pd.notna(row['13C_num']) else '13C=-'
            h = f"1H={row['1H_shift']:.2f}" if pd.notna(row['1H_shift']) else '1H=-'
            print(f"    {str(row['No']):>10}  {c}  {row['Type_norm']:>4}  {h}")

    if len(unmatched_ref) > 0:
        print(f"\n  --- 真实结构中未匹配 ({len(unmatched_ref)} 行) ---")
        for _, row in unmatched_ref.iterrows():
            c = f"13C={row['13C_num']:.1f}" if pd.notna(row['13C_num']) else '13C=-'
            h = f"1H={row['1H_shift']:.2f}" if pd.notna(row['1H_shift']) else '1H=-'
            print(f"    {str(row['No']):>10}  {c}  {row['Type_norm']:>4}  {h}")

    # ---- 综合结论 ----
    print_section('【3】综合结论')

    max_rows = max(len(df_prog), len(df_ref))
    match_rate = len(matches) / max_rows * 100 if max_rows > 0 else 0
    c13_both = df_prog['13C_num'].notna().sum() + df_ref['13C_num'].notna().sum()
    c13_matched_with_c = sum(1 for _, _, cd, _, _ in matches if cd is not None)

    print(f'\n  总行匹配率: {len(matches)}/{max_rows} = {match_rate:.1f}%')
    print(f'  13C 匹配对数: {c13_matched_with_c} (程序 {n_prog_c} vs 真实 {n_ref_c})')
    print(f'  碳类型不一致: {len(type_mismatches)}/{len(matches)}')
    print()

    if match_rate >= 80 and len(type_mismatches) == 0:
        verdict = '一致 ✓ — 两个归属表代表同一化合物'
    elif match_rate >= 60:
        verdict = '部分一致 △ — 可能同一化合物但归属有差异'
    else:
        verdict = '不一致 ✗ — 两个归属表代表不同化合物'

    print(f'  判定结果: {verdict}')
    print()

    # ---- 写文件 (默认行为) ----
    if args.no_write:
        return

    # 构造报告行 + prog_to_ref_no 映射 (复用 main 中已有结果)
    lines, _, _, prog_to_ref_no = build_report_lines(
        file_prog=file_prog, file_ref=file_ref,
        df_prog=df_prog, df_ref=df_ref,
        smiles_prog=smiles_prog, smiles_ref=smiles_ref,
        known_smiles=known_smiles, use_smiles=args.smiles,
        matches=matches, unmatched_prog=unmatched_prog,
        unmatched_ref=unmatched_ref, type_mismatches=type_mismatches)

    out_dir = file_prog.parent
    if args.in_place:
        dst_xlsx = file_prog
    else:
        dst_xlsx = out_dir / (file_prog.stem + '_matched.xlsx')
    report_path = out_dir / 'comparison_report.txt'

    try:
        write_joint_xlsx(
            prog_path=file_prog,
            ref_path=file_ref,
            df_prog=df_prog, df_ref=df_ref,
            matches=matches,
            out_path=dst_xlsx,
        )
        write_report(report_path, lines)
        print('=' * 80)
        print(f'  已写入拼接 xlsx (双sheet): {dst_xlsx}')
        print(f'    Sheet 1: 按标准表为准 (以标准为基, 行 {len(df_ref)} 行)')
        print(f'    Sheet 2: 按程序表为准 (以程序为基, 行 {len(df_prog)} 行)')
        print(f'  已写入报告:   {report_path}')
        print('=' * 80)
    except Exception as e:
        print(f'写入文件失败: {e}')


if __name__ == '__main__':
    main()
