# -*- coding: utf-8 -*-
import time
import numpy as np
import pandas as pd
from openpyxl import load_workbook

t0 = time.perf_counter()

# 不动原始"类型"列和 ppm 列，不删除任何行，只新增 type 列存放转换结果：
# 类型列中的 s/d/t/q 转换为 C/CH/CH2/CH3；OCH3 行（通过位置列中的 "OCH3" 标记识别）标记为 CH3。
# 用 openpyxl 在第一个 Sheet 上就地修改，原字体/列宽等格式原样保留；其余 Sheet 不动。
SRC = r'Y:\LXH数据处理\新化合物处理\程序解析\clxh81解析\clxh81参考.xlsx'

wb = load_workbook(SRC)
ws = wb.worksheets[0]

header = {str(c.value).strip().lower(): c.column for c in ws[1] if c.value is not None}

def find_col(candidates):
    for n in candidates:
        col = header.get(n.lower())
        if col is not None:
            return col
    raise KeyError(f"未找到列，候选: {candidates}，现有表头: {list(header)}")

pos_col = find_col(('position', 'pos.', 'pos', 'no.', 'no'))
type_col = find_col(('类型', 'type', 'types'))
ppm_col = find_col(('ppm', '13C NMR', 'δc'))
type_new_col = header.get('type', (max(header.values()) + 1))
ws.cell(row=1, column=type_new_col, value='type')

# 数据区读入 DataFrame，df 第 i 行对应工作表第 i+2 行（保持对齐，便于向量化计算）
data = pd.DataFrame(ws.iter_rows(min_row=2, values_only=True))

m = {'s': 'C', 'd': 'CH', 't': 'CH2', 'q': 'CH3'}
ts = data[type_col - 1].astype(str).str.strip().str.lower()
letters = ts.map(m)
invalid = letters.isna() & data[type_col - 1].notna()
if invalid.any():
    raise ValueError(f"类型列存在非法取值: {sorted(set(data[type_col - 1][invalid]))}")
pos_upper = data[pos_col - 1].astype(str).str.upper()
letters = letters.mask(pos_upper.str.contains('OCH3', na=False), 'CH3')
letters = letters.where(data[type_col - 1].notna(), None)  # 类型列为空的行，type 列留空

# 直接访问内部单元格字典 ws._cells，绕过逐格代理查找；只写有变化的格
cells = ws._cells
for i in range(len(letters)):
    r = i + 2
    tcell = cells.get((r, type_new_col))
    if tcell is None:
        tcell = ws.cell(row=r, column=type_new_col)
    if tcell.value != letters.iat[i]:
        tcell.value = letters.iat[i]

wb.save(SRC)
print(f'openpyxl 总耗时 {time.perf_counter() - t0:.2f} s')
