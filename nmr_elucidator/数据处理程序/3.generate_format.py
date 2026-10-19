# -*- coding: utf-8 -*-
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font

SRC = r'D:\01资料撰写\博士论文撰写\化合物数据归属\2026.10.04_已知化合物结构解析OK\CP13.xlsx'

# 输出文件：不能直接覆盖 SRC。openpyxl 保存时会丢弃公式的缓存计算结果，
# 一旦覆盖 SRC，ppm 列（公式列）下次读出来就全是 NaN
DST = SRC.replace('.xlsx', '_格式.xlsx')

df = pd.read_excel(SRC)

def find_col(candidates):
    # 按候选顺序优先匹配，列名不区分大小写
    for n in candidates:
        col = next((c for c in df.columns if str(c).strip().lower() == n.lower()), None)
        if col is not None:
            return col
    raise KeyError(f"未找到列，候选: {candidates}，现有列名: {df.columns.tolist()}")

pos_col = find_col(('position', 'pos.', 'pos', 'no.', 'no'))

# 撇号统一为 ASCII（′->', ″->'', ‴->''', ⁗->''''），与比对文件格式一致
pos = df[pos_col].astype(str).str.strip()
pos = pos.apply(lambda s: ''.join("'" if ch in '′″‴⁗' else ch for ch in s))

# 删除最后一个数字之前的所有字符串（OCH3 特殊值保留原样），与 1.sort_position.py 一致
is_special = pos.str.contains('OCH3', case=False)
pos = pos.where(is_special, pos.str.extract(r"(\d+[^\d]*)$")[0])

# 删除最后一个数字之后的 a/b/α/β，数字后的所有′保留：1a->1、6a'->6'、6''''a->6''''（OCH3 中无这些字母，不受影响）
pos = pos.str.replace(r'[aAbBαβ]', '', regex=True)

# 类型列：'type' 首选（不区分大小写），'types'、'类型' 次选
type_col = find_col(('type', 'types', '类型'))

types = df[type_col].astype(str).str.strip()
bad = sorted(set(types.str.lower()) - {'s', 'd', 't', 'q'})
if bad:
    raise ValueError(f"类型列存在非法取值 {bad}，只允许 s/d/t/q（不区分大小写）")

# ppm列：'ppm' 首选（不区分大小写），'您的碳谱'、'13C NMR' 次选
ppm_col = find_col(('ppm', '您的碳谱', '13C NMR', 'δC'))

# ppm列可能是文本格式，先整列转为浮点数，再由 f'{v:.1f}' 保留一位小数
ppm = pd.to_numeric(df[ppm_col], errors='raise')

fmt = ppm.map(lambda v: f'{v:.1f}') + ' (' + types + ', C-' + pos + ')'

# 参考比对文件 word格式 样式：将格式列各项按行顺序用", "拼接
word = ', '.join(fmt)

# 用 openpyxl 在 SRC 的基础上追加“格式”“word格式”列（已存在则覆盖），
# 结果写入 DST，SRC 保持不动（否则会清除 ppm 列公式的缓存值）
wb = load_workbook(SRC)
ws = wb.active
header = [c.value for c in ws[1]]
hdrs = [str(h).strip() if h is not None else '' for h in header]

if '格式' in hdrs:
    col_idx = hdrs.index('格式') + 1
else:
    col_idx = len(header) + 1
    ws.cell(row=1, column=col_idx, value='格式').font = Font(bold=True)
for i, v in enumerate(fmt, start=2):
    ws.cell(row=i, column=col_idx, value=v).font = Font(name='Arial')

if 'word格式' in hdrs:
    col_idx2 = hdrs.index('word格式') + 1
else:
    col_idx2 = len(header) + 2
    ws.cell(row=1, column=col_idx2, value='word格式').font = Font(bold=True)
ws.cell(row=2, column=col_idx2, value=word)

wb.save(DST)

print(df.assign(格式=fmt).head(10).to_string())

