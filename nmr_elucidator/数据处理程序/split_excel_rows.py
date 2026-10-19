# -*- coding: utf-8 -*-
"""把 Excel 中一个单元格含两行 δH 的行拆分为两行, 并补全其余列.

补全规则参考 assign_from_mnova.load_xlsx:
1. 单元格内换行: δH 含 '\\n' → 拆成 2 行, No./δc 向下填充 (补全行, 其余列继承原行);
2. 续行写法: No./δc 留空但 δH 有值 → 是上一行的第 2 个 H, 向下填充补全 No./δc.
'一' 表示该碳无氢 (季碳), 原样保留.
原始文件只读, 结果另存 <同名>_拆分.xlsx (补全版 + 置空版两个 sheet).
"""
import os
import sys

import pandas as pd

SRC = r'Y:\LXH数据处理\新化合物处理\程序解析\clxh122解析\clxh122_H参考.xlsx'


def main():
    df = pd.read_excel(SRC, header=0)
    df.columns = [str(c).strip() for c in df.columns]
    for col in ('No.', 'δ c', 'δH'):
        if col not in df.columns:
            raise SystemExit(f'ERROR: xlsx 缺少必需列 {col}: {SRC}')
    df = df.copy()
    # 统一 δH 为字符串, 便于按换行拆分
    df['δH'] = df['δH'].map(
        lambda v: '' if pd.isna(v) else str(v).replace('\r\n', '\n').strip())

    # ---- 续行归并: No./δc 留空但 δH 有值 → 上一行的第 2 个 H ----
    # 仅当 δH 以化学位移数值开头 (含小数点) 时才算续行, 排除 '3-Api' 之类文字行
    _is_shift = df['δH'].str.match(r'^\d+\.\d+')
    cont = df['No.'].isna() & df['δ c'].isna() & df['δH'].ne('') & _is_shift
    df['_cont'] = cont
    if cont.any():
        for c in ('No.', 'δ c'):
            df.loc[cont, c] = df[c].ffill()[cont]

    # ---- 单元格内换行拆分: δH 含 '\n' → 每个子行一行, No./δc 仅首行保留 ----
    rows = []
    for _, r in df.iterrows():
        parts = [p.strip() for p in str(r['δH']).split('\n') if p.strip()]
        if len(parts) <= 1:
            row = dict(r)
            row['δH'] = parts[0] if parts else ''
            rows.append(row)
            continue
        for i, p in enumerate(parts):
            row = dict(r)
            row['δH'] = p
            if i > 0:
                # 拆分出的第 2 行: 补全 No./δc, 标记为补全行 (置空版输出时置空)
                row['_cont'] = True
            rows.append(row)
    out = pd.DataFrame(rows)
    out = out.dropna(subset=['δ c']).reset_index(drop=True)

    # 补全版 (ppm/δc 完整) 与 置空版 (补全行 No./δc/δH-续行 置空, 参考 xlsx 续行格式)
    out_blank = out.copy()
    m_cont = out_blank['_cont'].fillna(False).astype(bool)
    out_blank['No.'] = out_blank['No.'].astype(object)
    out_blank['δ c'] = out_blank['δ c'].astype(object)
    out_blank.loc[m_cont, ['No.', 'δ c']] = ''
    out_full = out.drop(columns='_cont')
    out_blank = out_blank.drop(columns='_cont')

    dst = os.path.splitext(SRC)[0] + '_拆分.xlsx'
    with pd.ExcelWriter(dst, engine='openpyxl') as xw:
        out_full.to_excel(xw, sheet_name='拆分补全版', index=False)
        out_blank.to_excel(xw, sheet_name='拆分置空版', index=False)
        from openpyxl.styles import Alignment
        for si in (0, 1):
            ws = xw.book.worksheets[si]
            for col_name in ('δ c',):
                ci = list(out.columns).index(col_name) + 1
                for row in ws.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '0.0'
            ci = list(out.columns).index('δH') + 1
            for row in ws.iter_rows(min_row=2, min_col=ci, max_col=ci):
                for cell in row:
                    cell.number_format = '@'
                    cell.alignment = Alignment(horizontal='left')
    print(f'OK: {len(df)} 行 -> {len(out)} 行, 输出: {dst}')


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
