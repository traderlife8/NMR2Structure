# -*- coding: utf-8 -*-
"""独立函数检查: 归属表中所有氢是否归属完整.

提示词：首先要统计δH_ppm中相同的字符串出现了几次，正常一行算一次（例如CH和CH2中每行对应一个H），
但是特殊情况例如CH3如果没有分为三行，则一行对应三个H，相应的δH_ppm字符串在该行出现的次数要*3，
最后汇总到1H_NMR表中的δH_nH列，表示1H_NMR表中的δH_ppm每行在归属表中出现的次数，
并与nH列向量化比较每行是否相等，不相等的行数字标红。

原理: 统计归属表各行 δH_ppm 字符串占用的 H 数 (正常一行 1 个 H;
特殊情况 CH3 若未拆成多行, 一行代表 3 个 H, 该行计数 ×3),
汇总到 1H_NMR sheet 的 δH_ppm → δH_nH 列 (该信号在归属表中被占用的 H 数),
与 nH 列逐行向量化比较, 不相等的行 (δH_nH/nH 数值) 标红.

用法:
    python check_h_assignment.py <归属表.xlsx>
    python check_h_assignment.py            # 默认 clxh126A 输出表
"""
import os
import sys

import pandas as pd


def check_h_assignment(rows_df: pd.DataFrame, sig_df: pd.DataFrame) -> pd.DataFrame:
    """统计归属表各 δH_ppm 信号被占用的 H 数, 与积分表 nH 比对.

    rows_df: 归属表 (需含 ppm / Types / δH_ppm 列; δH_ppm 空 = 该行未归属 H, 不计数);
    sig_df: 1H_NMR 积分表 (需含 δH_ppm / nH 列).
    权重规则: 同 (ppm, Types) 槽仅 1 行且 Types == CH3 → 该行代表 3H (权重 3);
    其余每行权重 1 (CH/CH2 每行一个 H; CH3 拆成多行时每行也各算 1).
    返回: sig_df 副本附加 δH_nH 列 (该信号在归属表中被占用的 H 数, 未出现为 0).
    """
    from collections import Counter
    # 槽内行数: 同 (ppm, Types) 组的行数 (决定 CH3 权重)
    slot_rows: dict = {}
    counts: Counter = Counter()
    for _, r in rows_df.iterrows():
        s = '' if pd.isna(r.get('δH_ppm')) else str(r['δH_ppm']).strip()
        if not s or s.lower() == 'nan':
            continue  # 未归属行 (Cq / 待补 H) 不计数
        key = (round(float(r['ppm']), 2), str(r['Types']).strip())
        slot_rows[key] = slot_rows.get(key, 0) + 1
    for _, r in rows_df.iterrows():
        s = '' if pd.isna(r.get('δH_ppm')) else str(r['δH_ppm']).strip()
        if not s or s.lower() == 'nan':
            continue
        typ = str(r['Types']).strip()
        key = (round(float(r['ppm']), 2), typ)
        w = 3 if (typ == 'CH3' and slot_rows.get(key, 0) == 1) else 1
        counts[s] += w
    sig_df = sig_df.copy()
    sig_df['δH_nH'] = sig_df['δH_ppm'].map(lambda s: counts.get(str(s).strip(), 0))
    return sig_df


def report(xlsx_path: str) -> pd.DataFrame:
    """读取归属表输出 Excel, 执行检查并打印报告; 返回带 δH_nH 的 1H_NMR 表."""
    xl = pd.ExcelFile(xlsx_path)
    sheet_a = xl.sheet_names[0]   # 归属表 (补全版)
    sheet_s = next(s for s in xl.sheet_names if s.endswith('1H_NMR'))
    rows_df = pd.read_excel(xl, sheet_name=sheet_a)
    sig_df = pd.read_excel(xl, sheet_name=sheet_s)
    sig_df = check_h_assignment(rows_df, sig_df)
    m_bad = sig_df['δH_nH'] != sig_df['nH']
    print(f'[检查] {xlsx_path}')
    print(f'[检查] 积分信号 {len(sig_df)} 条, 归属完整 {int((~m_bad).sum())} 条, '
          f'不匹配 {int(m_bad.sum())} 条:')
    if m_bad.any():
        print(sig_df[m_bad][['δH', 'δH_ppm', 'nH', 'δH_nH']].to_string(index=False))
    return sig_df


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    path = sys.argv[1] if len(sys.argv) >= 2 else \
        r'Y:\LXH数据处理\新化合物处理\程序解析\clxh126A解析\clxh126A_归属表.xlsx'
    report(path)
