# -*- coding: utf-8 -*-
"""归属表对比脚本: 判断两个归属表代表的化合物结构是否一致.

核心逻辑:
  1. 以 (13C, 1H, 碳类型) 联合约束做 Hungarian 最优一对一匹配 (默认)
  2. 可选 --smiles 启用 RDKit SMILES 结构判等

默认仅比较 NMR 归属表, 不进行化合物结构判等.
行序可以不同, 列数/元信息行数/表头位置均可自适应.

用法:
    python compare_tables.py
    python compare_tables.py --prog /path/to/output.xlsx
    python compare_tables.py --ref /path/to/reference.xlsx
    python compare_tables.py --smiles    # 额外启用 SMILES 结构判等
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


def load_table(filepath: str | Path) -> tuple[pd.DataFrame, str | None]:
    """加载归属表 xlsx, 返回 (数据 DataFrame, SMILES 或 None).

    自适应: 表头位置 / 列数 / 列序 / 元信息行数.
    """
    raw = pd.read_excel(filepath, header=None)
    header_idx = _find_header_row(raw)
    col_map = _build_col_map(raw.iloc[header_idx])
    smiles = _extract_smiles_above_header(raw, header_idx)

    if 'No' not in col_map or '13C' not in col_map:
        raise ValueError(f'表头缺少必要列 (No/13C), 实际识别到: {list(col_map.keys())}')

    # 取表头下方的数据行
    data = raw.iloc[header_idx + 1:].reset_index(drop=True)

    # 按映射提取列
    result = pd.DataFrame()
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

    return result.reset_index(drop=True), smiles


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


def main():
    import argparse
    parser = argparse.ArgumentParser(description='归属表化合物结构一致性对比')
    parser.add_argument('--prog', default=None, help='程序输出 xlsx 路径')
    parser.add_argument('--ref', default=None, help='真实结构 xlsx 路径')
    parser.add_argument('--smiles', action='store_true',
                        help='启用 SMILES 结构判等 (默认仅做 NMR 归属表比较)')
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
    df_prog, smiles_prog = load_table(file_prog)
    df_ref, smiles_ref = load_table(file_ref)

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


if __name__ == '__main__':
    main()
