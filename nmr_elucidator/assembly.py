# -*- coding: utf-8 -*-
"""Part 3 Phase 7: 归属表组装.

从 2D 解谱中间结果 (C-H 对 + COSY/HMBC/ROESY 交叉峰) 组装化学家格式归属表.
不依赖打分公式, 纯粹从谱图数据推导每行内容.

输入:
  - ch_pairs: HSQC C-H 对 (from interpret_2d)
  - carbon_shifts_1d: 1D 13C 峰表
  - proton_shifts_1d: 1D 1H 峰表
  - cosy_pairs: COSY 交叉峰 (dH1, dH2)
  - hmbc_pairs: HMBC 交叉峰 (dH, dC)
  - roesy_pairs: ROESY 交叉峰 (dH1, dH2)

输出:
  - AssignmentRow 列表 (对标化学家 Excel 8 列)
  - CSV / XLSX 文件

依赖: 标准库 + numpy (可选 openpyxl 写 xlsx)
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---- 数据结构 ----

@dataclass
class AssignmentRow:
    """归属表一行 (对标化学家 Excel).

    列: No. / 13C NMR / Types / 1H NMR (mult., J in Hz) / COSY / HMBC / ROESY
    """

    label: str                      # 原子标签 (如 "1a", "Glc-1'") 或匿名 (如 "C01-Ha")
    c_ppm: Optional[float] = None   # 13C 化学位移
    carbon_type: str = ''           # CH/CH2/CH3/C
    h_ppm: Optional[float] = None   # 1H 化学位移
    h_mult: str = ''                # 多重性 (s/d/t/q/m/dd/...)
    h_j: List[float] = field(default_factory=list)  # J 耦合值 (Hz)
    cosy_partners: List[str] = field(default_factory=list)   # COSY 相关 H 标签
    hmbc_targets: List[str] = field(default_factory=list)    # HMBC 远程 C 标签
    roesy_contacts: List[str] = field(default_factory=list)  # ROESY 空间 H 标签
    is_ch2_stereo: bool = False     # 是否为 CH2 的立体氢 (a/b)
    is_quaternary: bool = False     # 是否为季碳行


# ---- 核心组装逻辑 ----

def assemble_assignment_table(
    ch_pairs: list,
    carbon_shifts_1d: List[float],
    proton_shifts_1d: List[float],
    cosy_pairs: List[Tuple[float, float]] = None,
    hmbc_pairs: List[Tuple[float, float]] = None,
    roesy_pairs: List[Tuple[float, float]] = None,
    h_tol: float = 0.08,
    c_tol: float = 2.0,
    atom_labels: Optional[Dict[int, str]] = None,
) -> List[AssignmentRow]:
    """组装完整归属表.

    Args:
        ch_pairs: CHPair 列表 (from interpret_2d, 含 h_ppm, c_ppm, multiplicity)
        carbon_shifts_1d: 1D 13C 观测位移列表
        proton_shifts_1d: 1D 1H 观测位移列表
        cosy_pairs: COSY 交叉峰 [(dH1, dH2), ...]
        hmbc_pairs: HMBC 交叉峰 [(dH, dC), ...]
        roesy_pairs: ROESY 交叉峰 [(dH1, dH2), ...]
        h_tol: 1H 匹配容差 (ppm)
        c_tol: 13C 匹配容差 (ppm)
        atom_labels: 可选原子编号映射 {ch_pair_index: label}

    Returns:
        归属表行列表 (按 dC 降序)
    """
    cosy_pairs = cosy_pairs or []
    hmbc_pairs = hmbc_pairs or []
    roesy_pairs = roesy_pairs or []

    rows: List[AssignmentRow] = []

    # 1. 按 dC 分组 C-H 对 (同一 dC 可能有多个 H, 如 CH2)
    c_groups: Dict[float, list] = {}
    for i, pair in enumerate(ch_pairs):
        c_key = round(pair.c_ppm, 1)
        if c_key not in c_groups:
            c_groups[c_key] = []
        c_groups[c_key].append((i, pair))

    # 2. 为每个 C-H 对生成行
    pair_to_label: Dict[int, str] = {}  # ch_pair index -> H label
    row_idx = 0

    # 按 dC 降序处理
    sorted_c_keys = sorted(c_groups.keys(), reverse=True)

    for c_key in sorted_c_keys:
        group = c_groups[c_key]
        c_ppm = group[0][1].c_ppm
        mult = group[0][1].multiplicity

        if len(group) == 1:
            # 单个 H (CH 或 CH3)
            idx, pair = group[0]
            row_idx += 1
            label = atom_labels.get(idx, f'C{row_idx:02d}') if atom_labels else f'C{row_idx:02d}'
            pair_to_label[idx] = f'H-{label}'

            # 1H 多重性
            h_mult_str = _format_multiplicity(pair)

            rows.append(AssignmentRow(
                label=label,
                c_ppm=round(c_ppm, 1),
                carbon_type=mult if mult != 'unknown' else 'CH',
                h_ppm=round(pair.h_ppm, 2),
                h_mult=h_mult_str,
                h_j=list(pair.j_couplings) if hasattr(pair, 'j_couplings') else [],
            ))
        else:
            # 多个 H (CH2 的 a/b)
            # 按 dH 降序排列 (高场在前)
            group_sorted = sorted(group, key=lambda x: x[1].h_ppm, reverse=True)
            suffixes = ['a', 'b', 'c', 'd']
            for si, (idx, pair) in enumerate(group_sorted):
                row_idx += 1
                suffix = suffixes[si] if si < len(suffixes) else str(si + 1)
                label = atom_labels.get(idx, f'C{row_idx:02d}{suffix}') if atom_labels else f'C{row_idx:02d}{suffix}'
                pair_to_label[idx] = f'H-{label}'

                h_mult_str = _format_multiplicity(pair)

                rows.append(AssignmentRow(
                    label=label,
                    c_ppm=round(c_ppm, 1),
                    carbon_type=mult,
                    h_ppm=round(pair.h_ppm, 2),
                    h_mult=h_mult_str,
                    h_j=list(pair.j_couplings) if hasattr(pair, 'j_couplings') else [],
                    is_ch2_stereo=True,
                ))

    # 3. 季碳行 (13C 峰表中无 HSQC 对应的碳)
    hsqc_c_set = set(round(p.c_ppm, 1) for p in ch_pairs)
    for c_ppm in sorted(carbon_shifts_1d, reverse=True):
        c_key = round(c_ppm, 1)
        if c_key not in hsqc_c_set:
            # 检查是否有近似匹配 (±c_tol)
            has_hsqc = any(abs(c_ppm - hc) <= c_tol for hc in [p.c_ppm for p in ch_pairs])
            if not has_hsqc:
                row_idx += 1
                rows.append(AssignmentRow(
                    label=f'Cq{row_idx:02d}',
                    c_ppm=round(c_ppm, 1),
                    carbon_type='C',
                    is_quaternary=True,
                ))

    # 4. 回填 COSY/HMBC/ROESY 列
    _backfill_cosy(rows, ch_pairs, cosy_pairs, pair_to_label, h_tol)
    _backfill_hmbc(rows, ch_pairs, hmbc_pairs, carbon_shifts_1d, pair_to_label, h_tol, c_tol)
    _backfill_roesy(rows, ch_pairs, roesy_pairs, pair_to_label, h_tol)

    # 5. 按 dC 降序排列 (季碳也参与排序)
    rows.sort(key=lambda r: -(r.c_ppm or 0))

    return rows


# ---- 辅助函数 ----

def _format_multiplicity(pair) -> str:
    """格式化 1H 多重性字符串."""
    j_vals = getattr(pair, 'j_couplings', []) or []
    mult = getattr(pair, 'multiplicity', 'unknown')

    # 过滤不合理的 J 值 (1H-1H 耦合通常 0-18 Hz)
    j_vals = [j for j in j_vals if 0.5 < abs(j) < 20.0]

    if not j_vals:
        # 无 J 数据, 根据碳级数猜测
        if mult == 'CH3':
            return 's'  # 甲基常为单峰 (无邻位 H)
        return 'm'  # 默认多重峰

    # 有 J 数据
    n_j = len(j_vals)
    if n_j == 1:
        return f'd ({j_vals[0]:.1f})'
    elif n_j == 2:
        return f'dd ({j_vals[0]:.1f}, {j_vals[1]:.1f})'
    else:
        j_str = ', '.join(f'{j:.1f}' for j in j_vals[:3])
        return f'm ({j_str})'


def _find_label_for_h(h_ppm: float, ch_pairs: list,
                      pair_to_label: Dict[int, str],
                      h_tol: float) -> Optional[str]:
    """根据 dH 找到对应的 H 标签."""
    for i, pair in enumerate(ch_pairs):
        if abs(pair.h_ppm - h_ppm) <= h_tol:
            return pair_to_label.get(i, f'H({h_ppm:.2f})')
    return f'H({h_ppm:.2f})'


def _find_c_label(c_ppm: float, carbon_shifts_1d: List[float],
                  ch_pairs: list, c_tol: float) -> str:
    """根据 dC 找到对应的 C 标签."""
    # 先找 C-H 对中的标签
    for pair in ch_pairs:
        if abs(pair.c_ppm - c_ppm) <= c_tol:
            return f'C({pair.c_ppm:.1f})'
    # 季碳
    for c1d in carbon_shifts_1d:
        if abs(c1d - c_ppm) <= c_tol:
            return f'C({c1d:.1f})'
    return f'C({c_ppm:.1f})'


def _backfill_cosy(rows: List[AssignmentRow], ch_pairs: list,
                   cosy_pairs: List[Tuple[float, float]],
                   pair_to_label: Dict[int, str],
                   h_tol: float) -> None:
    """回填 COSY 列: 对每个 C-H 对, 找 COSY 中与其 H 相关的其他 H."""
    if not cosy_pairs:
        return

    for i, pair in enumerate(ch_pairs):
        h_ppm = pair.h_ppm
        partners = []
        for h1, h2 in cosy_pairs:
            other_h = None
            if abs(h1 - h_ppm) <= h_tol:
                other_h = h2
            elif abs(h2 - h_ppm) <= h_tol:
                other_h = h1
            if other_h is not None and abs(other_h - h_ppm) > 0.02:
                label = _find_label_for_h(other_h, ch_pairs, pair_to_label, h_tol)
                if label and label not in partners:
                    partners.append(label)

        # 写入对应的 row
        if partners:
            # 找到对应的 row (按 h_ppm 匹配)
            for row in rows:
                if row.h_ppm is not None and abs(row.h_ppm - h_ppm) <= h_tol:
                    row.cosy_partners = partners[:5]  # 最多 5 个
                    break


def _backfill_hmbc(rows: List[AssignmentRow], ch_pairs: list,
                   hmbc_pairs: List[Tuple[float, float]],
                   carbon_shifts_1d: List[float],
                   pair_to_label: Dict[int, str],
                   h_tol: float, c_tol: float) -> None:
    """回填 HMBC 列: 对每个 C-H 对, 找 HMBC 中从该 H 出发的远程 C."""
    if not hmbc_pairs:
        return

    for i, pair in enumerate(ch_pairs):
        h_ppm = pair.h_ppm
        own_c = pair.c_ppm
        targets = []
        for h, c_remote in hmbc_pairs:
            if abs(h - h_ppm) <= h_tol:
                # 排除 1J 残留 (与自身 C 匹配)
                if abs(c_remote - own_c) <= c_tol:
                    continue
                c_label = _find_c_label(c_remote, carbon_shifts_1d, ch_pairs, c_tol)
                if c_label and c_label not in targets:
                    targets.append(c_label)

        if targets:
            for row in rows:
                if row.h_ppm is not None and abs(row.h_ppm - h_ppm) <= h_tol:
                    row.hmbc_targets = targets[:6]
                    break


def _backfill_roesy(rows: List[AssignmentRow], ch_pairs: list,
                    roesy_pairs: List[Tuple[float, float]],
                    pair_to_label: Dict[int, str],
                    h_tol: float) -> None:
    """回填 ROESY 列: 对每个 C-H 对, 找 ROESY 中与其 H 空间接近的其他 H."""
    if not roesy_pairs:
        return

    for i, pair in enumerate(ch_pairs):
        h_ppm = pair.h_ppm
        contacts = []
        for h1, h2 in roesy_pairs:
            other_h = None
            if abs(h1 - h_ppm) <= h_tol:
                other_h = h2
            elif abs(h2 - h_ppm) <= h_tol:
                other_h = h1
            if other_h is not None and abs(other_h - h_ppm) > 0.02:
                label = _find_label_for_h(other_h, ch_pairs, pair_to_label, h_tol)
                if label and label not in contacts:
                    contacts.append(label)

        if contacts:
            for row in rows:
                if row.h_ppm is not None and abs(row.h_ppm - h_ppm) <= h_tol:
                    row.roesy_contacts = contacts[:5]
                    break


# ---- 输出 ----

def write_assignment_csv(rows: List[AssignmentRow], path: str) -> None:
    """写出归属表 CSV (通用中间产物)."""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['No.', '13C NMR', 'Types', '1H NMR (mult., J in Hz)',
                    '1H-1H COSY', 'HMBC', 'ROESY'])
        for row in rows:
            # 1H 列格式: "dH, mult (J)"
            h_str = ''
            if row.h_ppm is not None:
                h_str = f'{row.h_ppm:.2f}, {row.h_mult}' if row.h_mult else f'{row.h_ppm:.2f}'

            w.writerow([
                row.label,
                f'{row.c_ppm:.1f}' if row.c_ppm is not None else '',
                row.carbon_type,
                h_str,
                '/'.join(row.cosy_partners),
                '/'.join(row.hmbc_targets),
                '/'.join(row.roesy_contacts),
            ])


def write_assignment_xlsx(rows: List[AssignmentRow], path: str) -> None:
    """写出归属表 XLSX (化学家 Excel 格式: 表头 + 数据行).

    **不写任何顶部元信息行** (Solvent / Formula / HRESIMS 等另行写到 meta.txt).
    与化学家手工 Excel 格式完全一致 (第 1 行表头, 后续数据行).
    """
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, Border, Side
    except ImportError:
        # openpyxl 未装, 回退到 CSV
        write_assignment_csv(rows, path.replace('.xlsx', '.csv'))
        return

    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    wb = Workbook()
    ws = wb.active
    ws.title = 'NMR Assignment'

    thin = Side(style='thin')
    border = Border(top=thin, bottom=thin, left=thin, right=thin)
    bold = Font(bold=True)

    # 表头 (第 1 行, 化学家 Excel 格式)
    headers = ['No.', '13C NMR', 'Types', '1H NMR (mult., J in Hz)',
               '1H-1H COSY', 'HMBC', 'ROESY']
    for col, h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=col, value=h)
        cell.font = bold
        cell.border = border

    # 数据行 (从第 2 行开始)
    row_num = 2
    for arow in rows:
        h_str = ''
        if arow.h_ppm is not None:
            h_str = f'{arow.h_ppm:.2f}, {arow.h_mult}' if arow.h_mult else f'{arow.h_ppm:.2f}'

        values = [
            arow.label,
            f'{arow.c_ppm:.1f}' if arow.c_ppm is not None else '',
            arow.carbon_type,
            h_str,
            '/'.join(arow.cosy_partners),
            '/'.join(arow.hmbc_targets),
            '/'.join(arow.roesy_contacts),
        ]
        for col, val in enumerate(values, 1):
            cell = ws.cell(row=row_num, column=col, value=val)
            cell.border = border
        row_num += 1

    # 列宽
    widths = [10, 10, 8, 25, 20, 25, 20]
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[chr(64 + i)].width = w

    wb.save(path)


def write_assignment_meta(path: str,
                          meta: Optional[Dict[str, str]] = None) -> None:
    """写出归属表元信息到独立 .txt 文件 (UTF-8).

    用于存放 Solvent / Formula / HRESIMS / [α] / IR / Field_MHz 等元信息.
    与化学家 Excel 表分开存放, 便于手工修改 + 程序读取.

    Args:
        path: meta.txt 输出路径
        meta: 字典, key 为字段名 (Solvent / Formula / HRESIMS 等), value 为字符串值.
              字段顺序按字典插入顺序输出 (Python 3.7+ 字典有序保证).
    """
    if meta is None:
        meta = {}
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    lines = []
    lines.append('# NMR 归属表元信息 (与 assignment_table.xlsx 分离存放)')
    lines.append('# 生成时间: ' + __import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    lines.append('')
    if not meta:
        lines.append('(无元信息)')
    else:
        for key, val in meta.items():
            lines.append(f'{key}: {val}')
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
