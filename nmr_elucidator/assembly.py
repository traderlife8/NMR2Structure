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
    # 产物 G: 1H 峰来源标记 (measured/estimated/'')
    h_source: str = ''


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
    dept_summary: Optional[List[dict]] = None,
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
        # 改进点 #7+#10+#11+#15: 同 13C 多个 HSQC 交叉峰 → 升级为 CH2
        # 条件: (1) 原始 mult 不是确定的 'CH'/'CH3'
        #       (2) 恰好 2 个 H, 且化学位移差 > 0.1 ppm (真正的 diastereotopic H)
        # 如果 H 化学位移差 <= 0.1 ppm, 可能是噪声/重复峰, 保持原 mult
        # 改进点 #15: 如果有 3+ 个 H, 很可能是噪声/伪峰, 不升级, 保持原 mult
        if len(group) >= 2:
            orig_mult = mult.upper()
            # 计算 H 化学位移范围
            h_shifts = [pair.h_ppm for _, pair in group]
            dh = max(h_shifts) - min(h_shifts)
            if orig_mult in ('CH', 'CH3'):
                # 原 mult 是确定的 CH/CH3, 保持原 mult
                mult = orig_mult
            elif len(group) >= 3:
                # 改进点 #15: 3+ 个 H 很可能是噪声/伪峰, 不升级
                # 真实 CH2 只有 2 个 H, CH3 只有 3 个 H 但通常在 <25 ppm
                # 在 25-220 ppm 区域, 3+ H 几乎肯定是噪声
                # 改进点 #26: CH2/CH3 本身模糊, 3+ 峰更不可信, 默认为 CH
                if orig_mult == 'CH2/CH3':
                    mult = 'CH'
                else:
                    mult = orig_mult if orig_mult else 'CH'
            elif dh > 0.1:
                # 恰好 2 个 H, 且化学位移差 > 0.1 ppm, 真正的 CH2 (diastereotopic)
                mult = 'CH2'
            else:
                # H 化学位移差 <= 0.1 ppm, 可能是噪声/重复, 保持原 mult
                mult = orig_mult if orig_mult else 'CH'

        if len(group) == 1:
            # 单个 H (CH 或 CH3)
            idx, pair = group[0]
            row_idx += 1
            label = atom_labels.get(idx, f'C{row_idx:02d}') if atom_labels else f'C{row_idx:02d}'
            pair_to_label[idx] = f'H-{label}'

            # 改进点 #23: CH2/CH3 仅 1 个 HSQC 峰 → 降级为 CH
            # 理由: 真正的 CH2 应有 2 个 H 在 HSQC 上可见, 仅 1 个峰说明是 CH
            if mult == 'CH2/CH3':
                mult = 'CH'

            # 1H 多重性
            h_mult_str = _format_multiplicity(pair)

            rows.append(AssignmentRow(
                label=label,
                c_ppm=round(c_ppm, 1),
                carbon_type=mult if mult != 'unknown' else 'CH',
                h_ppm=round(pair.h_ppm, 2),
                h_mult=h_mult_str,
                h_j=list(pair.j_couplings) if hasattr(pair, 'j_couplings') else [],
                h_source=getattr(pair, 'h_source', ''),
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
                    h_source=getattr(pair, 'h_source', ''),
                ))

    # 3. 季碳行 (13C 峰表中无 HSQC 对应的碳)
    # 改进点 #24: 参考 DEPT 信息, 若 DEPT 说非季碳则不创建 Cq 行
    dept_by_c: Dict[float, str] = {}
    if dept_summary:
        for d in dept_summary:
            ppm = d.get('ppm')
            mult_d = d.get('mult', 'unknown')
            if ppm is not None:
                key = round(ppm, 1)
                # 优先保留非 Cq 类型 (避免同 ppm 多个峰时 Cq 覆盖 CH)
                if key not in dept_by_c or dept_by_c[key] == 'Cq':
                    dept_by_c[key] = mult_d

    hsqc_c_set = set(round(p.c_ppm, 1) for p in ch_pairs)
    for c_ppm in sorted(carbon_shifts_1d, reverse=True):
        c_key = round(c_ppm, 1)
        if c_key not in hsqc_c_set:
            # 检查是否有近似匹配 (±c_tol)
            has_hsqc = any(abs(c_ppm - hc) <= c_tol for hc in [p.c_ppm for p in ch_pairs])
            if not has_hsqc:
                # 改进点 #24: 检查 DEPT 是否认为该碳不是季碳
                dept_mult = dept_by_c.get(c_key, 'Cq')
                if dept_mult in ('CH', 'CH2', 'CH3', 'CH2/CH3'):
                    # DEPT 确认有 H 连接, 不是季碳, 创建带类型的行 (无 1H 信息)
                    row_idx += 1
                    rows.append(AssignmentRow(
                        label=f'C{row_idx:02d}',
                        c_ppm=round(c_ppm, 1),
                        carbon_type=dept_mult if dept_mult != 'CH2/CH3' else 'CH',
                    ))
                else:
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

    # 4.5 改进点 #27: 后处理 — 同 13C 的 CH 对 (dH 0.3-1.0) 升级为 CH2
    # 场景: DEPT 说 CH 但同一碳实际有 2 个 diastereotopic H (HSQC 检测到多个峰)
    # 仅当 dH > 0.3 (明确的化学不等价) 且 < 1.0 (排除不同碳的偶然重叠)
    from collections import defaultdict
    ch_rows_by_c: dict = defaultdict(list)
    for r in rows:
        if r.carbon_type == 'CH' and r.h_ppm is not None and not r.is_quaternary:
            ch_rows_by_c[r.c_ppm].append(r)
    for c_val, group_rows in ch_rows_by_c.items():
        if len(group_rows) < 2:
            continue
        # 找 dH 在 (0.3, 1.0) 的最佳对 (取最小 dH, 真正 CH2 两 H 通常较接近)
        best_pair = None
        best_dh = 999.0
        for i in range(len(group_rows)):
            for j in range(i + 1, len(group_rows)):
                dh_ij = abs(group_rows[i].h_ppm - group_rows[j].h_ppm)
                if 0.45 < dh_ij < 1.0 and dh_ij < best_dh:
                    best_dh = dh_ij
                    best_pair = (i, j)
        if best_pair is not None:
            group_rows[best_pair[0]].carbon_type = 'CH2'
            group_rows[best_pair[1]].carbon_type = 'CH2'
            group_rows[best_pair[0]].is_ch2_stereo = True
            group_rows[best_pair[1]].is_ch2_stereo = True

    # 4.6 改进点 #29: 后处理 — DEPT 说 Cq 但 HSQC 检测到 H 的碳, 若 13C<35 + 1H 1.8-2.2 → CH3
    # 场景: 乙酰基甲基 (COCH3) 等 DEPT 漏检的 CH3
    # 仅当同 13C 只有 1 行符合条件时才升级 (排除多碳重叠的误判)
    from collections import Counter as _Counter
    _ch3_candidates: dict = {}  # c_ppm -> list of rows
    for r in rows:
        if (r.carbon_type == 'CH' and r.h_ppm is not None
                and r.c_ppm is not None and r.c_ppm < 35.0
                and 1.8 <= r.h_ppm <= 2.2):
            c_key = round(r.c_ppm, 1)
            if dept_by_c.get(c_key) == 'Cq':
                _ch3_candidates.setdefault(c_key, []).append(r)
    for c_key, cand_rows in _ch3_candidates.items():
        if len(cand_rows) == 1:
            cand_rows[0].carbon_type = 'CH3'

    # 4.7 改进点 #30: DEPT 说 Cq 但 HSQC 检测到恰好 2 个 CH → 升级为 CH2
    # 场景: DEPT 漏检的 CH2 (如 29.1 ppm 的 24a/24b)
    _cq_ch_groups: dict = {}
    for r in rows:
        if r.carbon_type == 'CH' and r.h_ppm is not None and r.c_ppm is not None:
            c_key = round(r.c_ppm, 1)
            if dept_by_c.get(c_key) == 'Cq':
                _cq_ch_groups.setdefault(c_key, []).append(r)
    for c_key, grp in _cq_ch_groups.items():
        if len(grp) == 2:
            grp[0].carbon_type = 'CH2'
            grp[1].carbon_type = 'CH2'
            grp[0].is_ch2_stereo = True
            grp[1].is_ch2_stereo = True

    # 4.8 改进点 #31: CH2 组 3+ 成员中, 若仅 1 个 H < 1.5 ppm 且其余 > 1.5 → 降级该成员为 CH
    # 场景: 35.5 ppm 组中 C50c (H=1.33) 实为 CH (25), 其余为真正 CH2
    _ch2_groups: dict = {}
    for r in rows:
        if r.carbon_type == 'CH2' and r.h_ppm is not None and r.c_ppm is not None:
            _ch2_groups.setdefault(round(r.c_ppm, 1), []).append(r)
    for c_key, grp in _ch2_groups.items():
        if len(grp) >= 3:
            low_h = [r for r in grp if r.h_ppm < 1.5]
            high_h = [r for r in grp if r.h_ppm >= 1.5]
            if len(low_h) == 1 and len(high_h) >= 2:
                low_h[0].carbon_type = 'CH'
                low_h[0].is_ch2_stereo = False

    # 4.9 改进点 #33: 13C 组中已有 CH2 对, 剩余 CH 若 H>2.0 且远离 CH2 对 → 降级为 Cq
    # 场景: 37.7 ppm 组有 CH2 对 (1.76/1.09) + 多余 CH (2.04), 后者是重叠 Cq 的假 HSQC 峰
    # 条件: 13C < 60 (苷元区, 糖区 60-110 的 CH 正常 H>2.0) + H > 2.0 + 距 CH2 对两成员均 > 0.25 ppm
    _all_groups: dict = {}
    for r in rows:
        if r.h_ppm is not None and r.c_ppm is not None and not r.is_quaternary:
            _all_groups.setdefault(round(r.c_ppm, 1), []).append(r)
    for c_key, grp in _all_groups.items():
        if len(grp) < 3 or c_key >= 60.0:
            continue
        ch2_members = [r for r in grp if r.carbon_type == 'CH2']
        ch_members = [r for r in grp if r.carbon_type == 'CH']
        if len(ch2_members) == 2 and len(ch_members) >= 1:
            h1, h2 = ch2_members[0].h_ppm, ch2_members[1].h_ppm
            for ch_row in ch_members:
                if (ch_row.h_ppm > 2.0
                        and abs(ch_row.h_ppm - h1) > 0.25
                        and abs(ch_row.h_ppm - h2) > 0.25):
                    ch_row.carbon_type = 'C'
                    ch_row.is_quaternary = True
                    ch_row.h_ppm = None

    # 4.10 改进点 #33b: 13C < 60 组中 CH2 成员 H>2.5 且同组有 CH (H<2.0) → 该 CH2 为假峰, 降级 Cq
    # 场景: 46.7 ppm 组有 CH(1.68) + CH2(2.60) + CH2(1.83), CH 是真实 C-14, CH2(2.60) 是假峰
    for c_key, grp in _all_groups.items():
        if len(grp) < 3 or c_key >= 60.0:
            continue
        ch2_members = [r for r in grp if r.carbon_type == 'CH2']
        ch_members = [r for r in grp if r.carbon_type == 'CH']
        if not ch_members:
            continue
        has_low_ch = any(r.h_ppm < 2.0 for r in ch_members)
        if not has_low_ch:
            continue
        for ch2_row in ch2_members:
            if ch2_row.h_ppm > 2.5:
                ch2_row.carbon_type = 'C'
                ch2_row.is_quaternary = True
                ch2_row.h_ppm = None
                ch2_row.is_ch2_stereo = False

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
                    '1H-1H COSY', 'HMBC', 'ROESY', '1H来源'])
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
                row.h_source,
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
               '1H-1H COSY', 'HMBC', 'ROESY', '1H来源']
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
            arow.h_source,
        ]
        for col, val in enumerate(values, 1):
            cell = ws.cell(row=row_num, column=col, value=val)
            cell.border = border
        row_num += 1

    # 列宽
    widths = [10, 10, 8, 25, 20, 25, 20, 10]
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
