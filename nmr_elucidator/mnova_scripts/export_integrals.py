"""从 .mnova 导出的 1H 谱 JSON 矩阵重建积分区段, 与 氢谱积分.txt 比对核对.

背景: .mnova 文档内未保存积分区间表 (spec.integrals() 均为空), 氢谱积分.txt 是
化学家的外部报告. 本脚本用已导出的 1D 1H 实数矩阵做分段积分:
  1. 噪声阈值 (基线 MAD) 分割信号区段
  2. 每段面积, 以最大单 H 段 (nH=1 的 d 峰) 为基准归一化
  3. 与 txt 的 (δ, nH) 列表逐段比对: δ 覆盖 + nH 偏差
直接运行: python export_integrals.py [化合物名]   (默认 _DEFAULT_COMPOUND)
"""
from __future__ import annotations

import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assign_from_mnova import _DEFAULT_COMPOUND, _DATA_DIR, _MNOVA_OUT_ROOT, parse_1h_report

# 溶剂/杂质参考峰 (CDCl3): 7.26 残留 CHCl3; 水峰 ~1.5; TMS 0
_SOLVENT_RANGES = [(7.19, 7.33), (1.49, 1.60)]  # 宽容差覆盖


def load_1h_json(json_dir: str):
    """找 1H 1D 谱 JSON (取峰数最多), 返回 (y, x_ppm)."""
    best = None
    for fn in sorted(os.listdir(json_dir)):
        if not fn.endswith('.json') or fn.startswith('_'):
            continue
        with open(os.path.join(json_dir, fn), encoding='utf-8') as f:
            j = json.load(f)
        if j.get('dim') != 1 or j.get('nucleus') != '1H' or 'data' not in j:
            continue
        if best is None or len(j['peaks']) > len(best['peaks']):
            best = j
    if best is None:
        raise SystemExit(f'ERROR: 未找到 1H 1D JSON: {json_dir}')
    y = np.array(best['data'], dtype=float)
    x = np.linspace(float(best['firstXppm']), float(best['lastXppm']), len(y))
    return y, x


def segment_integrals(y, x):
    """噪声阈值分割 → [(hi_ppm, lo_ppm, area)] 面积降序."""
    from scipy import ndimage
    sm = ndimage.uniform_filter1d(y, 9)
    thr = y.max() * 0.02
    lab, n = ndimage.label(sm > thr)
    segs = []
    for k in range(1, n + 1):
        idx = np.where(lab == k)[0]
        if len(idx) < 5:
            continue
        area = float(np.sum(np.clip(y[idx], 0, None)))
        segs.append({'hi': float(x[idx].max()), 'lo': float(x[idx].min()),
                     'area': area})
    segs.sort(key=lambda s: -s['area'])
    return segs


def is_solvent(lo: float, hi: float) -> bool:
    c = (lo + hi) / 2
    return any(a <= c <= b for a, b in _SOLVENT_RANGES)


def main():
    compound = sys.argv[1] if len(sys.argv) >= 2 else _DEFAULT_COMPOUND
    json_dir = os.path.join(_MNOVA_OUT_ROOT, compound)
    txt_path = os.path.join(_DATA_DIR, f'{compound}氢谱积分.txt')
    if not os.path.isdir(json_dir):
        raise SystemExit(f'ERROR: JSON 目录不存在: {json_dir} (先运行 assign_from_mnova.py)')
    if not os.path.isfile(txt_path):
        raise SystemExit(f'ERROR: 氢谱积分.txt 不存在: {txt_path}')

    y, x = load_1h_json(json_dir)
    segs = segment_integrals(y, x)
    sig_segs = [s for s in segs if not is_solvent(s['lo'], s['hi'])]
    sol_segs = [s for s in segs if is_solvent(s['lo'], s['hi'])]
    print(f'[积分重建] {len(segs)} 段 (溶剂/水 {len(sol_segs)} 段排除, 剩 {len(sig_segs)})')

    truth = parse_1h_report(txt_path)
    # 归一化: 用"未重叠且 nH=1 的最强 d/dd/t 尖峰段"做每 H 基准
    # 候选: txt 中 nH=1 且 δ±0.08 内只有一段的信号 → 段面积/nH
    unit_cands = []
    for t in truth:
        if t['nH'] != 1:
            continue
        grp = [r for r in sig_segs if abs((r['hi'] + r['lo']) / 2 - t['ppm']) <= 0.08]
        if len(grp) == 1 and grp[0]['area'] > 0:
            unit_cands.append(grp[0]['area'] / t['nH'])
    if unit_cands:
        # 取中位数基准 (抗单峰异常)
        unit = float(np.median(unit_cands))
        print(f'[归一化] 基准峰 {len(unit_cands)} 个 (nH=1 孤立段) → 每H {unit:.0f}')
    else:
        truth_total = sum(s['nH'] for s in truth)
        seg_total = sum(s['area'] for s in sig_segs)
        unit = seg_total / truth_total
        print(f'[归一化] 回退总面积比: 每H {unit:.0f}')

    rows = []
    for s in sig_segs:
        nh = s['area'] / unit
        rows.append({**s, 'nH_est': nh})
    rows.sort(key=lambda r: -r['hi'])

    # 与 txt 比对: 信号中心 (或区间) ±0.08 内的段求和
    print(f"\n{'txt δH':>8} {'mult':>6} {'nH_txt':>6} | {'段范围':>16} {'nH_积分':>7}")
    matched = 0
    for t in truth:
        grp = [r for r in rows if abs((r['hi'] + r['lo']) / 2 - t['ppm']) <= 0.08]
        if not grp:
            print(f"{t['ppm']:8.2f} {t['mult']:>6} {t['nH']:6.0f} | {'<未找到段>':>16}")
            continue
        area_sum = sum(g['area'] for g in grp)
        nh_est = area_sum / unit
        rng = f"{grp[0]['hi']:.2f}~{grp[-1]['lo']:.2f}"
        flag = 'OK' if abs(nh_est - t['nH']) <= max(0.6, t['nH'] * 0.25) else '偏差'
        if flag == 'OK':
            matched += 1
        print(f"{t['ppm']:8.2f} {t['mult']:>6} {t['nH']:6.0f} | {rng:>16} {nh_est:7.2f} {flag}")
    print(f'\n[比对] {matched}/{len(truth)} 信号 nH 一致 (±max(0.6, 25%))')


if __name__ == '__main__':
    main()
