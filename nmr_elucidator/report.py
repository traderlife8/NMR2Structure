"""报告输出: JSON + Markdown + HTML.

设计原则:
- JSON 含完整证据链, 便于程序读取 (含 dept_evidence + h_count_evidence)
- Markdown 给人阅读备查, 强调未匹配峰
- HTML 给人类阅读: 含解谱过程/结果/结构图/未匹配峰/DEPT 推断
- HTML 中所有图 base64 内嵌, 单文件可移植
"""
from __future__ import annotations

import base64
import json
import os
from dataclasses import asdict
from datetime import datetime
from io import BytesIO
from typing import Dict, List, Optional

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .evidence import CandidateEvidence, DeptEvidence, HCountEvidence
from .candidate_table import (
    dedupe_by_inchikey,
    evidences_to_dataframe,
    to_csv as candidates_to_csv,
    to_markdown as candidates_to_markdown,
)

# seaborn 可视化 (可选). 未装时函数返回 None, 报告自动跳过.
try:
    from .viz import (
        plot_candidate_score_heatmap,
        plot_constraint_match_matrix,
        plot_ppm_comparison_scatter,
    )
    _HAS_SEABORN = True
except ImportError:
    _HAS_SEABORN = False


def write_report_json(path: str, evidences: List[CandidateEvidence], meta: Dict) -> None:
    """写出 JSON 报告."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    obj = {
        'meta': {
            **meta,
            'generated_at': datetime.now().isoformat(timespec='seconds'),
            'n_candidates': len(evidences),
        },
        'candidates': [asdict(e) for e in evidences],
    }
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def write_report_markdown(path: str, evidences: List[CandidateEvidence], meta: Dict,
                          top_n: int = 10) -> None:
    """写出 Markdown 报告."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    lines: List[str] = []
    lines.append(f"# NMR 反推结构报告")
    lines.append('')
    lines.append(f"- 生成时间: {datetime.now().isoformat(timespec='seconds')}")
    for k, v in meta.items():
        if k in ('observed_carbon_shifts', 'observed_proton_shifts', 'obs_dept_summary'):
            continue  # 略过大数组
        lines.append(f"- {k}: {v}")
    lines.append(f"- 候选总数: **{len(evidences)}**")
    lines.append('')

    # DEPT 推断
    if meta.get('dept_counts'):
        lines.append('## DEPT 推断碳类型')
        lines.append('')
        for k, v in meta['dept_counts'].items():
            if v > 0:
                lines.append(f"- {k}: **{v}**")
        lines.append('')
        if meta.get('obs_dept_summary'):
            lines.append('| 13C (ppm) | 碳类型 | 方向 |')
            lines.append('|---:|---|:-:|')
            for d in meta['obs_dept_summary']:
                sign_icon = {'positive': '↑', 'negative': '↓', 'unknown': '?'}.get(d.get('sign', '?'), '?')
                lines.append(f"| {d['ppm']:.3f} | **{d['mult']}** | {sign_icon} |")
            lines.append('')

    # 2D NMR 关联验证
    if meta.get('used_2d') and meta.get('2d_experiments'):
        lines.append('## 2D NMR 关联验证')
        lines.append('')
        lines.append(f"检测到的 2D 实验: {', '.join(meta['2d_experiments'])}")
        lines.append('')
        corr_2d = meta.get('correlation_2d', {})
        if corr_2d:
            lines.append('| Rank | SMILES | HSQC | HMBC | COSY | 综合 |')
            lines.append('|---:|---|---:|---:|---:|---:|')
            for e in evidences[:top_n]:
                s2d = corr_2d.get(e.smiles)
                if s2d:
                    lines.append(
                        f"| {e.rank} | `{e.smiles[:30]}` "
                        f"| {s2d['hsqc_match_rate']:.0%} "
                        f"| {s2d['hmbc_match_rate']:.0%} "
                        f"| {s2d['cosy_match_rate']:.0%} "
                        f"| {s2d['combined_2d_score']:.0%} |"
                    )
            lines.append('')

    # 2D NMR 逐步解读过程 (化学家工作流)
    interp = meta.get('interpretation_2d')
    if interp:
        lines.append('## 2D NMR 解谱过程 (化学家工作流)')
        lines.append('')
        # HSQC C-H 对
        ch_pairs = interp.get('ch_pairs', [])
        if ch_pairs:
            lines.append(f'### 步骤 1: HSQC C-H 直接连接对 ({len(ch_pairs)} 个)')
            lines.append('')
            lines.append('| 编号 | δH (ppm) | δC (ppm) | 碳级数 |')
            lines.append('|---|---:|---:|---|')
            for p in ch_pairs:
                lines.append(f"| {p.get('c_label', '-')} "
                             f"| {p['h_ppm']:.3f} | {p['c_ppm']:.3f} "
                             f"| {p.get('multiplicity', '?')} |")
            lines.append('')
        # COSY 自旋系统
        spin_systems = interp.get('spin_systems', [])
        if spin_systems:
            lines.append(f'### 步骤 2: COSY 自旋系统 ({len(spin_systems)} 个)')
            lines.append('')
            for ss in spin_systems:
                h_str = ' → '.join(f"{h:.3f}" for h in ss.get('h_ppms', []))
                lines.append(f"- **{ss.get('label', '?')}**: {h_str}")
            lines.append('')
        # HMBC 远程连接
        hmbc_conns = interp.get('hmbc_connections', [])
        if hmbc_conns:
            lines.append(f'### 步骤 3: HMBC 远程连接 ({len(hmbc_conns)} 个)')
            lines.append('')
            lines.append('| δH (ppm) | → δC (ppm) | 键类型 |')
            lines.append('|---:|---:|---|')
            for conn in hmbc_conns:
                lines.append(f"| {conn['h_ppm']:.3f} | {conn['c_ppm']:.3f} "
                             f"| {conn.get('bond_type', '?')} |")
            lines.append('')
        # 季碳推断
        cq_shifts = interp.get('cq_shifts', [])
        if cq_shifts:
            lines.append(f'### 步骤 4: 推断季碳 ({len(cq_shifts)} 个)')
            lines.append('')
            lines.append('δC: ' + ', '.join(f'{c:.2f}' for c in cq_shifts))
            lines.append('')
        # 重叠峰分解 (步骤 1b)
        overlap_regions = interp.get('overlap_regions', [])
        if overlap_regions:
            n_resolved = sum(1 for r in overlap_regions if r.get('status') == 'resolved')
            lines.append(f'### 步骤 1b: 重叠峰分解 ({len(overlap_regions)} 个区域, '
                        f'已分解 {n_resolved})')
            lines.append('')
            lines.append('| 区域 | δH | 总积分 | ¹³C切片峰数 | ¹³C峰位置 | 分解H/峰 | 状态 |')
            lines.append('|---|---:|---:|---:|---|---:|---|')
            for i, r in enumerate(overlap_regions, 1):
                c_peaks_str = ', '.join(f"{c:.1f}" for c in r.get('c_peaks', []))
                n_pk = r.get('n_peaks_in_13c_slice', 0)
                h_per = f"{r['integral_total']/n_pk:.1f}" if n_pk > 0 else '—'
                status_map = {'resolved': '✔ 已分解', 'partial': '△ 部分',
                              'unresolvable': '✘ 未解'}
                st = status_map.get(r.get('status', ''), r.get('status', ''))
                lines.append(f"| {i} | {r['h_ppm_center']:.3f} "
                            f"| {r['integral_total']:.1f} | {n_pk} "
                            f"| {c_peaks_str} | {h_per} | {st} |")
            lines.append('')

    # 全归属表 (Top-1)
    fa = meta.get('full_assignment')
    if fa:
        lines.append('## Top-1 完整 NMR 归属表')
        lines.append('')
        lines.append(f"结构: `{fa.get('structure_smiles', '?')}`  "
                     f"置信度: {fa.get('confidence', 0):.0%}")
        lines.append('')
        lines.append('| 位置 | δC | 级数 | δH | HSQC | COSY | HMBC |')
        lines.append('|---|---:|---|---:|:-:|---|---|')
        for r in fa.get('rows', []):
            hsqc_icon = '✔' if r.get('hsqc_match') else '—'
            h_str = f"{r['h_ppm']:.3f}" if r.get('h_ppm') is not None else '—'
            cosy_str = ', '.join(r.get('cosy_partners', [])) or '—'
            hmbc_str = ', '.join(r.get('hmbc_to_c', [])) or '—'
            lines.append(f"| {r['atom_label']} | {r['c_ppm']:.2f} "
                         f"| {r['multiplicity']} | {h_str} | {hsqc_icon} "
                         f"| {cosy_str} | {hmbc_str} |")
        lines.append('')
        unassigned_c = fa.get('unassigned_obs_c', [])
        unassigned_h = fa.get('unassigned_obs_h', [])
        if unassigned_c or unassigned_h:
            lines.append(f"未归属观测峰: 13C={unassigned_c}, 1H={unassigned_h}")
            lines.append('')

    # 未匹配峰
    lines.append('## 重点关注: 未匹配的实测峰 (可能是新结构片段)')
    lines.append('')
    all_unmatched_c = sorted({round(s, 2)
                              for e in evidences[:top_n]
                              for s in e.unmatched_carbon_shifts})
    all_unmatched_h = sorted({round(s, 2)
                              for e in evidences[:top_n]
                              for s in e.unmatched_proton_shifts})
    if all_unmatched_c:
        lines.append(f'### 未匹配 13C 峰 (ppm)')
        lines.append(', '.join(f'{s:.2f}' for s in all_unmatched_c))
        lines.append('')
    else:
        lines.append('### 未匹配 13C 峰')
        lines.append('无 (top-N 候选都已覆盖所有实测 13C 峰)')
        lines.append('')
    if all_unmatched_h:
        lines.append(f'### 未匹配 1H 峰 (ppm)')
        lines.append(', '.join(f'{s:.2f}' for s in all_unmatched_h))
        lines.append('')
    else:
        lines.append('### 未匹配 1H 峰')
        lines.append('无 (top-N 候选都已覆盖所有实测 1H 峰)')
        lines.append('')

    # Top-N 候选
    lines.append('## Top-N 候选')
    lines.append('')
    lines.append('| Rank | Score | SMILES | nmrshiftdb2 ID | 13C cov | 1H cov | DEPT full | H 一致 | 冲突 |')
    lines.append('|---:|---:|---|---:|---:|---:|---:|---:|---|')
    for e in evidences[:top_n]:
        dept_full = ''
        if e.dept_evidence:
            dept_full = '%d/%d' % (
                sum(1 for d in e.dept_evidence if d.match_type == 'full'),
                len(e.dept_evidence))
        h_count_match = ''
        if e.h_count_evidence:
            h_count_match = '%d/%d' % (e.h_count_evidence.lib_heavy_h, e.h_count_evidence.expected)
        lines.append(
            f'| {e.rank} | {e.score:.1f} | `{e.smiles[:40]}` | {e.nmrshiftdb2_id or "-"} '
            f'| {sum(1 for p in e.carbon_evidence if p.matched)}/{len(e.carbon_evidence)} '
            f'| {sum(1 for p in e.proton_evidence if p.matched)}/{len(e.proton_evidence)} '
            f'| {dept_full} | {h_count_match} '
            f'| {"; ".join(e.conflicts) if e.conflicts else "-"} |'
        )
    lines.append('')

    # 详细
    for e in evidences[:top_n]:
        lines.append(f'### Rank {e.rank}: `{e.smiles}`')
        lines.append('')
        lines.append(f'- 分数: {e.score:.1f}')
        lines.append(f'- 分子式: {e.formula}')
        lines.append(f'- nmrshiftdb2 ID: {e.nmrshiftdb2_id}')
        lines.append(f'- InChIKey: {e.inchikey}')
        lines.append(f'- 库溶剂: {e.solvent_lib}')
        lines.append(f'- 冲突: {"; ".join(e.conflicts) if e.conflicts else "无"}')
        if e.h_count_evidence:
            hce = e.h_count_evidence
            lines.append(f'- H 总数: 实测={hce.expected} 候选={hce.lib_heavy_h} Δ={hce.delta:+d} ({hce.match_score:.0%})')
            for n in hce.notes:
                lines.append(f'  - {n}')
        lines.append('')
        if e.dept_evidence:
            lines.append('**DEPT 碳类型匹配**')
            lines.append('')
            lines.append('| 候选 C (ppm) | 实测 C (ppm) | Δ | 实测 mult | 候选 expected | 一致 |')
            lines.append('|---:|---:|---:|---|---|:-:|')
            for d in e.dept_evidence:
                l = '' if d.observed_carbon_ppm is None else f'{d.observed_carbon_ppm:.2f}'
                dt = '' if d.delta_ppm is None else f'{d.delta_ppm:+.2f}'
                ok = 'OK' if d.mult_match else '**不匹配**'
                lines.append(f'| {d.candidate_carbon_ppm:.2f} | {l} | {dt} | {d.observed_mult or "-"} | {d.expected_mult} | {ok} |')
            lines.append('')
        if e.carbon_evidence:
            lines.append('**13C 匹配**')
            lines.append('')
            lines.append('| Obs (ppm) | Lib (ppm) | Δ | 状态 |')
            lines.append('|---:|---:|---:|---|')
            for p in e.carbon_evidence:
                l = '' if p.closest_lib_shift is None else f'{p.closest_lib_shift:.2f}'
                d = '' if p.delta_ppm is None else f'{p.delta_ppm:+.2f}'
                status = 'OK' if p.matched else '**未匹配**'
                lines.append(f'| {p.observed_shift:.2f} | {l} | {d} | {status} |')
            lines.append('')
        if e.proton_evidence:
            lines.append('**1H 匹配**')
            lines.append('')
            lines.append('| Obs (ppm) | Lib (ppm) | Δ | 状态 |')
            lines.append('|---:|---:|---:|---|')
            for p in e.proton_evidence:
                l = '' if p.closest_lib_shift is None else f'{p.closest_lib_shift:.2f}'
                d = '' if p.delta_ppm is None else f'{p.delta_ppm:+.2f}'
                status = 'OK' if p.matched else '**未匹配**'
                lines.append(f'| {p.observed_shift:.2f} | {l} | {d} | {status} |')
            lines.append('')
        if e.notes:
            lines.append('**备注**')
            lines.append('')
            for n in e.notes:
                lines.append(f'- {n}')
            lines.append('')

    # pandas 候选表 (按 InChIKey 去重, 便于化学家筛选/排序; tabulate 未装时静默退回手写表)
    if evidences[:top_n]:
        try:
            df = evidences_to_dataframe(evidences[:top_n])
            df = dedupe_by_inchikey(df, keep='first')
            if not df.empty:
                lines.append('## pandas 候选表 (按 InChIKey 去重)')
                lines.append('')
                md_table = candidates_to_markdown(df, columns=[
                    'rank', 'score', 'formula', 'inchikey', 'solvent_lib',
                    'heavy_atoms', 'c_match', 'h_match', 'dept_full',
                    'h_count', 'conflicts'])
                lines.append(md_table)
                lines.append('')
        except Exception:
            pass  # pandas/tabulate 未装或出锺, 现有手工表已够用

    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


# ===================== HTML 报告 =====================

# 内联 CSS (中文友好, 表格清晰, 移动端可读)
_HTML_CSS = """
* { box-sizing: border-box; }
body {
  font-family: "Microsoft YaHei", "PingFang SC", "Hiragino Sans GB",
               "Source Han Sans CN", "Noto Sans CJK SC", "WenQuanYi Micro Hei",
               -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
  margin: 0; padding: 24px;
  background: #f5f7fa; color: #1f2937; line-height: 1.6;
}
h1 { color: #1e40af; font-size: 28px; margin: 0 0 8px 0; border-bottom: 3px solid #1e40af; padding-bottom: 12px; }
h2 { color: #1e3a8a; font-size: 22px; margin: 32px 0 16px 0;
     border-left: 5px solid #3b82f6; padding-left: 12px; }
h3 { color: #1e40af; font-size: 18px; margin: 24px 0 12px 0; }
h4 { color: #374151; font-size: 16px; margin: 16px 0 8px 0; }
.container { max-width: 1200px; margin: 0 auto; background: #fff;
             padding: 32px 40px; border-radius: 8px;
             box-shadow: 0 2px 8px rgba(0,0,0,0.05); }
.meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
             gap: 12px 24px; background: #f0f4f8; padding: 16px 20px;
             border-radius: 6px; margin: 16px 0; }
.meta-item { font-size: 14px; }
.meta-item .label { color: #6b7280; font-weight: 500; }
.meta-item .value { color: #111827; font-weight: 600; margin-left: 6px; }
.alert { background: #fef3c7; border-left: 4px solid #f59e0b;
         padding: 12px 16px; margin: 16px 0; border-radius: 4px; }
.alert-danger { background: #fee2e2; border-left: 4px solid #dc2626; }
.alert-info { background: #dbeafe; border-left: 4px solid #3b82f6; }
table { border-collapse: collapse; width: 100%; margin: 12px 0;
        background: #fff; font-size: 14px; }
th { background: #1e40af; color: #fff; padding: 10px 12px;
     text-align: left; font-weight: 600; }
td { padding: 8px 12px; border-bottom: 1px solid #e5e7eb; }
tr:hover td { background: #f9fafb; }
tr.unmatched td { background: #fef2f2; }
tr.matched td { background: #f0fdf4; }
tr.partial td { background: #fffbeb; }
.mol-img { max-width: 100%; height: auto; display: block;
           border: 1px solid #e5e7eb; border-radius: 4px; }
.mol-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 12px; margin: 16px 0; }
.mol-card { background: #f9fafb; border: 1px solid #e5e7eb;
            border-radius: 6px; padding: 12px; text-align: center; }
.mol-card .rank { display: inline-block; background: #1e40af; color: #fff;
                  width: 28px; height: 28px; line-height: 28px; border-radius: 50%;
                  font-weight: 600; margin-bottom: 8px; }
.mol-card .smiles { font-family: "Courier New", monospace; font-size: 12px;
                    color: #374151; word-break: break-all; margin-top: 8px; }
.mol-card .score { color: #059669; font-weight: 600; margin-top: 4px; }
.flow { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin: 16px 0; }
.flow-step { background: #dbeafe; color: #1e40af; padding: 8px 16px;
             border-radius: 20px; font-size: 14px; font-weight: 500; }
.flow-arrow { color: #6b7280; font-size: 18px; }
.dept-bar { display: inline-block; padding: 3px 8px; border-radius: 3px;
            font-weight: 600; font-size: 12px; }
.dept-CH { background: #dbeafe; color: #1e40af; }
.dept-CH2 { background: #fed7aa; color: #9a3412; }
.dept-CH3 { background: #d1fae5; color: #065f46; }
.dept-Cq { background: #e9d5ff; color: #6b21a8; }
.dept-unknown { background: #f3f4f6; color: #6b7280; }
.dept-CH2/CH3 { background: #fde68a; color: #92400e; }
.dept-mismatch { background: #fee2e2; color: #991b1b; font-weight: 700; }
.section-card { background: #fafbfc; border: 1px solid #e5e7eb;
               border-radius: 6px; padding: 16px 20px; margin: 16px 0; }
footer { margin-top: 48px; padding-top: 16px; border-top: 1px solid #e5e7eb;
         color: #6b7280; font-size: 12px; text-align: center; }
"""


def _img_to_base64(img_path: str) -> Optional[str]:
    """读图片为 base64 data URI (HTML 内嵌用)."""
    if not img_path or not os.path.exists(img_path):
        return None
    try:
        with open(img_path, 'rb') as f:
            data = f.read()
        b64 = base64.b64encode(data).decode('ascii')
        ext = os.path.splitext(img_path)[1].lstrip('.').lower()
        mime = 'image/png' if ext == 'png' else f'image/{ext}'
        return f'data:{mime};base64,{b64}'
    except Exception:
        return None


def _render_html_escaped(s: str) -> str:
    """HTML 转义."""
    return (s.replace('&', '&amp;')
             .replace('<', '&lt;')
             .replace('>', '&gt;')
             .replace('"', '&quot;'))


def _render_dept_pill(mult: str, sign: str = 'unknown') -> str:
    """DEPT 标签 pill."""
    cls = f'dept-{mult}'
    icon = {'positive': '↑', 'negative': '↓', 'unknown': '?'}.get(sign, '?')
    return f'<span class="dept-bar {cls}">{_render_html_escaped(mult)} {icon}</span>'


def _render_overview_chart(observed_c: List[float], observed_h: List[float],
                           obs_dept_summary: List[dict], out_path: str) -> Optional[str]:
    """渲染一张谱图概览图 (实测 13C + 1H + DEPT 标注) 到 PNG, 返回 PNG 路径."""
    try:
        fig, axes = plt.subplots(2, 1, figsize=(12, 5), gridspec_kw={'height_ratios': [1, 1]})
        # 13C
        ax_c = axes[0]
        if observed_c:
            for c in observed_c:
                ax_c.axvline(c, color='#1e40af', alpha=0.5, linewidth=0.8)
        # DEPT 颜色覆盖
        if obs_dept_summary:
            for d in obs_dept_summary:
                color_map = {'CH': '#3b82f6', 'CH2': '#f59e0b', 'CH3': '#10b981',
                             'Cq': '#a855f7', 'CH2/CH3': '#fbbf24'}
                c = color_map.get(d.get('mult', '?'), '#6b7280')
                ax_c.axvline(d['ppm'], color=c, alpha=0.9, linewidth=2.0)
        ax_c.set_xlim(max(observed_c) + 20 if observed_c else 250, min(observed_c) - 20 if observed_c else 0)
        ax_c.set_xlabel('13C chemical shift (ppm)')
        ax_c.set_title('13C + DEPT inferred carbon type (color)', fontsize=12)
        ax_c.grid(True, alpha=0.3)
        # 1H
        ax_h = axes[1]
        if observed_h:
            for h in observed_h:
                ax_h.axvline(h, color='#059669', alpha=0.5, linewidth=0.8)
        ax_h.set_xlim(max(observed_h) + 2 if observed_h else 12, min(observed_h) - 2 if observed_h else -1)
        ax_h.set_xlabel('1H chemical shift (ppm)')
        ax_h.set_title('1H peaks', fontsize=12)
        ax_h.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_path, dpi=120, bbox_inches='tight')
        plt.close(fig)
        return out_path
    except Exception:
        return None


def write_report_html(
    path: str,
    evidences: List[CandidateEvidence],
    meta: Dict,
    bridge_info: Optional[Dict] = None,
    top_n: int = 10,
) -> None:
    """写出 HTML 报告 (人类阅读, 含流程/结果/结构图/未匹配峰/DEPT 推断)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    out_dir = os.path.dirname(os.path.abspath(path))
    base = os.path.splitext(os.path.basename(path))[0]
    report_basename = base  # 供 2.5 节引用

    # 渲染额外的概览图
    overview_path = os.path.join(out_dir, f'{base}_overview.png')
    rendered_overview = _render_overview_chart(
        meta.get('observed_carbon_shifts', []),
        meta.get('observed_proton_shifts', []),
        meta.get('obs_dept_summary', []),
        overview_path,
    )
    overview_b64 = _img_to_base64(rendered_overview) if rendered_overview else None

    # 复用已有的 PNG: 推断路径
    topn_png = os.path.join(out_dir, f'{base}_topn.png')
    overlay_png = os.path.join(out_dir, f'{base}_overlay.png')
    topn_b64 = _img_to_base64(topn_png)
    overlay_b64 = _img_to_base64(overlay_png)
    rank_pngs: dict[int, str] = {}
    for e in evidences[:top_n]:
        rp = os.path.join(out_dir, f'{base}_rank{e.rank}.png')
        b = _img_to_base64(rp)
        if b:
            rank_pngs[e.rank] = b

    # seaborn 可视化 (可选, 未装时三项均为 None)
    heatmap_png = os.path.join(out_dir, f'{base}_heatmap.png')
    scatter_png = os.path.join(out_dir, f'{base}_scatter.png')
    matrix_png = os.path.join(out_dir, f'{base}_matrix.png')
    rendered_heatmap = plot_candidate_score_heatmap(evidences, heatmap_png, top_n) if _HAS_SEABORN else None
    rendered_scatter = plot_ppm_comparison_scatter(
        meta.get('observed_carbon_shifts', []),
        meta.get('observed_proton_shifts', []),
        evidences, scatter_png, top_n) if _HAS_SEABORN else None
    rendered_matrix = plot_constraint_match_matrix(evidences, matrix_png, top_n) if _HAS_SEABORN else None
    heatmap_b64 = _img_to_base64(rendered_heatmap) if rendered_heatmap else None
    scatter_b64 = _img_to_base64(rendered_scatter) if rendered_scatter else None
    matrix_b64 = _img_to_base64(rendered_matrix) if rendered_matrix else None

    # 未匹配峰汇总
    all_unmatched_c = sorted({round(s, 2) for e in evidences[:top_n] for s in e.unmatched_carbon_shifts})
    all_unmatched_h = sorted({round(s, 2) for e in evidences[:top_n] for s in e.unmatched_proton_shifts})

    # DEPT 推断表
    dept_counts = meta.get('dept_counts', {})
    obs_dept = meta.get('obs_dept_summary', [])

    # 渲染 HTML
    parts: List[str] = []
    parts.append('<!DOCTYPE html>')
    parts.append('<html lang="zh-CN">')
    parts.append('<head>')
    parts.append('<meta charset="UTF-8">')
    parts.append(f'<title>NMR 反推结构报告 - {os.path.basename(meta.get("folder", ""))}</title>')
    parts.append(f'<style>{_HTML_CSS}</style>')
    parts.append('</head>')
    parts.append('<body>')
    parts.append('<div class="container">')

    # 标题
    sample_name = os.path.basename(meta.get('folder', '')) or '未命名样本'
    parts.append(f'<h1>🔬 NMR 反推结构报告 — {sample_name}</h1>')
    parts.append(f'<div style="color:#6b7280;font-size:14px;">生成时间: {datetime.now().isoformat(timespec="seconds")}</div>')

    # 1. 元信息
    parts.append('<h2>1. 谱图概览</h2>')
    parts.append('<div class="meta-grid">')
    parts.append(f'<div class="meta-item"><span class="label">样本目录</span><span class="value">{_render_html_escaped(meta.get("folder", "-"))}</span></div>')
    parts.append(f'<div class="meta-item"><span class="label">溶剂</span><span class="value">{_render_html_escaped(str(meta.get("solvent_observed", "-")))}</span></div>')
    parts.append(f'<div class="meta-item"><span class="label">频率</span><span class="value">{meta.get("field_mhz", "-")} MHz</span></div>')
    parts.append(f'<div class="meta-item"><span class="label">13C 容差</span><span class="value">±{meta.get("c_tol", "-")} ppm</span></div>')
    parts.append(f'<div class="meta-item"><span class="label">1H 容差</span><span class="value">±{meta.get("h_tol", "-")} ppm</span></div>')
    parts.append(f'<div class="meta-item"><span class="label">候选库</span><span class="value">nmrshiftdb2</span></div>')
    if meta.get('formula_constraint'):
        parts.append(f'<div class="meta-item"><span class="label">分子式约束</span><span class="value">{_render_html_escaped(meta["formula_constraint"])}</span></div>')
    parts.append('</div>')

    # 概览图
    if overview_b64:
        parts.append(f'<img class="mol-img" src="{overview_b64}" alt="谱图概览">')
    else:
        parts.append('<div class="alert alert-info">谱图概览图未生成 (无数据)。</div>')

    # 2. 解谱过程
    parts.append('<h2>2. 解谱过程</h2>')
    parts.append('<div class="flow">')
    parts.append('<div class="flow-step">1. 读 Bruker 谱图 (peaklist.xml + parm.txt)</div>')
    parts.append('<span class="flow-arrow">→</span>')
    parts.append('<div class="flow-step">2. DEPT-90 / DEPT-135 推断碳类型</div>')
    parts.append('<span class="flow-arrow">→</span>')
    parts.append('<div class="flow-step">3. 强度过滤去噪 (top 20%)</div>')
    parts.append('<span class="flow-arrow">→</span>')
    parts.append('<div class="flow-step">4. nmrshiftdb2 候选检索 (58187 库)</div>')
    parts.append('<span class="flow-arrow">→</span>')
    parts.append('<div class="flow-step">5. 综合打分 (cov + RMS + DEPT + H_count)</div>')
    parts.append('<span class="flow-arrow">→</span>')
    parts.append('<div class="flow-step">6. Top-N 排序 + 证据链</div>')
    parts.append('</div>')

    # 2.5 2D NMR 解谱详细过程 (化学家工作流)
    interp = meta.get('interpretation_2d')
    if interp:
        parts.append('<h2>2.5  2D NMR 解谱过程 (化学家工作流)</h2>')
        parts.append('<div class="alert alert-info">按 MestReNova 5.4 节化学家顺序: '
                     'HSQC → COSY → HMBC → 季碳推断。详细表格与等高线图见各自中间文件 '
                     '(本节仅给路径, 不嵌入)。</div>')

        # 汇总数据 (后续不再嵌入单独表格/图片)
        ch_pairs = interp.get('ch_pairs', [])
        overlap_regions = interp.get('overlap_regions', [])
        spin_systems = interp.get('spin_systems', [])
        hmbc_conns = interp.get('hmbc_connections', [])
        cq_shifts = interp.get('cq_shifts', [])

        n_with_int = sum(1 for p in ch_pairs if p.get('h_integral') is not None)
        n_overlap = sum(1 for p in ch_pairs
                        if p.get('h_integral') is not None and p['h_integral'] >= 4.0)
        if overlap_regions:
            n_resolved = sum(1 for r in overlap_regions if r.get('status') == 'resolved')
            n_partial = sum(1 for r in overlap_regions if r.get('status') == 'partial')
            overlap_summary = (f'{len(overlap_regions)} 区 '
                               f'(已分解 {n_resolved}, 部分 {n_partial})')
        else:
            overlap_summary = '—'

        # 路径索引表
        parts.append('<table>')
        parts.append('<thead><tr><th>步骤</th><th>内容</th>'
                     '<th>产出数</th><th>中间结果文件/目录</th></tr></thead>')
        parts.append('<tbody>')

        # 步骤 1
        parts.append('<tr>'
                     '<td><strong>步骤 1</strong></td>'
                     '<td>HSQC C-H 直接连接对</td>'
                     f'<td>{len(ch_pairs)} 个 '
                     f'(含积分 {n_with_int}, 重叠 {n_overlap})</td>'
                     '<td>'
                     '<code>step1_hsqc_ch_pairs.csv</code><br>'
                     '<code>step1_hsqc_ch_pairs/normal/</code> (非重叠峰等高线)<br>'
                     '<code>step1_hsqc_ch_pairs/overlap/</code> (重叠峰等高线)'
                     '</td></tr>')

        # 步骤 1b
        parts.append('<tr>'
                     '<td><strong>步骤 1b</strong></td>'
                     '<td>重叠峰分解 (5.9 节流程)</td>'
                     f'<td>{overlap_summary}</td>'
                     '<td>'
                     '<code>step1b_overlap_analysis.csv</code> '
                     '(13C 切片检峰详情)<br>'
                     '<code>step1_hsqc_ch_pairs/overlap/</code> '
                     '(局部放大含切片)'
                     '</td></tr>')

        # 步骤 2 / 3 / 4
        parts.append('<tr>'
                     '<td><strong>步骤 2</strong></td>'
                     '<td>COSY 自旋系统</td>'
                     f'<td>{len(spin_systems)} 个系统</td>'
                     '<td><code>step2_cosy_spin_systems.csv</code></td></tr>')
        parts.append('<tr>'
                     '<td><strong>步骤 3</strong></td>'
                     '<td>HMBC 远程连接</td>'
                     f'<td>{len(hmbc_conns)} 个连接</td>'
                     '<td><code>step3_hmbc_connections.csv</code></td></tr>')
        parts.append('<tr>'
                     '<td><strong>步骤 4</strong></td>'
                     '<td>季碳推断 (含 DEPT 证据)</td>'
                     f'<td>{len(cq_shifts)} 个</td>'
                     '<td><code>step4_quaternary_carbons.csv</code></td></tr>')

        # 步骤 5: NOE 空间近邻 (若有)
        noe_contacts = interp.get('noe_contacts', [])
        if noe_contacts:
            parts.append('<tr>'
                         '<td><strong>步骤 5</strong></td>'
                         '<td>NOE 空间近邻 H-H 对</td>'
                         f'<td>{len(noe_contacts)} 对</td>'
                         '<td><code>step5_noe_spatial_contacts.csv</code></td></tr>')

        # 步骤 6: NOE 3D 验证 (若有)
        noe_verif = meta.get('noe_verification', [])
        if noe_verif:
            parts.append('<tr>'
                         '<td><strong>步骤 6</strong></td>'
                         '<td>NOE 3D 构象验证 (Top-N)</td>'
                         f'<td>{len(noe_verif)} 个候选</td>'
                         '<td><code>step6_noe_3d_verification.csv</code></td></tr>')

        # Part 2.C: 步骤 6b: Karplus J 3D 二面角验证 (若有)
        j_verif = meta.get('j_verification', [])
        if j_verif:
            parts.append('<tr>'
                         '<td><strong>步骤 6b</strong></td>'
                         '<td>Karplus J 3D 二面角验证 (Top-N)</td>'
                         f'<td>{len(j_verif)} 个候选</td>'
                         '<td><code>step6b_j_karplus_3d_verification.csv</code></td></tr>')

        # 步骤 7: 归属完整性核对
        av = interp.get('assignment_validation')
        if av:
            av_status = 'OK' if av.get('is_complete') else '有警告'
            parts.append('<tr>'
                         '<td><strong>步骤 7</strong></td>'
                         '<td>归属完整性核对</td>'
                         f'<td>C: {av.get("assigned_c", 0)}/{av.get("total_observed_c", 0)}; '
                         f'H: {av.get("assigned_h", 0)}/{av.get("total_observed_h", 0)}; '
                         f'{av_status}</td>'
                         '<td><code>step7_assignment_validation.csv</code></td></tr>')

        # Top-1 完整归属表
        fa = meta.get('full_assignment')
        # 统计生成的归属表 CSV (论文 Table 1 格式) 数量
        assignment_dir = os.path.join(out_dir, '完整NMR归属表(论文Table1)')
        n_assignment_csv = 0
        if os.path.isdir(assignment_dir):
            n_assignment_csv = len([f for f in os.listdir(assignment_dir) if f.endswith('.csv')])

        if fa:
            fa_info = (f'<code>{_render_html_escaped(fa.get("structure_smiles", "?"))}</code> '
                       f'置信度 <strong>{fa.get("confidence", 0):.0%}</strong> '
                       f'({len(fa.get("rows", []))} 行)')
        else:
            fa_info = '—'
        parts.append('<tr>'
                     '<td><strong>Top-1</strong></td>'
                     '<td>完整 NMR 归属表 (论文 Table 1)</td>'
                     f'<td>{fa_info}</td>'
                     f'<td><code>完整NMR归属表(论文Table1)/</code> '
                     f'子目录 ({n_assignment_csv} 个 CSV, '
                     f'每个 Top-1 CSV 含 位置/δC/级数/δH/HSQC/COSY/HMBC 列)</td></tr>')

        parts.append('</tbody></table>')
        parts.append('<div class="alert alert-info" style="margin-top:6px;">'
                     '<small>所有中间结果 (CSV/PNG) 与 HTML 报告同级目录;'
                     ' 化学家可在 Excel 或图片查看器中分别打开复核。</small></div>')
    # 3. DEPT 推断
    parts.append('<h2>3. DEPT 推断 (实测碳类型)</h2>')
    if not obs_dept:
        parts.append('<div class="alert">无 DEPT 数据, 跳过碳类型推断。</div>')
    else:
        # 统计
        parts.append('<div class="section-card">')
        parts.append('<h4>碳类型分布</h4>')
        parts.append('<div class="flow">')
        for k in ('CH', 'CH2', 'CH3', 'CH2/CH3', 'Cq'):
            v = dept_counts.get(k, 0)
            if v > 0:
                parts.append(f'<div>共 <strong>{v}</strong> 个 {_render_dept_pill(k)}</div>')
        parts.append('</div>')
        parts.append('</div>')
        # 详细表
        parts.append('<table>')
        parts.append('<thead><tr><th>13C (ppm)</th><th>碳类型</th><th>方向</th><th>说明</th></tr></thead>')
        parts.append('<tbody>')
        for d in obs_dept:
            sign = d.get('sign', 'unknown')
            mult = d.get('mult', '?')
            note_map = {
                'CH': '三级碳 (一氢), 芳/烯 CH',
                'CH2': '亚甲基 (二氢), DEPT-135 负峰',
                'CH3': '甲基 (三氢), DEPT-135 正峰 (化学位移 < 35 ppm)',
                'Cq': '季碳 (无氢), 排除 CH/CH2/CH3',
                'CH2/CH3': '亚甲基/甲基 不能直接区分 (需查 1r 方向)',
            }
            note = note_map.get(mult, '')
            parts.append(f'<tr><td><strong>{d["ppm"]:.3f}</strong></td>'
                         f'<td>{_render_dept_pill(mult, sign)}</td>'
                         f'<td>{sign}</td><td>{_render_html_escaped(note)}</td></tr>')
        parts.append('</tbody></table>')
        # H 总数预期
        h_expected = sum({'CH': 1, 'CH2': 2, 'CH3': 3}.get(d.get('mult'), 0) for d in obs_dept)
        parts.append(f'<div class="alert alert-info" style="margin-top:12px;">'
                     f'<strong>从 DEPT 推断的 H 总数 (非活泼 H):</strong> {h_expected}'
                     f'<br><small>公式: CH×1 + CH2×2 + CH3×3; Cq 计入 0; CH2/CH3 占位不计入</small></div>')

    # 3.5 2D NMR 关联验证
    if meta.get('used_2d') and meta.get('2d_experiments'):
        parts.append('<h2>3.5 2D NMR 关联验证</h2>')
        parts.append(f'<div class="alert alert-info">检测到的 2D 实验: '
                     f'<strong>{", ".join(meta["2d_experiments"])}</strong></div>')
        corr_2d = meta.get('correlation_2d', {})
        if corr_2d:
            parts.append('<table>')
            parts.append('<thead><tr><th>Rank</th><th>SMILES</th>'
                         '<th>HSQC</th><th>HMBC</th><th>COSY</th><th>综合</th></tr></thead>')
            parts.append('<tbody>')
            for e in evidences[:top_n]:
                s2d = corr_2d.get(e.smiles)
                if s2d:
                    parts.append(
                        f'<tr><td><strong>#{e.rank}</strong></td>'
                        f'<td><code>{_render_html_escaped(e.smiles[:40])}</code></td>'
                        f'<td>{s2d["hsqc_match_rate"]:.0%}</td>'
                        f'<td>{s2d["hmbc_match_rate"]:.0%}</td>'
                        f'<td>{s2d["cosy_match_rate"]:.0%}</td>'
                        f'<td><strong>{s2d["combined_2d_score"]:.0%}</strong></td></tr>'
                    )
            parts.append('</tbody></table>')
            parts.append('<div class="alert alert-info" style="margin-top:8px;">'
                         '<small>HSQC: ¹H-¹³C 直接连接匹配率; '
                         'HMBC: ¹H-¹³C 远程 (2-3键) 匹配率; '
                         'COSY: ¹H-¹H 邻接匹配率</small></div>')

    # 4. Top-N 结构图
    parts.append('<h2>4. Top-{0} 候选结构</h2>'.format(min(top_n, len(evidences))))
    if topn_b64:
        parts.append(f'<img class="mol-img" src="{topn_b64}" alt="Top-N 候选网格">')
    else:
        parts.append('<div class="alert">候选网格图未生成 (结构图绘制失败)。</div>')

    # 5. 候选详细
    parts.append('<h2>5. 候选详细证据</h2>')
    parts.append('<table>')
    parts.append('<thead><tr><th>Rank</th><th>Score</th><th>SMILES</th><th>13C cov</th><th>1H cov</th><th>DEPT</th><th>H</th><th>冲突</th></tr></thead>')
    parts.append('<tbody>')
    for e in evidences[:top_n]:
        c_cov = '%d/%d' % (sum(1 for p in e.carbon_evidence if p.matched), len(e.carbon_evidence))
        h_cov = '%d/%d' % (sum(1 for p in e.proton_evidence if p.matched), len(e.proton_evidence))
        dept_full = '%d/%d' % (
            sum(1 for d in e.dept_evidence if d.match_type == 'full'),
            len(e.dept_evidence)) if e.dept_evidence else '-'
        h_match = '%d/%d' % (e.h_count_evidence.lib_heavy_h, e.h_count_evidence.expected) if e.h_count_evidence else '-'
        conflicts = '; '.join(e.conflicts) if e.conflicts else '-'
        parts.append(f'<tr><td><strong>#{e.rank}</strong></td><td>{e.score:.1f}</td>'
                     f'<td><code>{_render_html_escaped(e.smiles[:50])}</code></td>'
                     f'<td>{c_cov}</td><td>{h_cov}</td><td>{dept_full}</td><td>{h_match}</td>'
                     f'<td>{_render_html_escaped(conflicts)}</td></tr>')
    parts.append('</tbody></table>')

    # 6. 每候选详细
    parts.append('<h2>6. 每候选完整证据</h2>')
    for e in evidences[:top_n]:
        parts.append(f'<div class="section-card">')
        parts.append(f'<h3>#{e.rank} { _render_html_escaped(e.smiles[:60])}</h3>')
        # 基本信息
        parts.append('<div class="meta-grid">')
        parts.append(f'<div class="meta-item"><span class="label">分数</span><span class="value">{e.score:.2f}</span></div>')
        parts.append(f'<div class="meta-item"><span class="label">分子式</span><span class="value">{_render_html_escaped(e.formula or "-")}</span></div>')
        parts.append(f'<div class="meta-item"><span class="label">InChIKey</span><span class="value" style="font-family:monospace;font-size:11px;">{_render_html_escaped(e.inchikey or "-")}</span></div>')
        parts.append(f'<div class="meta-item"><span class="label">重原子</span><span class="value">{e.heavy_atoms}</span></div>')
        parts.append(f'<div class="meta-item"><span class="label">库溶剂</span><span class="value">{_render_html_escaped(e.solvent_lib or "-")}</span></div>')
        parts.append(f'<div class="meta-item"><span class="label">nmrshiftdb2 ID</span><span class="value">{e.nmrshiftdb2_id or "-"}</span></div>')
        parts.append('</div>')
        # 结构图
        if e.rank in rank_pngs:
            parts.append(f'<img class="mol-img" src="{rank_pngs[e.rank]}" alt="Rank {e.rank}">')
        # H 总数
        if e.h_count_evidence:
            hce = e.h_count_evidence
            cls = 'alert-danger' if abs(hce.delta) > 2 else 'alert-info'
            parts.append(f'<div class="alert {cls}"><strong>H 总数一致性:</strong> '
                         f'实测 (DEPT 推断) = {hce.expected}; 候选 (RDKit) = {hce.lib_heavy_h}; '
                         f'Δ = {hce.delta:+d}; 一致性 = {hce.match_score:.0%}'
                         f'<ul style="margin:6px 0 0 0;">')
            for n in hce.notes:
                parts.append(f'<li>{_render_html_escaped(n)}</li>')
            parts.append('</ul></div>')
        # DEPT 匹配
        if e.dept_evidence:
            full_n = sum(1 for d in e.dept_evidence if d.match_type == 'full')
            partial_n = sum(1 for d in e.dept_evidence if d.match_type == 'partial')
            unmatched_n = sum(1 for d in e.dept_evidence if d.match_type == 'unmatched')
            parts.append(f'<h4>DEPT 碳类型匹配 (full={full_n}, partial={partial_n}, unmatched={unmatched_n})</h4>')
            parts.append('<table>')
            parts.append('<thead><tr><th>候选 C (ppm)</th><th>实测 C (ppm)</th><th>Δ</th><th>实测 mult</th><th>候选 expected</th><th>状态</th></tr></thead>')
            parts.append('<tbody>')
            for d in e.dept_evidence:
                l = '' if d.observed_carbon_ppm is None else f'{d.observed_carbon_ppm:.2f}'
                dt = '' if d.delta_ppm is None else f'{d.delta_ppm:+.2f}'
                row_cls = 'matched' if d.match_type == 'full' else ('unmatched' if d.match_type == 'unmatched' else 'partial')
                status = 'OK' if d.mult_match else '不匹配'
                parts.append(f'<tr class="{row_cls}"><td>{d.candidate_carbon_ppm:.2f}</td>'
                             f'<td>{l}</td><td>{dt}</td>'
                             f'<td>{_render_html_escaped(d.observed_mult or "-")}</td>'
                             f'<td>{_render_html_escaped(d.expected_mult)}</td>'
                             f'<td>{status}</td></tr>')
            parts.append('</tbody></table>')
        # 13C 匹配
        if e.carbon_evidence:
            parts.append('<h4>13C 峰匹配</h4>')
            parts.append('<table>')
            parts.append('<thead><tr><th>Obs (ppm)</th><th>Lib (ppm)</th><th>Δ</th><th>状态</th></tr></thead>')
            parts.append('<tbody>')
            for p in e.carbon_evidence:
                l = '' if p.closest_lib_shift is None else f'{p.closest_lib_shift:.2f}'
                dt = '' if p.delta_ppm is None else f'{p.delta_ppm:+.2f}'
                row_cls = '' if p.matched else 'unmatched'
                status = 'OK' if p.matched else '**未匹配**'
                parts.append(f'<tr class="{row_cls}"><td>{p.observed_shift:.2f}</td>'
                             f'<td>{l}</td><td>{dt}</td><td>{status}</td></tr>')
            parts.append('</tbody></table>')
        # 1H 匹配
        if e.proton_evidence:
            parts.append('<h4>1H 峰匹配</h4>')
            parts.append('<table>')
            parts.append('<thead><tr><th>Obs (ppm)</th><th>Lib (ppm)</th><th>Δ</th><th>状态</th></tr></thead>')
            parts.append('<tbody>')
            for p in e.proton_evidence:
                l = '' if p.closest_lib_shift is None else f'{p.closest_lib_shift:.2f}'
                dt = '' if p.delta_ppm is None else f'{p.delta_ppm:+.2f}'
                row_cls = '' if p.matched else 'unmatched'
                status = 'OK' if p.matched else '**未匹配**'
                parts.append(f'<tr class="{row_cls}"><td>{p.observed_shift:.2f}</td>'
                             f'<td>{l}</td><td>{dt}</td><td>{status}</td></tr>')
            parts.append('</tbody></table>')
        # 备注
        if e.notes:
            parts.append('<h4>备注</h4><ul>')
            for n in e.notes:
                parts.append(f'<li>{_render_html_escaped(n)}</li>')
            parts.append('</ul>')
        parts.append('</div>')

    # 7. 未匹配峰 (重点关注)
    parts.append('<h2>7. ⚠️ 重点关注: 未匹配峰 (可能的新结构片段)</h2>')
    parts.append('<div class="alert alert-danger">'
                 '以下实测峰在所有 Top-N 候选中均未被解释。需要 ML 模型推断或化学家人工复核。'
                 '</div>')
    if all_unmatched_c:
        parts.append(f'<h4>未匹配 13C 峰 (ppm): {", ".join(f"{s:.2f}" for s in all_unmatched_c)}</h4>')
    else:
        parts.append('<h4>未匹配 13C 峰: 无 (top-N 候选已覆盖所有实测 13C 峰)</h4>')
    if all_unmatched_h:
        parts.append(f'<h4>未匹配 1H 峰 (ppm): {", ".join(f"{s:.2f}" for s in all_unmatched_h)}</h4>')
    else:
        parts.append('<h4>未匹配 1H 峰: 无 (top-N 候选已覆盖所有实测 1H 峰)</h4>')

    # 8. seaborn 可视化 (候选评分热力图 / ppm 对比散点 / 约束矩阵)
    parts.append('<h2>8. 📊 候选评估可视化 (seaborn)</h2>')
    if not any((heatmap_b64, scatter_b64, matrix_b64)):
        parts.append('<div class="alert alert-info">seaborn 未装或数据不足, 跳过统计可视化。</div>')
    else:
        if heatmap_b64:
            parts.append(f'<h3>候选评分热力图</h3>')
            parts.append(f'<img class="mol-img" src="{heatmap_b64}" alt="评分热力图">')
        if scatter_b64:
            parts.append(f'<h3>实测 vs 候选 ppm 对比散点图</h3>')
            parts.append(f'<img class="mol-img" src="{scatter_b64}" alt="ppm 对比">')
        if matrix_b64:
            parts.append(f'<h3>约束匹配矩阵</h3>')
            parts.append(f'<img class="mol-img" src="{matrix_b64}" alt="约束匹配矩阵">')

    # overlay 对比图
    if overlay_b64:
        parts.append('<h3>实测 vs 候选对比</h3>')
        parts.append(f'<img class="mol-img" src="{overlay_b64}" alt="实测 vs 候选对比">')

    # footer
    parts.append('<footer>')
    parts.append('本报告由 <strong>nmr_elucidator</strong> 自动生成 | '
                 '数据来源: Bruker 原始谱图 + nmrshiftdb2 数据库 + RDKit 化学信息学')
    parts.append('</footer>')
    parts.append('</div></body></html>')

    html = '\n'.join(parts)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(html)
