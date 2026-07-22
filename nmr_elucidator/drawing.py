"""候选结构图与谱图对比可视化.

功能:
1. draw_single_mol: 单分子 2D 结构图 (RDKit 分子骨架 + matplotlib 中文标签)
2. draw_topn_grid: Top-N 候选网格图 (合成方案, 中文安全)
3. draw_unmatched_overlay: 实测 vs Top-1 候选化学位移对比图 (matplotlib)

中文字体兼容策略:
- 分子骨架: RDKit Draw 生成 PNG (无中文, 无乱码风险)
- 中文标签 / 标题 / 图例: matplotlib rcParams = SimHei
- 两阶段合成: PIL 把 RDKit 分子图贴在 matplotlib canvas 下方, 中文字渲染在上面
"""
from __future__ import annotations

import io
import os
import warnings
from typing import List, Optional

import matplotlib

matplotlib.use('Agg')

# ---- matplotlib 中文字体配置 (必须在 import pyplot 之前) ----
import matplotlib.pyplot as plt
from matplotlib import font_manager

_ZH_CANDIDATES = [
    ('SimHei', r'C:\Windows\Fonts\simhei.ttf'),
    ('Microsoft YaHei', r'C:\Windows\Fonts\msyh.ttc'),
    ('SimSun', r'C:\Windows\Fonts\simsun.ttc'),
    ('DejaVu Sans', None),  # 兜底
]
_chinese_font_name: Optional[str] = None
for _name, _path in _ZH_CANDIDATES:
    if _path and os.path.exists(_path):
        try:
            font_manager.fontManager.addfont(_path)
            _chinese_font_name = _name
            break
        except Exception:
            continue

if _chinese_font_name:
    plt.rcParams['font.sans-serif'] = [_chinese_font_name, 'DejaVu Sans']
else:
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 修复负号显示为方块

import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

from .evidence import CandidateEvidence


def _draw_mol_image(mol, size_px: int = 420) -> Optional[Image.Image]:
    """用 RDKit Draw 渲染无文字的分子 PNG, 返回 PIL.Image. None 表示 SMILES 解析失败."""
    if mol is None:
        return None
    try:
        # legend 留空, 由我们在 matplotlib 阶段叠中文
        return Draw.MolToImage(mol, size=(size_px, size_px), legend='')
    except Exception:
        return None


def draw_single_mol(smiles: str, out_path: str, legend: str = '') -> None:
    """单分子 2D 结构图 (中文 legend).

    smiles: 输入 SMILES
    out_path: PNG 文件路径
    legend: 图例文本 (中文安全)
    """
    mol = Chem.MolFromSmiles(smiles)
    mol_img = _draw_mol_image(mol, size_px=420)
    if mol_img is None:
        return

    fig_h = 5.0 if legend else 4.4
    fig, ax = plt.subplots(figsize=(4.5, fig_h))
    ax.axis('off')

    # 把 RDKit 分子图 (PIL Image) 显示在 axes 里
    ax_img = fig.add_axes([0.05, 0.15, 0.9, 0.8]) if legend else fig.add_axes([0.05, 0.05, 0.9, 0.9])
    ax_img.axis('off')
    ax_img.imshow(np.asarray(mol_img))

    if legend:
        ax.text(0.5, 0.04, legend, ha='center', va='bottom',
                fontsize=12, color='black', transform=fig.transFigure)

    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def draw_topn_grid(evidences: List[CandidateEvidence], out_path: str, top_n: int = 10) -> None:
    """Top-N 候选网格图.

    每行 2 个分子. RDKit 渲染结构, matplotlib 排版, 中文安全.
    """
    evidences = evidences[:top_n]
    items = []
    for e in evidences:
        mol = Chem.MolFromSmiles(e.smiles)
        if mol is None:
            continue
        items.append((mol, e))

    if not items:
        return

    n = len(items)
    cols = 2
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(9, 4.5 * rows))
    if rows == 1:
        axes = np.array([axes])
    axes_flat = axes.flatten()

    for i, (mol, e) in enumerate(items):
        ax = axes_flat[i]
        ax.axis('off')
        mol_img = _draw_mol_image(mol, size_px=320)
        if mol_img is not None:
            ax.imshow(np.asarray(mol_img))
        # 中文: rank + score + smiles 摘要
        title = f"#{e.rank}  score={e.score:.1f}\n{e.smiles[:34]}"
        ax.set_title(title, fontsize=11, pad=8)

    # 隐藏多余的 axes
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis('off')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def draw_unmatched_overlay(
    observed_c: List[float],
    observed_h: List[float],
    evidences: List[CandidateEvidence],
    out_path: str,
    top_n: int = 5,
) -> None:
    """实测谱 vs Top-N 候选预测谱的对比图 (中文安全).

    上下两个子图: 13C (上) 和 1H (下).
    """
    evidences = evidences[:top_n]

    fig, (ax_c, ax_h) = plt.subplots(2, 1, figsize=(12, 8))

    # 13C
    if observed_c:
        for shift in observed_c:
            ax_c.axvline(shift, color='red', alpha=0.4, linewidth=0.8)
        ax_c.set_xticks(sorted(set(round(s, 1) for s in observed_c)))
        ax_c.set_xticklabels([f'{s:.1f}' for s in sorted(set(round(s, 1) for s in observed_c))],
                              rotation=45, fontsize=7)
    colors = plt.cm.tab10.colors
    for i, e in enumerate(evidences):
        shifts = e.lib_13c_shifts
        if shifts:
            ys = [-(i + 1)] * len(shifts)
            ax_c.scatter(shifts, ys, color=colors[i % len(colors)], s=20,
                          alpha=0.6, label=f'#{e.rank}')
    ax_c.set_ylabel('Top-N 候选 (向下偏移)')
    ax_c.set_title('13C 化学位移对比 (红=实测, 散点=候选预测)')
    ax_c.legend(fontsize=7, loc='best')
    ax_c.set_yticks([])
    ax_c.invert_xaxis()
    ax_c.grid(True, alpha=0.3)

    # 1H
    if observed_h:
        for shift in observed_h:
            ax_h.axvline(shift, color='red', alpha=0.4, linewidth=0.8)
        ax_h.set_xticks(sorted(set(round(s, 2) for s in observed_h)))
        ax_h.set_xticklabels([f'{s:.2f}' for s in sorted(set(round(s, 2) for s in observed_h))],
                              rotation=45, fontsize=7)
    for i, e in enumerate(evidences):
        shifts = e.lib_1h_shifts
        if shifts:
            ys = [-(i + 1)] * len(shifts)
            ax_h.scatter(shifts, ys, color=colors[i % len(colors)], s=20,
                          alpha=0.6)
    ax_h.set_ylabel('Top-N 候选 (向下偏移)')
    ax_h.set_title('1H 化学位移对比 (红=实测, 散点=候选预测)')
    ax_h.set_yticks([])
    ax_h.invert_xaxis()
    ax_h.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def draw_real_spectrum_overlay(
    specs_processed: dict,
    observed_by_exp: dict,
    out_path: str,
) -> None:
    """真实谱曲线图 (来自 pdata/1/1r) + 各实验 peaklist 峰位标注.

    specs_processed: {exp_no: ProcessedSpectrum} (io_processed.read_processed_folder 的结果),
                     鸭子类型使用 .ppm / .intensity / .nucleus / .exp_no / .sf_mhz.
    observed_by_exp: {exp_no: [峰位 ppm, ...]} 来自 peaklist.xml, 叠加为红色竖线.
    out_path: PNG 路径.

    每个实验一个子图, 画真实连续谱线 (完整线形, 非竹签), 叠加 peaklist 检出峰位.
    这是"用 1r 最终谱图"的核心可视化: 直观展示真实谱线形 + 我们的峰检测落点是否吻合.
    """
    items = [specs_processed[k] for k in sorted(specs_processed.keys())]
    if not items:
        return

    n = len(items)
    fig, axes = plt.subplots(n, 1, figsize=(12, 2.6 * n), squeeze=False)
    for i, ps in enumerate(items):
        ax = axes[i][0]
        ax.plot(ps.ppm, ps.intensity, 'b-', linewidth=0.5,
                label=f'真实谱 1r (n={len(ps.ppm)})')
        obs = observed_by_exp.get(ps.exp_no, [])
        if obs:
            ax.vlines(obs, float(ps.intensity.min()), float(ps.intensity.max()),
                      colors='r', alpha=0.45, linewidth=0.8,
                      label=f'peaklist 检出 ({len(obs)} 峰)')
        ax.set_title(f'exp {ps.exp_no} - {ps.nucleus} - SF={ps.sf_mhz:.2f} MHz')
        ax.set_xlabel('化学位移 (ppm)')
        ax.invert_xaxis()  # ppm 惯例: 高场在右
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def write_unmatched_summary_csv(
    evidences: List[CandidateEvidence],
    out_path: str,
    top_n: int = 10,
) -> None:
    """未匹配峰汇总 CSV (Excel 友好, UTF-8 BOM).

    列: rank, smiles, unmatched_carbon_ppm, unmatched_proton_ppm
    """
    import csv
    evidences = evidences[:top_n]
    with open(out_path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'smiles', 'unmatched_13C_ppm', 'unmatched_1H_ppm',
                    'note'])
        for e in evidences:
            c_str = '; '.join(f'{s:.2f}' for s in e.unmatched_carbon_shifts)
            h_str = '; '.join(f'{s:.2f}' for s in e.unmatched_proton_shifts)
            note = '可能是新结构片段, 留待 ML 模型 + 化学家复核' if (e.unmatched_carbon_shifts or e.unmatched_proton_shifts) else ''
            w.writerow([e.rank, e.smiles, c_str, h_str, note])
