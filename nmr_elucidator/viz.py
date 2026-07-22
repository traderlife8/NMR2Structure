"""seaborn 增强可视化 (候选评分热力图 / ppm 对比散点图 / 约束矩阵).

设计目标:
1. 与现有 drawing.py 的 matplotlib 风格互补 (seaborn 提供统计聚合视图)
2. 每图独立 fig, ax, 保存后关闭
3. 中文字体复用 drawing.py 的 SimHei/Microsoft YaHei 兜底配置
4. 错误显式: seaborn 未装 / 数据为空 静默返回 None, 调用方跳过该图

依赖: seaborn>=0.12 (本项目环境已装 0.13.2).
"""
from __future__ import annotations

import os
import warnings
from typing import List, Optional

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 复用 drawing.py 的中文字体探测 (避免重新走一遍字体加载)
try:
    from matplotlib import font_manager
    _ZH_CANDIDATES = [
        ('SimHei', r'C:\Windows\Fonts\simhei.ttf'),
        ('Microsoft YaHei', r'C:\Windows\Fonts\msyh.ttc'),
        ('SimSun', r'C:\Windows\Fonts\simsun.ttc'),
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
    plt.rcParams['axes.unicode_minus'] = False
except Exception:
    pass

# seaborn 是可选依赖 (本项目已装, 但下游环境未必)
try:
    import seaborn as sns
    _HAS_SEABORN = True
except ImportError:
    sns = None
    _HAS_SEABORN = False

from .evidence import CandidateEvidence


# ---- 候选评分热力图 ----

def plot_candidate_score_heatmap(
    evidences: List[CandidateEvidence],
    out_path: str,
    top_n: int = 10,
) -> Optional[str]:
    """绘制候选评分热力图: 行=候选, 列=维度 (score / 13C_cov / 1H_cov / DEPT_full / H_match).

    seaborn.heatmap 把分数矩阵可视化, 一眼看出哪条候选在哪个维度最强.
    返回 PNG 路径或 None (seaborn 未装 / 数据空).
    """
    if not _HAS_SEABORN or not evidences:
        return None
    rows: List[dict] = []
    for e in evidences[:top_n]:
        c_cov = sum(1 for p in e.carbon_evidence if p.matched) / max(1, len(e.carbon_evidence))
        h_cov = sum(1 for p in e.proton_evidence if p.matched) / max(1, len(e.proton_evidence))
        dept_full = sum(1 for d in e.dept_evidence if d.match_type == 'full') / max(1, len(e.dept_evidence))
        h_match = (e.h_count_evidence.match_score if e.h_count_evidence else 0.0)
        rows.append({
            'rank': e.rank,
            'score': e.score,
            '13C cov': c_cov,
            '1H cov': h_cov,
            'DEPT full': dept_full,
            'H 一致': h_match,
        })
    if not rows:
        return None
    try:
        import pandas as pd
        df = pd.DataFrame(rows).set_index('rank')
    except ImportError:
        return None
    # 各列归一化到 0-1 (避免 score 100+ 与 cov 0-1 量纲差)
    norm = df.copy()
    for col in norm.columns:
        vmax = float(norm[col].max())
        if vmax > 0:
            norm[col] = norm[col] / vmax
    fig, ax = plt.subplots(figsize=(9, 0.55 * len(norm) + 1.5))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sns.heatmap(
            norm, annot=df.round(2), fmt='.2f',
            cmap='RdYlGn', vmin=0, vmax=1, cbar_kws={'label': '归一化分数'},
            linewidths=0.5, ax=ax,
        )
    ax.set_title('候选评分热力图 (各列独立归一化到 [0, 1])', fontsize=12, pad=10)
    ax.set_xlabel('维度')
    ax.set_ylabel('Rank')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out_path


# ---- 实测 vs 候选 ppm 散点对比 ----

def plot_ppm_comparison_scatter(
    observed_c: List[float],
    observed_h: List[float],
    evidences: List[CandidateEvidence],
    out_path: str,
    top_n: int = 5,
) -> Optional[str]:
    """实测 vs Top-N 候选 ppm 对比散点图 (seaborn 风格).

    两个子图: 13C (上) 与 1H (下). 实测为垂直线, 候选为散点.
    颜色按 Rank 编码, 让差异可见.
    """
    if not _HAS_SEABORN or not evidences:
        return None
    items = list(evidences[:top_n])
    if not items:
        return None

    try:
        import pandas as pd
        c_rows = []
        h_rows = []
        for e in items:
            for s in e.lib_13c_shifts:
                c_rows.append({'rank': e.rank, 'shift': float(s)})
            for s in e.lib_1h_shifts:
                h_rows.append({'rank': e.rank, 'shift': float(s)})
        c_df = pd.DataFrame(c_rows) if c_rows else pd.DataFrame(columns=['rank', 'shift'])
        h_df = pd.DataFrame(h_rows) if h_rows else pd.DataFrame(columns=['rank', 'shift'])
    except ImportError:
        return None

    fig, (ax_c, ax_h) = plt.subplots(2, 1, figsize=(12, 7))
    # 13C
    if observed_c:
        for s in observed_c:
            ax_c.axvline(s, color='red', alpha=0.4, linewidth=0.7)
    if not c_df.empty:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sns.stripplot(data=c_df, x='shift', y='rank', orient='h',
                         palette='tab10', size=5, alpha=0.7, ax=ax_c, jitter=0.15)
    ax_c.set_xlabel('13C 化学位移 (ppm)')
    ax_c.set_ylabel('候选 Rank')
    ax_c.set_title('实测 (红) vs 候选 (散点) 13C 对比', fontsize=11)
    ax_c.invert_xaxis()
    ax_c.grid(True, alpha=0.3)

    # 1H
    if observed_h:
        for s in observed_h:
            ax_h.axvline(s, color='red', alpha=0.4, linewidth=0.7)
    if not h_df.empty:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sns.stripplot(data=h_df, x='shift', y='rank', orient='h',
                         palette='tab10', size=5, alpha=0.7, ax=ax_h, jitter=0.15)
    ax_h.set_xlabel('1H 化学位移 (ppm)')
    ax_h.set_ylabel('候选 Rank')
    ax_h.set_title('实测 (红) vs 候选 (散点) 1H 对比', fontsize=11)
    ax_h.invert_xaxis()
    ax_h.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out_path


# ---- 约束匹配矩阵 ----

def plot_constraint_match_matrix(
    evidences: List[CandidateEvidence],
    out_path: str,
    top_n: int = 10,
) -> Optional[str]:
    """约束匹配矩阵: 行=候选, 列=约束维度 (13C / 1H / DEPT / H_count / 分子式 / 溶剂).

    每格: 1=满足, 0.5=部分满足, 0=不满足. seaborn.heatmap 二值化矩阵.
    """
    if not _HAS_SEABORN or not evidences:
        return None
    items = list(evidences[:top_n])
    if not items:
        return None
    rows: List[dict] = []
    for e in items:
        c_cov = sum(1 for p in e.carbon_evidence if p.matched) / max(1, len(e.carbon_evidence))
        h_cov = sum(1 for p in e.proton_evidence if p.matched) / max(1, len(e.proton_evidence))
        dept_full = sum(1 for d in e.dept_evidence if d.match_type == 'full') / max(1, len(e.dept_evidence))
        h_match = (e.h_count_evidence.match_score if e.h_count_evidence else 0.5)
        formula_ok = 0.0 if any('分子式' in c for c in e.conflicts) else 1.0
        solvent_ok = 0.0 if any('溶剂' in c for c in e.conflicts) else 1.0
        rows.append({
            'rank': e.rank,
            '13C': c_cov,
            '1H': h_cov,
            'DEPT': dept_full,
            'H_count': h_match,
            '分子式': formula_ok,
            '溶剂': solvent_ok,
        })
    try:
        import pandas as pd
        df = pd.DataFrame(rows).set_index('rank')
    except ImportError:
        return None
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(df) + 1.5))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sns.heatmap(
            df, annot=True, fmt='.2f',
            cmap='RdYlGn', vmin=0, vmax=1,
            cbar_kws={'label': '满足度'},
            linewidths=0.5, ax=ax,
        )
    ax.set_title('约束匹配矩阵 (1=满足, 0=不满足, 0.5=部分/未知)', fontsize=12, pad=10)
    ax.set_xlabel('约束维度')
    ax.set_ylabel('Rank')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out_path


__all__ = [
    'plot_candidate_score_heatmap',
    'plot_ppm_comparison_scatter',
    'plot_constraint_match_matrix',
]