"""候选结构表 (pandas DataFrame 表示).

把 nmr_elucidator.evidence.CandidateEvidence 列表转 pandas DataFrame,
便于:
1. 按 InChIKey 去重相似候选 (同结构不同盐/同位素)
2. 按多维度排序 (score / coverage / dept_full)
3. 直接输出 Markdown 表格嵌入 report.py
4. CSV 导出供化学家在 Excel 中排序/筛选

依赖: pandas>=2.0 (本项目环境已装).
"""
from __future__ import annotations

from typing import Iterable, List, Optional

import pandas as pd

from .evidence import CandidateEvidence, HCountEvidence


def _dept_full_n(e: CandidateEvidence) -> int:
    return sum(1 for d in e.dept_evidence if d.match_type == 'full')


def _dept_total_n(e: CandidateEvidence) -> int:
    return len(e.dept_evidence)


def _c_match_n(e: CandidateEvidence) -> int:
    return sum(1 for p in e.carbon_evidence if p.matched)


def _h_match_n(e: CandidateEvidence) -> int:
    return sum(1 for p in e.proton_evidence if p.matched)


def _h_count_str(e: CandidateEvidence) -> str:
    hce = e.h_count_evidence
    if hce is None:
        return '-'
    return f'{hce.lib_heavy_h}/{hce.expected}'


def _conflict_str(e: CandidateEvidence) -> str:
    return '; '.join(e.conflicts) if e.conflicts else '-'


def evidences_to_dataframe(evidences: Iterable[CandidateEvidence]) -> pd.DataFrame:
    """把 CandidateEvidence 列表转为 DataFrame.

    列: rank, score, smiles, nmrshiftdb2_id, formula, inchikey, solvent_lib,
        heavy_atoms, c_match, c_total, h_match, h_total, dept_full, dept_total,
        h_count, conflicts.

    不匹配的峰信息保留在 CandidateEvidence 原始对象 (DataFrame 适合排序/筛选, 不适合存数组).
    """
    rows: List[dict] = []
    for e in evidences:
        rows.append({
            'rank': e.rank,
            'score': round(float(e.score), 2),
            'smiles': e.smiles,
            'nmrshiftdb2_id': e.nmrshiftdb2_id,
            'formula': e.formula,
            'inchikey': e.inchikey,
            'solvent_lib': e.solvent_lib,
            'heavy_atoms': e.heavy_atoms,
            'c_match': _c_match_n(e),
            'c_total': len(e.carbon_evidence),
            'h_match': _h_match_n(e),
            'h_total': len(e.proton_evidence),
            'dept_full': _dept_full_n(e),
            'dept_total': _dept_total_n(e),
            'h_count': _h_count_str(e),
            'conflicts': _conflict_str(e),
        })
    if not rows:
        return pd.DataFrame(columns=[
            'rank', 'score', 'smiles', 'nmrshiftdb2_id', 'formula',
            'inchikey', 'solvent_lib', 'heavy_atoms',
            'c_match', 'c_total', 'h_match', 'h_total',
            'dept_full', 'dept_total', 'h_count', 'conflicts'])
    df = pd.DataFrame(rows)
    # 类型规整 (避免 object 列影响排序)
    df['score'] = df['score'].astype('float64')
    df['heavy_atoms'] = df['heavy_atoms'].astype('Int64')
    df['c_match'] = df['c_match'].astype('Int64')
    df['c_total'] = df['c_total'].astype('Int64')
    df['h_match'] = df['h_match'].astype('Int64')
    df['h_total'] = df['h_total'].astype('Int64')
    df['dept_full'] = df['dept_full'].astype('Int64')
    df['dept_total'] = df['dept_total'].astype('Int64')
    return df


def dedupe_by_inchikey(df: pd.DataFrame, keep: str = 'first') -> pd.DataFrame:
    """按 InChIKey 去重相似候选 (同结构不同盐形式或同位素标记).

    keep: 'first' (保留首次出现) / 'last' (保留末次) / 'max_score' (保留分数最高).
          'first' 与 pandas.drop_duplicates 默认行为一致.
    去重后 rank 重新从 1 开始.

    注意: InChIKey 缺失 (None) 的候选不被合并 (视为独立候选).
    """
    if df.empty or 'inchikey' not in df.columns:
        return df.reset_index(drop=True)
    if keep == 'max_score':
        # 按 inchikey 分组, 保留 score 最大行
        df_sorted = df.sort_values('score', ascending=False, na_position='last')
        df = df_sorted.drop_duplicates(subset=['inchikey'], keep='first')
    else:
        df = df.drop_duplicates(subset=['inchikey'], keep=keep)
    df = df.reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)
    return df


def to_markdown(df: pd.DataFrame, columns: Optional[List[str]] = None) -> str:
    """把 DataFrame 渲染为 Markdown 表格.

    columns: 指定要展示的列. None 用全部列.
    返回: Markdown 字符串 (含表头与对齐说明行).
    """
    if df.empty:
        return '_(空候选表)_'
    cols = columns if columns else list(df.columns)
    # to_markdown 在空 DataFrame 上会抛, 已在外层处理
    try:
        return df[cols].to_markdown(index=False)
    except ImportError:
        # tabulate 未装时退化为手写 Markdown
        return _manual_markdown(df[cols])


def _manual_markdown(df: pd.DataFrame) -> str:
    """不依赖 tabulate 的 Markdown 表格渲染 (兜底)."""
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(str(c) for c in cols) + ' |']
    lines.append('|' + '|'.join('---' for _ in cols) + '|')
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if v is None or (isinstance(v, float) and pd.isna(v)):
                cells.append('-')
            else:
                cells.append(str(v))
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def to_csv(df: pd.DataFrame, path: str) -> None:
    """写 CSV (UTF-8 with BOM, Excel 友好)."""
    import os
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    df.to_csv(path, index=False, encoding='utf-8-sig')