"""结构推断置信度评估.

融合多源证据, 输出置信度分级结论:
  - determined (确定): Top-1 综合分 > 0.85 且与 Top-2 差距 > 0.15
  - probable (很可能): Top-1 综合分 > 0.70 且与 Top-2 差距 > 0.08
  - ambiguous (不确定): 多个候选分数接近, 需 2D NMR 进一步确认
  - no_match (无匹配): 所有候选验证分数 < 0.50
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from .evidence import CandidateEvidence
from .nmr_predict import VerificationScore


# ---------- 置信度等级 ----------

CONFIDENCE_DETERMINED = 'determined'
CONFIDENCE_PROBABLE = 'probable'
CONFIDENCE_AMBIGUOUS = 'ambiguous'
CONFIDENCE_NO_MATCH = 'no_match'

_CONFIDENCE_LABELS = {
    CONFIDENCE_DETERMINED: '确定',
    CONFIDENCE_PROBABLE: '很可能',
    CONFIDENCE_AMBIGUOUS: '不确定 (需 2D NMR 确认)',
    CONFIDENCE_NO_MATCH: '无匹配',
}


# ---------- 数据结构 ----------

@dataclass
class DeterminationResult:
    """结构推断结论."""

    conclusion: str  # determined / probable / ambiguous / no_match
    conclusion_label: str  # 中文标签
    confidence: float  # 0-1 综合置信度
    top_candidate_smiles: Optional[str] = None
    top_candidate_score: float = 0.0
    gap_to_second: float = 0.0
    # 多源分数明细
    score_breakdown: dict = field(default_factory=dict)
    # 证据链摘要
    evidence_summary: List[str] = field(default_factory=list)
    # 建议
    recommendations: List[str] = field(default_factory=list)


# ---------- 融合打分 ----------

def compute_unified_score(
    evidence: CandidateEvidence,
    hose_verification: Optional[VerificationScore] = None,
    weights: Optional[dict] = None,
) -> float:
    """计算单个候选的统一验证分 (0-1).

    融合维度:
      - 13C 覆盖率 (来自 evidence)
      - 1H 覆盖率 (来自 evidence)
      - HOSE 13C 预测 RMS (来自 hose_verification)
      - DEPT 碳级数一致性 (来自 evidence)
      - 分子式/DoU 一致性 (来自 evidence)

    Args:
        evidence: 候选证据链
        hose_verification: HOSE 预测验证结果 (可选)
        weights: 各维度权重 (可选, 有默认值)

    Returns:
        统一分数 (0-1)
    """
    if weights is None:
        weights = {
            'c_coverage': 0.20,
            'h_coverage': 0.15,
            'hose_rms': 0.30,
            'dept_match': 0.15,
            'formula_match': 0.10,
            'peak_count': 0.10,
        }

    scores = {}

    # 13C 覆盖率
    c_matched = sum(1 for e in evidence.carbon_evidence if e.matched)
    c_total = len(evidence.carbon_evidence) if evidence.carbon_evidence else 1
    scores['c_coverage'] = c_matched / c_total if c_total > 0 else 0.0

    # 1H 覆盖率
    h_matched = sum(1 for e in evidence.proton_evidence if e.matched)
    h_total = len(evidence.proton_evidence) if evidence.proton_evidence else 1
    scores['h_coverage'] = h_matched / h_total if h_total > 0 else 0.0

    # HOSE 预测验证
    if hose_verification is not None:
        # RMS 分: 1 - rms/10 (rms=0 → 1.0, rms>=10 → 0.0)
        rms_score = max(0.0, 1.0 - hose_verification.hose_rms / 10.0)
        # 匹配率
        match_score = hose_verification.hose_match_rate
        # 综合: 0.6 * match_rate + 0.4 * rms_score
        scores['hose_rms'] = 0.6 * match_score + 0.4 * rms_score
    else:
        scores['hose_rms'] = 0.0  # 无 HOSE 数据时不加分

    # DEPT 碳级数一致性
    if evidence.dept_evidence:
        dept_full = sum(1 for d in evidence.dept_evidence if d.match_type == 'full')
        dept_total = len(evidence.dept_evidence)
        scores['dept_match'] = dept_full / dept_total if dept_total > 0 else 0.0
    else:
        scores['dept_match'] = 0.5  # 无 DEPT 数据时给中性分

    # 分子式一致性
    if evidence.conflicts:
        # 有冲突 → 扣分
        formula_conflicts = [c for c in evidence.conflicts if '分子式' in c]
        scores['formula_match'] = 0.0 if formula_conflicts else 0.7
    else:
        scores['formula_match'] = 1.0

    # 峰数差异 (未匹配峰越少越好)
    n_unmatched = len(evidence.unmatched_carbon_shifts) + len(evidence.unmatched_proton_shifts)
    n_total_peaks = c_total + h_total
    if n_total_peaks > 0:
        scores['peak_count'] = max(0.0, 1.0 - n_unmatched / n_total_peaks)
    else:
        scores['peak_count'] = 0.5

    # 加权求和
    total = sum(weights.get(k, 0.0) * v for k, v in scores.items())
    return round(total, 4)


# ---------- 置信度评估 ----------

def assess_determination(
    evidences: List[CandidateEvidence],
    hose_verifications: Optional[List[VerificationScore]] = None,
    weights: Optional[dict] = None,
) -> DeterminationResult:
    """对所有候选进行置信度评估, 输出结构推断结论.

    Args:
        evidences: 候选证据链列表 (已按 rank 排序)
        hose_verifications: 对应的 HOSE 验证结果 (可选, 与 evidences 等长)
        weights: 融合权重

    Returns:
        DeterminationResult 结构推断结论
    """
    if not evidences:
        return DeterminationResult(
            conclusion=CONFIDENCE_NO_MATCH,
            conclusion_label=_CONFIDENCE_LABELS[CONFIDENCE_NO_MATCH],
            confidence=0.0,
            evidence_summary=['无候选结构'],
            recommendations=['检查输入数据质量', '放宽容差参数重试'],
        )

    # 计算每个候选的统一分数
    unified_scores: List[float] = []
    for i, ev in enumerate(evidences):
        hv = hose_verifications[i] if (hose_verifications and i < len(hose_verifications)) else None
        unified_scores.append(compute_unified_score(ev, hv, weights))

    # 排序 (降序)
    indexed = sorted(enumerate(unified_scores), key=lambda x: x[1], reverse=True)
    top_idx, top_score = indexed[0]
    second_score = indexed[1][1] if len(indexed) > 1 else 0.0
    gap = top_score - second_score

    top_ev = evidences[top_idx]

    # 构建分数明细
    breakdown = {
        'top_smiles': top_ev.smiles,
        'top_unified_score': top_score,
        'second_unified_score': second_score,
        'gap': gap,
        'n_candidates': len(evidences),
    }

    # 证据链摘要
    summary = [
        f'Top-1: {top_ev.smiles} (统一分={top_score:.3f})',
        f'Top-1 原始打分: {top_ev.score:.1f} (rank #{top_ev.rank})',
    ]
    if len(evidences) > 1:
        second_ev = evidences[indexed[1][0]]
        summary.append(f'Top-2: {second_ev.smiles} (统一分={second_score:.3f})')
    summary.append(f'Top-1 与 Top-2 差距: {gap:.3f}')

    # 置信度分级
    recommendations: List[str] = []

    if top_score >= 0.85 and gap >= 0.15:
        conclusion = CONFIDENCE_DETERMINED
        confidence = min(1.0, top_score)
        recommendations.append('结构已确定, 可直接用于后续研究')
    elif top_score >= 0.70 and gap >= 0.08:
        conclusion = CONFIDENCE_PROBABLE
        confidence = top_score
        recommendations.append('结构很可能正确, 建议用 2D NMR (HSQC/HMBC) 进一步确认')
    elif top_score >= 0.50:
        conclusion = CONFIDENCE_AMBIGUOUS
        confidence = top_score * 0.7  # 不确定时降低置信度
        recommendations.append('多个候选分数接近, 需要 2D NMR 数据区分')
        recommendations.append(f'建议优先获取 HSQC 和 HMBC 谱图')
        # 列出前 3 候选
        for rank_i, (idx, sc) in enumerate(indexed[:3], 1):
            recommendations.append(f'  候选 {rank_i}: {evidences[idx].smiles} (分={sc:.3f})')
    else:
        conclusion = CONFIDENCE_NO_MATCH
        confidence = top_score
        recommendations.append('所有候选验证分数过低, 可能是新骨架化合物')
        recommendations.append('建议: 1) 检查分子式是否正确; 2) 获取 2D NMR; 3) 考虑片段组装')

    return DeterminationResult(
        conclusion=conclusion,
        conclusion_label=_CONFIDENCE_LABELS[conclusion],
        confidence=round(confidence, 4),
        top_candidate_smiles=top_ev.smiles,
        top_candidate_score=top_score,
        gap_to_second=round(gap, 4),
        score_breakdown=breakdown,
        evidence_summary=summary,
        recommendations=recommendations,
    )
