"""单元测试: score_candidates 打分逻辑."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from nmr_elucidator.evidence import score_candidates, ScoredCandidate
from nmr_elucidator.retrieve import RetrievalHit


def _hit(c_match, c_obs, c_lib, c_rms, h_match=0, h_obs=0, h_lib=0, h_rms=0):
    """构造 mock RetrievalHit."""
    from nmr_elucidator.io_sdf import NmrShiftDb2Entry
    e = NmrShiftDb2Entry(rec_id=1, smiles='C', canonical_smiles='C', mol_block='')
    return RetrievalHit(
        entry=e,
        c_match_count=c_match, c_total_obs=c_obs, c_total_lib=c_lib,
        c_coverage=c_match / max(c_obs, c_lib) if max(c_obs, c_lib) else 0,
        c_rms=c_rms,
        h_match_count=h_match, h_total_obs=h_obs, h_total_lib=h_lib,
        h_coverage=h_match / max(h_obs, h_lib) if max(h_obs, h_lib) else 0,
        h_rms=h_rms,
    )


def test_score_orders_by_coverage():
    """高覆盖率应排前面."""
    hits = [
        _hit(2, 5, 2, 0.5),   # cov=0.4
        _hit(4, 5, 5, 0.5),   # cov=0.8  ← 应该排第一
        _hit(1, 5, 5, 1.0),   # cov=0.2
    ]
    scored = score_candidates(hits)
    assert scored[0].rank == 1
    assert scored[0].score > scored[1].score > scored[2].score


def test_score_penalizes_rms():
    """RMS 越大分越低."""
    good = _hit(3, 5, 3, 0.0)
    bad = _hit(3, 5, 3, 5.0)
    s = score_candidates([bad, good])
    assert s[0].hit.entry.smiles == 'C'  # 顺序不影响
    # good 的分应高于 bad
    s2 = score_candidates([good, bad])
    assert s2[0].score > s2[1].score


def test_score_penalizes_count_diff():
    """峰数差异越大分越低."""
    matched = _hit(3, 3, 3, 0.5)        # lib = obs
    big_diff = _hit(3, 3, 10, 0.5)      # lib = 10, obs = 3, diff = 7
    s = score_candidates([big_diff, matched])
    # matched 的分应更高
    s2 = score_candidates([matched, big_diff])
    assert s2[0].score > s2[1].score


def test_score_ranks_assigned():
    """rank 字段自动赋值."""
    hits = [_hit(1, 5, 5, 0.5), _hit(2, 5, 5, 0.5)]
    scored = score_candidates(hits)
    assert scored[0].rank == 1
    assert scored[1].rank == 2


def test_score_empty():
    """空输入."""
    assert score_candidates([]) == []