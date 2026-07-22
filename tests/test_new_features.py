"""测试新增的 DEPT-135 sign 读取 + evidence dept 字段 + HTML 报告."""
import os
import sys
import tempfile

sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')
from nmr_elucidator.io_bruker import _parse_peaklist_xml, Peak, _infer_peak_sign
from nmr_elucidator.bridge_nmr2struct import assign_carbon_multiplicity, DEPT_CH, DEPT_CH2, DEPT_CH3, DEPT_CQ
from nmr_elucidator.evidence import (
    score_candidates, build_evidence, build_dept_evidence, compute_h_count_evidence,
    DeptEvidence, HCountEvidence,
)
from nmr_elucidator.retrieve import RetrievalHit
from nmr_elucidator.io_sdf import NmrShiftDb2Entry
from nmr_elucidator.report import write_report_html


# ---- _infer_peak_sign ----

def test_infer_peak_sign_positive():
    assert _infer_peak_sign(2.5) == 'positive'
    assert _infer_peak_sign(0.01) == 'positive'
    print('PASS test_infer_peak_sign_positive')


def test_infer_peak_sign_negative():
    assert _infer_peak_sign(-2.5) == 'negative'
    assert _infer_peak_sign(-0.01) == 'negative'
    print('PASS test_infer_peak_sign_negative')


def test_infer_peak_sign_zero():
    assert _infer_peak_sign(0.0) == 'unknown'
    print('PASS test_infer_peak_sign_zero')


# ---- assign_carbon_multiplicity 用 sign ----

def test_assign_with_signs_split_ch2_ch3():
    """DEPT-135 sign=positive + DEPT-90 无 → CH3; sign=negative → CH2."""
    main = [20.0, 50.0]
    d90 = []  # 无 DEPT-90
    d135 = [20.0, 50.0]
    d135_signs = ['positive', 'negative']
    out = assign_carbon_multiplicity(main, d90, d135, d135_signs)
    # 降序: [50, 20]
    by_shift = {round(o.shift, 3): o for o in out}
    assert by_shift[20.0].multiplicity == DEPT_CH3
    assert by_shift[20.0].sign == 'positive'
    assert by_shift[50.0].multiplicity == DEPT_CH2
    assert by_shift[50.0].sign == 'negative'
    print('PASS test_assign_with_signs_split_ch2_ch3')


def test_assign_dept90_overrides_positive_in_alkene_region():
    """DEPT-90 命中 → 总是 CH (芳/烯/三级碳), 即使 DEPT-135 sign=positive."""
    main = [130.0]  # 芳/烯区
    d90 = [130.0]
    d135 = [130.0]
    d135_signs = ['positive']
    out = assign_carbon_multiplicity(main, d90, d135, d135_signs)
    assert out[0].multiplicity == DEPT_CH
    print('PASS test_assign_dept90_overrides_positive_in_alkene_region')


# ---- score_candidates 加分 ----

def _hit(c_match=5, c_obs=10, c_lib=10, c_rms=0.5, h_match=5, h_obs=10, h_lib=10, h_rms=0.3,
         dept_total=10, dept_coverage=0.7, h_count_expected=20, h_count_lib=20):
    e = NmrShiftDb2Entry(rec_id=1, smiles='C', canonical_smiles='C', mol_block='')
    return RetrievalHit(
        entry=e,
        c_match_count=c_match, c_total_obs=c_obs, c_total_lib=c_lib,
        c_coverage=c_match / max(c_obs, c_lib), c_rms=c_rms,
        h_match_count=h_match, h_total_obs=h_obs, h_total_lib=h_lib,
        h_coverage=h_match / max(h_obs, h_lib), h_rms=h_rms,
        heavy_atoms=10,
        dept_total=dept_total, dept_coverage=dept_coverage,
        h_count_expected=h_count_expected, h_count_lib=h_count_lib,
    )


def test_score_dept_bonus_changes_ranking():
    """高 DEPT 覆盖率应排前面."""
    high_dept = _hit(dept_total=10, dept_coverage=0.9)
    low_dept = _hit(dept_total=10, dept_coverage=0.1)
    s = score_candidates([low_dept, high_dept])
    assert s[0].hit.dept_coverage == 0.9
    print('PASS test_score_dept_bonus_changes_ranking')


def test_score_hcount_match_bonus():
    """H 总数一致应加分."""
    perfect = _hit(h_count_expected=20, h_count_lib=20)
    mismatch = _hit(h_count_expected=20, h_count_lib=5)
    s = score_candidates([mismatch, perfect])
    assert s[0].hit.h_count_lib == 20  # perfect 排前
    assert s[0].score > s[1].score
    print('PASS test_score_hcount_match_bonus')


# ---- HTML 报告 ----

def test_write_report_html_basic():
    """HTML 报告能写出且含关键内容."""
    e = NmrShiftDb2Entry(rec_id=1, smiles='CCO', canonical_smiles='CCO', mol_block='')
    from nmr_elucidator.evidence import PeakAssignment
    hit = RetrievalHit(
        entry=e,
        c_match_count=1, c_total_obs=2, c_total_lib=2,
        c_coverage=0.5, c_rms=0.3,
        h_match_count=1, h_total_obs=2, h_total_lib=2,
        h_coverage=0.5, h_rms=0.1,
        heavy_atoms=2,
        dept_total=2, dept_coverage=1.0,
        h_count_expected=3, h_count_lib=3,
    )
    from nmr_elucidator.evidence import ScoredCandidate
    scored = ScoredCandidate(hit=hit, score=70.0, rank=1)
    obs_dept = [
        {'ppm': 50.0, 'mult': 'CH2', 'sign': 'negative'},
        {'ppm': 20.0, 'mult': 'CH3', 'sign': 'positive'},
    ]
    evid = build_evidence(
        scored, observed_c=[50.0, 20.0], observed_h=[1.2, 3.5],
        c_tol=2.0, h_tol=0.3, formula_constraint=None,
        obs_dept_summary=obs_dept,
    )
    meta = {
        'folder': 'test_folder', 'sdf_path': '/tmp/sdf.sd',
        'formula_constraint': None, 'solvent_observed': 'CDCl3',
        'field_mhz': 500.0, 'c_tol': 2.0, 'h_tol': 0.3,
        'top_n': 1, 'used_nmrsim': False,
        'dept_counts': {'CH': 0, 'CH2': 1, 'CH3': 1, 'CH2/CH3': 0, 'Cq': 0, '?': 0},
        'observed_carbon_shifts': [50.0, 20.0],
        'observed_proton_shifts': [1.2, 3.5],
        'obs_dept_summary': obs_dept,
    }
    out_path = r'E:\NMR2Structure\data\test_html.html'
    write_report_html(out_path, [evid], meta, top_n=1)
    assert os.path.exists(out_path)
    with open(out_path, 'r', encoding='utf-8') as f:
        html = f.read()
    # 关键内容
    assert 'DEPT 推断' in html
    assert 'CCO' in html
    assert 'CH2' in html
    assert 'CH3' in html
    assert 'H 总数' in html
    assert 'Top-1 候选' in html or 'Top-1' in html
    assert '解谱过程' in html
    assert '未匹配峰' in html
    print('PASS test_write_report_html_basic (size=%d bytes)' % os.path.getsize(out_path))
    # 清理
    os.remove(out_path)


if __name__ == '__main__':
    test_infer_peak_sign_positive()
    test_infer_peak_sign_negative()
    test_infer_peak_sign_zero()
    test_assign_with_signs_split_ch2_ch3()
    test_assign_dept90_overrides_positive_in_alkene_region()
    test_score_dept_bonus_changes_ranking()
    test_score_hcount_match_bonus()
    test_write_report_html_basic()
    print('\nALL PASSED')
