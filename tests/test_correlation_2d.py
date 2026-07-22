"""2D 关联验证测试.

验证:
  - HSQC/HMBC/COSY 关联提取
  - 候选分子 2D 打分
  - 无 2D 数据时的回退行为
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from nmr_elucidator.io_bruker_2d import read_bruker_2d_folder
from nmr_elucidator.correlation_2d import (
    Correlation2DScore,
    extract_cosy_pairs,
    extract_hmbc_pairs,
    extract_hsqc_pairs,
    score_candidate_2d,
)

CLXH100_FOLDER = r'E:\NMR2Structure\NMR\clxh100'


def test_extract_hsqc_pairs():
    """测试 HSQC 关联提取."""
    if not os.path.isdir(CLXH100_FOLDER):
        print(f'  SKIP: {CLXH100_FOLDER} 不存在')
        return
    spectra = read_bruker_2d_folder(CLXH100_FOLDER, quiet=True)
    if 'HSQC' not in spectra:
        print('  SKIP: 无 HSQC 数据')
        return
    pairs = extract_hsqc_pairs(spectra['HSQC'], min_snr=20.0)
    assert len(pairs) > 0, 'HSQC 应提取到关联对'
    # 检查格式: (1H_ppm, 13C_ppm)
    for h, c in pairs[:5]:
        assert -1 < h < 14, f'1H 位移超范围: {h}'
        assert 0 < c < 220, f'13C 位移超范围: {c}'
    print(f'  test_extract_hsqc_pairs: OK ({len(pairs)} pairs)')


def test_extract_hmbc_pairs():
    """测试 HMBC 关联提取."""
    if not os.path.isdir(CLXH100_FOLDER):
        print(f'  SKIP: {CLXH100_FOLDER} 不存在')
        return
    spectra = read_bruker_2d_folder(CLXH100_FOLDER, quiet=True)
    if 'HMBC' not in spectra:
        print('  SKIP: 无 HMBC 数据')
        return
    pairs = extract_hmbc_pairs(spectra['HMBC'], min_snr=10.0)
    assert len(pairs) > 0, 'HMBC 应提取到关联对'
    print(f'  test_extract_hmbc_pairs: OK ({len(pairs)} pairs)')


def test_extract_cosy_pairs():
    """测试 COSY 关联提取 (去对角线 + 去重复)."""
    if not os.path.isdir(CLXH100_FOLDER):
        print(f'  SKIP: {CLXH100_FOLDER} 不存在')
        return
    spectra = read_bruker_2d_folder(CLXH100_FOLDER, quiet=True)
    if 'COSY' not in spectra:
        print('  SKIP: 无 COSY 数据')
        return
    pairs = extract_cosy_pairs(spectra['COSY'], min_snr=15.0)
    assert len(pairs) > 0, 'COSY 应提取到关联对'
    # 检查: 所有对都是 (a, b) 且 a < b (去重复)
    for a, b in pairs:
        assert a < b, f'COSY 对应 a < b, 实得 ({a}, {b})'
        assert abs(a - b) >= 0.1, f'对角线峰未去除: ({a}, {b})'
    print(f'  test_extract_cosy_pairs: OK ({len(pairs)} pairs)')


def test_score_candidate_2d_no_data():
    """测试无 2D 数据时的打分 (应返回零分)."""
    from nmr_elucidator.io_sdf import NmrShiftDb2Entry
    # 构造一个空 entry
    entry = NmrShiftDb2Entry(
        rec_id=1,
        smiles='C',
        canonical_smiles='C',
        mol_block='',
    )
    score = score_candidate_2d(entry, {})
    assert score.combined_2d_score == 0.0
    assert score.n_hsqc_obs == 0
    print('  test_score_candidate_2d_no_data: OK')


def test_score_candidate_2d_with_data():
    """测试有 2D 数据时的打分 (需要 nmrshiftdb2 entry 有 per-atom 归属)."""
    if not os.path.isdir(CLXH100_FOLDER):
        print(f'  SKIP: {CLXH100_FOLDER} 不存在')
        return
    spectra = read_bruker_2d_folder(CLXH100_FOLDER, quiet=True)
    if not spectra:
        print('  SKIP: 无 2D 数据')
        return

    # 用一个简单的 entry 测试 (无 per-atom 归属时应返回零分)
    from nmr_elucidator.io_sdf import NmrShiftDb2Entry, SpectrumSignal
    entry = NmrShiftDb2Entry(
        rec_id=1,
        smiles='CCO',
        canonical_smiles='CCO',
        mol_block='',  # 无 mol_block 时无法构建预期关联
        spectrum_13c=[SpectrumSignal(shift=18.0, multiplicity='Q', atom_index=0)],
        spectrum_1h=[SpectrumSignal(shift=1.2, multiplicity='T', atom_index=1)],
    )
    score = score_candidate_2d(entry, spectra)
    # 无 mol_block 时无法构建预期关联, 应返回零分
    assert score.combined_2d_score == 0.0
    print('  test_score_candidate_2d_with_data: OK (无 mol_block 时零分)')


def test_correlation_2d_score_dataclass():
    """测试 Correlation2DScore 数据结构."""
    score = Correlation2DScore(
        hsqc_match_rate=0.8,
        hmbc_match_rate=0.6,
        cosy_match_rate=0.7,
        combined_2d_score=0.72,
        n_hsqc_obs=50,
        n_hmbc_obs=100,
        n_cosy_obs=80,
    )
    assert score.hsqc_match_rate == 0.8
    assert score.combined_2d_score == 0.72
    print('  test_correlation_2d_score_dataclass: OK')


if __name__ == '__main__':
    print('=== test_correlation_2d ===')
    test_extract_hsqc_pairs()
    test_extract_hmbc_pairs()
    test_extract_cosy_pairs()
    test_score_candidate_2d_no_data()
    test_score_candidate_2d_with_data()
    test_correlation_2d_score_dataclass()
    print('=== ALL PASSED ===')
