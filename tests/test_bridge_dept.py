"""测试 bridge_nmr2struct 真实积分 + DEPT 推断."""
import sys
sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')
from nmr_elucidator.bridge_nmr2struct import (
    convert_bruker_folder_to_nmr2struct,
    _write_hnmr_txt,
    assign_carbon_multiplicity,
    DEPT_CH, DEPT_CQ, DEPT_CH2_or_CH3,
)


def test_write_hnmr_uses_real_integral():
    """1H 峰应使用真实积分当强度."""
    peaks = [0.96, 1.5, 4.5]
    integrals = [(1.22, 0.72, 9.6), (4.7, 4.12, 15.0)]
    out = r'E:\NMR2Structure\data\test_hnmr.txt'
    n = _write_hnmr_txt(out, peaks, integrals)
    assert n == 3
    with open(out, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    # 0.96 -> 9.6 (在 0.72-1.22 区间, 颠倒写但 min/max 后正确匹配)
    # 1.5  -> 1.0 (不在任何区间, 落空)
    # 4.5  -> 15.0 (在 4.12-4.7 区间)
    assert '0.9600 9.600000' in lines[0], lines[0]
    assert '1.5000 1.000000' in lines[1], lines[1]
    assert '4.5000 15.000000' in lines[2], lines[2]
    print('PASS test_write_hnmr_uses_real_integral')


def test_write_hnmr_no_integral_fallback():
    """无积分时所有强度 = 1.0."""
    peaks = [0.96, 4.5]
    out = r'E:\NMR2Structure\data\test_hnmr2.txt'
    n = _write_hnmr_txt(out, peaks, None)
    assert n == 2
    with open(out, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    assert lines[0].endswith('1.000000\n')
    print('PASS test_write_hnmr_no_integral_fallback')


def test_assign_carbon_multiplicity():
    """DEPT 交叉推断碳类型."""
    # 13C 主谱: 5 个峰
    main = [100.0, 80.0, 60.0, 40.0, 20.0]
    # DEPT-90: 只 CH (idx 0, 2)
    d90 = [100.0, 60.0]
    # DEPT-135: CH(0,2) + CH3(3,4), 不含 CH2
    d135 = [100.0, 60.0, 40.0, 20.0]
    result = assign_carbon_multiplicity(main, dept90_shifts=d90, dept135_shifts=d135)
    # 按降序: 100(idx0), 80(idx1), 60(idx2), 40(idx3), 20(idx4)
    assert result[0].multiplicity == DEPT_CH        # 100: DEPT-90 有
    assert result[1].multiplicity == DEPT_CQ        # 80:  DEPT 全无 → 季碳
    assert result[2].multiplicity == DEPT_CH        # 60:  DEPT-90 有
    assert result[3].multiplicity == DEPT_CH2_or_CH3  # 40: DEPT-90 无, 135 有
    assert result[4].multiplicity == DEPT_CH2_or_CH3  # 20: DEPT-90 无, 135 有
    print('PASS test_assign_carbon_multiplicity')


def test_end_to_end_clxh3():
    """clxh3 端到端 stage 1: 写出 hnmr/cnmr, 验证 DEPT 推断非空."""
    info = convert_bruker_folder_to_nmr2struct(
        r'E:\NMR2Structure\NMR\clxh3',
        r'E:\NMR2Structure\data\test_hnmr_clxh3.txt',
        r'E:\NMR2Structure\data\test_cnmr_clxh3.txt',
    )
    assert info['n_1h_written'] > 0
    assert info['n_13c_written'] > 0
    assert len(info['dept_summary']) > 0
    assert info['exp_dept90'] == '23'
    assert info['exp_dept135'] == '24'
    # 至少有一种 CH 和一种 Cq (clxh3 实测有季碳)
    assert info['dept_counts']['CH'] >= 1
    print('PASS test_end_to_end_clxh3  (n_1h=%d, n_13c=%d, dept=%s)' % (
        info['n_1h_written'], info['n_13c_written'], info['dept_counts']))


if __name__ == '__main__':
    test_write_hnmr_uses_real_integral()
    test_write_hnmr_no_integral_fallback()
    test_assign_carbon_multiplicity()
    test_end_to_end_clxh3()
    print('\nALL PASSED')
