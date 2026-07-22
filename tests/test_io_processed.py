"""单元测试: io_processed (1r 读取 + resample_to_grid 重采样).

覆盖重点 (AGENTS.md 规则 13): 边界数值, 非法输入, 错误阶段, 真实数据保真.
resample_to_grid 是纯逻辑, 用构造数据严格验证插值/裁剪/归一化.
真实 1r 读取依赖 clxh3 原始数据, 缺失时跳过 (不制造脆弱失败).
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')

from nmr_elucidator.io_processed import (
    ProcessedReadError,
    read_processed_folder,
    read_processed_spectrum,
    resample_to_grid,
)

_CLXH3 = r'E:\NMR2Structure\NMR\clxh3'


# ---- resample_to_grid: 插值正确性 ----

def test_resample_exact_grid_points():
    """网格点与采样点重合时, 应精确还原强度 (normalize=False)."""
    ppm = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    inten = np.array([0.0, 10.0, 20.0, 10.0, 0.0])
    out = resample_to_grid(ppm, inten, 0.0, 4.0, 5, normalize=False)
    assert np.allclose(out, [0, 10, 20, 10, 0]), out
    print('PASS test_resample_exact_grid_points')


def test_resample_handles_descending_ppm():
    """NMR ppm 轴通常降序 (高场在右); argsort 后仍应正确插值."""
    ppm = np.array([4.0, 3.0, 2.0, 1.0, 0.0])  # 降序
    inten = np.array([0.0, 10.0, 20.0, 10.0, 0.0])  # 对应 4->0, 2->20
    out = resample_to_grid(ppm, inten, 0.0, 4.0, 5, normalize=False)
    # grid[2]=2.0 对应强度 20
    assert abs(out[2] - 20.0) < 1e-6, out
    assert abs(out[0] - 0.0) < 1e-6
    print('PASS test_resample_handles_descending_ppm')


def test_resample_out_of_range_filled_zero():
    """网格落在采样范围之外时填 0, 不外插造成假信号."""
    ppm = np.array([10.0, 11.0, 12.0])
    inten = np.array([5.0, 5.0, 5.0])
    out = resample_to_grid(ppm, inten, 0.0, 4.0, 5, normalize=False)
    assert np.all(out == 0.0), out
    print('PASS test_resample_out_of_range_filled_zero')


def test_resample_negative_clipped_to_zero():
    """负强度 (基线/相位残留) 应裁剪为 0."""
    ppm = np.array([0.0, 1.0, 2.0])
    inten = np.array([-5.0, 10.0, -3.0])
    out = resample_to_grid(ppm, inten, 0.0, 2.0, 3, normalize=False)
    assert out[0] == 0.0 and out[2] == 0.0, out
    assert abs(out[1] - 10.0) < 1e-6
    print('PASS test_resample_negative_clipped_to_zero')


def test_resample_normalize_to_max_one():
    """normalize=True 时峰值归一到 1.0."""
    ppm = np.array([0.0, 1.0, 2.0])
    inten = np.array([0.0, 20.0, 5.0])
    out = resample_to_grid(ppm, inten, 0.0, 2.0, 3, normalize=True)
    assert abs(out.max() - 1.0) < 1e-6, out
    assert abs(out[1] - 1.0) < 1e-6
    print('PASS test_resample_normalize_to_max_one')


def test_resample_length_and_dtype():
    """输出长度 == n_points, dtype == float32."""
    ppm = np.linspace(0, 10, 100)
    inten = np.abs(np.sin(ppm))
    out = resample_to_grid(ppm, inten, 0.0, 10.0, 80)
    assert out.shape == (80,), out.shape
    assert out.dtype == np.float32, out.dtype
    print('PASS test_resample_length_and_dtype')


def test_resample_all_zero_intensity_no_div0():
    """全零强度归一化不应除零 (返回全 0, 不是 nan)."""
    ppm = np.array([0.0, 1.0, 2.0])
    inten = np.array([0.0, 0.0, 0.0])
    out = resample_to_grid(ppm, inten, 0.0, 2.0, 3, normalize=True)
    assert np.all(out == 0.0), out
    assert not np.any(np.isnan(out))
    print('PASS test_resample_all_zero_intensity_no_div0')


# ---- resample_to_grid: 错误路径 ----

def test_resample_length_mismatch_raises():
    ppm = np.array([0.0, 1.0, 2.0])
    inten = np.array([1.0, 2.0])
    try:
        resample_to_grid(ppm, inten, 0.0, 2.0, 3)
        assert False, '应抛 ValueError'
    except ValueError as e:
        assert '长度不一致' in str(e)
    print('PASS test_resample_length_mismatch_raises')


def test_resample_bad_npoints_raises():
    ppm = np.array([0.0, 1.0])
    inten = np.array([1.0, 2.0])
    try:
        resample_to_grid(ppm, inten, 0.0, 1.0, 0)
        assert False, '应抛 ValueError'
    except ValueError as e:
        assert 'n_points' in str(e)
    print('PASS test_resample_bad_npoints_raises')


# ---- read 错误路径 (无需真实数据) ----

def test_read_folder_nonexistent_raises():
    try:
        read_processed_folder(r'E:\不存在的目录_xyz123')
        assert False, '应抛 ProcessedReadError'
    except ProcessedReadError as e:
        assert e.stage == 'check-folder', e.stage
    print('PASS test_read_folder_nonexistent_raises')


def test_read_spectrum_missing_1r_raises():
    """pdata 目录存在但无 1r 文件 -> stage=check-1r."""
    with tempfile.TemporaryDirectory() as d:
        try:
            read_processed_spectrum(d, exp_no=99, nucleus='1H')
            assert False, '应抛 ProcessedReadError'
        except ProcessedReadError as e:
            assert e.stage == 'check-1r', e.stage
            assert e.path == d
    print('PASS test_read_spectrum_missing_1r_raises')


# ---- 真实数据读取 (clxh3 缺失则跳过) ----

def test_read_real_clxh3_if_available():
    if not os.path.isdir(_CLXH3):
        print('SKIP test_read_real_clxh3_if_available (clxh3 不存在)')
        return
    from nmr_elucidator.io_bruker import read_bruker_folder
    specs = read_bruker_folder(_CLXH3)
    out = read_processed_folder(_CLXH3, specs=specs)
    # clxh3 有 4 个实验 21/22/23/24
    assert set(out.keys()) == {21, 22, 23, 24}, sorted(out.keys())
    # 核种标注应从 specs 正确传播 (不是 '?')
    assert out[21].nucleus == '1H', out[21].nucleus
    assert out[22].nucleus == '13C', out[22].nucleus
    # ppm 轴与强度等长, 频率 > 0
    for exp_no, ps in out.items():
        assert ps.ppm.shape[0] == ps.intensity.shape[0]
        assert ps.sf_mhz > 0.0, (exp_no, ps.sf_mhz)
    # 1H 频率 ~500, 13C 频率 ~125
    assert 480 < out[21].sf_mhz < 520, out[21].sf_mhz
    assert 110 < out[22].sf_mhz < 140, out[22].sf_mhz
    # 13C ppm 范围应覆盖芳香/羰基区 (粗略 0..234)
    lo, hi = out[22].ppm_range()
    assert lo < 5 and hi > 200, (lo, hi)
    print('PASS test_read_real_clxh3_if_available')


if __name__ == '__main__':
    test_resample_exact_grid_points()
    test_resample_handles_descending_ppm()
    test_resample_out_of_range_filled_zero()
    test_resample_negative_clipped_to_zero()
    test_resample_normalize_to_max_one()
    test_resample_length_and_dtype()
    test_resample_all_zero_intensity_no_div0()
    test_resample_length_mismatch_raises()
    test_resample_bad_npoints_raises()
    test_read_folder_nonexistent_raises()
    test_read_spectrum_missing_1r_raises()
    test_read_real_clxh3_if_available()
    print('\nALL PASSED')
