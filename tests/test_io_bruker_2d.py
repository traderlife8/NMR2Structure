"""2D Bruker 读取测试.

验证:
  - clxh100 的 2D 数据读取 (COSY/HSQC/HMBC/ROESY)
  - 实验类型检测
  - 峰检测
  - clxh3 (无 2D) 回退到空 dict
"""
import os
import sys

# 项目路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from nmr_elucidator.io_bruker_2d import (
    Bruker2DReadError,
    detect_2d_experiment_type,
    is_2d_pulse_program,
    pick_peaks_2d,
    read_bruker_2d_experiment,
    read_bruker_2d_folder,
)

# 测试数据路径
CLXH100_FOLDER = r'E:\NMR2Structure\NMR\clxh100'
CLXH3_FOLDER = r'E:\NMR2Structure\NMR\clxh3'


def test_detect_experiment_type():
    """测试脉冲程序 -> 实验类型检测."""
    assert detect_2d_experiment_type('cosygpmfqf', '1H', '1H') == 'COSY'
    assert detect_2d_experiment_type('hsqcetgpsi2', '13C', '1H') == 'HSQC'
    assert detect_2d_experiment_type('hmbcgplpndqf', '13C', '1H') == 'HMBC'
    assert detect_2d_experiment_type('croesyph', '1H', '1H') == 'ROESY'
    assert detect_2d_experiment_type('noesygppr1d', '1H', '1H') == 'NOESY'
    assert detect_2d_experiment_type('hmqcetgpsi', '13C', '1H') == 'HSQC'
    assert detect_2d_experiment_type('dipsi2ph', '1H', '1H') == 'TOCSY'
    # inv4 系列
    assert detect_2d_experiment_type('inv4gplplrndqf', '13C', '1H') == 'HMBC'
    assert detect_2d_experiment_type('inv4gpsi', '13C', '1H') == 'HSQC'
    print('  test_detect_experiment_type: OK')


def test_is_2d_pulse_program():
    """测试 2D 脉冲程序白名单."""
    assert is_2d_pulse_program('cosygpmfqf') is True
    assert is_2d_pulse_program('hsqcetgpsi2') is True
    assert is_2d_pulse_program('hmbcgplpndqf') is True
    assert is_2d_pulse_program('croesyph') is True
    assert is_2d_pulse_program('zg') is False
    assert is_2d_pulse_program('zgdc') is False
    assert is_2d_pulse_program('deptsp135') is False
    print('  test_is_2d_pulse_program: OK')


def test_read_clxh100_folder():
    """测试 clxh100 的 2D 数据读取."""
    if not os.path.isdir(CLXH100_FOLDER):
        print(f'  SKIP: {CLXH100_FOLDER} 不存在')
        return
    result = read_bruker_2d_folder(CLXH100_FOLDER, quiet=True)
    assert len(result) == 4, f'期望 4 个 2D 实验, 实得 {len(result)}'
    assert 'COSY' in result
    assert 'HSQC' in result
    assert 'HMBC' in result
    assert 'ROESY' in result

    # HSQC 检查
    hsqc = result['HSQC']
    assert hsqc.f1_nucleus == '13C'
    assert hsqc.f2_nucleus == '1H'
    assert hsqc.exp_no == '26'
    assert len(hsqc.peaks) > 0
    assert hsqc.data.ndim == 2

    # COSY 检查
    cosy = result['COSY']
    assert cosy.f1_nucleus == '1H'
    assert cosy.f2_nucleus == '1H'
    assert cosy.exp_no == '25'

    print(f'  test_read_clxh100_folder: OK '
          f'(HSQC peaks={len(hsqc.peaks)}, COSY peaks={len(cosy.peaks)})')


def test_read_clxh3_folder_fallback():
    """测试 clxh3 (无 2D 数据) 回退到空 dict."""
    if not os.path.isdir(CLXH3_FOLDER):
        print(f'  SKIP: {CLXH3_FOLDER} 不存在')
        return
    result = read_bruker_2d_folder(CLXH3_FOLDER, quiet=True)
    assert len(result) == 0, f'clxh3 无 2D 数据, 期望空 dict, 实得 {list(result.keys())}'
    print('  test_read_clxh3_folder_fallback: OK')


def test_read_single_experiment():
    """测试单个 2D 实验读取."""
    exp_dir = os.path.join(CLXH100_FOLDER, '26')  # HSQC
    if not os.path.isdir(exp_dir):
        print(f'  SKIP: {exp_dir} 不存在')
        return
    spec = read_bruker_2d_experiment(exp_dir)
    assert spec.experiment_type == 'HSQC'
    assert spec.f1_nucleus == '13C'
    assert spec.f2_nucleus == '1H'
    assert spec.pulse_program == 'hsqcetgpsi2'
    assert spec.solvent == 'Pyr'
    assert spec.field_mhz > 400  # 500 MHz 谱仪
    assert len(spec.peaks) > 0
    # 检查 top 峰的化学合理性 (HSQC 应有 C-H 对)
    top_peak = max(spec.peaks, key=lambda p: p.intensity)
    assert 0 < top_peak.f1_ppm < 220  # 13C 范围
    assert -1 < top_peak.f2_ppm < 14  # 1H 范围
    print(f'  test_read_single_experiment: OK (top peak: 13C={top_peak.f1_ppm:.1f}, 1H={top_peak.f2_ppm:.1f})')


def test_nonexistent_folder():
    """测试不存在的目录."""
    result = read_bruker_2d_folder(r'E:\NMR2Structure\NMR\nonexistent', quiet=True)
    assert result == {}
    print('  test_nonexistent_folder: OK')


if __name__ == '__main__':
    print('=== test_io_bruker_2d ===')
    test_detect_experiment_type()
    test_is_2d_pulse_program()
    test_read_clxh100_folder()
    test_read_clxh3_folder_fallback()
    test_read_single_experiment()
    test_nonexistent_folder()
    print('=== ALL PASSED ===')
