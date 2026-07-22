"""单元测试: retrieve._greedy_match 算法."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from nmr_elucidator.retrieve import _greedy_match


def test_perfect_match():
    """实测与库完全一致."""
    obs = [17.5, 18.0, 100.0]
    lib = [17.5, 18.0, 100.0]
    m, rms = _greedy_match(obs, lib, 2.0)
    assert m == 3
    assert rms == 0.0


def test_within_tolerance():
    """差 0.5 ppm, 容差 2.0."""
    obs = [17.5, 18.0]
    lib = [17.0, 18.5]
    m, rms = _greedy_match(obs, lib, 2.0)
    assert m == 2
    assert 0 < rms < 1.0


def test_outside_tolerance_excluded():
    """差 5 ppm, 容差 2.0 -> 不匹配."""
    obs = [17.5]
    lib = [22.5]
    m, _ = _greedy_match(obs, lib, 2.0)
    assert m == 0


def test_greedy_one_to_one():
    """贪心一对一: 同 lib 峰不能被两个 obs 复用."""
    obs = [17.0, 17.1]  # 都很接近 17.2
    lib = [17.2]
    m, _ = _greedy_match(obs, lib, 1.0)
    assert m == 1  # 只匹配一个


def test_empty():
    """空输入."""
    assert _greedy_match([], [], 2.0) == (0, 0.0)
    assert _greedy_match([17.0], [], 2.0) == (0, 0.0)
    assert _greedy_match([], [17.0], 2.0) == (0, 0.0)


def test_rms_computation():
    """RMS 计算正确."""
    obs = [17.0]
    lib = [19.0]
    m, rms = _greedy_match(obs, lib, 4.0)
    assert m == 1
    assert abs(rms - 2.0) < 1e-6  # (2.0)^2 = 4, sqrt(4/1) = 2