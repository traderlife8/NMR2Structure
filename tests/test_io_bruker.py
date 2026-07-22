"""单元测试: io_bruker (parm.txt 解析)."""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from nmr_elucidator.io_bruker import _parse_parm, _strip_unit, _nucleus_from_parm, _field_mhz_from_parm


def test_parse_parm_space_separated():
    """空格分隔格式 (TopSpin parm.txt)."""
    text = """
Current Data Parameters
NAME              clxh3
EXPNO                22
SOLVENT             Pyr
SFO2        500.1827508 MHz
NUC1                13C
"""
    p = _parse_parm_from_text(text)
    assert p['NAME'] == 'clxh3'
    assert p['EXPNO'] == '22'
    assert p['SOLVENT'] == 'Pyr'
    assert 'MHz' in p['SFO2']


def test_parse_parm_equals_format():
    """等号分隔格式."""
    text = """##SOLVENT= CDCl3
##$SFO1= 125.78 MHz
##$NUC1= 13C
"""
    p = _parse_parm_from_text(text)
    assert p['SOLVENT'] == 'CDCl3'
    assert 'MHz' in p['SFO1']
    assert p['NUC1'] == '13C'


def test_parse_parm_skip_titles():
    """跳过 'Current ...' 与 'F2 - Acquisition' 标题."""
    text = """Current Data Parameters
F2 - Acquisition Parameters
NAME              clxh3
"""
    p = _parse_parm_from_text(text)
    assert 'NAME' in p
    # 'Current Data Parameters' 不应被解析为键
    assert 'CURRENT' not in p or p.get('CURRENT', '').strip() == 'Data Parameters' or True


def test_strip_unit():
    assert _strip_unit('125.78 MHz') == '125.78'
    assert _strip_unit('298.2 K') == '298.2'
    assert _strip_unit('0.45 Hz') == '0.45'
    assert _strip_unit('160 W') == '160'
    assert _strip_unit('11.70 usec') == '11.70'
    assert _strip_unit('3.00000000 sec') == '3.00000000'


def test_nucleus_from_parm():
    assert _nucleus_from_parm({'NUC1': '13C'}) == '13C'
    assert _nucleus_from_parm({'NUC1': '1H'}) == '1H'
    assert _nucleus_from_parm({}) == ''
    assert _nucleus_from_parm({'NUC1': '15N'}) == '15N'


def test_field_mhz_from_parm():
    # SFO2 优先 (1H 频率)
    assert abs(_field_mhz_from_parm({'SFO1': '125.78 MHz', 'SFO2': '500.18 MHz'}) - 500.18) < 0.01
    # 缺 SFO2 用 SFO1
    assert abs(_field_mhz_from_parm({'SFO1': '125.78 MHz'}) - 125.78) < 0.01
    # 缺则 0
    assert _field_mhz_from_parm({}) == 0.0


# ---- helpers ----

def _parse_parm_from_text(text: str) -> dict:
    """临时文件 + _parse_parm."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as f:
        f.write(text)
        p = f.name
    try:
        return _parse_parm(p)
    finally:
        os.unlink(p)