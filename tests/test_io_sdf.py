"""单元测试: io_sdf 自定义 SDF 解析器."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from nmr_elucidator.io_sdf import _split_records, _parse_tags, _parse_spectrum_value


# 构造一个最小 SDF 字符串 (含 MOL block + tags)
SAMPLE_SDF = """\

  CDKD

  8 8  0  0  0  0  0  0  0  0999 V2000
   -2.4127   -1.4040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2376   -1.3998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4965   -2.1831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8315   -2.6713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1618   -2.1898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3666   -1.9706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4499   -2.5998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8331   -0.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  1  1  0  0  0  0
  1  8  1  0  0  0  0
M  END
> <INChI key>
TEST-KEY-123

> <INChI>
InChI=1S/C8H12O

> <nmrshiftdb2 ID>
12345

> <Spectrum 13C 0>
17.6;0.0Q;8|18.3;0.0T;0|22.6;0.0Q;12|26.5;0.0T;6|

> <Spectrum 1H 0>
1.43;0.0;0|2.34;0.0;1|1.63;0.0;2|

> <Solvent>
0:Chloroform-D1 (CDCl3)

> <Field Strength [MHz]>
0:50.328

"""


def test_split_records_with_separator():
    """两条记录用 $$$$\\r\\n 分隔."""
    raw = (SAMPLE_SDF + "$$$$\r\n" + SAMPLE_SDF + "$$$$\r\n").encode('utf-8')
    parts = list(_split_records(raw))
    assert len(parts) == 2


def test_parse_tags_extracts_fields():
    mol_block, tags = _parse_tags(SAMPLE_SDF)
    assert 'M  END' in mol_block
    assert tags.get('INChI key') == 'TEST-KEY-123'
    assert tags.get('nmrshiftdb2 ID') == '12345'
    assert tags.get('Solvent') == '0:Chloroform-D1 (CDCl3)'
    assert 'Spectrum 13C 0' in tags
    assert 'Spectrum 1H 0' in tags


def test_parse_spectrum_value_13c():
    """13C 谱字段解析."""
    val = "17.6;0.0Q;8|18.3;0.0T;0|22.6;0.0Q;12|"
    sigs = _parse_spectrum_value(val)
    assert len(sigs) == 3
    assert sigs[0].shift == 17.6
    assert sigs[0].multiplicity == 'Q'
    assert sigs[0].atom_index == 8
    assert sigs[1].shift == 18.3
    assert sigs[1].multiplicity == 'T'
    assert sigs[2].atom_index == 12


def test_parse_spectrum_value_1h():
    """1H 谱字段解析 (可能无 multiplicity)."""
    val = "1.43;0.0;0|2.34;0.0;1|"
    sigs = _parse_spectrum_value(val)
    assert len(sigs) == 2
    assert sigs[0].shift == 1.43
    assert sigs[0].atom_index == 0
    assert sigs[1].shift == 2.34


def test_parse_spectrum_value_empty():
    """空值."""
    assert _parse_spectrum_value('') == []