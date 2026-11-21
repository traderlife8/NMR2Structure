# -*- coding: utf-8 -*-
"""把 '数据格式处理.xlsx' 中各 '解析' sheet 的糖基块与苷元块合并,
处理为第一个 '归属表' sheet 的格式 (No. | δH, mult. (J ) | δC | type).

规则:
- 类型参考 2.replace_type_sdtq.py 映射为小写字母: C->s, CH->d, CH2->t, CH3->q,
  Pos. 含 OCH3 时强制 q; δC 保留一位小数.
- CH2 行的 1H-NMR 含 ';' 且不含 '2H' 时拆成两行:
  数值较大的部分放第一行, No. 后加 α, 保留 δC/type;
  数值较小的部分放第二行, No. 后加 β, δC 和 type 两列置空.
- No. 列参考 1.sort_position.py: 撇号归一 (″/全角' → 两个/一个半角 ' 后再转回
  示例表用的 ′ U+2032 单字符); 糖基连接位 (如 '4′′′-Rha-1′′′′' 之类完整字符串)
  原样保留不截断, 纯位置编号才删除最后一个数字之前的所有字符串 (OCH3 保留原样);
  按 (′数量, 数字, α<β) 升序排序, 糖基连接位字符串按 (′数量, 末尾数字) 参与排序.
- No./Pos. 列名均兼容读取 (不区分大小写).
- 字符格式与示例归属表 (clxh4归属数据表) 一致: Times New Roman 16 号;
  全表撇号统一为 ′ (U+2032), α/β 为 U+03B1/U+03B2.
- δH 格式统一为 "ppm, mult (J, J)" 风格 (与示例一致): 逗号后一个空格, J 值用
  括号列出且括号前有空格, 去掉 ' Hz, 1H' 等冗余; 仅 2H 的 CH2 保留一行并标注 ', 2H, mult'.
- 表头 rich text: δ 斜体, H/C 下标, J 斜体; 全表 Times New Roman 16 号;
  δC 列数字格式 '0.0_ ', A~D 列宽/行高、左对齐与示例表一致.
- 原始文件只读, 结果另存 <同名>_合并归属表.xlsx, 每个 '解析' sheet 输出一个归属表 sheet.
"""
import os
import re

import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.styles import Font, Alignment

IF_I = InlineFont(rFont='Times New Roman', sz=16, i=True, family=1)
IF_S = InlineFont(rFont='Times New Roman', sz=16, family=1, vertAlign='subscript')
IF_N = InlineFont(rFont='Times New Roman', sz=16, family=1)

SRC = r'Y:\LXH数据处理\新化合物处理\数据格式处理\数据格式修改.xlsx'
DST = os.path.splitext(SRC)[0] + '_合并归属表.xlsx'

TYPE_MAP = {'C': 's', 'CH': 'd', 'CH2': 't', 'CH3': 'q', 'OCH3': 'q'}
FONT = Font(name='Times New Roman', size=16)
ALIGN = Alignment(horizontal='left', vertical='center')

# 表头 rich text: δ 斜体 + H 下标 (仅 H, 逗号正体) + ' mult. (' + J 斜体 + ')'
HDR_D = CellRichText(TextBlock(IF_I, 'δ'), TextBlock(IF_S, 'H'),
                     TextBlock(IF_N, ', mult. ('), TextBlock(IF_I, 'J'),
                     TextBlock(IF_N, ')'))
HDR_C = CellRichText(TextBlock(IF_I, 'δ'), TextBlock(IF_S, 'C'))
HDR = ['No.', HDR_D, HDR_C, 'type']


def norm_no(v) -> str:
    """No. 归一 (参考 1.sort_position.py _norm_prime): 先把撇号统一为半角 '
    (双撇号字符→两个 '), 再把每个 ' 转回示例表使用的 ′ (U+2032).
    糖基连接位字符串 (含 '-' 前缀, 如 '4′′′-Rha-1′′′′') 原样保留不截断;
    纯位置编号才删除最后一个数字之前的所有字符串 (OCH3 特殊值保留原样)."""
    s = str(v).strip()
    s = (s.replace('\u2033', "''").replace('\uff02', "''").replace('\"', "''")
          .replace('\u2032', "'").replace('\uff07', "'"))
    s = s.replace("'", '\u2032')
    if 'OCH3' in s.upper() or '-' in s:         # 糖基连接位整体保留
        return s
    m = re.search(r"(\d+(?:\u2032)*(?:[\u03b1\u03b2])?)$", s)
    if m is None:
        raise ValueError(f'No. 无法解析: {v!r}')
    return m.group(1)


def no_sort_key(no: str):
    """排序键: (′数量, 末尾数字, α=0/β=1), 参考 1.sort_position.py 排序规则.
    糖基连接位字符串取最后一个数字及其后的 ′ 参与排序.
    无法解析的 No. (如 '16′-COOH') 返回 None, 由调用方统一排到最后."""
    m = re.search(r"(\d+)((?:\u2032)*)(?:[\u03b1\u03b2])?$", no)
    if m is None:
        return None
    lastnum, primes = m.group(1), m.group(2)
    suffix = 1 if no.endswith('\u03b2') else 0
    return (len(primes), int(lastnum), suffix)


def find_block_cols(block_header):
    """在一个块的表头行里定位 ppm / 类型 / Pos. / 1H 列 (块内相对索引)."""
    names = [str(h).strip() if h is not None else '' for h in block_header]
    def find(cands):
        for n in cands:
            for i, h in enumerate(names):
                if h.lower() == n.lower():
                    return i
        raise KeyError(f'未找到列, 候选: {cands}, 表头: {names}')
    ppm_i = find(('ppm',))
    type_i = find(('类型', 'types'))
    pos_i = find(('pos.', 'pos', 'position', 'no.', 'no'))
    h_i = find(('ppm（h, j）', '1h-nmr', 'ppm(h, j)'))
    return ppm_i, type_i, pos_i, h_i


def extract_num(s):
    """从 1H 文本开头提取化学位移数值, 无数值返回 None."""
    m = re.match(r'^\s*(\d+(?:\.\d+)?)', s)
    return float(m.group(1)) if m else None


def has_two_ppm(h):
    """H 文本是否含两个化学位移 (用 ';' 分隔, 或同一段写 '1.2, 1.4 (m)').
    同段写法要求第二个数值紧跟第一个且后面不是 J 值列表 (如 'd, 7.8' 不算).
    '1.46, m; 2H' 这类 '; 2H' 是氢数标记 (等同 '(2H, m)'), 不是第二个位移."""
    parts = [p.strip() for p in h.split(';') if p.strip()]
    if not parts:
        return False
    if len(parts) == 1:
        return re.match(r'^\s*\d+(?:\.\d+)?\s*(?:[(\u03b1\u03b2)][^,]*?)?\s*[,;]\s*'
                        r'\d+(?:\.\d+)?\s*(?:[(\u03b1\u03b2]|$)', h) is not None
    return (len(parts) == 2
            and all(extract_num(p) is not None for p in parts)
            and not re.fullmatch(r'\dH', parts[1]))                    # '; 2H' 不拆


def split_ab(h):
    """把 CH2 双氢文本拆成 (α 段, β 段).
    - 标注了 α/β (如 '2.69(β), m; 2.44(α) (d, ...)') 时按标注: α 段在前,
      另一段即 β; 只标注了一边的, 另一边为另一段.
    - 未标注时按化学位移大小: 大值为 α, 小值为 β.
    - 两个位移写在同一段 '1.2, 1.4 (m)' 时按位置拆开, 数值后的多重峰归后一个."""
    parts = [p.strip() for p in h.split(';') if p.strip()]
    if len(parts) != 2:
        m2 = re.match(r'^\s*(\d+(?:\.\d+)?\s*(?:[(\u03b1\u03b2)]*)?)\s*[,;]\s*'
                      r'(\d+(?:\.\d+)?)\s*(.*)$', h.strip())
        if not m2:
            return None
        # 原写法把多重峰放在第一个数值之后 ('1.4 (m), 1.2'), 归给哪个不明确:
        # 两个数值都没带多重峰时按末尾片段补上, 否则只保留数值
        rest = m2.group(3).strip()
        p1, p2 = m2.group(1).strip(), m2.group(2).strip()
        if rest and not re.search(r'[a-zA-Z]', p1 + p2):
            p2 = f'{p2}{rest if rest.startswith("(") else " " + rest}'
        parts = [p1, p2]
    tags = []
    for p in parts:
        t = re.search(r'[\u03b1\u03b2]', p)
        tags.append(t.group(0) if t else None)
    if '\u03b1' in tags:
        a_i = tags.index('\u03b1')
    elif '\u03b2' in tags:
        a_i = 1 - tags.index('\u03b2')
    else:
        a_i = 0 if (extract_num(parts[0]) or 0) >= (extract_num(parts[1]) or 0) else 1
    return parts[a_i], parts[1 - a_i]


def num_key(s):
    """排序键: 无化学位移数值的片段排到最后, 不与 float 混比."""
    v = extract_num(s)
    return (v is None, -(v if v is not None else 0.0))


def fmt_h(s):
    """δH 归一为示例风格 'ppm, mult (J, J)'.

    '4.85 (d, J = 7.8 Hz, 1H)' / '4.97, d (7.2)' / '1.87' / '1.46, (2H, m)'
    → '4.85, d (7.8)' / '4.97, d (7.2)' / '1.87' / '1.46, (2H, m)'.
    所有 'J = ' 字样去除, J 值直接并入括号; 仅 2H 标记保留 (写成 '(2H, mult)'),
    其余氢数标记 (1H/3H 等) 全部去除; 逗号后统一一个空格."""
    s = (s.replace('\uff08', '(').replace('\uff09', ')')
          .replace('\uff0c', ',').replace('\uff1b', ';')
          .replace('\u3000', ' '))
    # 化学位移后的括号可能只含 α/β 标记 (如 '2.69(β), m'), 去掉该括号但保留标记
    s = re.sub(r'^(\d+(?:\.\d+)?)\s*\(\s*([\u03b1\u03b2])\s*\)',
               r'\1\2', s.strip().rstrip(';').strip())
    # 去掉 ' Hz'; 氢数标记仅保留 2H, 其余 (1H/3H 等) 去除
    s = re.sub(r'\s*Hz\b', '', s)
    s = re.sub(r',\s*(?!2H\b)\dH\b', '', s)
    # 首段化学位移与其余部分 (多重峰/J 值) 用 ', ' 连接
    m = re.match(r'^(\d+(?:\.\d+)?)[\s,]*(.*)$', s)
    if not m:
        return s
    ppm, rest = m.group(1), m.group(2).strip()
    if not rest:
        return ppm
    # 原数据的多余标记: 化学位移后、多重峰字母前单独的 α/β 或孤立 'H'
    # ('1.42, α, m' → '1.42, m'; '2, H' → '2')
    rest = re.sub(r'^[\s,\u03b1\u03b2]*[\u03b1\u03b2][\s,]*', '', rest)
    rest = re.sub(r'^[\s,]*H(?=[\s,]|$)[\s,]*', '', rest)
    if not rest:
        return ppm
    # '(2H, m)' 与 'm; 2H' 含义相同: 分号统一成逗号, 两种写法都归一为 '(2H, m)'
    rest = rest.replace(';', ',')
    rest = rest.strip('()')                     # 去掉整体包裹的括号
    # 去掉所有 'J = ' 字样 (J 值留下)
    rest = re.sub(r'\bJ\s*=\s*', '', rest)
    # 多重峰字母后接数值列表时把数值并入括号: 'd, 7.8' → 'd (7.8)';
    # 'α (d, 7.8)' 也统一成 'α, d (7.8)'
    m2 = (re.match(r'^([a-zA-Z\u03b1\u03b2]+)\s*,\s*\(?([\d.,\s]+?)\)?$', rest)
          or re.match(r'^([\u03b1\u03b2])\s*\(\s*([a-zA-Z]+)\s*,\s*([\d.,\s]+?)\s*\)$', rest))
    if m2:
        g = m2.groups()
        mult, nums = (g[0], g[2]) if len(g) == 3 else (g[0], g[1])
        vals = ", ".join(j.strip() for j in nums.split(",") if j.strip())
        rest = f'{mult} ({vals})' if len(g) == 2 else f'{g[0]}, {mult} ({vals})'
    # 逗号后统一一个空格; 保证括号配对
    rest = re.sub(r'\s*,\s*', ', ', rest).strip(' ,')
    if rest.count('(') > rest.count(')'):
        rest += ')'
    # 含 2H 标记的行统一写成 '(2H, mult)': 2H 提到括号最前
    # ('m; 2H' / 'd, 6.6, 2H' → '(2H, m)' / '(2H, d, 6.6)')
    if (re.search(r'(?<![\d.])2H\b', rest)
            and not re.fullmatch(r'[\u03b1\u03b2]', rest.split(' (')[0])):
        rest = rest.replace('2H', '', 1).strip(' ,')
        rest = re.sub(r'\s*,\s*', ', ', rest).strip(' ,')
        rest = f'(2H, {rest})' if rest else '(2H)'
    return f'{ppm}, {rest}'


def no_rich(no: str):
    """No. 单元格值 rich text: 糖基连接位 'O-Glc' 中的 O 斜体;
    α/β 用斜体; 其余部分 (含 ′) 用正体. 无特殊部分时直接存普通字符串."""
    parts = re.split(r'(?<![A-Za-z0-9])(O)(?=-)', no)
    out = []
    for part in parts:                          # part == 'O' 的段为斜体 O
        if part == 'O':
            out.append(TextBlock(IF_I, part))
        elif part.endswith('\u03b2'):
            out += [TextBlock(IF_N, part[:-1]), TextBlock(IF_I, '\u03b2')]
        elif part.endswith('\u03b1'):
            out += [TextBlock(IF_N, part[:-1]), TextBlock(IF_I, '\u03b1')]
        elif part:
            out.append(TextBlock(IF_N, part))
    if len(out) == 1 and not out[0].font.i:
        return no                               # 无斜体部分, 存普通字符串
    return CellRichText(*out)


def block_rows(df, c0, c1):
    """把一个块 (列区间 [c0, c1)) 转成归属表行; 返回 list[dict].
    仅 CH2 (t) 类型保留 δH 中的 2H 标记; 其它类型 (CH3 等) 的 nH 标记一律去除."""
    header = df.iloc[1, c0:c1].tolist()
    ppm_i, type_i, pos_i, h_i = find_block_cols(header)
    body = df.iloc[2:, c0:c1]
    rows = []
    for _, r in body.iterrows():
        pos = r.iloc[pos_i]
        if pd.isna(pos) or not str(pos).strip():
            continue                      # 块尾部空行
        pos = str(pos).strip()
        typ_raw = '' if pd.isna(r.iloc[type_i]) else str(r.iloc[type_i]).strip()
        if typ_raw in TYPE_MAP:
            typ = TYPE_MAP[typ_raw]
        elif 'OCH3' in pos.upper():
            typ = 'q'
        else:
            raise ValueError(f'类型列存在非法取值: {typ_raw!r} (Pos.={pos})')
        ppm = round(float(r.iloc[ppm_i]), 1)
        h = '' if pd.isna(r.iloc[h_i]) else str(r.iloc[h_i]).strip()

        # 注意: 预处理阶段不能先删 nH 标记, 否则 '2.69(β), m; 2.44(α) (d, J=8.2 Hz, 1H)'
        # 会丢掉末尾 ')', 影响下面按 ';' 拆分成 α/β 两行. 统一交给 fmt_h 处理, 且每个
        # 半行单独调用 fmt_h (各自去掉 nH 标记).
        if typ == 't' and has_two_ppm(h):
            # CH2 双氢拆分: α 段放第一行 (保留 δC/type), β 段放第二行 (δC/type 置空).
            # 标注了 α/β 则按标注, 未标注则大值作 α、小值作 β (见 split_ab).
            a_part, b_part = split_ab(h)
            rows.append({'No.': pos + 'α', 'δH': fmt_h(a_part),
                         'δC': ppm, 'type': typ})
            rows.append({'No.': pos + 'β', 'δH': fmt_h(b_part),
                         'δC': None, 'type': None})
        else:
            rows.append({'No.': pos, 'δH': fmt_h(h), 'δC': ppm, 'type': typ})
    return rows


def main():
    wb_in = load_workbook(SRC, read_only=True)
    parse_sheets = [n for n in wb_in.sheetnames if n.endswith('解析')]
    if not parse_sheets:
        raise SystemExit(f'ERROR: 未找到 "解析" sheet: {wb_in.sheetnames}')
    wb_in.close()

    dfs = pd.read_excel(SRC, sheet_name=parse_sheets, header=None)

    wb_out = Workbook()
    wb_out.remove(wb_out.active)
    log = []                              # 无法排序的 No. 说明, 另存到源文件同目录
    for name, df in dfs.items():
        # 块边界: 第 1 行标题含 '糖基数据' / '排除糖基后数据' 的列起点
        titles = df.iloc[0].fillna('').astype(str)
        starts = [i for i, t in enumerate(titles)
                  if '糖基数据' in t or '排除糖基后数据' in t]
        if len(starts) != 2:
            raise SystemExit(f'ERROR: {name} 未找到糖基/苷元两个块: {starts}')
        bounds = starts + [df.shape[1]]
        rows = []
        for c0, c1 in zip(bounds[:-1], bounds[1:]):
            rows += block_rows(df, c0, c1)

        # No. 归一并排序 (糖基块与苷元块整体一起排, 参考 1.sort_position.py 规则)
        for r in rows:
            r['No.'] = no_rich(norm_no(r['No.']))
        # 无法排序的 No. 排到最后 (Python 排序稳定, 保持原有先后顺序)
        for r in rows:
            r['_plain'] = str(r['No.'])   # CellRichText 的 str() 即纯文本

        def sort_key(r):
            k = no_sort_key(r['_plain'])
            return (1, (0, 0, 0)) if k is None else (0, k)
        bad = [r for r in rows if no_sort_key(r['_plain']) is None]
        rows.sort(key=sort_key)
        if bad:
            msg = (f'{name}: "No." 无法排序的 {len(bad)} 行已排到最后: '
                   + '; '.join(r['_plain'] for r in bad))
            log.append(msg)
            print(msg)
            for r in bad:
                print(f'  - {r["_plain"]}  δH={r["δH"]}  δC={r["δC"]}  '
                      f'type={r["type"]}')

        ws = wb_out.create_sheet(f'{name}归属表')
        ws.append(HDR)
        for r in rows:
            ws.append([r['No.'], no_rich(r['δH']), r['δC'], r['type']])
        # 格式与示例归属表一致: 字体/对齐/列宽/行高/δC 数字格式
        for row in ws.iter_rows():
            for c in row:
                c.font = FONT
                c.alignment = ALIGN
        for c in ('A', 'B', 'C', 'D'):
            ws.column_dimensions[c].width = {'A': 20.78, 'B': 33.0,
                                             'C': 9.55, 'D': 7.0}[c]
        ws.row_dimensions[1].height = 22.95
        for row in ws.iter_rows(min_row=2, min_col=3, max_col=3):
            row[0].number_format = '0.0_ '
        print(f'{name}: {len(rows)} 行')

    wb_out.save(DST)
    print(f'OK: {DST}')
    # 无法排序的 No. 说明: 同时输出到终端并另存 txt (与源文件同目录)
    LOG = os.path.splitext(SRC)[0] + '_无法排序No说明.txt'
    with open(LOG, 'w', encoding='utf-8') as f:
        if log:
            f.write('以下 No. 无法按位置编号排序, 已排到各自归属表的最后一行:\n')
            f.write('\n'.join(log) + '\n')
        else:
            f.write('所有 No. 均可排序, 无需特殊处理.\n')
    print(f'说明文件: {LOG}')


if __name__ == '__main__':
    main()
