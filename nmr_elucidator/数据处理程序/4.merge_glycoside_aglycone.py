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

SRC = r'Y:\LXH数据处理\新化合物处理\数据格式处理\2026.11.16_数据格式修改（E环芳构化）.xlsx'
DST = os.path.splitext(SRC)[0] + '_合并归属表.xlsx'

TYPE_MAP = {'C': 's', 'CH': 'd', 'CH2': 't', 'CH3': 'q', 'OCH3': 'q'}
CH2_TAIL = re.compile(r'[ab\u03b1\u03b2]$')      # CH2 行的 α/β 或 a/b 构型后缀
# CH2 两行 H 数值相同时的输出方式 (可由命令行参数 --single-2h / --two-rows 覆盖):
#   False (默认): 输出两行相同的 H, 不标 2H;
#   True: 输出一行并标注 2H (如 '1.67, (2H, m)'), 与原程序一致.
SINGLE_2H = False
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
    # 末尾兼容 CH2 的构型字母 a/b (如 '5′a'), 统一在 norm_no_ch2 中转成 α/β
    m = re.search(r"(\d+(?:\u2032)*(?:[\u03b1\u03b2ab])?)$", s)
    if m is None:
        raise ValueError(f'No. 无法解析: {v!r}')
    return m.group(1)


def norm_no_ch2(v) -> str:
    """norm_no 的扩展: CH2 行的 a/b 构型写法统一为 α/β (a→α, b→β),
    使同位点的两行按 α<β 规则排序 (a 行第一、b 行第二)."""
    s = norm_no(v)
    if CH2_TAIL.search(s):
        s = s[:-1] + {'a': '\u03b1', 'b': '\u03b2'}.get(s[-1], s[-1])
    return s


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
        t = re.search(r'[\u03b1\u03b2]|(?<![A-Za-z])[ab](?![A-Za-z])', p)
        tags.append(t.group(0) if t else None)
    if '\u03b1' in tags or 'a' in tags:
        a_i = tags.index('\u03b1') if '\u03b1' in tags else tags.index('a')
    elif '\u03b2' in tags or 'b' in tags:
        a_i = 1 - (tags.index('\u03b2') if '\u03b2' in tags else tags.index('b'))
    else:
        a_i = 0 if (extract_num(parts[0]) or 0) >= (extract_num(parts[1]) or 0) else 1
    return parts[a_i], parts[1 - a_i]


def num_key(s):
    """排序键: 无化学位移数值的片段排到最后, 不与 float 混比."""
    v = extract_num(s)
    return (v is None, -(v if v is not None else 0.0))


def fmt_h(s, keep_2h=False):
    """δH 归一为示例风格 'ppm, mult (J, J)'.

    '4.85 (d, J = 7.8 Hz, 1H)' / '4.97, d (7.2)' / '1.87' / '1.46, (2H, m)'
    → '4.85, d (7.8)' / '4.97, d (7.2)' / '1.87' / '1.46, (2H, m)'.
    所有 'J = ' 字样去除, J 值直接并入括号; 其余氢数标记 (1H/3H 等) 全部去除;
    逗号后统一一个空格. 2H 标记 (写成 '(2H, mult)') 默认不保留, 仅
    keep_2h=True (单行 2H 输出模式) 时保留."""
    s = (s.replace('\uff08', '(').replace('\uff09', ')')
          .replace('\uff0c', ',').replace('\uff1b', ';')
          .replace('\u3000', ' '))
    # 氢数标记统一去掉: '(2H, m)' / 'm; 2H' / 's, 3H' / 'Hz, 1H' 等 → 去掉 nH 本身;
    # 仅在 keep_2h 模式 (单行 2H 输出) 下 2H 标记先保留, 最后统一改写成 '(2H, mult)'
    if keep_2h:
        s = re.sub(r',?\s*(?!2H\b)\dH\b', '', s)
    else:
        s = re.sub(r',?\s*\(?\dH[,)]?', '', s)
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
    # ('1.42, α, m' → '1.42, m'; '2, H' → '2'); 剩余多重峰列表里的也清掉
    rest = re.sub(r'(?<![A-Za-z])[\u03b1\u03b2](?![A-Za-z])', '', rest)
    rest = re.sub(r'\s*,\s*', ', ', rest).strip(' ,')
    rest = re.sub(r'^[\s,]*H(?=[\s,]|$)[\s,]*', '', rest)
    if not rest:
        return ppm
    if not keep_2h:
        rest = re.sub(r',?\s*\(?\dH[,)]?', '', rest).strip(' ,')   # 去掉 nH 标记
        if not rest:
            return ppm
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
    # ('m; 2H' / 'd, 6.6, 2H' → '(2H, m)' / '(2H, d, 6.6)'); 仅 keep_2h 模式保留
    if (keep_2h and re.search(r'(?<![\d.])2H\b', rest)
            and not re.fullmatch(r'[\u03b1\u03b2]', rest.split(' (')[0])):
        rest = rest.replace('2H', '', 1).strip(' ,')
        rest = re.sub(r'\s*,\s*', ', ', rest).strip(' ,')
        rest = f'(2H, {rest})' if rest else '(2H)'
    return f'{ppm}, {rest}'


def no_rich(no: str):
    """No. 单元格值 rich text: 糖基连接位 'O-Glc' 中的 O 斜体;
    α/β 用斜体, CH2 的 a/b 构型后缀用正体 (不斜体); 其余部分 (含 ′) 用正体.
    无特殊部分时直接存普通字符串."""
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
            out.append(TextBlock(IF_N, part))   # a/b 后缀走这里, 正体
    if len(out) == 1 and not out[0].font.i:
        return no                               # 无斜体部分, 存普通字符串
    return CellRichText(*out)


def _is_blank(v) -> bool:
    """单元格是否为空 (NaN / None / 空白串 / 字面量 'nan')."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return True
    s = str(v).strip()
    return s == '' or s.lower() == 'nan'


def block_rows(df, c0, c1, sheet=''):
    """把一个块 (列区间 [c0, c1)) 转成归属表行; 返回 list[dict].
    仅 CH2 (t) 类型保留 δH 中的 2H 标记; 其它类型 (CH3 等) 的 nH 标记一律去除.
    sheet 仅用于报错时指明来源, 便于人工排查.

    预处理 (参考 assign_from_mnova.py 的续行处理):
    - ppm 为空但 1H-NMR 有内容的行是上一行的续行 (老格式 CH2 拆两行), ppm/类型/Pos.
      用上一行补齐; Pos. 重复的续行再依次补 α/a、β/b, 避免两行 Pos. 相同被合并;
    - 若有 '构型' 列, 该行 Pos. 与该行构型字符串相加 (如 '6′' + 'b' → '6′b');
    - 同一 CH2 (同 Pos. 基名) 的两行分为一组: α/a 行在前、β/b 行在后,
      δC/type 只保留第一行."""
    header = df.iloc[1, c0:c1].tolist()
    ppm_i, type_i, pos_i, h_i = find_block_cols(header)
    names = [str(h).strip() if h is not None else '' for h in header]
    cfg_i = next((i for i, h in enumerate(names) if h == '构型'), None)
    body = df.iloc[2:, c0:c1]

    # ---- 1) 续行补齐: 空缺的 ppm/类型/Pos. (及构型) 用上一行数据填充 ----
    recs = []
    last = None
    for idx, r in body.iterrows():
        ppm_v = r.iloc[ppm_i]
        typ_v = r.iloc[type_i]
        pos_v = r.iloc[pos_i]
        h = '' if _is_blank(r.iloc[h_i]) else str(r.iloc[h_i]).strip()
        cfg = ('' if cfg_i is None or _is_blank(r.iloc[cfg_i])
               else str(r.iloc[cfg_i]).strip())
        # 各列是否为空 (统一按空串/NaN 判断, 不依赖 pandas NaN 取值形态)
        m_ppm, m_typ = _is_blank(ppm_v), _is_blank(typ_v)
        m_pos = _is_blank(pos_v)
        blank_pos = m_pos                # 原表 Pos. 是否为空 (是否为"缺位续行")
        blank = m_ppm or m_typ or m_pos
        if blank:
            if last is None:
                continue                  # 块开头空行
            if not h:
                continue                  # 尾部/中间空行
            # 续行: 用上一行补齐 (构型列每个续行自带, 不补齐)
            if m_ppm:
                ppm_v = last['ppm']
            if m_typ:
                typ_v = last['typ_raw']
            if m_pos:
                pos_v = str(last['pos']).strip()   # 只取基名, 后缀由组装步骤决定
                pos_v = CH2_TAIL.sub('', pos_v)
        pos = str(pos_v).strip()
        typ_raw = str(typ_v).strip()
        if typ_raw in TYPE_MAP:
            typ = TYPE_MAP[typ_raw]
        elif 'OCH3' in pos.upper():
            typ = 'q'
        else:
            # Excel 行号 = 表头占 2 行 + pandas 数据行序号
            raise ValueError(f'sheet "{sheet}" 第 {idx + 1} 行: '
                             f'Type列存在非法取值: {typ_raw!r} (Pos.={pos}, 1H-NMR={h!r})')
        rc = {'idx': idx, 'ppm': round(float(ppm_v), 1), 'typ_raw': typ_raw,
              'typ': typ, 'pos': pos, 'h': h, 'cfg': cfg, 'cont': blank,
              'raw_pos': blank_pos}
        # 补齐后 Pos. 与上一行相同 (老格式把 Pos. 写了两遍, 或上一行 Pos. 是空的):
        # 按同 Pos. 出现次数依次补 α、β, 否则两行同 Pos. 会被当成同一行合并/覆盖.
        # 只有原 Pos. 确实写了 (与上一行重复) 才补; 原本就缺 Pos. 的续行交给组装步骤
        # 无条件并入上一碳, 保持原表的 α/β 先后
        if (blank and typ == 't' and last is not None and blank_pos
                and pos == CH2_TAIL.sub('', str(last['pos']).strip())):
            n_same = sum(1 for x in recs if x['pos'] == pos)
            rc['pos'] = CH2_TAIL.sub('', pos) + ('\u03b1' if n_same == 0 else '\u03b2')
        recs.append(rc)
        last = rc

        # ---- 构型列: 原表写了自己的 Pos. 且带 '构型' 时, Pos. 与构型相拼
        #      (如 '6′' + 'b' → '6′b'); 缺 Pos. 的续行与已补后缀的行不再重复拼 ----
        if cfg and not blank_pos and not CH2_TAIL.search(pos):
            pos = pos + cfg
            rc['pos'] = pos

    # ---- 3) 组装为归属表行; 同一 CH2 的两行合成一组 ----
    # 按原表行序输出: 每个 CH2 组在其首行位置展开成两行; 其余行原位输出.
    items = []                            # (行序号, 子序号, 行/组)
    groups = {}                           # pos 基名 → CH2 组
    last_group = None                     # 最近一个 CH2 组 (续行并入它)
    for i, rc in enumerate(recs):
        pos, h, typ = rc['pos'], rc['h'], rc['typ']
        base = CH2_TAIL.sub('', pos)
        # 行号后缀: 有 '构型' 列用 a/b 或 α/β (按该行构型), 无则默认 α/β
        style = _ch2_style(rc['cfg']) if rc['cfg'] else '\u03b1\u03b2'
        if typ == 't' and has_two_ppm(h):
            # 同一段里写两个位移 ('2.69(β), m; 2.44(α) (d, ...)'): 拆成 α/β 两行
            a_part, b_part = split_ab(h)
            items.append((i, 0, {'No.': pos + '\u03b1', 'δH': fmt_h(a_part),
                                 'δC': rc['ppm'], 'type': typ}))
            items.append((i, 1, {'No.': pos + '\u03b2', 'δH': fmt_h(b_part),
                                 'δC': None, 'type': None}))
        elif typ == 't' and rc['cont'] and last_group is not None \
                and base == last_group['pos']:
            # 缺 ppm/类型的续行: 无条件并入上一碳 (作为第二行), 不另开新组
            letter = rc['cfg']
            if letter in ('\u03b2', 'b') and last_group['\u03b1'] is None:
                last_group['\u03b1'] = h      # 原表先写 β 后写 α
            elif last_group['\u03b2'] is None:
                last_group['\u03b2'] = h
            else:
                last_group['\u03b1'] = last_group['\u03b1'] or h
            last_group['_styles'].add(style)
        elif typ == 't':
            # CH2 首行: 带 a/b/α/β 后缀或 '构型' 的按该字母归为 α 或 β
            letter = rc['cfg'] or pos[len(base):]
            slot = '\u03b2' if letter in ('b', '\u03b2') else '\u03b1'
            g = groups.get(base)
            if g is None:
                g = groups[base] = {'pos': base, 'ppm': rc['ppm'],
                                    '\u03b1': None, '\u03b2': None,
                                    '_styles': set()}
                items.append((i, 0, g))
            if g['ppm'] is None:
                g['ppm'] = rc['ppm']
            g['_styles'].add(style)
            if g[slot] is None:
                g[slot] = h
            last_group = g
        else:
            items.append((i, 0, {'No.': pos, 'δH': fmt_h(h),
                                 'δC': rc['ppm'], 'type': typ}))

    rows = []
    for _, _, val in items:
        if '_styles' in val:              # CH2 组: 展开成 1~2 行
            rows += _emit_ch2(val, SINGLE_2H)
        else:
            rows.append(val)
    return rows


def _ch2_style(letter: str) -> str:
    """构型字母对应的行号后缀样式: a/b → 'ab', α/β (或无) → 'αβ'."""
    return 'ab' if letter in ('a', 'b') else '\u03b1\u03b2'


def _emit_ch2(g, single_2h: bool = False):
    """输出一个 CH2 组.
    - 组内两行 H 不同: 输出两行, 第一行 (α/a) 保留 δC/type, 第二行 (β/b) 留空;
    - 组内两行 H 相同 (或只有一行):
        single_2h=True  输出一行并标注 2H (如 '2.63, (2H, m)');
        否则输出两行相同的 H, 且都不带 2H 标记 (如 '2.63, m'), 第二行 δC/type 留空;
    - 行号后缀: 组内出现过 a/b 时用 a/b, 否则用 α/β (α/a 行在前, β/b 行在后)."""
    styles = g.get('_styles') or set()
    style = 'ab' if styles and styles <= {'ab'} else '\u03b1\u03b2'
    alpha = fmt_h(g['\u03b1'] or '')
    same = g['\u03b2'] is None or fmt_h(g['\u03b2']) == alpha
    if same and single_2h:
        # 只输出一行: '2.63, m' → '2.63, (2H, m)'
        return [{'No.': g['pos'], 'δH': _as_2h(fmt_h(g['\u03b1'] or '', keep_2h=True)),
                 'δC': g['ppm'], 'type': 't'}]
    beta = alpha if same else fmt_h(g['\u03b2'])
    return [
        {'No.': g['pos'] + style[0], 'δH': alpha,
         'δC': g['ppm'], 'type': 't'},
        {'No.': g['pos'] + style[1], 'δH': beta,
         'δC': None, 'type': None},
    ]


def _as_2h(h: str) -> str:
    """单个 H 文本标注为 2H: '1.67, m' → '1.67, (2H, m)'; 本就含 2H 的原样返回."""
    m = re.match(r'^(\d+(?:\.\d+)?)(.*)$', h)
    if not m or re.search(r'(?<![\d.])2H\b', h):
        return h
    rest = m.group(2).strip().lstrip(',').strip()
    return f'{m.group(1)},{f" (2H, {rest})" if rest else " (2H)"}'


def blank_ch2_second(rows):
    """最终输出前: 同一 CH2 (同 No. 基名) 的两行, 第二行 (β/b) 的 δC 与 type 置空.
    参考 assign_from_mnova.py 的 _blank_dup_carbon: 仅对 CH2 生效,
    组内行序即显示行序, 除首行外全部置空."""
    groups = {}
    for r in rows:
        if str(r['type']).strip() != 't':
            continue
        groups.setdefault(CH2_TAIL.sub('', r['_plain']), []).append(r)
    for grp in groups.values():
        for r in grp[1:]:
            r['δC'] = None
            r['type'] = None


def check_input(df, bounds, name):
    """输入核查 (只读, 不改动 df / 输入文档): 用原始表的 类型/Pos./构型 列判断输入是否有误.

    规则: 只有 CH2 允许在 Pos. 末尾带构型后缀 (α/β 或 a/b);
    其余类型 (C/CH/CH3/OCH3) 带构型后缀即为输入错误.
    续行 (ppm/类型/Pos. 为空) 的构型是正常写法, 不检查."""
    issues = []
    header = df.iloc[1].tolist()
    for c0, c1 in zip(bounds[:-1], bounds[1:]):
        names = [str(h).strip() if h is not None else '' for h in header[c0:c1]]
        try:
            _, type_i, pos_i, _ = find_block_cols(header[c0:c1])
        except KeyError:
            continue
        cfg_i = next((i for i, h in enumerate(names) if h == '构型'), None)
        for idx, r in df.iloc[2:, c0:c1].iterrows():
            if _is_blank(r.iloc[type_i]) or _is_blank(r.iloc[pos_i]):
                continue                     # 续行: 构型由上一行继承, 属 CH2 正常写法
            typ = str(r.iloc[type_i]).strip()
            pos = str(r.iloc[pos_i]).strip()
            cfg = ('' if cfg_i is None or _is_blank(r.iloc[cfg_i])
                   else str(r.iloc[cfg_i]).strip())
            has_tail = pos[-1:] in ('\u03b1', '\u03b2', 'a', 'b')
            if typ == 'CH2' or not (has_tail or cfg):
                continue                     # CH2 允许构型后缀; 未带后缀正常
            issues.append(f'第 {idx + 1} 行: 类型 {typ} 的 Pos. 不应带构型后缀/构型值 '
                          f'(Pos.={pos!r}, 构型={cfg!r})')
    return issues


def check_rows(rows):
    """输出前逐行检查待输出内容是否符合期待, 用来发现"输入数据"的问题.
    只读检查, 不修改 rows / 输入文档 / 程序逻辑; 返回问题列表 (str), 写入 LOG 说明.

    检查项 (参考示例归属表 clxh4 的约定):
    - No.: 必须能解析出位置编号;
    - δC: 有值时应为 0~250 的数值, 且同一 No. 不应重复出现 δC;
    - δH: 有值时应以化学位移数值开头 (0~20), 括号必须配对, 位置编号后应有空白;
    - type: ∈ {s, d, t, q} 或置空;
    - CH2: 两行写法时 α/a 行须有 δC+type、β/b 行两者都应置空; 一行写法须同时有 δC+type;
    - 连续性: 排序后每位点编号应连续递增 (按 ′ 数量分组, 无 ′ 与 1 个 ′ 视为同一组),
      缺号即为输入遗漏; 组内编号也不应重复 (如 '7′α' 与 '7β' 混用)."""
    issues = []
    ppm_of = {}                            # 去后缀的 No. → 首个 δC, 用于查重
    ch2 = {}                               # 去后缀的 No. → [α 行, β 行...]
    # 位置编号序列: 去后缀、去糖基前缀后的 (′ 数量, 数字); α/β 两行属同一位置, 只记一次.
    # 按 ′ 数量分组 [苷元(无'), 1 个 ', 2 个 ', ...] 逐组检查连续性
    seq = {}                               # ′数量 → {数字: 首次出现的 No.}
    for r in rows:
        no = r.get('_plain', '')
        dh = str(r.get('δH') or '')
        dc = r.get('δC')
        ty = r.get('type')
        ty_s = '' if ty is None else str(ty).strip()
        if not no:
            issues.append('No. 为空')
            continue
        if no_sort_key(no) is None:
            issues.append(f'{no}: No. 无法解析位置编号 (该行已排到最后)')
        else:
            # ---- 连续性采集: 去掉 α/β/a/b 后缀, 数字 + 末尾 ′ 数量 ----
            stripped = CH2_TAIL.sub('', no)
            m = re.search(r'(\d+)(\u2032*)$', stripped)
            if m:
                primes = len(m.group(2))
                seq.setdefault(primes, {}).setdefault(
                    int(m.group(1)), no)
        # ---- δC ----
        if dc is not None:
            if not isinstance(dc, (int, float)) or not 0 <= float(dc) <= 250:
                issues.append(f'{no}: δC 数值不合理 ({dc!r})')
            key = CH2_TAIL.sub('', no)
            if key in ppm_of:
                issues.append(f'{no}: No. 重复 (δC={dc}, 首条 δC={ppm_of[key]})')
            else:
                ppm_of[key] = dc
        # ---- δH ----
        if dh:
            v = extract_num(dh)
            if v is None:
                issues.append(f'{no}: δH 不以化学位移数值开头 ({dh!r})')
            elif not 0 <= v <= 20:
                issues.append(f'{no}: δH 数值超出 0~20 ({dh!r})')
            if dh.count('(') != dh.count(')'):
                issues.append(f'{no}: δH 括号不配对 ({dh!r})')
            if not re.match(r'^\d+(?:\.\d+)?', dh):
                issues.append(f'{no}: δH 缺少化学位移后空格/格式 ({dh!r})')
        elif ty_s and ty_s != 's':
            issues.append(f'{no}: type={ty_s} 但 δH 为空 (应注意核对积分/归属)')
        # ---- type ----
        if ty_s and ty_s not in ('s', 'd', 't', 'q'):
            issues.append(f'{no}: type 取值异常 ({ty_s!r})')
        if ty_s == 't':
            ch2.setdefault(CH2_TAIL.sub('', no), []).append(r)
    # ---- CH2 两行: 第一行保留 δC/type, 后续行必须置空 ----
    for key, grp in ch2.items():
        if not grp:
            continue
        if len(grp) > 2:
            issues.append(f'{key}: CH2 出现 {len(grp)} 行 (应为 1~2 行)')
        first = grp[0]
        if first.get('δC') is None or not str(first.get('type') or '').strip():
            issues.append(f'{key}: CH2 首行缺 δC 或 type')
        for r in grp[1:]:
            if r.get('δC') is not None or r.get('type') is not None:
                issues.append(f'{key}: CH2 第二行 δC/type 未置空 '
                              f'(δC={r.get("δC")!r}, type={r.get("type")!r})')
        if len(grp) == 2:
            h1 = str(grp[0].get('δH') or '')
            h2 = str(grp[1].get('δH') or '')
            if not h1 or not h2:
                issues.append(f'{key}: CH2 两行中有一行 δH 为空')
    # ---- 编号连续性: 按 ′ 数量分组 ([苷元无 ', 1 个 ', 2 个 ', ...]),
    #      组内编号必须从 1 起连续递增; 缺号即为输入遗漏 ----
    for primes in sorted(seq):
        nums = seq[primes]
        if not nums:
            continue
        lo, hi = min(nums), max(nums)
        missing = sorted(set(range(lo, hi + 1)) - set(nums))
        if missing:
            tag = {0: ' (苷元, 无 ′)', 1: ' (1 个 ′)'}.get(
                primes, f' ({primes} 个 ′)')
            issues.append(f'位点编号不连续{tag}: 缺少 '
                          + ', '.join(str(x) for x in missing)
                          + f' (现有 {lo}~{hi})')
    return issues


def sort_key(r):
    """排序键: 可排序的按 (′数量, 数字, α/β) 升序; 无法解析的排到最后."""
    k = no_sort_key(r['_plain'])
    return (1, (0, 0, 0)) if k is None else (0, k)


def main():
    global SINGLE_2H
    import argparse
    ap = argparse.ArgumentParser(description='合并糖基/苷元块为归属表')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--single-2h', action='store_true',
                   help="CH2 两行 H 相同时输出一行并标注 2H (如 '1.67, (2H, m)')")
    g.add_argument('--two-rows', action='store_true',
                   help='CH2 两行 H 相同时输出两行相同的 H (默认行为)')
    args = ap.parse_args()
    if args.single_2h:
        SINGLE_2H = True
    print(f'CH2 两行 H 相同: {"单行 + 2H 标注" if SINGLE_2H else "两行相同 H"}')

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
            rows += block_rows(df, c0, c1, name)
        # 输入核查 (在改写 No. 之前): 只有 CH2 允许带构型后缀; 其余类型带后缀说明输入有误
        raw_issues = check_input(df, bounds, name)

        # No. 归一并排序 (糖基块与苷元块整体一起排, 参考 1.sort_position.py 规则).
        # 同一位点 a/b 与 α/β 并存时, 构型字母不影响排序键 (α/a 行在前, β/b 行在后),
        # 靠稳定排序保持组装时的先后 (见 _emit_ch2).
        for r in rows:
            r['_pos'] = norm_no_ch2(r['No.'])     # 先留纯文本供排序, 再转 rich text
            r['_plain'] = r['_pos']
        rows.sort(key=sort_key)
        for r in rows:
            r['No.'] = no_rich(r['_pos'])
        bad = [r for r in rows if no_sort_key(r['_plain']) is None]
        if bad:
            msg = (f'{name}: "No." 无法排序的 {len(bad)} 行已排到最后: '
                   + '; '.join(r['_plain'] for r in bad))
            log.append(msg)
            print(msg)

        ws = wb_out.create_sheet(f'{name}归属表')
        ws.append(HDR)
        # 输出前最后一步: 同一位点 CH2 两行中的第二行 (β/b 行) 的 ppm/类型列置空
        blank_ch2_second(rows)
        # 输出前逐行检查 (只读, 不修改任何内容): 发现的问题只写进 LOG 说明
        issues = raw_issues + check_rows(rows)
        if issues:
            header = f'{name}: 待输出内容有 {len(issues)} 处异常 (疑似输入数据问题)'
            log.append(header)
            log += [f'  - {x}' for x in issues]
            print(header)
            for x in issues:
                print(f'  - {x}')
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
    # 说明文件: 同时输出到终端并另存 txt (与源文件同目录)
    LOG = os.path.splitext(SRC)[0] + '_说明.txt'
    with open(LOG, 'w', encoding='utf-8') as f:
        f.write('输出前逐行检查结果 (仅记录, 未修改程序与输入文档)\n')
        f.write('=' * 60 + '\n')
        if log:
            f.write('\n'.join(log) + '\n')
        else:
            f.write('所有 sheet 的待输出内容均符合期待, 未发现输入数据问题.\n')
    print(f'说明文件: {LOG}')


if __name__ == '__main__':
    main()
