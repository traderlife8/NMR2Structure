"""从 MNova JSON 导出物 + 氢谱积分.txt + 数据归属.xlsx (ppm/Types) 推导 1H-NMR 归属列.

泛化算法 (不依赖特定化合物):
  1. 解析氢谱积分.txt: "δ (mult, nH), ..." → 精细 1H 信号单 [{ppm, mult, J, nH}]
     (化学家报告值, 比仪器 GSD 峰单更可靠 — 融合的多重峰只占一行)
  2. 碳槽: xlsx 同 (ppm, Types) 合并; 需求 H 数:
     - CH: 1   CH2: 2   CH3: 1 (3H 等价, 若 HSQC 实际分离出 >1 个 δH 则多输出)
     - 理论 CHn 输出 n 个不同 δH; HSQC 交叉峰不足 n 个时按实际数输出
  3. HSQC 候选: |ΔC| <= tol; 强度 <= 1.0 的 GSD 弱峰剔除
  4. 全局一对一匹配: (|ΔC|, -强度) 升序贪心; 1D 信号占用互斥;
     CH3 额外加成: 与 1D CH3 信号 (nH=3) 的匹配距离罚减
  5. 匹配后: HSQC δH 对齐到最近的 1D 1H 信号 (±tol) 报告;
     CH2 若 HSQC 只给出 1 个 δH 但 1D 信号 nH>=2, 用 HMBC/COSY 核对补充第 2 个 H
     (fallback: 同碳 HSQC 峰簇内次强峰)

用法:
  直接运行 (默认处理脚本顶部 _DEFAULT_COMPOUND 指定的化合物):
    python assign_from_mnova.py
  换化合物: 手工修改顶部 _DEFAULT_COMPOUND; 或命令行传入化合物名:
    python assign_from_mnova.py clxh101
  输入约定 (相对 _DATA_DIR / mnova 导出目录):
    <_DATA_DIR>/<化合物名>数据归属.xlsx
    <_DATA_DIR>/<化合物名>氢谱积分.txt
    <_DATA_DIR>/mnova_export/<化合物名>/   (mnovaBatchExportJson 导出的 JSON)
"""
from __future__ import annotations

import json
import os
import re
import sys

import pandas as pd

# ---- 默认输入配置 (直接运行即用; 换化合物时手工修改这三行) ----
_DEFAULT_COMPOUND = 'clxh122'   # 化合物名
_DATA_DIR = rf'Y:\LXH数据处理\新化合物处理\程序解析\{_DEFAULT_COMPOUND}解析'
_MNOVA_EXE = r'D:\Program Files\Mestrelab Research S.L\MestReNova\MestReNova.exe'
# row_match 差异容差 (ppm): |δH_程序 - δH_xlsx| <= 此值时判 '一致 (ΔδH=...)', 保留差异值
_ROW_MATCH_TOL = 0.015
# mnova JSON 输出根 = _DATA_DIR/mnova_export (qs 从文档路径动态推导同一目录, 两端一致)
_MNOVA_OUT_ROOT = os.path.join(_DATA_DIR, 'mnova_export')


def _deploy_qs(qs_source: str, qs_deployed: str) -> None:
    """部署导出脚本到 Mnova scripts 目录 (每次运行都调用).

    qs 输出根从文档自身路径动态推导 (<.mnova 所在目录>/mnova_export), 脚本中不写死任何路径,
    因此部署即原样复制; 不存在或内容与源脚本不一致时重新部署.
    """
    with open(qs_source, 'rb') as f:
        src_data = f.read()
    try:
        need_deploy = True
        if os.path.isfile(qs_deployed):
            with open(qs_deployed, 'rb') as f:
                need_deploy = f.read() != src_data
        if need_deploy:
            with open(qs_deployed, 'wb') as f:
                f.write(src_data)
            print(f'[mnova] 已部署导出脚本 → {qs_deployed}')
    except OSError as e:
        raise SystemExit(f'ERROR: 导出脚本部署失败 (需管理员权限): {e}')


def _export_mnova_json(compound: str) -> str:
    """调用 MestReNova.exe 无界面导出 .mnova → JSON 目录 (_DATA_DIR/mnova_export/<compound>). 返回 JSON 目录.

    qs 脚本由本文件同目录的 mnovaBatchExportJson.qs 部署 (每次运行校验), 输出根由 qs 从
    文档路径动态推导为 <.mnova 所在目录>/mnova_export, 即 _DATA_DIR/mnova_export.
    """
    mnova_path = os.path.join(_DATA_DIR, f'{compound}.mnova')
    if not os.path.isfile(mnova_path):
        raise SystemExit(f'ERROR: .mnova 文件不存在: {mnova_path}')
    if not os.path.isfile(_MNOVA_EXE):
        raise SystemExit(f'ERROR: 未找到 MestReNova.exe: {_MNOVA_EXE}')
    qs_source = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mnovaBatchExportJson.qs')
    qs_deployed = os.path.join(os.path.dirname(_MNOVA_EXE), 'scripts', 'mnovaBatchExportJson.qs')
    if not os.path.isfile(qs_source):
        raise SystemExit(f'ERROR: 导出脚本缺失: {qs_source}')
    _deploy_qs(qs_source, qs_deployed)
    # 每次导出全新, 避免旧数据混入
    out_dir = os.path.join(_MNOVA_OUT_ROOT, compound)
    import shutil
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    import subprocess
    cmd = [_MNOVA_EXE, '--nogui', '-w', mnova_path, '--sf', 'mnovaBatchExportJson']
    print(f'[mnova] 无界面导出 JSON: {os.path.basename(mnova_path)} ...')
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise SystemExit(f'ERROR: MestReNova 导出失败 exit={proc.returncode}\n'
                         f'stdout: {proc.stdout[-2000:] if proc.stdout else ""}\n'
                         f'stderr: {proc.stderr[-2000:] if proc.stderr else ""}')
    if not os.path.isdir(out_dir) or not any(
            fn.endswith('.json') for fn in os.listdir(out_dir)):
        raise SystemExit(f'ERROR: 导出完成但未发现 JSON 产物: {out_dir}')
    n = sum(1 for fn in os.listdir(out_dir) if fn.endswith('.json') and not fn.startswith('_'))
    print(f'[mnova] 导出 {n} 条谱 JSON → {out_dir}')
    return out_dir

# ---- 容差 (ppm) ----
HSQC_C_TOL = 1.5     # xlsx 13C 与 HSQC F1 的溶剂/校准偏差 (CDCl3 13C vs Pyr HSQC)
H_ALIGN_TOL = 0.10   # HSQC δH 对齐到 1D 1H 信号容差
DUP_TOL = 0.04       # 同碳槽内两 H 视为相同的最小间距
WEAK_INT = 1.0       # HSQC GSD 弱峰强度阈值 (低于视为噪声)
CH3_BONUS = 0.8      # CH3 匹配 1D nH=3 信号的距离罚减 (ppm)


def _typical_h(c_ppm: float, typ: str) -> float:
    """按 δC 经验区间估计该碳的典型 δH (借用时的距离基准)."""
    if c_ppm >= 110:      # 芳/烯碳
        return 7.0
    if 90 <= c_ppm < 110:  # 异头碳
        return 5.0
    if 60 <= c_ppm < 90:   # 糖/氧代碳
        return 4.0
    return 1.5             # 脂肪碳


def _round_h_in_str(s):
    """把 '4.645, m; 3.88, d (11.2)' 类字符串中的 δH (每项开头的数值) 四舍五入到 2 位小数.

    括号内的 J 值 (如 'dd (9.9, 7.7)') 是 txt 原文, 不改小数位数.
    """
    if pd.isna(s) or not str(s).strip():
        return s

    def _sub(m):
        try:
            return f'{round(float(m.group(1)), 2):.2f}'
        except ValueError:
            return m.group(0)
    # 逐项处理: 按 '; ' 拆分, 仅替换每项开头 (δH 位置) 的数字, 括号内不动
    parts = str(s).split('; ')
    fixed = []
    for part in parts:
        # δH 是项内第一个数值 (可带区间 '-'); 只替换第一个匹配
        fixed.append(re.sub(r'^(\d+\.\d+)(-\d+\.?\d*)?', _sub, part, count=1))
    return '; '.join(fixed)


# ---- 积分.txt 解析 ----

_RE_ENTRY = re.compile(
    r'(\d+\.\d+)\s*\((?:([a-z]+(?:,\s*[a-z]+)*)?),?\s*(?:J\s*=\s*([^)]+))?\s*(\d+(?:\.\d+)?)H\)')


def parse_1h_report(txt_path: str):
    """解析 1H NMR 报告行 → [{'ppm','mult','nH','J'}, ...]."""
    with open(txt_path, 'r', encoding='utf-8') as f:
        text = f.read()
    # 截取 δ 之后的部分
    m = re.search(r'δ\s*(.+)$', text, re.S)
    if m:
        text = m.group(1)
    # 去掉区间记法 "0.92 – 0.97" 保留每个端点 (端点都算信号中心候选)
    out = []
    # 逐项: "0.60 (s, 3H)" / "0.92 – 0.97 (m, 2H)" / "1.17 (t, J = 8.2 Hz, 1H)"
    for m in re.finditer(
            r'(\d+\.\d+)(?:\s*[–-]\s*(\d+\.\d+))?\s*'      # δ 或 δ-δ 区间
            r'\(([^)]*)\)', text):                          # (mult, nH) 括号
        ppm1 = float(m.group(1))
        ppm2 = float(m.group(2)) if m.group(2) else None
        inner = m.group(3)
        nH = 1
        # nH 是括号内最后一个 "数字H"; 注意勿把 "8.2 Hz" 的 J 值误当 nH
        mn = re.search(r'(\d+(?:\.\d+)?)\s*H(?![a-z])', inner.split(',')[-1].strip()) \
            if inner else None
        if mn:
            nH = float(mn.group(1))
        mult = inner.split(',')[0].strip() if inner else ''
        # J 值列表: "t, J = 8.8 Hz, 1H" → ['8.8']; "tt, J = 7.9, 3.2 Hz, 1H" → ['7.9','3.2']
        js: list = []
        jm = re.search(r'J\s*=\s*([\d.,\s]+?)\s*Hz', inner)
        if jm:
            js = [x.strip() for x in jm.group(1).split(',') if x.strip()]
        center = (ppm1 + ppm2) / 2 if ppm2 is not None else ppm1
        # 区间原样文本 (如 '1.64-1.88') 与排序键 (区间取大端, 单值取自身)
        if ppm2 is not None:
            hi, lo = max(ppm1, ppm2), min(ppm1, ppm2)
            ppm_raw = f'{hi:.2f}-{lo:.2f}'
            sort_key = hi
        else:
            ppm_raw = f'{ppm1:.2f}'
            sort_key = ppm1
        out.append({'ppm': round(center, 3), 'mult': mult, 'nH': nH, 'J': js,
                    'ppm_raw': ppm_raw, 'sort_key': sort_key})
    return out


def export_1h_nmr(txt_path: str, out_dir: str = '', compound: str = '') -> pd.DataFrame:
    """氢谱积分.txt → 积分表 DataFrame (δH, δH_ppm, type, J, nH; 按 δH 降序).

    δH = 数值列 (区间取中点, 四舍五入 2 位小数); δH_ppm = txt 原文 (区间保留 'a-b');
    J 列为 '(a, b)' 文本. 不单独落盘, 由归属表写入 <compound>_归属表.xlsx 的
    '<compound>_1H_NMR' sheet (若该文件已含此 sheet 则直接读取复用).
    返回 DataFrame (J 还原为 list).
    """
    xlsx_path = os.path.join(out_dir or os.path.dirname(txt_path),
                             f'{compound}_归属表.xlsx' if compound else '归属表.xlsx')
    sheet_name = f'{compound}_1H_NMR' if compound else '1H_NMR'
    if os.path.isfile(xlsx_path):
        try:
            xl = pd.ExcelFile(xlsx_path)
            if sheet_name in xl.sheet_names:
                df = pd.read_excel(xl, sheet_name=sheet_name, dtype={'J': str})
                # 复用: 还原内部格式 (J 列字符串 → list)
                df['J'] = df['J'].fillna('').map(
                    lambda s: [x.strip() for x in str(s).strip('()').split(',') if x.strip()])
                df['nH'] = df['nH'].astype(float)
                return df
        except (OSError, ValueError) as e:
            print(f'[积分表][WARN] 复用 {xlsx_path} 失败 ({e}), 重新解析 txt')
    sigs = parse_1h_report(txt_path)
    rows = [{'δH': round(s['ppm'], 2), 'δH_ppm': s['ppm_raw'], 'type': s['mult'],
             'J': f"({', '.join(s['J'])})" if s['J'] else '',
             'nH': int(s['nH']), '_k': s['sort_key']} for s in sigs]
    df = pd.DataFrame(rows).sort_values('_k', ascending=False).drop(columns='_k').reset_index(drop=True)
    return df.assign(J=df['J'].map(
        lambda s: [x.strip() for x in str(s).strip('()').split(',') if x.strip()]))


def _pick_2d(json_dir: str, etype: str):
    """按 experimentType 精确取 2D JSON (COSY/HMBC), 无则 None."""
    return pick_docs(json_dir, 2, etype)


def _cosy_partners(json_dir: str, h_v: float):
    """COSY: δH 的邻位 H 列表 [(h2, I)] 强度降序 (仅 I>=1.0 强峰)."""
    cosy = _pick_2d(json_dir, 'COSY')
    if cosy is None:
        return []
    agg: dict = {}
    for pk in cosy.get('peaks', []):
        a, b = pk.get('ppm1'), pk.get('ppm2')
        if a is None or b is None:
            continue
        if abs(a - h_v) <= 0.03 and abs(b - h_v) > 0.05:
            key = round(float(b), 2)
            agg[key] = max(agg.get(key, 0.0), float(pk.get('intensity', 0) or 0))
    strong = [(b, i) for b, i in agg.items() if i >= 1.0]
    strong.sort(key=lambda x: -x[1])
    return strong


def _hmbc_carbons(json_dir: str, h_v: float):
    """HMBC: δH 的 2-3 键碳列表 [(δC, I)] 强度降序."""
    hmbc = _pick_2d(json_dir, 'HMBC')
    if hmbc is None:
        return []
    cands = []
    for pk in hmbc.get('peaks', []):
        a, b = pk.get('ppm1'), pk.get('ppm2')
        if a is None or b is None:
            continue
        if abs(b - h_v) <= 0.04:
            cands.append((float(a), abs(float(pk.get('intensity', 0) or 0))))
    cands.sort(key=lambda x: -x[1])
    return cands


def third_pass_cosy_hmbc(out: pd.DataFrame, sig_df: pd.DataFrame,
                         h1d_signals, json_dir: str) -> pd.DataFrame:
    """第三步: COSY/HMBC 复核全部 H 归属 + 推断无 HSQC 交叉峰的信号.

    全表复核规则 (对每个已归属行):
    - COSY 交叉峰 = 2-3 键 H-H 耦合 → 相邻碳上的 H;
    - HMBC 交叉峰 = 2-3 键 C-H → 活泼氢 (OH/NH, 无 HSQC 峰) 也可见, 常为弱峰;
    复核项:
    1. HSQC (第二阶段已做): C-H 直接键连;
    2. COSY 一致性: 本 H 的 COSY 邻位 H 应归属到相邻碳 (δC 差 <= 3.5 ppm 内的碳槽);
       本 H 邻位列表为空且多重峰非 s/m → 标 'COSY 无邻位?';
    3. HMBC 一致性: 本 H 的 HMBC 碳应包含本碳 (2J 本身常不出现) 及邻碳 (3J);
       HMBC 碳列表含本碳 (±0.5) → 'HMBC含本碳'; 含邻位 H 的碳 → 'HMBC含邻碳'.
    未归属信号 (无 HSQC 交叉峰) 沿用推断: COSY 邻位传播 / HMBC 碳列表 / 活泼氢判据.
    返回待附加 sheet 的 DataFrame (每行 = 1 个 H 信号/未归属信号).
    """
    rows = []
    # ---- 部分 A: 全部已归属行的 COSY/HMBC 复核 ----
    for _, r in out.iterrows():
        if pd.isna(r.get('assigned_1H')) or not str(r['assigned_1H']).strip():
            continue
        try:
            c_v, h_v = float(r['ppm']), float(str(r['assigned_1H']))
        except ValueError:
            continue
        entry = {'ppm': c_v, 'Types': r['Types'], 'assigned_1H': str(r['assigned_1H']),
                 'COSY邻位': '', 'HMBC碳': '', '复核结论': ''}
        # COSY 一致性
        partners = _cosy_partners(json_dir, h_v)
        mult = str(r.get('1H-NMR', '')).split(',')[-1].strip() if pd.notna(r.get('1H-NMR')) else ''
        if partners:
            adj_txt = []
            for b, i in partners[:3]:
                c_of_b = None
                m_row = out[(out['assigned_1H'].notna())
                            & (out['assigned_1H'].astype(str).str.strip() != '')]
                for _, rr in m_row.iterrows():
                    try:
                        if abs(float(str(rr['assigned_1H'])) - b) <= 0.02:
                            c_of_b = float(rr['ppm'])
                            break
                    except ValueError:
                        continue
                # 邻位碳应与本碳相邻 (δC 差 <= 3.5 为合理邻位范围)
                adjacent_ok = c_of_b is not None and abs(c_of_b - c_v) <= 3.5
                adj_txt.append(f'{b:.2f}' + (f'(C{c_of_b:.2f})' if c_of_b else '')
                               + ('✓' if adjacent_ok else '?'))
            entry['COSY邻位'] = '; '.join(adj_txt)
        elif mult and mult not in ('s', 'm'):
            entry['复核结论'] = 'COSY 无邻位? '
        # HMBC 一致性
        hmbc_cs = _hmbc_carbons(json_dir, h_v)
        if hmbc_cs:
            entry['HMBC碳'] = '; '.join(f'{c:.2f}(I={i:.2f})' for c, i in hmbc_cs[:6])
            has_own = any(abs(c - c_v) <= 0.5 for c, i in hmbc_cs)
            entry['复核结论'] += ('HMBC含本碳' if has_own else 'HMBC不含本碳(常见, 2J未必出现)')
        else:
            entry['复核结论'] += 'HMBC无相关'
        rows.append(entry)
    # ---- 部分 B: 未归属 (无 HSQC 交叉峰) 信号的推断 ----
    used_h = [float(str(r['assigned_1H'])) for _, r in out.iterrows()
              if pd.notna(r.get('assigned_1H')) and str(r['assigned_1H']).strip()]
    for _, s in sig_df.iterrows():
        h_v = float(s['δH'])
        if any(abs(h_v - u) <= 0.025 for u in used_h):
            continue
        entry = {'ppm': '', 'Types': '(无HSQC峰)', 'assigned_1H': s['δH_ppm'],
                 'COSY邻位': '', 'HMBC碳': '', '复核结论': ''}
        partners = _cosy_partners(json_dir, h_v)
        if partners:
            adj_txt = []
            for b, i in partners[:3]:
                c_of_b = None
                m_row = out[(out['assigned_1H'].notna())
                            & (out['assigned_1H'].astype(str).str.strip() != '')]
                for _, rr in m_row.iterrows():
                    try:
                        if abs(float(str(rr['assigned_1H'])) - b) <= 0.02:
                            c_of_b = float(rr['ppm'])
                            break
                    except ValueError:
                        continue
                adj_txt.append(f'{b:.2f}' + (f'(C{c_of_b:.2f})' if c_of_b else ''))
            entry['COSY邻位'] = '; '.join(adj_txt)
        hmbc_cs = _hmbc_carbons(json_dir, h_v)
        if hmbc_cs:
            entry['HMBC碳'] = '; '.join(f'{c:.2f}(I={i:.2f})' for c, i in hmbc_cs[:6])
        if not partners and (not hmbc_cs or all(i < 0.3 for c, i in hmbc_cs)):
            entry['复核结论'] = '疑似活泼氢/交换氢 (无 HSQC 峰, COSY 无邻位, HMBC 弱)'
        elif partners and hmbc_cs:
            entry['复核结论'] = 'COSY 邻位传播: 该 H 位于邻位 H 所在碳的相邻碳上; 结合 HMBC 碳列表人工确认'
        elif hmbc_cs:
            entry['复核结论'] = '仅 HMBC 相关 (2-3 键), 参考 HMBC 碳列表推断'
        else:
            entry['复核结论'] = 'COSY/HMBC 均无相关 (无法推断, 人工判定)'
        rows.append(entry)
    return pd.DataFrame(rows)


def _write_excel_formatted(df: pd.DataFrame, path: str,
                           float2dp=(), text_cols=()):
    """写 xlsx: float2dp 列固定 2 位小数显示, text_cols 列强制文本格式."""
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        df.to_excel(xw, index=False)
        ws = xw.book.worksheets[0]
        from openpyxl.styles import Alignment
        for col_name in float2dp:
            if col_name in df.columns:
                ci = df.columns.get_loc(col_name) + 1
                for row in ws.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '0.00'
        for col_name in text_cols:
            if col_name in df.columns:
                ci = df.columns.get_loc(col_name) + 1
                for row in ws.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '@'
                        cell.alignment = Alignment(horizontal='left')


# ---- 数据加载 ----

def load_xlsx(path: str) -> pd.DataFrame:
    """读取归属 xlsx: 必需列 ppm / Types(兼容 '类型'); 1H-NMR 可选 (无则不做行级比较).

    其它所有列原样保留 (输出归属表时一并带回). CH2 行补全两种化学家写法, 保证每个
    CH2 在返回结果中恰好占 2 行 (每行=1 个 H, 由程序归属); 补全行其它列继承原行:
    1. 续行写法 (clxh98 sheet1): 第 2 个 H 写在下一行, ppm/Types 留空 → 向下填充归并;
    2. 单行写法 (clxh119): 只有 1 行 → 自动复制该行并在其下增加相同一行 (1H-NMR 留空).
    """
    df = pd.read_excel(path, header=0)
    df.columns = [str(c).strip() for c in df.columns]
    if 'ppm' not in df.columns:
        raise SystemExit(f'ERROR: xlsx 缺少必需列 ppm: {path}')
    type_col = 'Types' if 'Types' in df.columns else ('类型' if '类型' in df.columns else None)
    if type_col is None:
        raise SystemExit(f"ERROR: xlsx 缺少必需列 Types (或 '类型'): {path}")
    df = df.rename(columns={type_col: 'Types'})
    if '1H-NMR' not in df.columns:
        df['1H-NMR'] = None
    df = df.copy()
    # ---- 续行归并: ppm 留空但 1H-NMR 有值 → 是上一行的第 2 个 H ----
    # ppm/Types 向下填充; 其它列 (No. 等) 跟随原行一起填充补全
    cont = df['ppm'].isna() & df['1H-NMR'].notna() & df['1H-NMR'].astype(str).str.strip().ne('')
    if cont.any():
        for c in df.columns:
            if c == '_cont':
                continue
            df.loc[cont, c] = df[c].ffill()[cont]
    df['_cont'] = cont  # 补全/续行标记: 归属表"置空版" sheet 输出时置空 ppm/Types
    df = df.dropna(subset=['ppm']).reset_index(drop=True)
    # ---- CH2 单行补全: 同 (ppm, CH2) 只出现 1 次时, 紧跟原行插入相同一行 (1H-NMR 留空, 由程序补) ----
    ch2_count: dict = {}
    for _, r in df.iterrows():
        if str(r['Types']).strip() == 'CH2':
            key = round(float(r['ppm']), 2)
            ch2_count[key] = ch2_count.get(key, 0) + 1
    rows = []
    for _, r in df.iterrows():
        rows.append(dict(r))  # 原行 (含全部其它列)
        if str(r['Types']).strip() == 'CH2':
            key = round(float(r['ppm']), 2)
            if ch2_count[key] == 1:
                dup = dict(r)
                dup['1H-NMR'] = None
                dup['_cont'] = True  # 补全行: "置空版" sheet 输出时置空 ppm/Types
                rows.append(dup)
    df = pd.DataFrame(rows)
    df = df.dropna(subset=['ppm']).dropna(subset=['Types']).reset_index(drop=True)
    return df


def load_json(path: str) -> dict:
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def pick_docs(json_dir: str, dim: int, name_filter='', exclude=''):
    """按 experimentType 选 JSON (同类型取峰数最多)."""
    best = None
    for fn in sorted(os.listdir(json_dir)):
        if not fn.endswith('.json') or fn.startswith('_'):
            continue
        doc = load_json(os.path.join(json_dir, fn))
        if doc.get('dim') != dim:
            continue
        et = (doc.get('experimentType') or '').upper()
        if name_filter and name_filter not in et:
            continue
        if exclude and exclude in et:
            continue
        if best is None or len(doc.get('peaks', [])) > len(best[1].get('peaks', [])):
            best = (fn, doc)
    return best[1] if best else None


def hsqc_pairs(hsqc_doc: dict):
    """HSQC 峰 → [(c_ppm, h_ppm, intensity)] 强度降序, 剔除弱噪声."""
    out = []
    for pk in hsqc_doc.get('peaks', []):
        c, h = pk.get('ppm1'), pk.get('ppm2')
        if c is None or h is None:
            continue
        if not (-5 <= c <= 240 and -5 <= h <= 15):
            continue
        inten = float(pk.get('intensity', 0.0) or 0.0)
        if abs(inten) < WEAK_INT:
            continue  # GSD 弱峰/噪声 (真实交叉峰强度一般 > 2)
        out.append((float(c), float(h), inten))
    out.sort(key=lambda x: -x[2])
    return out


# ---- 第二阶段: 矩阵复核 (模拟化学家看 HSQC 等高线图) ----
# 参考: 氢谱积分 1D 信号的 nH 是硬约束; 若某 δH 附近只有一团等高线但 nH>=2,
# 说明是斜向融合椭圆 → 在矩阵上做单/双 2D 高斯拟合拆分出多个 (δC, δH) 中心;
# 逐行验证初步归属: 每行 (C, H) 必须落在某个交叉峰中心附近.

# 复核容差 (ppm)
REV_C_TOL = 0.75      # 交叉峰中心与归属 δC 的最大偏差
REV_H_TOL = 0.04      # 交叉峰中心与归属 δH 的最大偏差
REV_CLUSTER_SIG = 5.0  # 峰检测阈值: 噪声 MAD 的倍数
REV_SPLIT_RESID = 0.35  # 双高斯判定: 单高斯残差 RMS / 峰高 > 此值 → 尝试拆分


def load_hsqc_matrix(json_dir: str):
    """读 HSQC JSON 的 2D 矩阵. 返回 (matrix, c_axis, h_axis) 或 None.

    矩阵布局: 行 = F1 碳维 (firstXppm→lastXppm), 列 = F2 氢维 (firstYppm→lastYppm).
    """
    doc = pick_docs(json_dir, 2, 'HSQC')
    if doc is None:
        return None
    m = doc.get('data2d')
    if not m or not m[0]:
        return None
    if doc.get('firstXppm') is None or doc.get('firstYppm') is None:
        return None
    try:
        import numpy as np
        mat = np.array(m, dtype=float)
        c_ax = np.linspace(float(doc['firstXppm']), float(doc['lastXppm']), mat.shape[0])
        h_ax = np.linspace(float(doc['firstYppm']), float(doc['lastYppm']), mat.shape[1])
        return mat, c_ax, h_ax
    except (ValueError, TypeError) as e:
        print(f'[复核][WARN] HSQC 矩阵加载失败: {e}')
        return None


def _cluster_blobs(mat, c_ax, h_ax):
    """连通域分割: 阈值以上像素聚成峰团. 返回 [{rows, cols, peak}] 列表 (按峰高降序)."""
    import numpy as np
    from scipy import ndimage
    noise = float(np.median(np.abs(mat - np.median(mat)))) * 1.4826  # MAD→σ
    thr = noise * REV_CLUSTER_SIG
    mask = mat > max(thr, mat.max() * 0.05)
    lab, n = ndimage.label(mask)
    blobs = []
    for k in range(1, n + 1):
        ys, xs = np.where(lab == k)
        pk = mat[ys, xs].max()
        if len(ys) < 4 or pk < thr * 2:
            continue
        blobs.append({'rows': ys, 'cols': xs, 'peak': float(pk)})
    blobs.sort(key=lambda b: -b['peak'])
    return blobs


def _fit_gaussians(mat, c_ax, h_ax, blob, n_gauss):
    """在峰团窗口内拟合 n 个 2D 高斯 (共享 σ). 返回 (centers, resid_ratio) 或 None."""
    import numpy as np
    from scipy.optimize import curve_fit
    ys, xs = blob['rows'], blob['cols']
    pad_c = max(3, int(0.06 / max(abs(c_ax[1] - c_ax[0]), 1e-9)))
    pad_h = max(3, int(0.03 / max(abs(h_ax[1] - h_ax[0]), 1e-9)))
    r0, r1 = max(ys.min() - pad_c, 0), min(ys.max() + pad_c + 1, mat.shape[0])
    c0, c1 = max(xs.min() - pad_h, 0), min(xs.max() + pad_h + 1, mat.shape[1])
    sub = np.clip(mat[r0:r1, c0:c1], 0, None)
    if sub.size < n_gauss * 12 or sub.shape[0] == 0 or sub.shape[1] == 0:
        return None
    C, H = np.meshgrid(c_ax[r0:r1], h_ax[c0:c1], indexing='ij')
    pk = float(sub.max())
    dc = abs(c_ax[1] - c_ax[0])
    dh = abs(h_ax[1] - h_ax[0])
    # 初始中心: 极大点; 多高斯时第 2 个放在峰团对角另一端
    imax = np.unravel_index(sub.argmax(), sub.shape)
    centers0 = [(c_ax[r0:r1][imax[0]], h_ax[c0:c1][imax[1]])]
    if n_gauss >= 2:
        # 第 2 中心: 远离第 1 中心的次极大
        tmp = sub.copy()
        r_off = max(1, int(0.5 / dc))
        c_off = max(1, int(0.02 / dh))
        tmp[max(imax[0]-r_off,0):imax[0]+r_off+1, max(imax[1]-c_off,0):imax[1]+c_off+1] = 0
        i2 = np.unravel_index(tmp.argmax(), tmp.shape)
        centers0.append((c_ax[r0:r1][i2[0]], h_ax[c0:c1][i2[1]]))

    def model(XH, *args):
        C_, H_ = XH
        out = np.full(C_.shape, args[-1])
        for g in range(n_gauss):
            a, c, h = args[g * 3:(g + 1) * 3]
            out = out + a * np.exp(-((C_ - c) ** 2 / (2 * args[-3] ** 2) +
                                     (H_ - h) ** 2 / (2 * args[-2] ** 2)))
        return out.ravel()

    p0, lo, hi_ = [], [], []
    dc0 = min(1.5, max(0.2, (c_ax[r0:r1][-1] - c_ax[r0:r1][0]) / (4 * n_gauss)))
    dh0 = min(0.05, max(0.004, (h_ax[c0:c1][-1] - h_ax[c0:c1][0]) / (4 * n_gauss)))
    # 轴可能降序 (F1 从高 ppm 到低 ppm), 边界必须按 min/max 排序
    c_lo, c_hi2 = min(c_ax[r0], c_ax[r1 - 1]), max(c_ax[r0], c_ax[r1 - 1])
    h_lo, h_hi2 = min(h_ax[c0], h_ax[c1 - 1]), max(h_ax[c0], h_ax[c1 - 1])
    for (c0g, h0g) in centers0:
        p0 += [pk, c0g, h0g]
        lo += [pk * 0.1, c_lo, h_lo]
        hi_ += [pk * 3, c_hi2, h_hi2]
    p0 += [dc0, dh0, 0.0]
    lo += [dc0 * 0.2, dh0 * 0.2, -pk * 0.2]
    hi_ += [abs(c_hi2 - c_lo) * 2, abs(h_hi2 - h_lo) * 2, pk * 0.2]
    try:
        p, _ = curve_fit(model, (C, H), sub.ravel(), p0=p0,
                         bounds=(lo, hi_), maxfev=30000)
    except (RuntimeError, ValueError):
        return None
    # 残差只在峰团主体上算 (背景 0 会稀释 RMS, 导致融合峰漏检)
    fit = model((C, H), *p).reshape(sub.shape)
    body = sub > pk * 0.15
    resid = float(np.std((sub - fit)[body])) if body.any() else float(np.std(sub.ravel()))
    centers = [(float(p[g * 3 + 1]), float(p[g * 3 + 2]), float(p[g * 3]))
               for g in range(n_gauss)]
    return centers, resid / pk if pk > 0 else 9.9


def split_blobs(mat, c_ax, h_ax, want_splits: dict):
    """峰团拆分: 对峰团做单高斯; 残差大且该区域期望 H 数>=2 时做双高斯.

    want_splits: {(round(c,1), round(h,1)): n_h} — 1D 信号 nH 约束 (提示拆分).
    返回 centers: [(c_ppm, h_ppm, amp)], 及 note 列表.
    """
    import numpy as np
    blobs = _cluster_blobs(mat, c_ax, h_ax)
    centers, notes = [], []
    for b in blobs:
        ys, xs = b['rows'], b['cols']
        cc = float(np.median(c_ax[ys]))
        hh = float(np.median(h_ax[xs]))
        want = 0
        for (kc, kh), nh in want_splits.items():
            if abs(kc - cc) <= 0.35 and abs(kh - hh) <= 0.06:
                want = max(want, nh)
        r1 = _fit_gaussians(mat, c_ax, h_ax, b, 1)
        if r1 is None:
            # 无拟合回退: 极大点位置
            i = np.argmax(mat[b['rows'], b['cols']])
            centers.append((float(c_ax[b['rows'][i]]), float(h_ax[b['cols'][i]]), b['peak']))
            continue
        (cen1, rr1) = r1
        do_split = want >= 2 and rr1 > REV_SPLIT_RESID
        if not do_split:
            # 即使无 1D 提示, 残差很大也尝试双高斯 (交叉验证)
            do_split = rr1 > REV_SPLIT_RESID * 2.2
        if do_split:
            r2 = _fit_gaussians(mat, c_ax, h_ax, b, 2)
            if r2 is not None:
                (cen2, rr2) = r2
                if rr2 < rr1 * 0.8 and abs(cen2[0][0] - cen2[1][0]) > 0.3 \
                        and abs(cen2[0][1] - cen2[1][1]) > 0.008:
                    centers += cen2
                    notes.append(f'融合峰 ({cc:.2f},{hh:.2f}) 拆分为 '
                                 f'{cen2[0][0]:.2f}/{cen2[0][1]:.2f} + {cen2[1][0]:.2f}/{cen2[1][1]:.2f}'
                                 f' (残差 {rr1:.2f}→{rr2:.2f})')
                    continue
        centers.append((cen1[0][0], cen1[0][1], cen1[0][2]))
    return centers, notes


def fmt_1h_nmr(h_str, sig):
    """拼 xlsx 1H-NMR 格式: '4.31, t (8.8)'; 无 1D 信号时仅 δ.
    3 个及以上 J 值的复杂裂分 (如 ddt) 按惯例简写为 m."""
    if not sig:
        return h_str
    mult = sig.get('mult') or 'm'
    js = sig.get('J') or []
    if len(js) >= 3:
        return f"{h_str}, m"
    if js:
        return f"{h_str}, {mult} ({', '.join(js)})"
    return f"{h_str}, {mult}"


def review_assignment(out: pd.DataFrame, truth: pd.DataFrame,
                      h1d_signals, json_dir: str, sig_df=None) -> pd.DataFrame:
    """第二阶段矩阵复核: 逐行验证初步归属, 融合峰拆分修正.

    返回带 'review' 列的 out (复核结论: OK / 偏移 / 复核修正 / 未命中).
    回填仅使用拆分中心 (融合峰拆出的第 2 中心), 且要求本行 δC 与中心 δC 匹配、
    中心 δH 对齐到某 1D 信号 (有真实氢谱峰支撑), 避免激进改写.
    """
    loaded = load_hsqc_matrix(json_dir)
    if loaded is None:
        out['review'] = '无HSQC矩阵'
        return out
    mat, c_ax, h_ax = loaded
    # 融合峰拆分 (1D nH 提示 + 残差/双高斯验证), 输出交叉峰中心集合
    centers, split_flags, notes = split_blobs_with_hints(mat, c_ax, h_ax, h1d_signals)
    for nmsg in notes:
        print(f'[复核] {nmsg}')
    # 逐行核对: (δC, δH) 是否落在某交叉峰中心 ±(REV_C_TOL, REV_H_TOL) 内
    center_used = [False] * len(centers)
    rows_idx = []
    for i, (_, r) in enumerate(out.iterrows()):
        if pd.isna(r['assigned_1H']) or not str(r['assigned_1H']).strip():
            continue
        rows_idx.append((float(r['ppm']), float(str(r['assigned_1H'])), i))
    rev = [''] * len(out)
    for c_v, h_v, i in rows_idx:
        hit, hit_k, best_d = None, None, 9e9
        for k, (cc, hh, amp) in enumerate(centers):
            d = abs(cc - c_v) / REV_C_TOL + abs(hh - h_v) / REV_H_TOL
            if abs(cc - c_v) <= REV_C_TOL and abs(hh - h_v) <= REV_H_TOL and d < best_d:
                hit, hit_k, best_d = (cc, hh), k, d
        if hit is None:
            # 回填 (仅拆分中心): δC 匹配 + 中心 δH 对齐到某 1D 信号 (有氢谱峰支撑)
            cand = []
            for k, (cc, hh, amp) in enumerate(centers):
                if not (split_flags[k] and not center_used[k]
                        and abs(cc - c_v) <= 0.35):
                    continue
                s1d = min(h1d_signals, key=lambda x: abs(x['ppm'] - hh),
                          default=None)
                if s1d and abs(s1d['ppm'] - hh) <= H_ALIGN_TOL:
                    cand.append((abs(cc - c_v), abs(hh - h_v), k))
            if cand:
                cand.sort()
                _, _, k = cand[0]
                cc, hh, _ = centers[k]
                center_used[k] = True
                rev[i] = f'复核修正→({cc:.2f},{hh:.2f})'
                h_aligned = f'{hh:.2f}'
                s1d = min(h1d_signals, key=lambda x: abs(x['ppm'] - hh),
                          default=None)
                if s1d and abs(s1d['ppm'] - hh) <= H_ALIGN_TOL:
                    h_aligned = f"{s1d['ppm']:.2f}"
                out.at[i, 'assigned_1H'] = h_aligned
                # 第二步修正后的格式化 1H-NMR (δ, mult (J)) 及 nH / δH_ppm
                out.at[i, '1H-NMR-v2'] = fmt_1h_nmr(h_aligned, s1d)
                if s1d and s1d.get('nH') is not None:
                    out.at[i, 'nH'] = int(s1d['nH'])
                if s1d is not None and 'δH_ppm' in out.columns:
                    m_sig = sig_df[pd.to_numeric(sig_df['δH'], errors='coerce')
                                   .sub(float(h_aligned)).abs().le(0.011)] \
                        if isinstance(sig_df, pd.DataFrame) and 'δH' in sig_df.columns \
                        else pd.DataFrame()
                    if not m_sig.empty:
                        out.at[i, 'δH_ppm'] = str(m_sig.iloc[0]['δH_ppm'])
                continue
            rev[i] = f'未命中交叉峰 ({c_v:.2f})'
            continue
        center_used[hit_k] = True
        if abs(hit[0] - c_v) <= 0.15 and abs(hit[1] - h_v) <= 0.02:
            rev[i] = f"OK ({hit[0]:.2f},{hit[1]:.2f})"
        else:
            rev[i] = f"偏移→({hit[0]:.2f},{hit[1]:.2f})"
    # review 为空的行 (Cq/未归属): 从 mnova 13C 谱峰表找本行 δC 附近碳峰补入 review
    # Cq 独立碳峰命中 → "OK (C)"; 无 H 对应且 13C 也找不到 → "未命中交叉峰 (C)"
    c13_peaks = _load_c13_peaks(json_dir)
    for i, (_, r) in enumerate(out.iterrows()):
        if rev[i] or pd.isna(r['ppm']):
            continue
        c_v = float(r['ppm'])
        near = [(abs(p - c_v), p) for p in c13_peaks if abs(p - c_v) <= 0.5]
        if near:
            near.sort()
            rev[i] = f'OK ({near[0][1]:.2f})'
        else:
            rev[i] = f'未命中交叉峰 ({c_v:.2f})'
    out['review'] = rev
    return out


def _load_c13_peaks(json_dir: str):
    """读 13C 1D 谱 JSON 峰表 (GSD), 返回 δC 列表 (强度降序去重)."""
    peaks = []
    for fn in sorted(os.listdir(json_dir)):
        if not fn.endswith('.json') or fn.startswith('_'):
            continue
        try:
            with open(os.path.join(json_dir, fn), encoding='utf-8') as f:
                doc = json.load(f)
        except (OSError, ValueError):
            continue
        if doc.get('dim') != 1 or doc.get('nucleus') != '13C':
            continue
        for pk in doc.get('peaks', []):
            c = pk.get('ppm1')
            if c is not None and -5 <= float(c) <= 240:
                peaks.append(float(c))
    # 去重 (0.01 内合并), 按值排序
    peaks.sort()
    out = []
    for p in peaks:
        if not out or abs(p - out[-1]) > 0.01:
            out.append(p)
    return out


def split_blobs_with_hints(mat, c_ax, h_ax, h1d_signals):
    """带 1D nH 提示的峰团拆分 (模拟化学家看融合椭圆).

    两遍处理: 第一遍所有团做单高斯并收集 σ_C 全局中位数;
    第二遍对 (1D nH>=2 提示) 且 (残差偏高 或 σ_C 显著宽于中位) 的团做双高斯拆分.
    """
    import numpy as np
    blobs = _cluster_blobs(mat, c_ax, h_ax)
    infos = []
    for b in blobs:
        ys, xs = b['rows'], b['cols']
        # 轴可能降序, 区间端点取 min/max
        h_lo, h_hi_ = min(h_ax[xs.min()], h_ax[xs.max()]), max(h_ax[xs.min()], h_ax[xs.max()])
        want = 0
        for s in h1d_signals:
            if s.get('nH', 0) >= 2 and h_lo - 0.01 <= s['ppm'] <= h_hi_ + 0.01:
                want = max(want, s['nH'])
        r1 = _fit_gaussians(mat, c_ax, h_ax, b, 1)
        infos.append({'blob': b, 'want': want, 'r1': r1,
                      'cc': float(np.median(c_ax[ys])), 'hh': float(np.median(h_ax[xs]))})
    # 全局残差中位数 (拟合正常团) — 残差显著高于中位的融合团更可疑
    sig_c = [inf['r1'][1] for inf in infos if inf['r1'] is not None]
    med_resid = float(np.median(sig_c)) if sig_c else 0.1
    centers, split_flags, notes = [], [], []
    for inf in infos:
        b, want, r1 = inf['blob'], inf['want'], inf['r1']
        if r1 is None:
            ys, xs = b['rows'], b['cols']
            i = int(np.argmax(mat[ys, xs]))
            centers.append((float(c_ax[ys[i]]), float(h_ax[xs[i]]), b['peak']))
            split_flags.append(False)
            continue
        cen1, rr1 = r1
        # 拆分条件: 1D nH>=2 提示 (该 δH 区域确实有多个 H) → 总是尝试双高斯,
        # 接受标准: 残差显著下降 且 两中心间距合理 (ΔC 0.3~3.5, ΔH > 0.008)
        do_split = want >= 2
        if not do_split:
            do_split = rr1 > REV_SPLIT_RESID * 2.2  # 残差极大: 无提示也拆 (交叉验证)
        if do_split:
            r2 = _fit_gaussians(mat, c_ax, h_ax, b, 2)
            if r2 is not None:
                cen2, rr2 = r2
                dc2 = abs(cen2[0][0] - cen2[1][0])
                if rr2 < rr1 * 0.85 and 0.3 <= dc2 <= 3.5 \
                        and abs(cen2[0][1] - cen2[1][1]) > 0.008:
                    centers += cen2
                    split_flags += [True, True]
                    notes.append(f'融合峰 ({inf["cc"]:.2f},{inf["hh"]:.2f}) 拆分为 '
                                 f'({cen2[0][0]:.2f},{cen2[0][1]:.2f}) + ({cen2[1][0]:.2f},{cen2[1][1]:.2f})'
                                 f' [ΔC={dc2:.2f}, 残差 {rr1:.3f}→{rr2:.3f}, 1D nH={want}]')
                    continue
        centers.append((cen1[0][0], cen1[0][1], cen1[0][2]))
        split_flags.append(False)
    return centers, split_flags, notes


# ---- 核心归属 ----

def assign(xlsx_path: str, h_report_path: str, json_dir: str, out_dir: str,
           do_step3: bool = False) -> pd.DataFrame:
    df = load_xlsx(xlsx_path)
    # 化合物名: xlsx 文件名去 "数据归属.xlsx" (与归属表输出同名前缀)
    base = re.sub(r'数据归属\.xlsx$', '', os.path.basename(xlsx_path)) or 'compound'
    # 两步归属共用积分表: txt → <base>_1H_NMR.csv (已存在则复用) → 内部信号 list
    # (δH 列已是数值: 区间存中点)
    sig_df = export_1h_nmr(h_report_path, out_dir, compound=base)
    h1d_signals = [{'ppm': float(r['δH']), 'mult': r['type'], 'nH': float(r['nH']),
                    'J': r['J']} for _, r in sig_df.iterrows()]
    hsqc = pick_docs(json_dir, 2, 'HSQC')
    if hsqc is None:
        raise SystemExit('ERROR: 未找到 HSQC JSON')
    pairs = hsqc_pairs(hsqc)

    # ---- 碳槽: 同 (ppm 精确到 0.01, Types) 合并; n_h = 槽内 xlsx 行数 ----
    # 化学家 xlsx 每行 = 一个可观测 H 信号: CH 两行(如 72.32)表示两个可分辨 δH;
    # CH3 一行 (3H 等价) → 1 个 δH; 同 ppm 的 CH2 两行 = 非对映 2H.
    # 注意: 不做 0.06 ppm 模糊合并 (79.12/79.09 是两个真实碳, 各 1H)
    slots = []
    for idx, row in df.iterrows():
        c = round(float(row['ppm']), 2)
        typ = str(row['Types']).strip()
        slot = next((s for s in slots if abs(s['ppm'] - c) < 0.005 and s['type'] == typ), None)
        if slot is None:
            slot = {'ppm': c, 'type': typ, 'n_h': 0, 'rows': [], 'hs': []}
            slots.append(slot)
        slot['rows'].append(idx)
        slot['n_h'] += 1
        # Cq 恒 0
        if typ == 'C':
            slot['n_h'] = 0

    # ---- 全局一对一匹配: 密集簇匈牙利最优 + 其余贪心 ----
    # 先识别密集簇: 2.5 ppm 窗口内 >=3 个槽且 >=3 个 HSQC 峰 → 簇内匈牙利全局最优
    # (代价 = |ΔC| + 0.5 * |ΔH_align|, ΔH 利用 1D 报告信号消歧), 簇外按贪心.
    def _align_h(h):
        if not h1d_signals:
            return h
        s = min(h1d_signals, key=lambda x: abs(x['ppm'] - h))
        return s['ppm'] if abs(s['ppm'] - h) <= H_ALIGN_TOL else h

    active_slots = [s for s in slots if s['n_h'] > 0]
    assigned = [False] * len(active_slots)  # 匈牙利已处理的槽
    used_pairs: set = set()

    try:
        from scipy.optimize import linear_sum_assignment as _lsa
        import numpy as _np
        used_cluster = set()
        for si, slot in enumerate(active_slots):
            if si in used_cluster:
                continue
            # 簇: si 与其后 2.5 ppm 内的所有槽
            members = [si]
            for sj in range(si + 1, len(active_slots)):
                if abs(active_slots[sj]['ppm'] - slot['ppm']) <= 2.5:
                    members.append(sj)
            # 簇内槽 >=3 且可用峰 >=3 才匈牙利
            cpeaks = [(pi, p) for pi, p in enumerate(pairs)
                      if any(abs(p[0] - active_slots[m]['ppm']) <= HSQC_C_TOL
                             for m in members)]
            if len(members) >= 3 and len(cpeaks) >= 3:
                for m in members:
                    used_cluster.add(m)
                # 代价矩阵: 行=槽(展开 H 位), 列=峰
                rows = []
                slot_of_row = []
                for m in members:
                    for _k in range(active_slots[m]['n_h']):
                        rows.append(m)
                        slot_of_row.append(m)
                n_rows, n_cols = len(rows), len(cpeaks)
                cost = _np.full((n_rows, n_cols), 9.9)
                for ri, m in enumerate(rows):
                    sl = active_slots[m]
                    for ci, (pi, p) in enumerate(cpeaks):
                        dc = abs(p[0] - sl['ppm'])
                        if dc > HSQC_C_TOL:
                            continue
                        dh = abs(_align_h(p[1]) - p[1])
                        # CH3 加成
                        adj = 0.0
                        if sl['type'] == 'CH3':
                            s_near = min(h1d_signals, key=lambda x: abs(x['ppm'] - p[1]),
                                         default=None)
                            if s_near and s_near['nH'] >= 3 \
                                    and abs(s_near['ppm'] - p[1]) <= H_ALIGN_TOL:
                                adj = -CH3_BONUS
                        cost[ri, ci] = dc + 0.5 * dh + adj
                ri_idx, ci_idx = _lsa(cost)
                for ri, ci in zip(ri_idx, ci_idx):
                    if cost[ri, ci] >= 9.9:
                        continue
                    sl = active_slots[slot_of_row[ri]]
                    h = cpeaks[ci][1][1]
                    if any(abs(h - ph) < DUP_TOL for ph in sl['hs']):
                        continue
                    if len(sl['hs']) < sl['n_h']:
                        sl['hs'].append(h)
                        used_pairs.add(cpeaks[ci][0])
    except ImportError:
        pass

    # 贪心兜底 (非簇槽 / 匈牙利未覆盖的 H 位)
    edges = []
    for pi, (c, h, inten) in enumerate(pairs):
        # 该 H 在 1D 报告中的信号 (用于 CH3 加成)
        h1d_near = min(h1d_signals, key=lambda s: abs(s['ppm'] - h)) if h1d_signals else None
        for si, slot in enumerate(slots):
            if slot['n_h'] == 0:
                continue
            d = abs(c - slot['ppm'])
            if d > HSQC_C_TOL:
                continue
            eff = d
            if slot['type'] == 'CH3' and h1d_near and h1d_near['nH'] >= 3 \
                    and abs(h1d_near['ppm'] - h) <= H_ALIGN_TOL:
                eff -= CH3_BONUS
            edges.append((eff, -inten, pi, si, c, h))
    edges.sort()
    for eff, negint, pi, si, c, h in edges:
        slot = slots[si]
        if pi in used_pairs or len(slot['hs']) >= slot['n_h']:
            continue
        if any(abs(h - ph) < DUP_TOL for ph in slot['hs']):
            continue
        slot['hs'].append(h)
        used_pairs.add(pi)

    # ---- 第二遍: 复用扫描 (重叠 H) ----
    # 化学家报告的多重峰 (如 '1.48-1.62 m 6H') 是多碳 H 融合, 同一 HSQC 交叉峰
    # 可同时属于相邻两个碳. 当某槽仍缺 H 且没有未用峰可配时, 允许复用
    # 已被邻近槽 (ΔC<=1.0) 占用的峰, 并在 note 中标记 '重叠复用'.
    for pass_i in range(2):  # 最多两遍复用 (CH2 可能缺 2 个)
        remaining = [s for s in slots if len(s['hs']) < s['n_h']]
        if not remaining:
            break
        progressed = False
        for slot in remaining:
            # 找可复用峰: ΔC 最近的槽已占用, 且该槽与本槽 ΔC<=1.0
            cands = []
            for pi, (c, h, inten) in enumerate(pairs):
                if abs(c - slot['ppm']) > HSQC_C_TOL:
                    continue
                if any(abs(h - ph) < DUP_TOL for ph in slot['hs']):
                    continue
                if pi in used_pairs:
                    # 找占用者: 其他槽含同 h 且碳距 <=1.2 (重叠/融合 H 场景,
                    # 融合多重峰如 '3.95-4.05 m 4H' 可横跨多个碳)
                    owner = None
                    for s2 in slots:
                        if s2 is slot:
                            continue
                        if any(abs(h - ph) < DUP_TOL for ph in s2['hs']) \
                                and abs(c - s2['ppm']) <= HSQC_C_TOL \
                                and abs(s2['ppm'] - slot['ppm']) <= 1.2:
                            owner = s2
                            break
                    if owner is None:
                        continue
                    cands.append((abs(c - slot['ppm']), -inten, pi, h))
                else:
                    cands.append((abs(c - slot['ppm']), -inten, pi, h))
            if cands:
                cands.sort()
                _, _, pi, h = cands[0]
                slot['hs'].append(h)
                used_pairs.add(pi)
                slot.setdefault('reused', []).append(h)
                progressed = True
        if not progressed:
            break

    # ---- 第三遍: 融合信号借用 (谱图分辨率极限场景) ----
    # 若某槽仍缺 H: 找一个 1D 融合信号 (nH >= 全局已认领该信号的 H 数 + 本槽缺数),
    # 其 δH 附近的碳槽都已配满 → 该槽的 H 必然淹没在此融合信号中, 直接借用其 δH.
    # 标注 '融合借用', 供人工复核.
    used_h_count: dict = {}
    for slot in slots:
        for h in slot['hs']:
            s = min(h1d_signals, key=lambda x: abs(x['ppm'] - h), default=None)
            if s and abs(s['ppm'] - h) <= H_ALIGN_TOL:
                used_h_count[id(s)] = used_h_count.get(id(s), 0) + 1
    for slot in slots:
        while len(slot['hs']) < slot['n_h']:
            cands = []
            # 借用的 1D 信号 δH 必须接近本槽已配 H (同一融合峰区域) 或本槽 δC 附近的典型 1H 区
            anchor_hs = slot['hs'] or []
            for s in h1d_signals:
                n_used = used_h_count.get(id(s), 0)
                # 融合信号: nH 明显多于已认领 → 有未解释的 H 淹没其中
                if s['nH'] - n_used >= slot['n_h'] - len(slot['hs']) \
                        and s['nH'] >= 2:
                    # δH 合理性: 接近本槽已配的某个 H (融合区一致) 距离 <=0.5
                    if anchor_hs and not any(abs(s['ppm'] - ah) <= 0.5 for ah in anchor_hs):
                        continue
                    # 该信号附近 1.0 ppm 内的槽都应已配满 (否则应走复用而非借用)
                    near_slots = [s2 for s2 in slots
                                  if abs(s2['ppm'] - slot['ppm']) <= 1.0 and s2 is not slot]
                    all_full = all(len(s2['hs']) >= s2['n_h'] for s2 in near_slots)
                    if all_full:
                        # 借用距离: 优先离本槽已配 H 近的融合信号; 无锚时离本槽 δC
                        # 常识范围 (脂肪 0.5-2.5 / 糖 3-5.5 / 芳 6-8.5) 内最近的
                        if anchor_hs:
                            dkey = min(abs(s['ppm'] - ah) for ah in anchor_hs)
                        else:
                            dkey = abs(s['ppm'] - _typical_h(slot['ppm'], slot['type']))
                        cands.append((dkey, id(s), s))
            if not cands:
                break
            s = cands[0][2]
            slot['hs'].append(s['ppm'])
            used_h_count[id(s)] = used_h_count.get(id(s), 0) + 1
            slot.setdefault('borrowed', []).append(s['ppm'])

    # ---- CH2 补第 2 个 H: 1D 信号 nH>=2 且未被占用时, 找同碳最近未用 HSQC 峰 ----
    used_h1d = set()
    for slot in slots:
        for h in slot['hs']:
            s = min(h1d_signals, key=lambda x: abs(x['ppm'] - h), default=None)
            if s and abs(s['ppm'] - h) <= H_ALIGN_TOL:
                used_h1d.add(id(s))
    for slot in slots:
        if slot['type'] != 'CH2' or len(slot['hs']) >= 2:
            continue
        # 已对齐的 1D 信号若 nH>=2, 其覆盖的其他 δH 可作为第 2 个 H 候选 (多重峰融合)
        cands = []
        for pi, (c, h, inten) in enumerate(pairs):
            if pi in used_pairs:
                continue
            if abs(c - slot['ppm']) > HSQC_C_TOL:
                continue
            s = min(h1d_signals, key=lambda x: abs(x['ppm'] - h), default=None)
            if s and abs(s['ppm'] - h) <= H_ALIGN_TOL and s['nH'] >= 2 \
                    and id(s) not in used_h1d:
                cands.append((abs(c - slot['ppm']), -inten, pi, h))
        if cands:
            cands.sort()
            _, _, pi, h = cands[0]
            if not any(abs(h - ph) < DUP_TOL for ph in slot['hs']):
                slot['hs'].append(h)
                used_pairs.add(pi)

    # ---- 输出 ----
    def align(h):
        """HSQC δH → (对齐文本, 对齐的 1D 信号 dict|None)."""
        if not h1d_signals:
            return f'{h:.2f}', None
        s = min(h1d_signals, key=lambda x: abs(x['ppm'] - h))
        if abs(s['ppm'] - h) <= H_ALIGN_TOL:
            return f"{s['ppm']:.2f}", s
        return f'{h:.2f}', None

    results = []
    for slot in slots:
        hs_sorted = sorted(slot['hs'], reverse=True)
        h_strs, flags, sigs = [], [], []
        for h in hs_sorted:
            s, sig = align(h)
            h_strs.append(s)
            sigs.append(sig)
            flags.append(sig is not None)
        # CH3: 3H 等价理论 1 行, 但 HSQC 若分离出多个 δH 全部报告 (供人工判断)
        row_outputs = [h_strs[k] if k < len(h_strs) else '' for k in range(len(slot['rows']))]
        for k, ridx in enumerate(slot['rows']):
            row = df.loc[ridx]
            # 1H-NMR 列: xlsx 报告格式 (δ, mult (J)); 每行取对应的那个 H
            if k < len(h_strs) and h_strs[k]:
                nmr_str = fmt_1h_nmr(h_strs[k], sigs[k])
                # nH: 对应 1D 信号积分括号里的 H 数; δH_ppm: 该信号 txt 原文
                nh_val = int(sigs[k]['nH']) if sigs[k] and sigs[k].get('nH') is not None else ''
                dh_ppm = ''
                if sigs[k] is not None:
                    m_sig = sig_df[pd.to_numeric(sig_df['δH'], errors='coerce')
                                   .sub(float(h_strs[k])).abs().le(0.011)] \
                        if 'δH' in sig_df.columns else pd.DataFrame()
                    if not m_sig.empty:
                        dh_ppm = str(m_sig.iloc[0]['δH_ppm'])
            else:
                nmr_str = ''
                nh_val = ''
                dh_ppm = ''
            # 原 xlsx 其它列原样保留 (除程序生成的列)
            other_cols = {c: row[c] for c in df.columns
                          if c not in ('ppm', 'Types', '1H-NMR', '_cont')}
            results.append({
                'ppm': float(row['ppm']), 'Types': slot['type'],
                **other_cols,  # 输入 xlsx 的其它列 (原样)
                '1H-NMR-v1': nmr_str,  # 第一步初匹配结果 (复核修正前)
                'nH': nh_val,  # 归属 1D 信号的积分 H 数
                'δH_ppm': dh_ppm,  # 归属 1D 信号在积分表中的 δH 原文 (区间保留)
                'assigned_1H': row_outputs[k],
                'all_assigned_H': '; '.join(h_strs) if k == 0 else '',
                'aligned_1D': ('; '.join(f'{h_strs[i]}{"✓" if flags[i] else "✗"}'
                                         for i in range(len(h_strs)))) if k == 0 else '',
                '_cont': bool(row.get('_cont', False)),
            })
    out = pd.DataFrame(results)
    # 第二阶段: HSQC 矩阵复核 (逐行验证交叉峰 + 融合峰拆分)
    # review 内部对修正行生成 1H-NMR-v2 及更新 δH_ppm/nH; 未修正行 v2 = v1
    out = review_assignment(out, df, h1d_signals, json_dir, sig_df=sig_df)
    if '1H-NMR-v2' not in out.columns:
        out['1H-NMR-v2'] = out['1H-NMR-v1']
    else:
        m_missing = out['1H-NMR-v2'].isna() | (out['1H-NMR-v2'].astype(str).str.strip() == '')
        out.loc[m_missing, '1H-NMR-v2'] = out.loc[m_missing, '1H-NMR-v1']
    # 第三步 (COSY/HMBC 全表复核) 不改写 H — v3 = v2; 复核证据在第三 sheet
    # 开关 do_step3 默认关闭: 不执行第三步, 输出到 1H-NMR-v2 为止;
    # 原依赖 v3 的列 (row_match / 1H-NMR_real 对齐基准等) 改用 v2.
    if do_step3:
        out['1H-NMR-v3'] = out['1H-NMR-v2']
    _ver = '1H-NMR-v3' if do_step3 else '1H-NMR-v2'  # 后续判分/对齐基准列
    # 排序: ppm 从大到小; 同 ppm 按 assigned_1H 从大到小 (Cq/空值排最后)
    # ascending=False: ppm 与 δH 均降序
    out['_h_sort'] = out['assigned_1H'].map(
        lambda x: float(x) if not pd.isna(x) and str(x).strip() else -9e9)
    out = out.sort_values(['ppm', '_h_sort'], ascending=False).drop(columns='_h_sort').reset_index(drop=True)
    os.makedirs(out_dir, exist_ok=True)
    # 输出文件名带化合物名前缀 (base 已在函数开头从 xlsx 文件名推导)
    out_xlsx = os.path.join(out_dir, f'{base}_归属表.xlsx')

    # ---- 对比判分 (有 1H-NMR 已知列时) ----
    truth = load_xlsx(xlsx_path)

    def parse_all_h(s):
        if pd.isna(s) or not str(s).strip():
            return []
        vals = []
        for part in str(s).split(';'):
            try:
                vals.append(round(float(part.split(',')[0].strip().split(' ')[0]), 2))
            except ValueError:
                pass
        return vals

    truth['truth_all_H'] = truth['1H-NMR'].map(parse_all_h)

    # 集合级判分: 按 (ppm, Types) 槽聚合 (truth 多行/程序多行 笛卡尔积会错位,
    # 槽级 "程序 H 集合 vs truth H 集合" 贪心二分匹配才公平)
    prog_slots: dict = {}
    for _, r in out.iterrows():
        key = (round(float(r['ppm']), 2), str(r['Types']).strip())
        v = None if pd.isna(r['assigned_1H']) or not str(r['assigned_1H']).strip() \
            else float(r['assigned_1H'])
        prog_slots.setdefault(key, []).append(v)
    truth_slots: dict = {}
    for _, r in truth.iterrows():
        key = (round(float(r['ppm']), 2), str(r['Types']).strip())
        truth_slots.setdefault(key, []).extend(
            r['truth_all_H'] if isinstance(r['truth_all_H'], list) else [])

    cmp_rows = []
    n_ok = n_bad = 0
    for key in sorted(set(prog_slots) | set(truth_slots), key=lambda k: -k[0]):
        tset = truth_slots.get(key, [])
        pset = [h for h in prog_slots.get(key, []) if h is not None]
        n_prog_rows = len(prog_slots.get(key, []))
        if not tset:
            judge = 'OK' if not pset else 'MULTI?'
        elif not pset:
            judge = 'MISS'
        else:
            used: set = set()
            matched = 0
            for p in pset:
                cand = [(abs(p - t), i) for i, t in enumerate(tset) if i not in used]
                if not cand:
                    continue
                d, i = min(cand)
                if d <= 0.06:
                    used.add(i)
                    matched += 1
            if matched == len(tset) and len(pset) <= len(tset):
                judge = 'OK'
            elif matched > 0:
                judge = f'PART({matched}/{len(tset)})'
            elif any(abs(p - t) <= 0.35 for p in pset for t in tset):
                judge = 'SWAP?'
            else:
                judge = 'WRONG'
        if judge.startswith('OK'):
            n_ok += 1
        else:
            n_bad += 1
        # 简并槽对标记: 同 Types 相邻槽 (Δppm<=0.15) 的 SWAP 物理上不可区分
        cmp_rows.append({'ppm': key[0], 'Types': key[1],
                         'assigned_H': '; '.join(str(h) for h in pset),
                         'truth_H': '; '.join(str(h) for h in tset),
                         'judge': judge})
    # 后处理: 简并槽 SWAP? → OK(DEGENERATE)
    cmp_df = pd.DataFrame(cmp_rows)
    for i, r in cmp_df.iterrows():
        if not r['judge'].startswith(('SWAP?', 'WRONG')):
            continue
        for j, r2 in cmp_df.iterrows():
            if i == j or r['Types'] != r2['Types']:
                continue
            if abs(r['ppm'] - r2['ppm']) <= 0.15 and r2['judge'].startswith(('SWAP?', 'WRONG')):
                # 互换验证: 交换后都能对上
                a = [float(x) for x in str(r['assigned_H']).split(';') if x]
                b = [float(x) for x in str(r2['assigned_H']).split(';') if x]
                ta = [float(x) for x in str(r['truth_H']).split(';') if x]
                tb = [float(x) for x in str(r2['truth_H']).split(';') if x]
                if a and b:
                    merged_set = a + b
                    tset2 = ta + tb
                    used = set()
                    hits = 0
                    for p in merged_set:
                        cand = [(abs(p - t), k) for k, t in enumerate(tset2) if k not in used]
                        if cand:
                            dd, k = min(cand)
                            if dd <= 0.06:
                                used.add(k)
                                hits += 1
                    if hits == len(tset2):
                        cmp_df.at[i, 'judge'] = 'OK(DEGEN)'
                        cmp_df.at[j, 'judge'] = 'OK(DEGEN)'
                    break

    # ---- 把槽级 judge 回写 到 out 的每一行 (同 ppm+Types 的行共用槽 judge) ----
    slot_judge = {(r['ppm'], r['Types']): r['judge'] for _, r in cmp_df.iterrows()}
    out['judge'] = out.apply(
        lambda r: slot_judge.get(
            (round(float(r['ppm']), 2), str(r['Types']).strip()), ''),
        axis=1)

    # ---- 行级比较 (compare_assignment 逻辑) 回写 row_match 列 ----
    def _h_key(s):
        if pd.isna(s) or not str(s).strip():
            return 9e9
        try:
            return -float(str(s).split(',')[0].strip().split(' ')[0])
        except ValueError:
            return 9e9

    truth_rows = truth[['ppm', 'Types', '1H-NMR']].dropna(subset=['ppm']).copy()
    truth_rows['_k'] = truth_rows['1H-NMR'].map(_h_key)
    truth_rows = truth_rows.sort_values(['ppm', '_k']).reset_index(drop=True)
    prog_rows = out[['ppm', 'Types', _ver]].copy()
    prog_rows['_k'] = prog_rows[_ver].map(_h_key)
    prog_rows = prog_rows.sort_values(['ppm', '_k']).reset_index(drop=True)

    def _norm(v):
        if pd.isna(v):
            return ''
        if isinstance(v, float):
            return f'{v:.2f}'
        return str(v).strip()

    def _first_num(v):
        """取字符串首个数值 (δH), 无数值返回 None."""
        m = re.search(r'-?\d+\.?\d*', str(v))
        return float(m.group()) if m else None

    row_match_map: dict = {}
    n = max(len(truth_rows), len(prog_rows))
    for i in range(n):
        rt = truth_rows.iloc[i] if i < len(truth_rows) else None
        rp = prog_rows.iloc[i] if i < len(prog_rows) else None
        ok_row = True
        for col in ('ppm', 'Types', _ver):
            vt = _norm(rt[col] if col != _ver else rt['1H-NMR']) \
                if rt is not None else '<缺失>'
            vp = _norm(rp[col]) if rp is not None else '<缺失>'
            if col == 'ppm':
                try:
                    if abs(float(vt) - float(vp)) <= 0.005:
                        vp = vt
                except ValueError:
                    pass
            if vt != vp:
                ok_row = False
        if rp is not None:
            # 定位回 out 原行: 按 ppm+Types+基准列 三元组; 差异时附上 δH 差值
            label = '一致'
            if not ok_row:
                label = '差异'
                nh_p, nh_t = _first_num(rp[_ver]), _first_num(rt['1H-NMR']) if rt is not None else None
                if nh_p is not None and nh_t is not None:
                    diff = nh_p - nh_t
                    if abs(diff) < 0.002:
                        # Δ≈0: 只输出 '一致' 二字
                        label = '一致'
                    elif abs(diff) < _ROW_MATCH_TOL:
                        # 容差内: 保留差异值但改判一致
                        label = f'一致 (ΔδH={diff:+.2f})'
                    else:
                        label += f' (ΔδH={diff:+.2f})'
            row_match_map[(round(float(rp['ppm']), 2), str(rp['Types']).strip(),
                           _norm(rp[_ver]))] = label
    out['row_match'] = out.apply(
        lambda r: row_match_map.get(
            (round(float(r['ppm']), 2), str(r['Types']).strip(), _norm(r[_ver])), ''),
        axis=1)

    # ---- C_match 列: ppm (xlsx 碳峰) vs review 中的碳数值, 规则仿 row_match ----
    # review 碳数值来源: "OK (C,H)" / "偏移→(C,H)" / "复核修正→(C,H)" / "13C峰 C"
    def _first_c(v):
        """取 review 字符串中碳数值; 跳过 '13C峰' 前缀里的 13."""
        s = re.sub(r'13C峰\s*', '', str(v))
        m = re.search(r'-?\d+\.?\d*', s)
        return float(m.group()) if m else None

    c_match = []
    for _, r in out.iterrows():
        if pd.isna(r['ppm']) or not str(r['ppm']).strip():
            c_match.append('')
            continue
        c_v = float(r['ppm'])
        cr = _first_c(r.get('review'))
        if cr is None:
            c_match.append('')
            continue
        diff = cr - c_v
        if abs(diff) < 0.002:
            # Δ≈0: 只输出 '一致' 二字
            c_match.append('一致')
        elif abs(diff) < _ROW_MATCH_TOL:
            c_match.append(f'一致 (ΔδC={diff:+.2f})')
        else:
            c_match.append(f'差异 (ΔδC={diff:+.2f})')
    out['C_match'] = c_match

    # 重新排序后落盘; ppm/δH 均降序
    out['_h_sort'] = out['assigned_1H'].map(
        lambda x: float(x) if not pd.isna(x) and str(x).strip() else -9e9)
    out = out.sort_values(['ppm', '_h_sort'], ascending=False).drop(columns='_h_sort').reset_index(drop=True)

    # 氢谱数值统一四舍五入保留 2 位小数
    out['ppm'] = out['ppm'].map(lambda x: round(float(x), 2))
    for col in ('1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v3', 'assigned_1H',
                'all_assigned_H', 'aligned_1D'):
        if col in out.columns:
            out[col] = out[col].map(_round_h_in_str)

    # xlsx 有 1H-NMR truth 时: 插入 1H-NMR_real (xlsx 按 out 相同规则排序后行对齐复制) 于基准列 (_ver) 后
    has_truth = truth['1H-NMR'].notna().any()
    if has_truth:
        # 与 out 落盘排序同一规则: ppm 降序, 同 ppm 按 1H-NMR 中的 δH 降序
        def _h_val(s):
            if pd.isna(s) or not str(s).strip():
                return -9e9
            try:
                return float(str(s).split(',')[0].strip().split(' ')[0])
            except ValueError:
                return -9e9
        truth_sorted = truth[['ppm', 'Types', '1H-NMR']].dropna(subset=['ppm']).copy()
        truth_sorted['_k'] = truth_sorted['1H-NMR'].map(_h_val)
        truth_sorted = truth_sorted.sort_values(['ppm', '_k'], ascending=False).reset_index(drop=True)
        if len(truth_sorted) == len(out):
            out['1H-NMR_real'] = truth_sorted['1H-NMR'].map(
                lambda v: '' if pd.isna(v) else str(v)).values
            col_order = list(out.columns)
            col_order.remove('1H-NMR_real')
            # 关第三步时插到 1H-NMR-v2 后; 开第三步时插到 1H-NMR-v3 后
            idx = col_order.index(_ver) + 1
            col_order.insert(idx, '1H-NMR_real')
            out = out[col_order]
    else:
        out = out.drop(columns=[c for c in ('1H-NMR_real', 'row_match') if c in out.columns])
    # ---- 列序: ppm, Types, 1H-NMR-v1, [1H-NMR_real], 1H-NMR-v2, [1H-NMR-v3], δH_ppm, nH, review, judge ----
    front = ['ppm', 'Types', '1H-NMR-v1', '1H-NMR_real', '1H-NMR-v2', '1H-NMR-v3',
             'δH_ppm', 'nH']
    rest = [c for c in out.columns if c not in front]
    ordered = [c for c in front if c in out.columns]
    for c in rest:
        if c in ('review', 'judge', 'C_match', 'row_match'):
            ordered.append(c)
    ordered += [c for c in rest if c not in ordered]
    out = out[ordered]
    # ---- 输出两个归属表 sheet: 补全版 (ppm/Types 完整) + 置空版 (补全/续行行置空, 参考 xlsx sheet1 续行格式) ----
    out_blank = out.copy()
    if '_cont' in out_blank.columns:
        m_cont = out_blank['_cont'].fillna(False).astype(bool)
        out_blank['ppm'] = out_blank['ppm'].astype(object)
        out_blank['Types'] = out_blank['Types'].astype(object)
        out_blank.loc[m_cont, ['ppm', 'Types']] = ''
        out_blank = out_blank.drop(columns='_cont')
    out_full = out.drop(columns='_cont', errors='ignore')
    # ppm (碳) 列固定 2 位小数显示; 判定列文本格式 (防 '-' 开头误转公式)
    text_cols = [c for c in ('1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v3', '1H-NMR_real',
                             'review', 'judge', 'C_match', 'row_match', 'assigned_1H',
                             'all_assigned_H', 'aligned_1D', 'Types') if c in out.columns]
    # 双 sheet: 归属表 + 1H_NMR 积分表 (sheet 名沿用文件名前缀)
    # J 列落盘还原为 '(a, b)' 文本 (内存中是 list)
    sig_out = sig_df.copy()
    if 'J' in sig_out.columns:
        sig_out['J'] = sig_out['J'].map(
            lambda js: f"({', '.join(js)})" if isinstance(js, list) and js else
                       (str(js) if not isinstance(js, list) else ''))
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as xw:
        out_blank.to_excel(xw, sheet_name=f'{base}_归属表', index=False)
        out_full.to_excel(xw, sheet_name=f'{base}_归属表补全', index=False)
        sig_out.to_excel(xw, sheet_name=f'{base}_1H_NMR', index=False)
        ws = xw.book.worksheets[2]
        from openpyxl.styles import Alignment, Font
        for col_name in ('δH',):
            if col_name in sig_df.columns:
                ci = list(sig_df.columns).index(col_name) + 1
                for row in ws.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '0.00'
        for col_name in ('J', 'type', 'δH_ppm'):
            if col_name in sig_df.columns:
                ci = list(sig_df.columns).index(col_name) + 1
                for row in ws.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '@'
                        cell.alignment = Alignment(horizontal='left')
        # 两个归属表 sheet 的格式 (补全版 sheet0 + 置空版 sheet1)
        red = Font(color='FF0000')
        col_v = [c for c in ('1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v3') if c in out.columns]
        for si, df_v in ((0, out_full), (1, out_blank)):
            wsi = xw.book.worksheets[si]
            if 'ppm' in df_v.columns:
                ci = list(df_v.columns).index('ppm') + 1
                for row in wsi.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '0.00'
            for col_name in text_cols:
                if col_name in df_v.columns:
                    ci = list(df_v.columns).index(col_name) + 1
                    for row in wsi.iter_rows(min_row=2, min_col=ci, max_col=ci):
                        for cell in row:
                            cell.number_format = '@'
                            cell.alignment = Alignment(horizontal='left')
            # 红色标示: 三步结果列中, 每一步与上一步不同的行 (该三列单元格整行红字)
            for idx_row in range(len(df_v)):
                changed = False
                for k in range(1, len(col_v)):
                    prev = df_v.iloc[idx_row][col_v[k - 1]]
                    curr = df_v.iloc[idx_row][col_v[k]]
                    ps = '' if pd.isna(prev) else str(prev).strip()
                    cs = '' if pd.isna(curr) else str(curr).strip()
                    if ps != cs:
                        changed = True
                        break
                if changed:
                    for cname in col_v:
                        ci = list(df_v.columns).index(cname) + 1
                        wsi.cell(row=idx_row + 2, column=ci).font = red
        # 第三 sheet: 第三步 COSY/HMBC 推断 (仅 do_step3 开启时; 无 HSQC 交叉峰的信号; 全部已归属时跳过)
        if do_step3:
            third_df = third_pass_cosy_hmbc(out, sig_df, h1d_signals, json_dir)
            if not third_df.empty:
                third_df.to_excel(xw, sheet_name=f'{base}_第三步_无HSQC信号', index=False)
                ws2 = xw.book.worksheets[3]
                if 'δH' in third_df.columns:
                    ci = list(third_df.columns).index('δH') + 1
                    for row in ws2.iter_rows(min_row=2, min_col=ci, max_col=ci):
                        for cell in row:
                            cell.number_format = '0.00'
                for col_name in ('δH_ppm', 'type', 'J', 'COSY邻位', 'HMBC碳', '推断结论'):
                    if col_name in third_df.columns:
                        ci = list(third_df.columns).index(col_name) + 1
                        for row in ws2.iter_rows(min_row=2, min_col=ci, max_col=ci):
                            for cell in row:
                                cell.number_format = '@'
                                cell.alignment = Alignment(horizontal='left')
                print(f'[第三步] {len(third_df)} 个无 HSQC 交叉峰信号 → sheet {base}_第三步_无HSQC信号')
    print(f'[判分] {len(cmp_df)} 碳槽: OK={n_ok} 非OK={n_bad}')
    print(f'输出: {out_xlsx}')
    return cmp_df


def main():
    # 可选命令行参数: <化合物名|xlsx路径> [json_dir] [--step3]
    # --step3: 开启第三步 COSY/HMBC 复核 (默认关闭, 输出到 1H-NMR-v2 为止)
    do_step3 = '--step3' in sys.argv
    sys.argv = [a for a in sys.argv if a != '--step3']
    compound = _DEFAULT_COMPOUND
    if len(sys.argv) >= 2:
        compound = sys.argv[1]
    if os.path.isfile(compound) and compound.endswith('.xlsx'):
        xlsx_path = compound
        base = re.sub(r'数据归属\.xlsx$', '', os.path.basename(compound))
    else:
        base = compound
        xlsx_path = os.path.join(_DATA_DIR, f'{base}数据归属.xlsx')
    h_report = os.path.join(_DATA_DIR, f'{base}氢谱积分.txt')
    json_dir = sys.argv[2] if len(sys.argv) >= 3 else _export_mnova_json(base)
    out_dir = os.path.dirname(os.path.abspath(xlsx_path))
    # if not os.path.isfile(xlsx_path):
    #     raise SystemExit(f'ERROR: 归属表不存在: {xlsx_path}')
    if not os.path.isfile(h_report):
        raise SystemExit(f'ERROR: 氢谱积分.txt 不存在: {h_report}')
    print(f'[输入] xlsx={xlsx_path}\n[输入] 积分={h_report}\n[输入] JSON={json_dir}')
    cmp_df = assign(xlsx_path, h_report, json_dir, out_dir, do_step3=do_step3)
    print(cmp_df[['ppm', 'Types', 'assigned_H', 'truth_H', 'judge']].to_string(index=False))


if __name__ == '__main__':
    main()
