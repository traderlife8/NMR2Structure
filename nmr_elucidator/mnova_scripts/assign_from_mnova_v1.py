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
    ~/mnova_export/<化合物名>/   (mnovaBatchExportJson 导出的 JSON)
"""
from __future__ import annotations

import json
import os
import re
import sys

import pandas as pd

# ---- 默认输入配置 (直接运行即用; 换化合物时手工修改这三行) ----
_DEFAULT_COMPOUND = 'clxh119'   # 化合物名
_DATA_DIR = r'X:\OCF投稿\2026-10-11新化合物'
_MNOVA_EXE = r'D:\Program Files\Mestrelab Research S.L\MestReNova\MestReNova.exe'
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
        raise SystemExit(f'ERROR: MestReNova 导出失败 exit={proc.returncode}')
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
    """把 '4.645, m; 3.88, d (11.2)' 类字符串中的氢谱数值四舍五入到 2 位小数."""
    if pd.isna(s) or not str(s).strip():
        return s

    def _sub(m):
        try:
            return f'{round(float(m.group(1)), 2):.2f}'
        except ValueError:
            return m.group(0)
    # 仅替换 δH 位置的数字 (分号/逗号/分号分隔项的开头数字); J 值已是 txt 原文不动
    return re.sub(r'(?<![\w.])(\d+\.\d+)', _sub, str(s))


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
        out.append({'ppm': round(center, 3), 'mult': mult, 'nH': nH, 'J': js})
    return out


# ---- 数据加载 ----

def load_xlsx(path: str) -> pd.DataFrame:
    """读取归属 xlsx: 必需列 ppm / Types(兼容 '类型'); 1H-NMR 可选 (无则不做行级比较).

    CH2 行补全两种化学家写法, 保证每个 CH2 在返回结果中恰好占 2 行 (每行=1 个 H, 由程序归属):
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
    df = df[['ppm', 'Types', '1H-NMR']].copy()
    # ---- 续行归并: ppm 留空但 1H-NMR 有值 → 是上一行的第 2 个 H, 填充 ppm/Types ----
    cont = df['ppm'].isna() & df['1H-NMR'].notna() & df['1H-NMR'].astype(str).str.strip().ne('')
    if cont.any():
        df.loc[cont, 'ppm'] = df['ppm'].ffill()[cont]
        df.loc[cont, 'Types'] = df['Types'].ffill()[cont]
    df['_cont'] = cont  # 续行标记: 最终 CSV 输出时置空 ppm/Types (仅保留 1H-NMR)
    df = df.dropna(subset=['ppm']).reset_index(drop=True)
    # ---- CH2 单行补全: 同 (ppm, CH2) 只出现 1 次时, 紧跟原行插入相同一行 (1H-NMR 留空, 由程序补) ----
    ch2_count: dict = {}
    for _, r in df.iterrows():
        if str(r['Types']).strip() == 'CH2':
            key = round(float(r['ppm']), 2)
            ch2_count[key] = ch2_count.get(key, 0) + 1
    rows = []
    for _, r in df.iterrows():
        rows.append({'ppm': r['ppm'], 'Types': r['Types'], '1H-NMR': r['1H-NMR'],
                     '_cont': bool(r['_cont'])})
        if str(r['Types']).strip() == 'CH2':
            key = round(float(r['ppm']), 2)
            if ch2_count[key] == 1:
                rows.append({'ppm': r['ppm'], 'Types': r['Types'], '1H-NMR': None,
                             '_cont': True})  # 补全行: 输出时与续行同样置空 ppm/Types
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


# ---- 核心归属 ----

def assign(xlsx_path: str, h_report_path: str, json_dir: str, out_dir: str) -> pd.DataFrame:
    df = load_xlsx(xlsx_path)
    h1d_signals = parse_1h_report(h_report_path)
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
            else:
                nmr_str = ''
            results.append({
                'ppm': float(row['ppm']), 'Types': slot['type'],
                '1H-NMR': nmr_str,
                'assigned_1H': row_outputs[k],
                'all_assigned_H': '; '.join(h_strs) if k == 0 else '',
                'aligned_1D': ('; '.join(f'{h_strs[i]}{"✓" if flags[i] else "✗"}'
                                         for i in range(len(h_strs)))) if k == 0 else '',
                'note': (f"槽需 {slot['n_h']}H, HSQC 命中 {len(slot['hs'])}"
                         if k == 0 else ''),
                '_cont': bool(row.get('_cont', False)),
            })
    out = pd.DataFrame(results)
    # 排序: ppm 从大到小; 同 ppm 按 assigned_1H 从大到小 (Cq/空值排最后)
    # ascending=False: ppm 与 δH 均降序
    out['_h_sort'] = out['assigned_1H'].map(
        lambda x: float(x) if not pd.isna(x) and str(x).strip() else -9e9)
    out = out.sort_values(['ppm', '_h_sort'], ascending=False).drop(columns='_h_sort').reset_index(drop=True)
    os.makedirs(out_dir, exist_ok=True)
    # 输出文件名带化合物名前缀 (取 xlsx 文件名去 "数据归属.xlsx"), 防止同目录多化合物覆盖
    base = os.path.basename(xlsx_path)
    base = re.sub(r'数据归属\.xlsx$', '', base) or 'compound'
    out_csv = os.path.join(out_dir, f'{base}_归属表.csv')
    out.to_csv(out_csv, index=False, encoding='utf-8-sig')

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
    prog_rows = out[['ppm', 'Types', '1H-NMR']].copy()
    prog_rows['_k'] = prog_rows['1H-NMR'].map(_h_key)
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
        for col in ('ppm', 'Types', '1H-NMR'):
            vt = _norm(rt[col]) if rt is not None else '<缺失>'
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
            # 定位回 out 原行: 按 ppm+Types+1H-NMR 三元组; 差异时附上 δH 差值
            label = '一致'
            if not ok_row:
                label = '差异'
                nh_p, nh_t = _first_num(rp['1H-NMR']), _first_num(rt['1H-NMR']) if rt is not None else None
                if nh_p is not None and nh_t is not None:
                    label += f' (ΔδH={nh_p - nh_t:+.2f})'
            row_match_map[(round(float(rp['ppm']), 2), str(rp['Types']).strip(),
                           _norm(rp['1H-NMR']))] = label
    out['row_match'] = out.apply(
        lambda r: row_match_map.get(
            (round(float(r['ppm']), 2), str(r['Types']).strip(), _norm(r['1H-NMR'])), ''),
        axis=1)

    # 重新排序后落盘; ppm/δH 均降序
    out['_h_sort'] = out['assigned_1H'].map(
        lambda x: float(x) if not pd.isna(x) and str(x).strip() else -9e9)
    out = out.sort_values(['ppm', '_h_sort'], ascending=False).drop(columns='_h_sort').reset_index(drop=True)

    # 氢谱数值统一四舍五入保留 2 位小数
    out['ppm'] = out['ppm'].map(lambda x: round(float(x), 2))
    for col in ('assigned_1H', 'all_assigned_H', 'aligned_1D'):
        if col in out.columns:
            out[col] = out[col].map(_round_h_in_str)

    # xlsx 有 1H-NMR truth 时: 插入 1H-NMR_real (xlsx 按 out 相同规则排序后行对齐复制) 于 1H-NMR 列后
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
            idx = col_order.index('1H-NMR') + 1
            col_order.insert(idx, '1H-NMR_real')
            out = out[col_order]
    else:
        out = out.drop(columns=[c for c in ('1H-NMR_real', 'row_match') if c in out.columns])
    # ---- 最终输出格式: 补全/续行 (CH2 第 2 行) 置空 ppm/Types, 仅保留 1H-NMR (参考 xlsx sheet1 续行格式) ----
    if '_cont' in out.columns:
        m_cont = out['_cont'].fillna(False).astype(bool)
        out['ppm'] = out['ppm'].astype(object)
        out['Types'] = out['Types'].astype(object)
        out.loc[m_cont, ['ppm', 'Types']] = ''
        out = out.drop(columns='_cont')
    out.to_csv(out_csv, index=False, encoding='utf-8-sig')
    cmp_csv = os.path.join(out_dir, f'{base}_归属表_对比.csv')
    cmp_df.to_csv(cmp_csv, index=False, encoding='utf-8-sig')
    print(f'[判分] {len(cmp_df)} 碳槽: OK={n_ok} 非OK={n_bad}')
    print(f'输出: {out_csv}\n对比: {cmp_csv}')
    return cmp_df


def main():
    # 可选命令行参数覆盖默认配置: <化合物名|xlsx路径> [json_dir]
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
    cmp_df = assign(xlsx_path, h_report, json_dir, out_dir)
    print(cmp_df[['ppm', 'Types', 'assigned_H', 'truth_H', 'judge']].to_string(index=False))


if __name__ == '__main__':
    main()
