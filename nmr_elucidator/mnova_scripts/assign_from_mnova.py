"""从 MNova JSON 导出物 + 氢谱积分.txt + 数据归属.xlsx (ppm/Types) 推导 1H-NMR 归属列.

泛化算法 (不依赖特定化合物):
  1. 解析氢谱积分.txt: "δ (mult, nH), ..." → 精细 1H 信号单 [{ppm, mult, J, nH}]
     (化学家报告值, 比仪器 GSD 峰单更可靠 — 融合的多重峰只占一行)
  2. 碳槽: xlsx 同 (ppm, Types) 合并; 需求 H 数:
     - CH: 1   CH2: 2   CH3: 1 (3H 等价, 若 HSQC 实际分离出 >1 个 δH 则多输出)
     - 理论 CHn 输出 n 个不同 δH; HSQC 交叉峰不足 n 个时按实际数输出
  3. HSQC 候选: |ΔC| <= tol; 强度 <= 1.0 的 GSD 弱峰剔除
  4. 全局一对一匹配: 密集簇 (2.5 ppm 内 >=3 槽 >=3 峰) 匈牙利全局最优,
     代价 = |ΔC| + 0.5*|ΔH对齐| (CH3 与 1D nH=3 信号距离罚减); 其余贪心兜底;
     第二遍复用扫描 (ΔC<=1.0 邻槽重叠融合 H) + 第三遍融合信号借用 (1D nH 余量)
     + CH2 补第 2 个 H (1D nH>=2 未占用信号)
  5. 1H-NMR-v1: HSQC δH 对齐到最近 1D 1H 信号 (±tol), 格式 'δ, mult (J)';
     同步填 nH (1D 积分 H 数) 与 δH_ppm (积分表原文); v2 = v1
  6. 归属完整性迭代复核 (结合 2D HSQC, 内存 df): 检查 1H_NMR 的 δH_nH vs nH,
     差异惩罚 = Σ|差| (nH==1 行权重 1000 不允许差异, 其余权重 1/nH — nH 越大容差越大);
     逐轮 C-H 对重匹配: 途径 A 参考驱动 (如有参考: 参考 δH + ΔH整体 偏离超 _REF_H_TOL
     的行, 2D 中心支撑下重归属; nH==1 必须有 2D 支撑) + 途径 B 2D C-H 对交换
     (缺 H 信号从占用过多信号移入行; 1D J 类型规则: s/t 取截面峰顶, d/q 取双峰
     中间谷; 匹配值允许落在 δH_ppm 区间内任意值); 轮换候选直到全对/无改进/达 _MAX_RETRY
  7. HSQC 矩阵复核 (review): 逐行验证 (δC, δH) 落在交叉峰中心 ±tol;
     融合峰双高斯拆分 (1D nH 提示 + 残差验证), 拆分中心可回填修正行;
     修正结果写新列 1H-NMR-v2.5 (未修正行 v2.5 = v2), 修正行 nH/δH_ppm 同步更新
  8. [--ref 参考xlsx] 参考化合物归属 → v2.6 (v2.5 之后执行; 无参考则跳过不输出 v2.6):
     a. 参考表拆分/补全 (δH 单元格换行 + 续行向下填充, ppm/Types/No. 三列补全);
     b. C 对应 (糖基/苷元分区): 未知糖基碳 (No. 带 ') 与参考糖基碳、苷元碳与参考
        苷元碳各自分区内 |Δppm(C)| 最小一对一匹配 (同 No. 直接对应); 新增 Δppm(C) 列
        (未知 ppm − δc参考);
     c. 整体偏差校准: C 对应行统计 C/H 整体偏差 (中位数), H 容差 = _REF_H_TOL + |ΔH整体|;
     d. H 归属: 参考 C 对应的 δH参考 + ΔH整体 在自身 1D 积分表 (±容差) 找最近
        未占用信号, 补 v2 空行 → v2.6 (未补行 v2.6 = v2.5);
        参考列 No.参考/δc参考/Δppm(C)/δH参考 全组行齐全 (CH2 第二行不缺)
  9. [--step3] 第三步 COSY/HMBC 全表复核 (默认关闭): 不改写 H (v3 = 基准列),
     复核证据独立输出到第三 sheet (COSY 邻位一致性 / HMBC 碳列表 / 无 HSQC 信号推断)
 10. 输出: <base>_归属表.xlsx (置空版 + 补全版 + 1H_NMR 积分表, + 糖基/苷元拆分 sheet:
     No./Pos. 列撇号归一后带 ' 的行 = 糖基, 不带 ' = 苷元); 判分列 (judge 仅在有
     1H-NMR truth 时输出), 每步结果变化的行红字标示, 1H_NMR 中 δH_nH != nH 行红字,
     冻结首行表头 + 隔行浅黄背景
基准列 (_ver): 有参考 = v2.6, 无参考 = v2.5; --step3 时 v3 = 基准列并作为 _ver

用法:
  直接运行 (默认处理脚本顶部 _DEFAULT_COMPOUND 指定的化合物):
    python assign_from_mnova.py
  换化合物: 手工修改顶部 _DEFAULT_COMPOUND; 或命令行传入化合物名:
    python assign_from_mnova.py clxh101
  可选开关: --step3 (开启第三步 COSY/HMBC 复核, 默认关闭)
            --ref <参考xlsx> (类似已知化合物 H 归属表, 按 No. 对齐补齐未归属 H)
  输入约定 (相对 _DATA_DIR / mnova 导出目录, 文件名兼容有无下划线):
    <_DATA_DIR>/<化合物名>数据归属.xlsx
    <_DATA_DIR>/<化合物名>氢谱积分.txt
    <_DATA_DIR>/<化合物名>_H参考.xlsx (可选, 兼容 H参考/_参考/参考 变体)
    <_DATA_DIR>/mnova_export/<化合物名>/   (mnovaBatchExportJson 导出的 JSON)
"""
from __future__ import annotations

import json
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from numpy.lib.stride_tricks import sliding_window_view

import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from scipy import ndimage
from scipy.optimize import curve_fit, linear_sum_assignment

# ---- 默认输入配置 (直接运行即用; 换化合物时手工修改 _DEFAULT_COMPOUND 一行,
#      其余路径/开关均随化合物名联动; 命令行参数可覆盖) ----
_DEFAULT_COMPOUND = 'clxh141C'   # 化合物名
_DATA_DIR = rf'Y:\LXH数据处理\新化合物处理\程序解析\{_DEFAULT_COMPOUND}解析'
_MNOVA_EXE = r'D:\Program Files\Mestrelab Research S.L\MestReNova\MestReNova.exe'
# row_match 差异容差 (ppm): |δH_程序 - δH_xlsx| <= 此值时判 '一致 (ΔδH=...)', 保留差异值
_ROW_MATCH_TOL = 0.015
# mnova JSON 输出根 = _DATA_DIR/mnova_export (qs 从文档路径动态推导同一目录, 两端一致)
_MNOVA_OUT_ROOT = os.path.join(_DATA_DIR, 'mnova_export')
# 参考化合物 H 补归属: 参考 δH 与自身 1D 信号匹配容差 (类似化合物 δH 允许偏移);
# 迭代复核重归属时的 HSQC 中心候选 δH 窗口同样用 δH ± _REF_H_TOL
_REF_H_TOL = 0.3
# v2.6 参考归属容差: 若 v1~v2.6 各版本中该行 H 数值不一致, 取与参考 δH差别最小的版本;
# 差异 < _REF_MATCH_TOL → 直接 OK 采用; 否则再在 1D 信号中找与参考 δH 差异更小的 H
_REF_MATCH_TOL = 0.05
# 归属完整性迭代复核最大尝试次数 (δH_nH 与 nH 不一致时轮换重归属, 达上限留人工复核)
_MAX_RETRY = 100
# 可选开关默认值 (命令行 --step3 / --ref 可覆盖):
# 第三步 COSY/HMBC 复核, 默认关闭
_DEFAULT_DO_STEP3 = False
# 参考化合物 H 归属表 (空 = 不使用): 默认 <_DATA_DIR>/<化合物名>_H参考.xlsx,
# 兼容有无下划线及 H 前缀 ('clxh126A_H参考' / 'clxh126AH参考' / 'clxh126A_参考' /
# 'clxh126A参考' 均可); 修改 _DEFAULT_COMPOUND 即自动切换到该化合物的参考表
_DEFAULT_REF = lambda: _find_data_file(
    _DATA_DIR, _DEFAULT_COMPOUND, '_H参考', extra_variants=('H参考', '_参考', '参考'))


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


def _find_data_file(data_dir: str, compound: str, suffix: str, ext: str = '.xlsx',
                    extra_variants: tuple = ()) -> str:
    """在 data_dir 中查找化合物数据文件, 兼容化合物名与后缀之间有无下划线.

    变体顺序: '<compound><suffix>' → '<compound>_<suffix>' (suffix 自带下划线时
    两种自然覆盖 'clxh126A参考' 与 'clxh126A_参考'); extra_variants 追加其它变体
    (如 H参考 与 参考: '<c>_H参考.xlsx' / '<c>H参考.xlsx' / '<c>_参考.xlsx' /
    '<c>参考.xlsx'). 找不到返回 ''.
    """
    cands = [f'{compound}{suffix}{ext}', f'{compound}_{suffix}{ext}']
    cands += [f'{compound}{v}{ext}' for v in extra_variants]
    for fn in cands:
        p = os.path.join(data_dir, fn)
        if os.path.isfile(p):
            return p
    return ''


def _export_mnova_json(compound: str, data_dir: str = '') -> str:
    """调用 MestReNova.exe 无界面导出 .mnova → JSON 目录 (<data_dir>/mnova_export/<compound>). 返回 JSON 目录.

    qs 脚本由本文件同目录的 mnovaBatchExportJson.qs 部署 (每次运行校验), 输出根由 qs 从
    文档路径动态推导为 <.mnova 所在目录>/mnova_export. data_dir 缺省用顶部 _DATA_DIR
    (命令行传化合物名时随化合物名联动).
    """
    data_dir = data_dir or _DATA_DIR
    mnova_path = os.path.join(data_dir, f'{compound}.mnova')
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
    out_dir = os.path.join(data_dir, 'mnova_export', compound)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
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
    # ---- 糖基 sheet 合并 (clxh129 等): 数据归属 xlsx 含第二个 sheet 且该 sheet
    # 有 Pos. 列时, 把糖基 Pos. 按 (ppm 取 2 位小数, 类型) 对应回填到主 sheet 的
    # No. 列, 然后写回原 Excel (主 sheet 落盘保存 Pos. 列, 糖基 sheet 保留不动).
    try:
        xl_all = pd.ExcelFile(path)
        gly_sheet = next((s for s in xl_all.sheet_names[1:]
                          if 'Pos' in str(pd.read_excel(xl_all, sheet_name=s, nrows=0).columns)
                          .replace('.', '')), None)
    except (OSError, ValueError) as e:
        print(f'[糖基合并][WARN] 读取糖基 sheet 失败 ({e}), 跳过合并: {path}')
        gly_sheet = None
    if gly_sheet is not None:
        gly = pd.read_excel(xl_all, sheet_name=gly_sheet)
        gly.columns = [str(c).strip() for c in gly.columns]
        c_pos = next((c for c in gly.columns
                      if str(c).strip().rstrip('.').lower() == 'pos'), None)
        if c_pos is not None:
            merged = 0
            gly_type_col = 'Types' if 'Types' in gly.columns else \
                ('类型' if '类型' in gly.columns else None)
            # 糖基 Pos. 按 (ppm 2位小数, 类型) 分组 (保序列表) — 原位写回用
            _fill_map: dict = {}
            for _, gr in gly.iterrows():
                if pd.isna(gr.get('ppm')) or pd.isna(gr[c_pos]):
                    continue
                typ_g = '' if gly_type_col is None or pd.isna(gr[gly_type_col]) \
                    else str(gr[gly_type_col]).strip()
                _fill_map.setdefault((round(float(gr['ppm']), 2), typ_g), []).append(
                    str(gr[c_pos]).strip())
            # df 内存标记 (供本次归属流程使用 No. 列); 原位写回在下方 openpyxl 执行
            if _fill_map:
                if 'No.' not in df.columns:
                    df['No.'] = ''
                    df['No.'] = df['No.'].astype(object)
                for (k_p, k_t), poss in _fill_map.items():
                    m = (df['ppm'].map(lambda v: not pd.isna(v) and round(float(v), 2) == k_p)
                         & df['Types'].astype(str).str.strip().eq(k_t))
                    n = 0
                    for oi in df.index[m & df['No.'].astype(str).str.strip().isin(['', 'nan'])]:
                        if n < len(poss):
                            df.at[oi, 'No.'] = poss[n]
                            n += 1
                            merged += 1
            if merged:
                # 写回原 Excel: openpyxl 原位只在主 sheet 追加 No. 列 (表头+新填的
                # 空位值), 不重写任何已有单元格 — 原始文档的格式/文本完全保留.
                # (此前用 pandas 整 sheet 重写会丢失原格式并覆盖原始文本, 已废弃)
                from openpyxl import load_workbook as _lw
                wb_w = _lw(path)
                ws_w = wb_w.worksheets[0]
                # No. 列已存在则定位, 否则追加到表尾
                no_ci = next((c for c, cell in enumerate(ws_w[1], 1)
                              if str(cell.value).strip() == 'No.'), None)
                if no_ci is None:
                    no_ci = ws_w.max_column + 1
                    ws_w.cell(row=1, column=no_ci, value='No.')
                # 逐数据行: 该行 (ppm, 类型) 在本次回填集合中且 No. 空 → 写 Pos. 原文
                tci = next((c for c, cell in enumerate(ws_w[1], 1)
                            if str(cell.value).strip() in ('Types', '类型')), 1)
                for r_w in range(2, ws_w.max_row + 1):
                    c_v = ws_w.cell(row=r_w, column=1).value
                    if c_v is None:
                        continue
                    k_w = (round(float(c_v), 2),
                           str(ws_w.cell(row=r_w, column=tci).value).strip())
                    # 糖基 (ppm, 类型) 匹配且该行 No. 空 → 按序取一个未写的 Pos.
                    if k_w in _fill_map:
                        cur = ws_w.cell(row=r_w, column=no_ci).value
                        if cur is None or str(cur).strip() == '':
                            ws_w.cell(row=r_w, column=no_ci,
                                      value=_fill_map[k_w].pop(0))
                            if not _fill_map[k_w]:
                                del _fill_map[k_w]
                wb_w.save(path)
                print(f'[糖基合并] {gly_sheet}: {merged} 个 Pos. 按 (ppm, 类型) 回填主 sheet No. 列 (原位写入, 原格式保留): {path}')
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


def _norm_prime(v) -> str:
    """撇号归一: 全角/半角撇号统一 — 单引号字符 (′ U+2032 / ' U+FF07 / ') → 一个
    半角 ' ; 双引号字符 (″ U+2033 / " U+FF02 / ") → 两个半角 '."""
    s = str(v)
    s = (s.replace('\u2033', "''").replace('\uff02', "''").replace('\"', "''")
          .replace('\u2032', "'").replace('\uff07', "'"))
    return s


def _norm_no_key(v) -> str:
    """No./Pos. 键归一: 撇号归一 (_norm_prime), 去 Api-/Rha- 前缀, 取数字+撇号部分.

    连接位写法兼容: '4'''-Rha-1'''' (异头碳前带糖名) 取**首个** 数字+'* (连接位 4''');
    普通写法 ('4''' / '6''a' / '6′α') 取末尾 数字+'* (去掉末尾 a/b/α/β 构型字母 —
    代表 CH2 的两个不同 H).
    """
    s = _norm_prime(str(v).strip())
    s = re.sub(r'^(Api|Rha)-', '', s)
    # 'N'-Rha-... / N-O-Glc-... 型: 连接位编号在最前, 取第一个 数字+'*
    m_first = re.match(r"^(\d+'*)[-\u2014]", s)
    if m_first:
        return m_first.group(1)
    m = re.search(r"(\d+)'*[abαβ]?$", s)
    if not m:
        return s
    # 去掉末尾构型字母 (a/b/α/β — CH2 两个不同 H), 同位置两行归一组
    key = m.group(0)
    return re.sub(r"[abαβ]$", '', key)


def load_ref_h(path: str) -> pd.DataFrame:
    """读取参考化合物 H 归属表 (必需列 δH; No./δc/Types 可选, 空缺自动补全).

    列名兼容: No./No, δ c/δc/ppm, δH/1H-NMR, Types/类型.
    参考表存在两种多行写法, 均拆分/补全为每行一个 H:
    1. δH 单元格内换行 (\\n) → 拆成多行;
    2. No./δc 留空的续行 → 向下填充归并 (仅 δH 以数值开头才算续行, 排除 '3-Api' 文字行).
    续行 / 单元格换行拆分出的行, 其 No./δc/Types 一并向下填充 (ppm/Types/No. 三列补全).
    返回列: No.参考 / no_key / δc参考 / δH参考 / Types参考 (可空).
    """
    df = pd.read_excel(path, header=0)
    df.columns = [str(c).strip() for c in df.columns]
    c_no = next((c for c in df.columns
                 if str(c).strip().rstrip('.').lower() in ('no', 'pos')), None)
    c_c = next((c for c in df.columns if c in ('δ c', 'δc', 'ppm')), None)
    c_h = next((c for c in df.columns if c in ('δH', '1H-NMR')), None)
    c_t = next((c for c in df.columns if c in ('Types', '类型')), None)
    if c_h is None:
        raise SystemExit(f'ERROR: 参考 xlsx 缺少 δH (或 1H-NMR) 列: {path}')
    ren = {c_h: 'dH'}
    if c_no is not None:
        ren[c_no] = 'No.'
    if c_c is not None:
        ren[c_c] = 'dC'
    if c_t is not None:
        ren[c_t] = 'dT'
    df = df.rename(columns=ren).copy()
    for c, default in (('No.', ''), ('dC', None), ('dT', '')):
        if c not in df.columns:
            df[c] = default
    # 脏数据清理: 数值与多重度/J 值中间的多余空格直接去掉 ('1. 76 m' → '1.76 m',
    # '2.29q(7. 2)' → '2.29q(7.2)'), 避免 δH 解析截断成 1.0 导致排序/匹配错乱
    df['dH'] = df['dH'].map(
        lambda v: '' if pd.isna(v) else
        re.sub(r'(\d)\.\s+(\d)', r'\1.\2',          # 小数点后空格: '1. 76' → '1.76'
               re.sub(r'([a-z]\()\s+', r'\1',      # 括号内开头空格: 'q( 7.2)' → 'q(7.2)'
                      re.sub(r'\s+(?=[\d.]+\))', '',  # 括号内结尾空格: '(7. 2 )' → '(7. 2)'→再由上一条处理
                             str(v).replace('\r\n', '\n').strip()))).strip())
    # 续行归并: No. 留空且 δH 以化学位移数值开头 → 上一行的第 2 个 H;
    # ppm (dC) / Types (dT) / No. 三列向下填充补全 (clxh126A 参考表 δc 也有续行留空)
    _is_shift = df['dH'].str.match(r'^\d+\.\d+')
    cont = df['No.'].isna() & df['dH'].ne('') & _is_shift
    df['_cont'] = cont
    if cont.any():
        for c in ('No.', 'dC', 'dT'):
            df.loc[cont, c] = df[c].ffill()[cont]
    # 单元格内换行拆分: 每个子行一行, No./dC/dT 仅首行保留 (拆分行同样三列补全)
    rows = []
    for _, r in df.iterrows():
        parts = [p.strip() for p in str(r['dH']).split('\n') if p.strip()]
        if len(parts) <= 1:
            row = dict(r)
            row['dH'] = parts[0] if parts else ''
            rows.append(row)
            continue
        for i, p in enumerate(parts):
            row = dict(r)
            row['dH'] = p
            if i > 0:
                row['_cont'] = True
            rows.append(row)
    out = pd.DataFrame(rows)
    out = out.dropna(subset=['dC']).reset_index(drop=True)
    out['no_key'] = out['No.'].map(_norm_no_key)
    out = out.rename(columns={'No.': 'No.参考', 'dC': 'δc参考', 'dH': 'δH参考', 'dT': 'Types参考'})
    return out[['No.参考', 'no_key', 'δc参考', 'δH参考', 'Types参考']]


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
        mat = np.array(m, dtype=float)
        c_ax = np.linspace(float(doc['firstXppm']), float(doc['lastXppm']), mat.shape[0])
        h_ax = np.linspace(float(doc['firstYppm']), float(doc['lastYppm']), mat.shape[1])
        return mat, c_ax, h_ax
    except (ValueError, TypeError) as e:
        print(f'[复核][WARN] HSQC 矩阵加载失败: {e}')
        return None


def _cluster_blobs(mat, c_ax, h_ax):
    """连通域分割: 阈值以上像素聚成峰团. 返回 [{rows, cols, peak}] 列表 (按峰高降序)."""
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


def check_h_counts(out: pd.DataFrame, sig_df: pd.DataFrame) -> pd.DataFrame:
    """归属完整性检查 (与 数据处理程序/check_h_assignment.py 同规则).

    统计归属表各行 δH_ppm 占用的 H 数: 同 (ppm, Types) 槽仅 1 行且 CH3 → 权重 3
    (CH3 未拆行时一行代表 3H), 其余每行 1; 未归属行 (δH_ppm 空) 不计数.
    返回 sig_df 副本附 δH_nH 列 (该信号在归属表中被占用的 H 数, 未出现为 0).
    """
    slot_rows: dict = {}
    counts: Counter = Counter()
    for _, r in out.iterrows():
        s = '' if pd.isna(r.get('δH_ppm')) else str(r['δH_ppm']).strip()
        if not s or s.lower() == 'nan':
            continue
        key = (round(float(r['ppm']), 2), str(r['Types']).strip())
        slot_rows[key] = slot_rows.get(key, 0) + 1
    for _, r in out.iterrows():
        s = '' if pd.isna(r.get('δH_ppm')) else str(r['δH_ppm']).strip()
        if not s or s.lower() == 'nan':
            continue
        typ = str(r['Types']).strip()
        key = (round(float(r['ppm']), 2), typ)
        w = 3 if (typ == 'CH3' and slot_rows.get(key, 0) == 1) else 1
        counts[s] += w
    sig_df = sig_df.copy()
    sig_df['δH_nH'] = sig_df['δH_ppm'].map(lambda s: counts.get(str(s).strip(), 0))
    return sig_df


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
                      h1d_signals, json_dir: str, sig_df=None,
                      centers_out: list = None) -> pd.DataFrame:
    """第二阶段矩阵复核: 逐行验证初步归属, 融合峰拆分修正.

    返回带 'review' 列的 out (复核结论: OK / 偏移 / 复核修正 / 未命中).
    修正结果写入新列 '1H-NMR-v2.5' (以 1H-NMR-v2 为基准, 参考 H 补归属之后执行);
    回填仅使用拆分中心 (融合峰拆出的第 2 中心), 且要求本行 δC 与中心 δC 匹配、
    中心 δH 对齐到某 1D 信号 (有真实氢谱峰支撑), 避免激进改写.
    centers_out: 传入 list 时, 拆分中心集合 [(c_ppm, h_ppm, amp)] 写入其中
    (供归属完整性迭代复核轮换候选用).
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
                out.at[i, '1H-NMR-v2.5'] = fmt_1h_nmr(h_aligned, s1d)
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
    if centers_out is not None:
        centers_out.clear()
        centers_out.extend(centers)
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

def fill_missing_h_from_ref(out: pd.DataFrame, sig_df: pd.DataFrame,
                            h1d_signals, ref_path: str,
                            ref_groups_out: dict = None) -> pd.DataFrame:
    """参考已知化合物归属表, 按 No. 的 ppm 差异重排 C 对应后补齐未归属 H.

    流程 (糖基/苷元分开处理):
    1. 参考表 load_ref_h 拆分/补全为每行一个 H;
    2. 新增 Δppm(C) 列: 同 No. 归一键下未知 ppm − δc参考 (参考无该 No. 则空);
    3. C 对应 (糖基/苷元分开): 未知化合物的糖基碳 (No. 归一键带 ') 与参考糖基碳
       按 |Δppm(C)| 最小一一对应 (同 No. 直接对应; 无同 No. 时在分区内按 δc 差最小
       匹配, 一对一互斥); 苷元同理. 对应结果回写三个参考列 (糖基/苷元各自区内对齐).
    4. 整体偏差: 用 C 对应成立的行统计 C/H 整体偏差 (中位数),
       H 匹配容差 = _REF_H_TOL + |ΔH整体|;
    5. H 归属: 参考 C 对应的参考 H (δH参考) + ΔH整体 在自身 1D 积分表 (±容差内)
       找最近未占用信号, 填 1H-NMR-v2 / δH_ppm / nH / assigned_1H (苷元同样规则).
    返回追加列后的 out (行数不变).
    """
    ref = load_ref_h(ref_path)
    out = out.copy()
    # 参考列先建为 object 空列 (避免 FutureWarning / 后续混合 dtype 写入报错)
    for col in ('No.参考', 'δc参考', 'δH参考'):
        out[col] = ''
        out[col] = out[col].astype(object)
    out['_no_key'] = out['No.'].map(_norm_no_key) if 'No.' in out.columns else ''
    # 糖基/苷元分区: 归一化原始 No. 中含 ' 即糖基 (6'a / 3-O-Glc-1' 等, 任何撇号)
    out['_is_gly'] = out['No.'].map(
        lambda v: isinstance(v, str) and "'" in _norm_prime(v)) \
        if 'No.' in out.columns else False
    ref['_is_gly'] = ref['No.参考'].map(
        lambda v: isinstance(v, str) and "'" in _norm_prime(v))
    # ---- C 对应 (糖基/苷元各自分区内, |Δppm| 最小一对一匹配) ----
    # 同 No. 直接对应; 其余未知碳在分区内找 δc 差最小的未用参考碳 (匈牙利式贪心,
    # 按差值升序全局占坑, 保证一对一互斥)
    for is_gly in (True, False):
        uk = out[out['_is_gly'] == is_gly]
        rf = ref[ref['_is_gly'] == is_gly]
        if rf.empty:
            continue
        # 同 No. 直接对应
        direct: dict = {}  # out idx -> 参考行
        used_ref: set = set()
        for oi in uk.index:
            k = out.at[oi, '_no_key']
            rows = rf[rf['no_key'] == k]
            if not rows.empty:
                # 同 No. 多行 (CH2 两行): 各取一行, 按 δH 降序已排
                n = sum(1 for v in direct.values() if v['no_key'] == k)
                if n < len(rows):
                    direct[oi] = rows.iloc[n]
                    used_ref.add(rows.index[n])
        # 同 No. 未覆盖的未知碳: 分区内 |δc 差| 最小的未用参考碳
        pend = [oi for oi in uk.index if oi not in direct]
        cands = []
        for oi in pend:
            try:
                c_uk = float(out.at[oi, 'ppm'])
            except (TypeError, ValueError):
                continue
            for ri, rr in rf.iterrows():
                if ri in used_ref:
                    continue
                try:
                    d = abs(c_uk - float(rr['δc参考']))
                except (TypeError, ValueError):
                    continue
                cands.append((d, oi, ri))
        cands.sort()
        for d, oi, ri in cands:
            if oi in direct or ri in used_ref:
                continue
            direct[oi] = rf.loc[ri]
            used_ref.add(ri)
        # 写参考列 + Δppm(C) 列
        if 'Δppm(C)' not in out.columns:
            out['Δppm(C)'] = ''
            out['Δppm(C)'] = out['Δppm(C)'].astype(object)
        for oi, rr in direct.items():
            out.at[oi, 'No.参考'] = rr['No.参考']
            out.at[oi, 'δc参考'] = rr['δc参考']
            out.at[oi, 'δH参考'] = rr['δH参考']
            try:
                out.at[oi, 'Δppm(C)'] = round(
                    float(out.at[oi, 'ppm']) - float(rr['δc参考']), 2)
            except (TypeError, ValueError):
                pass
    # ---- 偏差校准: 共同归属行的 C/H 整体偏差 (中位数), 叠加到信号匹配容差 ----
    # 共同行 = 参考列已配对 且 参考行 δH 为数值 且 目标行自身 v2 已归属 (有 δH)
    d_cs, d_hs = [], []
    for oi in out.index:
        rv = _ref_h_val(out.at[oi, 'δH参考'])
        if rv is None or str(out.at[oi, 'δc参考']).strip() in ('', 'nan'):
            continue
        pv = _ref_h_val(out.at[oi, '1H-NMR-v2'])
        if pv is None:
            continue
        d_hs.append(pv - rv)
        try:
            d_cs.append(float(out.at[oi, 'ppm']) - float(out.at[oi, 'δc参考']))
        except (TypeError, ValueError):
            pass
    d_h_med = statistics.median(d_hs) if d_hs else 0.0
    d_c_med = statistics.median(d_cs) if d_cs else 0.0
    h_tol_eff = _REF_H_TOL + abs(d_h_med)  # H 整体偏差叠加到匹配容差
    print(f'[参考补H] 整体偏差: ΔC中位={d_c_med:+.2f} (n={len(d_cs)}), '
          f'ΔH中位={d_h_med:+.2f} (n={len(d_hs)}); H 匹配容差={h_tol_eff:.2f}')
    # ---- H 归属: 参考 C 对应的参考 H (δH参考) + ΔH整体 对齐自身 1D 信号 ----
    # 糖基/苷元均已由 C 对应得到 δH参考, 统一处理; 同一 1D 信号不重复占用.
    used_sig: set = set()
    m_fill = out['1H-NMR-v2'].isna() | \
        (out['1H-NMR-v2'].astype(str).str.strip().isin(['', 'nan']))
    m_has_h = out['Types'].astype(str).str.strip().isin(['CH', 'CH2', 'CH3'])
    for oi in out.index[m_fill & m_has_h]:
        ref_h = _ref_h_val(out.at[oi, 'δH参考'])
        if ref_h is None:
            continue  # 参考 C 无对应 H (季碳) 或未对应 → 留待人工
        # 参考 H + H 整体偏差 ±叠加容差 内最近的未占用信号
        target_h = ref_h + d_h_med
        s_best, d_best = None, 9e9
        for s in h1d_signals:
            if id(s) in used_sig:
                continue
            d = abs(s['ppm'] - target_h)
            if d <= h_tol_eff and d < d_best:
                s_best, d_best = s, d
        if s_best is None:
            continue
        used_sig.add(id(s_best))
        h_str = f"{s_best['ppm']:.2f}"
        out.at[oi, '1H-NMR-v2'] = fmt_1h_nmr(h_str, s_best)
        out.at[oi, 'assigned_1H'] = h_str
        if s_best.get('nH') is not None:
            out.at[oi, 'nH'] = int(s_best['nH'])
        if 'δH_ppm' in out.columns and 'δH' in sig_df.columns:
            m_sig = sig_df[pd.to_numeric(sig_df['δH'], errors='coerce')
                           .sub(s_best['ppm']).abs().le(0.011)]
            if not m_sig.empty:
                out.at[oi, 'δH_ppm'] = str(m_sig.iloc[0]['δH_ppm'])
    return out.drop(columns=['_no_key', '_is_gly'], errors='ignore')


def _ref_groups_of(ref_path: str) -> dict:
    """读取参考表并按 no_key 分组 (各组 δH 降序). 供迭代复核使用参考数据."""
    ref = load_ref_h(ref_path)
    groups: dict = {}
    for _, r in ref.iterrows():
        groups.setdefault(r['no_key'], []).append(r)
    for k in groups:
        groups[k].sort(key=lambda r: -(_ref_h_val(r['δH参考']) or -9e9))
    return groups


def _ref_h_val(s):
    """参考 δH 字符串 → 首个数值 (无数值返回 None).

    兼容小数点后带空格的写法 ('1. 76 m' → 1.76, '7. 2' → 7.2), 否则会被
    截断成 1.0 导致组内降序排序错乱."""
    m = re.search(r'-?\d+\.?\d*', str(s))
    return float(m.group()) if m else None


def assign(xlsx_path: str, h_report_path: str, json_dir: str, out_dir: str,
           do_step3: bool = False, ref_path: str = '') -> pd.DataFrame:
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

    # ---- 全局一对一匹配: 密集簇匈牙利最优 (scipy 缺失时跳过, 由贪心兜底) ----
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
            cost = np.full((n_rows, n_cols), 9.9)
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
            ri_idx, ci_idx = linear_sum_assignment(cost)
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
    # 第一阶段 v2: 初匹配结果先落到 1H-NMR-v2 (此时 v2 = v1)
    if '1H-NMR-v2' not in out.columns:
        out['1H-NMR-v2'] = out['1H-NMR-v1']
    else:
        m_missing = out['1H-NMR-v2'].isna() | (out['1H-NMR-v2'].astype(str).str.strip() == '')
        out.loc[m_missing, '1H-NMR-v2'] = out.loc[m_missing, '1H-NMR-v1']
    # 参考化合物标记 (如有): 参考归属重排统一放到 v2.5 之后 (v2.6 步),
    # 此处只建 no_key 供迭代复核的参考途径按 δc参考 有无判断.
    ref_groups_iter = _ref_groups_of(ref_path) if ref_path and os.path.isfile(ref_path) else None
    out['_no_key_iter'] = out['No.'].map(_norm_no_key) if 'No.' in out.columns else ''
    # 第二阶段: HSQC 矩阵复核 (在参考补 H 之后): 逐行验证 + 融合峰拆分修正,
    # 修正结果写新列 1H-NMR-v2.5; 未修正行 v2.5 = v2; δH_ppm/nH 随 v2.5 修正更新
    _centers_cache: list = []  # 拆分中心缓存 (供归属完整性迭代复核轮换候选)
    out = review_assignment(out, df, h1d_signals, json_dir, sig_df=sig_df,
                            centers_out=_centers_cache)
    if '1H-NMR-v2.5' not in out.columns:
        out['1H-NMR-v2.5'] = out['1H-NMR-v2']
    else:
        m_missing = out['1H-NMR-v2.5'].isna() | (out['1H-NMR-v2.5'].astype(str).str.strip() == '')
        out.loc[m_missing, '1H-NMR-v2.5'] = out.loc[m_missing, '1H-NMR-v2']
    # 最终归属基准列 = 1H-NMR-v2.5 (δH_ppm/nH 等参考-v2 补齐的列同步以 v2.5 为准:
    # v2.5 变化但 nH/δH_ppm 为空的行, 按自身 v2.5 的 δH 回查积分表补齐)
    def _final_fill(i):
        v = out.at[i, '1H-NMR-v2.5']
        hv = _ref_h_val(v)
        if hv is None:
            return
        if pd.isna(out.at[i, 'nH']) or str(out.at[i, 'nH']).strip() == '':
            s = min(h1d_signals, key=lambda x: abs(x['ppm'] - hv), default=None)
            if s and abs(s['ppm'] - hv) <= H_ALIGN_TOL and s.get('nH') is not None:
                out.at[i, 'nH'] = int(s['nH'])
        if 'δH_ppm' in out.columns and (pd.isna(out.at[i, 'δH_ppm'])
                                        or str(out.at[i, 'δH_ppm']).strip() == ''):
            if isinstance(sig_df, pd.DataFrame) and 'δH' in sig_df.columns:
                m_sig = sig_df[pd.to_numeric(sig_df['δH'], errors='coerce')
                               .sub(hv).abs().le(0.011)]
                if not m_sig.empty:
                    out.at[i, 'δH_ppm'] = str(m_sig.iloc[0]['δH_ppm'])
    for i in out.index:
        _final_fill(i)
    # ---- 归属完整性迭代复核 (内存中间 df, 输出 Excel 之前) ----
    # 结合 2D HSQC 数据反复匹配归属表中的 C-H 对, 使 1H_NMR 的 δH_nH 与 nH
    # 整体差异尽量小: 目标函数 = Σ(差异惩罚), 其中 nH==1 的行不允许差异 (权重极大),
    # nH 越大允许的差异越大 (惩罚随 nH 增大而衰减 → 权重 1/nH).
    # 不取消 H 容差限制意义上的盲匹配 — 候选 C-H 对仍须来自 2D HSQC 交叉峰中心
    # (物理证据), H 端无 ppm 限制 (取消 _REF_H_TOL), 匹配自由度交给整体优化.
    def _h_assign_state():
        """当前归属状态: [(slot_key, δH_ppm_str,assigned δH list)] 及全局差异评估."""
        sig_chk = check_h_counts(out, sig_df)
        m_bad = sig_chk['δH_nH'] != sig_chk['nH']
        # 差异惩罚: nH==1 行不允许差异 (权重 1000); 其余按 nH 越大容差越大 (权重 1/nH)
        n_bad = int(m_bad.sum())
        penalty = 0.0
        for _, sr in sig_chk[m_bad].iterrows():
            w = 1000.0 if int(sr['nH']) == 1 else 1.0 / max(int(sr['nH']), 1)
            penalty += abs(int(sr['δH_nH']) - int(sr['nH'])) * w
        return sig_chk, n_bad, penalty

    def _rewrite_rows(m_rows, s1d):
        """把归属表中 δH_ppm == m_rows 的行改归属到 1D 信号 s1d."""
        h_new = f"{s1d['ppm']:.2f}"
        dh_new = ''
        if isinstance(sig_df, pd.DataFrame) and 'δH' in sig_df.columns:
            m_sig = sig_df[pd.to_numeric(sig_df['δH'], errors='coerce')
                           .sub(s1d['ppm']).abs().le(0.011)]
            if not m_sig.empty:
                dh_new = str(m_sig.iloc[0]['δH_ppm'])
        for oi in out.index[m_rows]:
            out.at[oi, 'assigned_1H'] = h_new
            out.at[oi, '1H-NMR-v2.5'] = fmt_1h_nmr(h_new, s1d)
            if s1d.get('nH') is not None:
                out.at[oi, 'nH'] = int(s1d['nH'])
            if dh_new:
                out.at[oi, 'δH_ppm'] = dh_new

    # 预计算: 每个 HSQC 中心 → 对齐的 1D 信号 (2D 物理证据 + 1D 支撑)
    _center_to_sig = []
    for cc, hh, amp in _centers_cache:
        s1d = min(h1d_signals, key=lambda x: abs(x['ppm'] - hh), default=None)
        _center_to_sig.append((cc, hh, s1d if s1d and abs(s1d['ppm'] - hh) <= H_ALIGN_TOL else None))

    def _hsqc_h_profile(c_v: float):
        """2D HSQC 在碳维 c_v 处的氢维截面 (h_axis, profile) 或 None.

        截面 = 碳维 c_v ±2 像素的最大投影 (抗碳维抖动), 归一化到峰高 1.
        """
        if _hsqc_mat is None:
            return None
        ci = int(np.argmin(np.abs(_hsqc_c_ax - c_v)))
        lo, hi = max(ci - 2, 0), min(ci + 3, _hsqc_mat.shape[0])
        prof = _hsqc_mat[lo:hi, :].max(axis=0)
        pk = prof.max()
        if pk <= 0:
            return None
        return _hsqc_h_ax, prof / pk

    def _hsqc_h_peaks(c_v: float, mult: str):
        """按 1D 多重峰类型从 2D HSQC 截面提取 δH 候选值 (结合 review 与迭代复核).

        规则 (化学家读谱方式):
        - s (单峰) / t (三峰) / brs 等对称峰: 取截面峰顶 (局部极大) 的 H 数值;
        - d (双峰) / q (四重峰) / dd/dt/dq 等两主峰型: 找两个主峰顶, 返回二者
          **中间谷** 的 H 数值 (积分/二维中心落在双峰中线);
        - 其余 (m 等): 返回截面最大点 H 数值.
        返回候选 δH 列表 (降序按截面强度), 空截面返回 [].
        """
        prof = _hsqc_h_profile(c_v)
        if prof is None:
            return []
        h_ax, p = prof
        # 平滑后找局部极大 (峰顶) 与局部极小 (谷)
        w = 3
        if len(p) < w * 2 + 1:
            return [float(h_ax[int(np.argmax(p))])]
        ps = np.convolve(p, np.ones(w) / w, mode='same')
        win = sliding_window_view(ps, 3)
        imax = [(i + 1) for i in range(len(win)) if win[i][1] >= win[i][0]
                and win[i][1] >= win[i][2] and win[i][1] > 0.15]
        imin = [(i + 1) for i in range(len(win)) if win[i][1] <= win[i][0]
                and win[i][1] <= win[i][2]]
        imax = [i for i in imax if 0 < i < len(p) - 1]
        mult_l = (mult or '').lower()
        two_peak = any(k in mult_l for k in ('d', 'q', 'dd', 'dt', 'dq', 'td'))
        if two_peak and len(imax) >= 2:
            # 两个主峰 (强度最高) 之间的谷 → 双峰/四重峰中线 H 数值
            i1, i2 = sorted(sorted(imax, key=lambda i: -p[i])[:2])
            valleys = [i for i in imin if i1 < i < i2]
            gi = min(valleys, key=lambda i: ps[i]) if valleys else (i1 + i2) // 2
            return [float(h_ax[gi])]
        if imax:
            imax.sort(key=lambda i: -p[i])
            return [float(h_ax[i]) for i in imax[:3]]
        return [float(h_ax[int(np.argmax(p))])]

    _hsqc_loaded = load_hsqc_matrix(json_dir)
    if _hsqc_loaded is not None:
        _hsqc_mat, _hsqc_c_ax, _hsqc_h_ax = _hsqc_loaded
    else:
        _hsqc_mat = None
    _retry = 0
    _best_penalty = None
    while _retry < _MAX_RETRY:
        sig_chk, n_bad, penalty = _h_assign_state()
        if n_bad == 0:
            break
        if _best_penalty is None or penalty < _best_penalty:
            _best_penalty = penalty
        _retry += 1
        # ---- 整体优化: 对每个 δH 缺 H (δH_nH < nH) 的信号, 找 HSQC 中心落在
        # 其碳附近 (±REV_C_TOL) 的候选行 (优先 δH_nH > nH 的信号行, 即占用过多者),
        # 按 2D C-H 对把该行改归属到缺 H 信号 → C-H 对整体重匹配 ----
        # 参考化合物约束 (如有): 大部分 No. 两化合物归属一致, 因此参考列有值的行
        # 优先参与 C-H 对交换 (其参考 δH + ΔH整体 指明该位置的真实 δH),
        # 且参考 δH 指向的 1D 信号若未被占用, 可直接补到该行 (nH==1 行优先).
        m_over = sig_chk[sig_chk['δH_nH'] > sig_chk['nH']]   # 占用过多 → 有行可移出
        m_under = sig_chk[sig_chk['δH_nH'] < sig_chk['nH']]  # 缺 H → 接收
        if not m_under.empty:
            # 轮换偏移: 同一状态多次尝试时换不同候选, 避免死循环
            _shift = _retry - 1
            progressed = False
            # ---- 途径 A (参考化合物驱动): δc参考 有值 (C 对应成立) + 该行归属 H 与
            # 参考 δH + ΔH整体 偏离超过 _REF_H_TOL → 按 2D 截面/参考 δH 重归属.
            # 参考化合物大部分 No. 归属一致 → 参考列是可靠的先验证据.
            if ref_groups_iter is not None:
                # ΔH 整体偏差即时统计 (v2.6 之前的迭代用 δc参考 缺失时跳过该行)
                _d_hs_it = []
                for oi2 in out.index:
                    rv2 = _ref_h_val(out.at[oi2, 'δH参考']) \
                        if 'δH参考' in out.columns else None
                    pv2 = _ref_h_val(out.at[oi2, '1H-NMR-v2'])
                    if rv2 is not None and pv2 is not None:
                        _d_hs_it.append(pv2 - rv2)
                _d_h_med_iter = statistics.median(_d_hs_it) if _d_hs_it else 0.0
                h_tol_eff = _REF_H_TOL + abs(_d_h_med_iter)
                for oi in out.index:
                    if progressed:
                        break
                    rk = out.at[oi, '_no_key_iter']
                    if rk not in ref_groups_iter:
                        continue
                    # 参考 δH 从参考组数据取 (参考列在 v2.6 步才生成, 此处不依赖列);
                    # 该行按显示顺序对应组内第 n 个参考行 (与最终参考列分配规则一致)
                    ref_rows_g = ref_groups_iter[rk]
                    # 同 no_key 的目标行序号 (组内)
                    same = [j for j in out.index if out.at[j, '_no_key_iter'] == rk]
                    gi_n = same.index(oi)
                    ref_h = _ref_h_val(ref_rows_g[gi_n]['δH参考']) \
                        if gi_n < len(ref_rows_g) else None
                    if ref_h is None:
                        continue
                    pv = _ref_h_val(out.at[oi, 'assigned_1H'])
                    if pv is not None and abs(pv - (ref_h + _d_h_med_iter)) <= _REF_H_TOL:
                        continue  # 已与参考一致
                    # 该位置参考 δH 指向的 1D 信号 (未占用或本行已占用)
                    target_h = ref_h + _d_h_med_iter
                    s_ref = min(h1d_signals, key=lambda x: abs(x['ppm'] - target_h),
                                default=None)
                    if s_ref is None or abs(s_ref['ppm'] - target_h) > h_tol_eff:
                        continue
                    cur = str(out.at[oi, 'assigned_1H'])
                    if cur == f"{s_ref['ppm']:.2f}":
                        continue
                    # 2D 证据: 该行碳附近存在对齐到 s_ref 的中心 (nH==1 必须匹配正确)
                    c_v = float(out.at[oi, 'ppm'])
                    ok2d = any(t[2] is not None and t[2]['ppm'] == s_ref['ppm']
                               and abs(t[0] - c_v) <= REV_C_TOL
                               for t in _center_to_sig)
                    if int(out.at[oi, 'nH'] if not pd.isna(out.at[oi, 'nH']) else 0) == 1 \
                            and not ok2d:
                        continue  # nH==1 行必须 2D 支撑
                    m_one = pd.Series(False, index=out.index)
                    m_one[oi] = True
                    _rewrite_rows(m_one, s_ref)
                    progressed = True
            # ---- 途径 B (2D C-H 对整体重匹配): 缺 H 信号从占用过多信号移入行 ----
            for _, urow in m_under.iterrows():
                dh_u = float(urow['δH'])
                dh_u_str = str(urow['δH_ppm']).strip()
                # 缺 H 信号的碳位置: 归属表中已引用该信号的行 (或其参考 δc)
                m_ref = out[out['δH_ppm'].astype(str).str.strip() == dh_u_str]
                c_u = float(m_ref['ppm'].iloc[0]) if not m_ref.empty else None
                # 候选 1D 信号 (δH ± 无限制 — 取消 H 容差): 与缺 H 信号 δH 最近的
                # 若干个, 依次尝试把"占用过多"信号的行转归属过来 (必须经 2D 中心)
                # 缺 H 信号的 1D 多重峰类型 → 2D 截面提取的 δH 真值候选:
                # s/t → 峰顶; d/q → 双峰中间谷; 匹配值允许落在 δH_ppm 区间内
                # 任意位置 (不必等于区间中点), nH==1 的 C-H 对必须匹配正确
                u_type = ''
                if isinstance(sig_df, pd.DataFrame) and 'type' in sig_df.columns:
                    m_t = sig_df[sig_df['δH_ppm'].astype(str).str.strip() == dh_u_str]
                    u_type = str(m_t.iloc[0]['type']) if not m_t.empty else ''
                h_true_cands = _hsqc_h_peaks(c_u if c_u is not None else dh_u, u_type)
                # 缺 H 信号自身碳维处无截面 (碳归属错位) 时用其 2D 中心 δH 兜底
                if not h_true_cands and cands_u:
                    h_true_cands = [float(cands_u[_shift % len(cands_u)][1])]
                for _, orow in m_over.iterrows():
                    dh_o_str = str(orow['δH_ppm']).strip()
                    m_o_rows = out['δH_ppm'].astype(str).str.strip() == dh_o_str
                    if not m_o_rows.any():
                        continue
                    c_o = float(out[m_o_rows]['ppm'].iloc[0])
                    # 2D HSQC 证据: 占用过多信号与缺 H 信号的碳之间必须存在交叉峰
                    # 中心网络 (本碳中心各自的 s1d 互为邻近), 用中心轮换建立映射
                    cands_o = [t for t in _center_to_sig
                               if t[2] is not None and abs(t[0] - c_o) <= REV_C_TOL]
                    cands_u = [t for t in _center_to_sig
                               if t[2] is not None and c_u is not None
                               and abs(t[0] - c_u) <= REV_C_TOL]
                    if not cands_o or not cands_u:
                        continue
                    # 从占用过多信号中取一行 (轮换), 改归属到缺 H 信号 δH 附近
                    # 的 1D 信号 (2D 中心对齐支撑, 无 H ppm 限制)
                    oi = list(out.index[m_o_rows])[_shift % len(out.index[m_o_rows])]
                    # 目标 1D 信号: 优先 2D 截面 δH 真值候选对齐的信号 (J 类型规则),
                    # 无则退回缺 H 信号 δH 对齐的 2D 中心对应信号
                    sig_u = None
                    for ht in h_true_cands:
                        s_c = min(h1d_signals, key=lambda x: abs(x['ppm'] - ht),
                                  default=None)
                        if s_c and abs(s_c['ppm'] - ht) <= H_ALIGN_TOL:
                            sig_u = s_c
                            break
                    if sig_u is None and cands_u:
                        sig_u = cands_u[_shift % len(cands_u)][2]
                    if sig_u is None:
                        continue
                    if abs(sig_u['ppm'] - dh_u) > _REF_H_TOL + 1.0:
                        continue  # 防御: 目标信号应接近缺 H 信号的 δH (宽裕度)
                    # 单行改归属: 构造与 out 同索引的布尔掩码
                    m_one = pd.Series(False, index=out.index)
                    m_one[oi] = True
                    _rewrite_rows(m_one, sig_u)
                    progressed = True
                    break  # 每轮只移动一行, 下轮重新评估
                if progressed:
                    break
            if not progressed:
                break  # 无可改进的移动 → 收敛退出
    out = out.drop(columns='_no_key_iter', errors='ignore')
    sig_chk, n_bad_final, _ = _h_assign_state()
    if _retry:
        print(f'[检查] 归属完整性迭代复核: 尝试 {_retry} 次, '
              f'最终不匹配 {n_bad_final}/{len(sig_chk)} 条'
              + (' (无 2D 支撑的移动, 剩余行待人工复核)'
                 if n_bad_final and _retry < _MAX_RETRY else
                 (' (达上限, 剩余行待人工复核)' if n_bad_final else '')))
    # ---- 参考化合物归属 (v2.6, 在 v2.5 之后): 有参考表时执行, 无则跳过不输出 v2.6 ----
    # C 对应 (糖基/苷元分区, Δppm(C) 列) + 整体偏差校准 + 参考 H 对齐 1D 信号.
    # v2.6 = v2.5 基础上按参考化合物重归属 (v2 空的行参考补 H, v2 已有但偏离参考的行
    # 由后续迭代复核的参考途径修正); 参考列 (No.参考/δc参考/Δppm(C)/δH参考) 在此追加.
    _ver = '1H-NMR-v2.5'
    if ref_path and os.path.isfile(ref_path):
        out = fill_missing_h_from_ref(out, sig_df, h1d_signals, ref_path)
        out['1H-NMR-v2.6'] = out['1H-NMR-v2.5']
        m_missing = out['1H-NMR-v2.6'].isna() | \
            (out['1H-NMR-v2.6'].astype(str).str.strip().isin(['', 'nan']))
        out.loc[m_missing, '1H-NMR-v2.6'] = out.loc[m_missing, '1H-NMR-v2']

        def _rewrite_rows_v26(df_v, m_rows, s1d, ver_cols):
            """把 v2.6 (及 assigned_1H/nH/δH_ppm) 改归属到 1D 信号 s1d (仅 v2.6 生效)."""
            h_new = f"{s1d['ppm']:.2f}"
            dh_new = ''
            if isinstance(sig_df, pd.DataFrame) and 'δH' in sig_df.columns:
                m_sig = sig_df[pd.to_numeric(sig_df['δH'], errors='coerce')
                               .sub(s1d['ppm']).abs().le(0.011)]
                if not m_sig.empty:
                    dh_new = str(m_sig.iloc[0]['δH_ppm'])
            for oj in df_v.index[m_rows]:
                df_v.at[oj, '1H-NMR-v2.6'] = fmt_1h_nmr(h_new, s1d)
                df_v.at[oj, 'assigned_1H'] = h_new
                if s1d.get('nH') is not None:
                    df_v.at[oj, 'nH'] = int(s1d['nH'])
                if dh_new:
                    df_v.at[oj, 'δH_ppm'] = dh_new
        # ΔH 整体偏差 (参考 δH vs v2.6, C 对应行) — 仲裁目标 = 参考 δH + ΔH整体
        _d_hs_ref = []
        for oi2 in out.index:
            rh2 = _ref_h_val(out.at[oi2, 'δH参考'])
            pv2 = _ref_h_val(out.at[oi2, '1H-NMR-v2.6'])
            if rh2 is not None and pv2 is not None:
                _d_hs_ref.append(pv2 - rh2)
        _d_h_med_iter_ref = statistics.median(_d_hs_ref) if _d_hs_ref else 0.0
        # ---- 版本一致性仲裁 (v2.6 参考容差 _REF_MATCH_TOL) ----
        # 若该行 v1~v2.6 各版本 H 数值不完全一致: 取与参考 δH (经 ΔH整体校准) 差别
        # 最小的版本值; 差异 < _REF_MATCH_TOL → OK 直接采用该值; 否则在 1D 信号中
        # 找与参考 δH 数值差异更小的 H (±_REF_H_TOL 内), 替换 v2.6.
        def _v_val(v):
            return _ref_h_val(v)
        ver_cols = [c for c in ('1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v2.5', '1H-NMR-v2.6')
                    if c in out.columns]
        n_arb_ok = n_arb_fixed = 0
        for oi in out.index:
            rh = _ref_h_val(out.at[oi, 'δH参考'])
            if rh is None:
                continue
            target_h = rh + _d_h_med_iter_ref
            # 该行各版本的 H 数值
            vals = {}
            for vc in ver_cols:
                v = _v_val(out.at[oi, vc])
                if v is not None:
                    vals[vc] = v
            if not vals:
                continue
            # 取与参考差别最小的版本值
            best_vc = min(vals, key=lambda c: abs(vals[c] - target_h))
            best_v = vals[best_vc]
            d_best = abs(best_v - target_h)
            uniq = set(round(v, 2) for v in vals.values())
            if len(uniq) == 1:
                n_arb_ok += 1  # 各版本一致, 无需仲裁
                continue
            if d_best < _REF_MATCH_TOL:
                # 差异在容差内 → OK, 采用最接近参考的版本值
                s_c = min(h1d_signals, key=lambda x: abs(x['ppm'] - best_v), default=None)
                if s_c and abs(s_c['ppm'] - best_v) <= H_ALIGN_TOL:
                    m_one = pd.Series(False, index=out.index)
                    m_one[oi] = True
                    _rewrite_rows_v26(out, m_one, s_c, ver_cols)
                n_arb_ok += 1
            else:
                # 差异超容差 → 在 1D 信号中找与参考 δH 差异更小的 H (±_REF_H_TOL)
                s_c, d_c = None, 9e9
                for s in h1d_signals:
                    d = abs(s['ppm'] - target_h)
                    if d < d_best and d < _REF_H_TOL and d < d_c:
                        s_c, d_c = s, d
                if s_c is not None:
                    m_one = pd.Series(False, index=out.index)
                    m_one[oi] = True
                    _rewrite_rows_v26(out, m_one, s_c, ver_cols)
                    n_arb_fixed += 1
        print(f'[参考归属] v2.6 仲裁: 版本一致 {n_arb_ok} 行, '
              f'按参考容差(_REF_MATCH_TOL={_REF_MATCH_TOL})修正 {n_arb_fixed} 行')
        _ver = '1H-NMR-v2.6'
        print(f'[参考归属] v2.6 完成: 参考列 + Δppm(C) 已追加, 未归属行经参考补 H')
    # 第三步 (COSY/HMBC 全表复核) 不改写 H — v3 = 基准列; 复核证据在第三 sheet
    # 开关 do_step3 默认关闭: 不执行第三步, 最终归属以基准列 (v2.6/v2.5) 为准;
    # 原依赖 v3 的列 (row_match / 1H-NMR_real 对齐基准等) 改用基准列.
    if do_step3:
        out['1H-NMR-v3'] = out[_ver]
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
    for col in ('1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v2.5', '1H-NMR-v2.6', '1H-NMR-v3',
                'all_assigned_H'):
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
            # 1H-NMR_real 插在 1H-NMR-v2.5 后 (v2.5 恒存在, 是最终归属基准列)
            idx = col_order.index('1H-NMR-v2.5') + 1
            col_order.insert(idx, '1H-NMR_real')
            out = out[col_order]
    else:
        out = out.drop(columns=[c for c in ('1H-NMR_real', 'row_match') if c in out.columns])
    # ---- 列序: ppm, Types, No., [No.参考, δc参考, δH参考], 1H-NMR-v1, 1H-NMR-v2, 1H-NMR-v2.5, [1H-NMR_real], [1H-NMR-v3], δH_ppm, nH, review, judge ----
    front = ['ppm', 'Types', 'No.', 'No.参考', 'δc参考', 'Δppm(C)', 'δH参考',
             '1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v2.5', '1H-NMR-v2.6', '1H-NMR_real', '1H-NMR-v3',
             'δH_ppm', 'nH']
    rest = [c for c in out.columns if c not in front]
    ordered = [c for c in front if c in out.columns]
    for c in rest:
        if c in ('review', 'judge', 'C_match', 'row_match'):
            ordered.append(c)
    ordered += [c for c in rest if c not in ordered]
    out = out[ordered]
    # ---- 输出前裁剪: assigned_1H / aligned_1D 为中间过程列, 不输出;
    #      judge 为对比判分列, 无真实归属表 (无 truth) 时也不输出 ----
    out = out.drop(columns=['assigned_1H', 'aligned_1D'], errors='ignore')
    if not truth['1H-NMR'].notna().any():
        out = out.drop(columns=['judge'], errors='ignore')
    # ---- 参考列重排与补全: 参考列在补 H 时按"补 H 配对序"写入, 但最终按 δH 降序排序后
    #      行序可能变化 (如 No.15 补全行 1.22 排到未补行前, 而其参考 δH 1.47 < 2.14),
    #      导致 δH参考 列在组内不是从大到小. 泛化修复: 组内按当前显示行序重新分配
    #      降序参考行 (参考行已按 δH 降序), 保证 δH参考 列从大到小.
    #      补全要求: 有参考化合物时, 参考组存在的 No. 其组内所有行 (含 CH2 第二行)
    #      参考三列都必须有值 — 参考行不足时复用该组最后一个参考行 (同位置 CH2
    #      第二 H 的 No./δc 相同, δH参考 也沿用同位置数据).
    #      注意: 重排以 C 对应 (fill_missing_h_from_ref 写入的参考列) 为准, 仅对
    #      参考列有值的行按组内显示顺序重排参考行; 同 no_key 的糖基异头碳与苷元
    #      碳 no_key 相同但分区不同, 盲按 no_key 覆盖会串区 (clxh141E 苷元 C3 被覆
    #      盖成糖基异头碳参考), 因此未匹配到参考的行保持 C 对应结果 (空).
    if 'δH参考' in out.columns and ref_path and os.path.isfile(ref_path):
        ref_groups_all = _ref_groups_of(ref_path)
        out['_no_key2'] = out['No.'].map(_norm_no_key) if 'No.' in out.columns else ''
        out['_has_ref2'] = out['δc参考'].astype(str).str.strip().ne('') & \
            out['δc参考'].astype(str).str.strip().ne('nan')
        # 按 (no_key, 原 δc参考) 分组 — 同 no_key 不同 δc参考 属于不同区 (糖基异头
        # vs 苷元同号碳), 不互相串; 组内按显示顺序重排降序参考行
        for (key, dc_ref), gi in out[out['_has_ref2']].groupby(
                ['_no_key2', 'δc参考'], sort=False):
            ref_rows = ref_groups_all.get(key, [])
            if not ref_rows:
                continue
            # 与本组 δc参考 一致的参考行 (异头碳连接位组可能含多条), 按显示顺序分配
            rows_matched = [rr for rr in ref_rows
                            if str(rr['δc参考']).strip() == str(dc_ref).strip()] or ref_rows
            for j, oi in enumerate(gi.index):
                rr = rows_matched[min(j, len(rows_matched) - 1)]
                out.at[oi, 'No.参考'] = rr['No.参考']
                out.at[oi, 'δc参考'] = rr['δc参考']
                out.at[oi, 'δH参考'] = rr['δH参考']
        # CH2 第二行参考补全: 同 (no_key, ppm, Types) 组内首行有参考而次行无 →
        # 复用首行参考 (同位置第二 H, No./δc/δH参考 相同)
        for (key, _p, _t), gi in out.groupby(['_no_key2', 'ppm', 'Types'], sort=False):
            rows = list(gi.index)
            first_ref = out.at[rows[0], 'No.参考'] if rows else ''
            if len(rows) > 1 and str(first_ref).strip() not in ('', 'nan'):
                for oi in rows[1:]:
                    if str(out.at[oi, 'No.参考']).strip() in ('', 'nan'):
                        out.at[oi, 'No.参考'] = out.at[rows[0], 'No.参考']
                        out.at[oi, 'δc参考'] = out.at[rows[0], 'δc参考']
                        out.at[oi, 'δH参考'] = out.at[rows[0], 'δH参考']
        out = out.drop(columns=['_no_key2', '_has_ref2'])
    # ---- 输出两个归属表 sheet: 补全版 (ppm/Types 完整) + 置空版 (组内非首行置空) ----
    # 置空规则泛化: 最终排序后, 同 (ppm, Types) 组内除第一行外全部置空 ppm/Types
    # (化学家续行格式: 碳标签只在组内首行). 不依赖 _cont 标记 — 补全行可能因 δH
    # 降序排到组内首行, 此时首行保留标签、原行置空, 与行序自然一致.
    # ---- No./Pos. 列撇号归一 (输出前, 所有 sheet 生效): 单引号→' ; 双引号→''.
    # 必须在置空版 out_blank 复制之前执行, 否则归属表 sheet 仍是全角撇号.
    def _pos_col(df_v):
        return next((c for c in df_v.columns
                     if str(c).strip().rstrip('.').lower() in ('no', 'pos')), None)
    _pos = _pos_col(out)
    # No./Pos. 列保留原始字符串 (含全角′″), 不做归一替换;
    # 归一 (_norm_prime) 仅在识别糖基 (_norm_no_key) 时使用.
    out_blank = out.copy()
    if '_cont' in out_blank.columns:
        m_blank = pd.Series(False, index=out_blank.index)
        for (_, _), gi in out_blank.groupby(['ppm', 'Types']).groups.items():
            rows = list(gi)
            m_blank.loc[rows[1:]] = True
        out_blank['ppm'] = out_blank['ppm'].astype(object)
        out_blank['Types'] = out_blank['Types'].astype(object)
        out_blank.loc[m_blank, ['ppm', 'Types']] = ''
        out_blank = out_blank.drop(columns='_cont')
    out_full = out.drop(columns='_cont', errors='ignore')
    # ppm (碳) 列固定 2 位小数显示; 判定列文本格式 (防 '-' 开头误转公式)
    text_cols = [c for c in ('No.', 'Pos.', '1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v2.5', '1H-NMR-v2.6', '1H-NMR-v3', '1H-NMR_real',
                             'review', 'C_match', 'row_match', 'all_assigned_H', 'Types') if c in out.columns]
    # 双 sheet: 归属表 + 1H_NMR 积分表 (sheet 名沿用文件名前缀)
    # J 列落盘还原为 '(a, b)' 文本 (内存中是 list)
    sig_out = sig_df.copy()
    # 新增 1H-NMR 列: 与归属表 1H-NMR 列同一格式规则 (fmt_1h_nmr):
    # 'δ, mult (J...)'; J 值 >=3 个 (tdd/ddd 等) 逗号后直接写 m
    if 'δH' in sig_out.columns and 'type' in sig_out.columns:
        sig_out['1H-NMR'] = sig_out.apply(
            lambda r: fmt_1h_nmr(f"{float(r['δH']):.2f}",
                                 {'mult': r['type'],
                                  'J': r['J'] if isinstance(r.get('J'), list) else []}),
            axis=1)
    if 'J' in sig_out.columns:
        sig_out['J'] = sig_out['J'].map(
            lambda js: f"({', '.join(js)})" if isinstance(js, list) and js else
                       (str(js) if not isinstance(js, list) else ''))
    # 1H_NMR sheet 列序固定: δH, type, J, 1H-NMR, δH_ppm, nH
    sig_cols = [c for c in ('δH', 'type', 'J', '1H-NMR', 'δH_ppm', 'nH')
                if c in sig_out.columns]
    sig_out = sig_out[sig_cols]
    # ---- 检查复核 (输出 Excel 之前): 归属完整性检查 ----
    # 统计归属表各行 δH_ppm 占用 H 数 (同 (ppm, Types) 槽仅 1 行且 CH3 → 权重 3,
    # 其余每行 1), 汇总到 1H_NMR 表的 δH_nH 列, 与 nH 向量化比较,
    # 不相等行整行标红 (与 数据处理程序/check_h_assignment.py 同规则)
    _slot_rows: dict = {}
    _counts = Counter()
    for _, _r in out.iterrows():
        _s = '' if pd.isna(_r.get('δH_ppm')) else str(_r['δH_ppm']).strip()
        if _s and _s.lower() != 'nan':
            _k = (round(float(_r['ppm']), 2), str(_r['Types']).strip())
            _slot_rows[_k] = _slot_rows.get(_k, 0) + 1
    for _, _r in out.iterrows():
        _s = '' if pd.isna(_r.get('δH_ppm')) else str(_r['δH_ppm']).strip()
        if _s and _s.lower() != 'nan':
            _typ = str(_r['Types']).strip()
            _k = (round(float(_r['ppm']), 2), _typ)
            _w = 3 if (_typ == 'CH3' and _slot_rows.get(_k, 0) == 1) else 1
            _counts[_s] += _w
    sig_out['δH_nH'] = sig_out['δH_ppm'].map(lambda s: _counts.get(str(s).strip(), 0))
    m_bad = sig_out['δH_nH'] != sig_out['nH']
    n_bad_h = int(m_bad.sum())
    print(f'[检查] 归属完整性: 信号 {len(sig_out)} 条, 完整 {len(sig_out) - n_bad_h} 条, '
          f'不匹配 {n_bad_h} 条 (1H_NMR sheet 标红)')
    # ---- 糖基/苷元拆分 sheet: No. 或 Pos. (空间位置, '.'可有可无) 列存在时,
    # No. 末尾带 ' 的糖基部分行单独输出一个 sheet, 不带 ' 的苷元行单独输出一个 sheet.
    # No.列本身可能被置空版置空 (组内非首行), 此处用补全版 out 的 No. 判断归属哪部分.
    def _pos_col(df_v):
        return next((c for c in df_v.columns
                     if str(c).strip().rstrip('.').lower() in ('no', 'pos')), None)
    _pos = _pos_col(out)
    if _pos is not None:
        # No. 保留原始字符串; 识别糖基 = 归一化原始字符串中含 ' (6'a / 6'β 等
        # CH2 构型后缀、3-O-Glc-1' 等连接位写法, 任何含撇号的均属糖基)
        _has_prime = out[_pos].map(
            lambda v: not pd.isna(v) and "'" in _norm_prime(v))
        gly_df = out[_has_prime].copy()
        agly_df = out[~_has_prime].copy()
    else:
        gly_df = agly_df = None
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as xw:
        out_blank.to_excel(xw, sheet_name=f'{base}_归属表', index=False)
        out_full.to_excel(xw, sheet_name=f'{base}_归属表补全', index=False)
        sig_out.to_excel(xw, sheet_name=f'{base}_1H_NMR', index=False)
        # 糖基 / 苷元 拆分 sheet (追加在 1H_NMR 之后)
        if gly_df is not None:
            gly_df.to_excel(xw, sheet_name=f'{base}_糖基', index=False)
            agly_df.to_excel(xw, sheet_name=f'{base}_苷元', index=False)
            print(f'[拆分] 糖基 {len(gly_df)} 行 (No. 带 \'), 苷元 {len(agly_df)} 行 '
                  f'→ sheet {base}_糖基 / {base}_苷元')
        ws = xw.book.worksheets[2]
        red = Font(color='FF0000')
        # 糖基/苷元 sheet 文本列格式 (与归属表 sheet 同规则)
        for si_extra, df_e in ((3, gly_df), (4, agly_df)):
            if df_e is None:
                continue
            wse = xw.book.worksheets[si_extra]
            if 'ppm' in df_e.columns:
                ci = list(df_e.columns).index('ppm') + 1
                for row in wse.iter_rows(min_row=2, min_col=ci, max_col=ci):
                    for cell in row:
                        cell.number_format = '0.00'
            for col_name in text_cols:
                if col_name in df_e.columns:
                    ci = list(df_e.columns).index(col_name) + 1
                    for row in wse.iter_rows(min_row=2, min_col=ci, max_col=ci):
                        for cell in row:
                            cell.number_format = '@'
                            cell.alignment = Alignment(horizontal='left')
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
        # 1H_NMR sheet: δH_nH != nH 的行整行红字
        if n_bad_h:
            for r_i in range(len(sig_out)):
                if bool(m_bad.iloc[r_i]):
                    for c_i in range(1, len(sig_out.columns) + 1):
                        ws.cell(row=r_i + 2, column=c_i).font = red
        # 两个归属表 sheet 的格式 (补全版 sheet0 + 置空版 sheet1)
        col_v = [c for c in ('1H-NMR-v1', '1H-NMR-v2', '1H-NMR-v2.5', '1H-NMR-v2.6', '1H-NMR-v3') if c in out.columns]
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
        # 每个 sheet 视图固定首行表头 (上下滑动时表头不动); 隔行黄色背景便于阅读
        for wsi in xw.book.worksheets:
            wsi.freeze_panes = 'A2'
            yellow = PatternFill(start_color='FFFF80', end_color='FFFF80',
                                 fill_type='solid')
            for r_i, row in enumerate(wsi.iter_rows(min_row=2), start=2):
                if (r_i % 2) == 0:  # 数据区隔行 (第 2/4/6... 行) 填黄
                    for cell in row:
                        cell.fill = yellow
    print(f'[判分] {len(cmp_df)} 碳槽: OK={n_ok} 非OK={n_bad}')
    print(f'输出: {out_xlsx}')
    return cmp_df


def main():
    # 可选命令行参数: <化合物名|xlsx路径> [json_dir] [--step3] [--ref <参考xlsx>]
    # 开关默认值取顶部配置 (_DEFAULT_DO_STEP3 / _DEFAULT_REF), 命令行可覆盖;
    # --ref 未指定时默认用 <_DATA_DIR>/<化合物名>_H参考.xlsx (不存在则不启用参考补 H)
    do_step3 = _DEFAULT_DO_STEP3 or '--step3' in sys.argv
    sys.argv = [a for a in sys.argv if a != '--step3']
    ref_path = _DEFAULT_REF()
    if '--ref' in sys.argv:
        ri = sys.argv.index('--ref')
        if ri + 1 < len(sys.argv):
            ref_path = sys.argv[ri + 1]
            del sys.argv[ri:ri + 2]
        else:
            raise SystemExit('ERROR: --ref 需要参考 xlsx 路径参数')
    compound = _DEFAULT_COMPOUND
    if len(sys.argv) >= 2:
        compound = sys.argv[1]
    if os.path.isfile(compound) and compound.endswith('.xlsx'):
        xlsx_path = compound
        base = re.sub(r'数据归属\.xlsx$', '', os.path.basename(compound))
        data_dir = os.path.dirname(xlsx_path)
    else:
        base = compound
        # 传入化合物名时 数据目录 / 参考表随化合物名联动 (文件名兼容有无下划线)
        data_dir = rf'Y:\LXH数据处理\新化合物处理\程序解析\{base}解析'
        xlsx_path = _find_data_file(data_dir, base, '数据归属')
        if not xlsx_path:
            xlsx_path = os.path.join(data_dir, f'{base}数据归属.xlsx')
        if not ref_path:
            ref_path = _find_data_file(
                data_dir, base, '_H参考', extra_variants=('H参考', '_参考', '参考'))
    h_report = _find_data_file(data_dir, base, '氢谱积分', ext='.txt')
    if not h_report:
        h_report = os.path.join(data_dir, f'{base}氢谱积分.txt')
    if len(sys.argv) >= 3:
        json_dir = sys.argv[2]
    else:
        json_dir = _export_mnova_json(base, data_dir)
    out_dir = os.path.dirname(os.path.abspath(xlsx_path))
    # if not os.path.isfile(xlsx_path):
    #     raise SystemExit(f'ERROR: 归属表不存在: {xlsx_path}')
    if not os.path.isfile(h_report):
        raise SystemExit(f'ERROR: 氢谱积分.txt 不存在: {h_report}')
    print(f'[输入] xlsx={xlsx_path}\n[输入] 积分={h_report}\n[输入] JSON={json_dir}')
    cmp_df = assign(xlsx_path, h_report, json_dir, out_dir,
                    do_step3=do_step3, ref_path=ref_path)
    # print(cmp_df[['ppm', 'Types', 'assigned_H', 'truth_H', 'judge']].to_string(index=False))


if __name__ == '__main__':
    main()
