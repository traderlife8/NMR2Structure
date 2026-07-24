"""CLI 入口: 端到端运行 NMR 反推结构.

用法:
    python -m nmr_elucidator --folder <原始数据/<sample>>
                             [--db <数据库目录>] [--formula <分子式>]
                             [--top-n 10] [--retrieve] [--no-derive] [--no-2d]

默认行为:
    - 启用 --derive (从谱图直接推导归属表) 和 --use-2d (自动检测 2D 数据)
    - 默认仅推导, 不进入候选检索/打分 (加 --retrieve 继续检索)
    - 不强制 --formula (适用于未知化合物)
    - --db 默认为 NMR2Structure 外层平级 data/ 目录
    - 程序内部强制 UTF-8 编码, 无需 python -X utf8
"""
from __future__ import annotations

import argparse
import io
import os
import re
import sys
from pathlib import Path

# 程序内部强制 UTF-8 编码 (替代 python -X utf8 命令行标志)
# Windows 默认编码为 GBK/cp936, 不配置则中文输出可能乱码
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (AttributeError, OSError):
        pass  # 非标准 stream (如 IDE 重定向) 时静默跳过

from .constraints import Constraints
from .bridge_nmr2struct import convert_bruker_folder_to_nmr2struct
from .drawing import (
    draw_real_spectrum_overlay,
    draw_single_mol,
    draw_topn_grid,
    draw_unmatched_overlay,
    write_unmatched_summary_csv,
)
from .evidence import build_evidence, score_candidates
from .io_bruker import read_bruker_folder
from .io_processed import ProcessedReadError, read_processed_folder
from .nmrsim_check import check_candidate_with_nmrsim
from .report import write_report_html, write_report_json, write_report_markdown
from .retrieve import ProblemPeaks, retrieve_candidates
from .standard_assignment import (
    discover_known_compound as _discover_known_compound,
    read_standard_assignment_table as _read_standard_assignment_table,
    compare_program_to_standard as _compare_program_to_standard,
    write_joint_xlsx as _write_joint_xlsx,
    write_comparison_report as _write_comparison_report,
)

# 默认数据库目录: NMR2Structure 包外层平级的 data 文件夹
# 即 <project_root>/data (与 NMR2Structure/ 同级)
_DEFAULT_DB_DIR = str(Path(__file__).resolve().parent.parent.parent / 'data')

# 默认已知化合物根目录: E:\\NMR2Structure\\真实结构
# 当 <real_structure_dir>/<compound_name>/ 下存在 *数据归属*.xlsx 时,
# 走"已知化合物流程" (derive 不用 SMILES, 复核阶段才用 SMILES 做结构判等).
_DEFAULT_REAL_STRUCTURE_DIR = r'E:\NMR2Structure\真实结构'

# 项目 logger (在需要时调用 setup_logging() 初始化)
from .log_setup import setup_logging as _setup_logging
from .log_setup import get_logger as _get_logger


def _discover_databases(db_path: str) -> dict:
    """从 --db 参数自动发现可用数据库.

    支持两种模式:
      - 目录: 自动扫描 nmrshiftdb2/*.sd 和 lucy-ng-derep_db/*.db
      - 文件: 按扩展名判断类型 (.sd/.sdf → nmrshiftdb2, .db → lucy-ng)

    Returns:
        {'sdf': str|None, 'lucy_ng': str|None}
    """
    result = {'sdf': None, 'lucy_ng': None}
    db_path = os.path.abspath(db_path)

    if os.path.isdir(db_path):
        # 目录模式: 自动扫描已知数据库布局
        # 1. nmrshiftdb2 SDF
        sdf_candidates = [
            os.path.join(db_path, 'nmrshiftdb2', 'nmrshiftdb2withsignals.sd'),
        ]
        # 也扫描根目录下的 .sd 文件
        try:
            for f in os.listdir(db_path):
                if f.endswith(('.sd', '.sdf')):
                    sdf_candidates.append(os.path.join(db_path, f))
        except OSError:
            pass
        for p in sdf_candidates:
            if os.path.isfile(p):
                result['sdf'] = p
                break
        # 2. lucy-ng SQLite
        lucy_candidates = [
            os.path.join(db_path, 'lucy-ng-derep_db', 'lucy-ng-derep.db'),
            os.path.join(db_path, 'lucy-ng-derep.db'),
        ]
        try:
            for f in os.listdir(db_path):
                if f.endswith('.db'):
                    lucy_candidates.append(os.path.join(db_path, f))
        except OSError:
            pass
        for p in lucy_candidates:
            if os.path.isfile(p):
                result['lucy_ng'] = p
                break
    elif os.path.isfile(db_path):
        # 文件模式: 按扩展名判断类型
        ext = os.path.splitext(db_path)[1].lower()
        if ext in ('.sd', '.sdf'):
            result['sdf'] = db_path
        elif ext == '.db':
            result['lucy_ng'] = db_path
        else:
            # 未知扩展名: 尝试当 SDF
            result['sdf'] = db_path

    return result


def _cleanup_embedded_pngs(out_prefix: str, top_n: int, evidences) -> None:
    """HTML 已内嵌所有图片为 base64, 删除输出目录中的独立 PNG 文件.

    保留输出目录简洁: 只留 .json / .md / .html / .csv.
    删除失败时静默跳过 (不影响主流程).
    """
    import glob as _glob
    patterns = [
        out_prefix + '_topn.png',
        out_prefix + '_overlay.png',
        out_prefix + '_real_spectrum.png',
        out_prefix + '_heatmap.png',
        out_prefix + '_scatter.png',
        out_prefix + '_matrix.png',
        out_prefix + '_overview.png',
    ]
    # 单分子图 (含历史残留的 rank PNG)
    patterns.extend(_glob.glob(out_prefix + '_rank*.png'))
    for p in patterns:
        try:
            if os.path.isfile(p):
                os.remove(p)
        except OSError:
            pass


def _filter_peaks_by_intensity(shifts, intensities, percentile=0.20):
    """保留强度 ≥ 第 (1-percentile) 百分位的峰 (去最弱噪声).

    shifts 和 intensities 等长. 如果峰数 ≤ 5 则不过滤.

    语义修正 (Part 3.A): 原实现取前 percentile 比例的强峰, 把 41 峰 13C 砍到 8 峰,
    丢失 sp3 区弱峰. 现在按 docstring 描述"去底部 percentile 比例的弱峰":
      percentile=0.20 → 去掉最弱的 20% (从排序数组尾部取 cut 处的强度作阈值).
    """
    if not shifts or not intensities or len(shifts) <= 5:
        return list(shifts)
    paired = sorted(zip(intensities, shifts), reverse=True)  # 强→弱
    n = len(paired)
    cut = max(1, int(n * percentile))
    # 取第 cut 个最弱峰的强度作为阈值 (>= 该阈值即保留)
    thresh = paired[-cut][0]
    return sorted(s for i, s in paired if i >= thresh)


def _build_problem_from_bruker(specs, keep_solvent: bool = False) -> ProblemPeaks:
    """从读到的 BrukerSpectrum 列表构造 ProblemPeaks.

    选择规则：
      - 第一个 1H 实验作为 1H 峰表
      - 第一个 13C 实验（非 DEPT）作为 13C 峰表
      - 溶剂取第一个非空
    """
    h_spec = next((s for s in specs if s.nucleus == '1H'), None)
    # 碳谱：优先非 DEPT
    c_spec = next((s for s in specs if s.nucleus == '13C' and 'DEPT' not in s.nucleus), None)
    if c_spec is None:
        c_spec = next((s for s in specs if s.nucleus == '13C' or s.nucleus == '13C-DEPT'), None)

    c_shifts = [p.shift for p in c_spec.peaks] if c_spec else []
    h_shifts = [p.shift for p in h_spec.peaks] if h_spec else []
    h_integrals = []
    if h_spec:
        for region in h_spec.integrals:
            center = (region.low + region.high) / 2.0
            h_integrals.append((center, region.value))

    # 强度过滤: 13C 仅去最弱 10% (Bruker peaklist 已筛选); 1H 去最弱 65% (含 t1 噪声/旋转边峰)
    # Part 3.A: 修复 percentile 语义后, 采用绝对按比例过滤, 不再一刀切.
    c_shifts = _filter_peaks_by_intensity(
        c_shifts, [p.intensity for p in c_spec.peaks], percentile=0.10,
    ) if c_spec else []
    h_shifts = _filter_peaks_by_intensity(
        h_shifts, [p.intensity for p in h_spec.peaks], percentile=0.65,
    ) if h_spec else []

    solvent = None
    for s in specs:
        if s.solvent:
            solvent = s.solvent
            break

    # Part 2.B: 溶剂残余峰过滤 (默认启用, 可通过 --keep-solvent 跳过)
    # 1D 与候选检索可见一个明显的 CDCl3 残余峰 7.26 ppm 时, 它会污染检索.
    # 这里采用过滤策略: 检测到可用溶剂名时, 简单调用溶剂表过滤 proton/c shifts.
    n_h_filtered = 0
    n_c_filtered = 0
    if solvent and not keep_solvent:
        from .solvent_filter import is_solvent_peak as _is_sol_peak
        # 先统计原始谱中的溶剂残余峰 (用于 meta)
        if h_spec:
            n_h_filtered = sum(
                1 for p in h_spec.peaks if _is_sol_peak(p.shift, '1H', solvent)
            )
        if c_spec:
            n_c_filtered = sum(
                1 for p in c_spec.peaks if _is_sol_peak(p.shift, '13C', solvent)
            )
        h_shifts = [v for v in h_shifts if not _is_sol_peak(v, '1H', solvent)]
        c_shifts = [v for v in c_shifts if not _is_sol_peak(v, '13C', solvent)]

    return ProblemPeaks(
        carbon_shifts=c_shifts,
        proton_shifts=h_shifts,
        proton_integrals=h_integrals,
        solvent=solvent,
        field_mhz=h_spec.field_mhz if h_spec else (c_spec.field_mhz if c_spec else None),
        n_h_solvent_filtered=n_h_filtered,
        n_c_solvent_filtered=n_c_filtered,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='nmr_elucidator',
        description='基于 nmrshiftdb2 数据库 + nmrglue Bruker 解析的 NMR 反推结构工具',
    )
    parser.add_argument('--folder', required=True,
                        help='原始数据样本目录, e.g. NMR/clxh3')
    parser.add_argument('--db', default=_DEFAULT_DB_DIR,
                        help='数据库路径: 可传入总目录 (自动发现 nmrshiftdb2 + lucy-ng), '
                             '或单个数据库文件路径 (.sd/.db). '
                             '默认: NMR2Structure 外层平级 data/ 目录')
    parser.add_argument('--formula', default=None,
                        help='可选：硬约束分子式, e.g. C8H8O2. '
                             '默认不传 (适用于未知化合物推导)')
    parser.add_argument('--top-n', type=int, default=10,
                        help='输出 Top-N 候选')
    parser.add_argument('--max-hits', type=int, default=None,
                        help='限制候选命中数（用于快速验证）')
    parser.add_argument('--c-tol', type=float, default=2.0,
                        help='13C 化学位移容差 (ppm)')
    parser.add_argument('--h-tol', type=float, default=0.3,
                        help='1H 化学位移容差 (ppm)')
    parser.add_argument('--output-root', default=r'E:\NMR2Structure\output',
                        help='输出根目录（每个化合物分子目录）')
    parser.add_argument('--compound-name', default=None,
                        help='化合物名（用作子目录名）。默认从 --folder 路径末段推断，'
                             '例: .../clxh3 推断为 clxh3')
    parser.add_argument('--use-nmrsim', action='store_true',
                        help='启用 nmrsim 反向校验 (cosine 相似度)')
    parser.add_argument('--validate-with-nmrpy', action='store_true',
                        help='用 nmrpy 从 fid 时域信号二次处理 (FT/相位/基线), '
                             '与 peaklist.xml 交叉验证, 产出 *_nmrpy_validation.{json,md,png}')
    parser.add_argument('--with-real-spectrum', action='store_true',
                        help='读取 pdata/1/1r 真实处理后谱, 绘制真实谱线形图 (叠加 peaklist 峰位), '
                             '产出 *_real_spectrum.png')
    parser.add_argument('--use-hose', action='store_true',
                        help='启用 HOSE code 13C 预测验证 (从 nmrshiftdb2 自建查找表)')
    parser.add_argument('--use-pubchem', action='store_true', default=None,
                        help='强制启用 PubChem 扩展检索 (即使本地库有命中也检索). '
                             '默认: 本地库无命中时自动触发')
    parser.add_argument('--no-pubchem', action='store_true',
                        help='禁用 PubChem 扩展检索 (即使本地库无命中也不联网)')
    parser.add_argument('--use-coconut', action='store_true', default=None,
                        help='强制启用 COCONUT (开放天然产物库, 69.5万+ NP 结构) 联网检索. '
                             '默认: 本地库无命中时自动触发 (天然产物类样本优先)')
    parser.add_argument('--no-coconut', action='store_true',
                        help='禁用 COCONUT 检索 (即使本地库无命中也不联网)')
    parser.add_argument('--use-fragment-assembly', action='store_true',
                        help='启用片段组装候选生成 (库外新骨架)')
    parser.add_argument('--use-2d', action='store_true', default=True,
                        help='启用 2D NMR 关联验证 (HSQC/HMBC/COSY); '
                             '无 2D 数据时自动回退到 1D 流程 (默认开启)')
    parser.add_argument('--no-2d', action='store_true', default=False,
                        help='禁用 2D NMR 检测, 强制使用 1D 流程')
    parser.add_argument('--derive', action='store_true', default=True,
                        help='从谱图数据直接推导归属表 '
                             '(不依赖候选检索/打分, 仅用 HSQC/COSY/HMBC/ROESY; 默认开启)')
    parser.add_argument('--no-derive', action='store_true', default=False,
                        help='禁用 derive 模式, 直接走候选检索+打分流程')
    parser.add_argument('--retrieve', action='store_true', default=False,
                        help='derive 完成后继续执行候选检索+打分 (默认仅推导不检索)')
    parser.add_argument('--determine', action='store_true',
                        help='输出结构推断结论 (含置信度分级)')
    parser.add_argument('--hose-cache', default=None,
                        help='HOSE 查找表缓存文件路径 (避免每次重建)')
    parser.add_argument('--hose-max-records', type=int, default=50000,
                        help='从 SDF 构建 HOSE 查找表的最大记录数 (默认 50000; '
                             '0=处理全部; 10000~50000 够用, 100000+ 需约 5 分钟)')
    parser.add_argument('--lucy-ng-db', default=None,
                        help='lucy-ng 数据库路径覆盖 (优先于 --db 自动发现). '
                             '默认从 --db 目录自动查找 lucy-ng-derep.db')
    # lucy-ng 检索参数 (当发现 lucy-ng 数据库时自动启用)
    parser.add_argument('--lucy-ng-carbon-tol', type=int, default=2,
                        help='lucy-ng 检索时碳数容差 (±, 默认 2)')
    parser.add_argument('--lucy-ng-min-matches', type=int, default=10,
                        help='lucy-ng 检索最少 13C 匹配数 (默认 10)')
    parser.add_argument('--lucy-ng-max-results', type=int, default=500,
                        help='lucy-ng 检索最大返回数 (默认 500)')
    # Part 3 Phase 9: 化学家元信息接口
    parser.add_argument('--hresims', default='',
                        help='HRESIMS 数据, 例: "[M+Na]+ 647.3407". 会写入 Excel 顶部元信息行.')
    parser.add_argument('--optical-rotation', default='',
                        help='比旋光度 [α], 例: "(calcd. for C33H52O11Na+, 647.3402)".')
    parser.add_argument('--ir', default='',
                        help='IR/DBE 数据, 例: "DBE=8".')
    parser.add_argument('--compound-class', default='',
                        help='化合物分类提示, 例: "diterpenoid glucoside". 仅作为额外 meta 字段输出.')
    parser.add_argument('--quiet', action='store_true',
                        help='不打印检索进度')
    parser.add_argument('--log-dir', default=None,
                        help='日志输出目录 (默认 E:\\NMR2Structure\\logs); '
                             '每次运行生成 NMR2Structure_YYYY-MM-DD_HH-MM-SS.log')
    parser.add_argument('--no-log', action='store_true',
                        help='禁用日志记录 (默认每次运行都会把屏幕内容写入日志文件)')
    # Part 2.B: 溶剂残余峰过滤开关
    parser.add_argument('--keep-solvent', action='store_true',
                        help='禁用溶剂残余峰过滤 (默认启用 -- 仪器记录的 CDCl3 7.26 '
                             '/ DMSO 2.50 等溶剂信号会被剔除, 避免污染检索)')
    # 已知化合物流程: 自动检测 + 复核对比
    parser.add_argument('--real-structure-dir', default=_DEFAULT_REAL_STRUCTURE_DIR,
                        help='已知化合物 (真实结构 + 数据归属表) 所在根目录. '
                             '当 --folder 路径末段化合物名的子目录中存在 *数据归属*.xlsx 时, '
                             '自动切换到已知化合物流程 (derive + 复核). '
                             '默认 E:\\NMR2Structure\\真实结构.')
    parser.add_argument('--known-c-tol', type=float, default=0.5,
                        help='已知化合物复核对比时的 13C 匹配容差 (ppm, 默认 0.5)')
    parser.add_argument('--known-h-tol', type=float, default=0.1,
                        help='已知化合物复核对比时的 1H 匹配容差 (ppm, 默认 0.1)')
    args = parser.parse_args(argv)

    # 处理 --no-2d / --no-derive 互斥开关
    if args.no_2d:
        args.use_2d = False
    if args.no_derive:
        args.derive = False

    # 初始化日志 (默认开启, --no-log 可禁用)
    if not args.no_log:
        log_file = _setup_logging(log_dir=args.log_dir)
        log = _get_logger('main')
        log.info(f'日志文件: {log_file}')

    folder = os.path.abspath(args.folder)
    if not os.path.isdir(folder):
        print(f"ERROR: 原始数据目录不存在: {folder}", file=sys.stderr)
        return 2

    # 推断化合物名 (优先 --compound-name, 否则取 folder 路径最后一段)
    compound_name = args.compound_name or os.path.basename(folder.rstrip(r'\/'))

    # Part 3.5: 检测已知化合物 (在 E:\NMR2Structure\真实结构/<compound_name>/ 下查找 *数据归属*.xlsx)
    # 已知化合物流程: derive 不用 SMILES, 复核阶段才用 SMILES 做结构判等
    known_table_path: str | None = None
    known_smiles: str = ''
    # 自动检测: real_structure_dir 不存在 / <compound> 子目录缺失 / 无 *数据归属*.xlsx
    # 时 _discover_known_compound 返回 None, 自动跳过已知化合物复核 (无需额外参数).
    try:
        known_result = _discover_known_compound(compound_name, args.real_structure_dir)
        if known_result:
            known_table_path, known_smiles = known_result
            print(f"[已知化合物] 检测到标准归属表: {known_table_path}", file=sys.stderr)
            if known_smiles:
                print(f"[已知化合物] 已知 SMILES (仅复核用): "
                      f"{known_smiles[:60]}{'...' if len(known_smiles) > 60 else ''}",
                      file=sys.stderr)
            else:
                print(f"[已知化合物] 未检测到同目录 smiles.txt (仅复核数据可用, "
                      f"无法做 SMILES 结构判等)", file=sys.stderr)
    except Exception as ex:
        print(f"WARN: 已知化合物检测失败: {ex}", file=sys.stderr)

    # 数据库自动发现 (--db 支持目录或单文件)
    dbs = _discover_databases(args.db)
    # --lucy-ng-db 显式指定时覆盖自动发现
    if args.lucy_ng_db:
        dbs['lucy_ng'] = os.path.abspath(args.lucy_ng_db)
    sdf_path = dbs['sdf']
    lucy_ng_path = dbs['lucy_ng']
    if not args.quiet:
        print(f"[数据库] 发现: SDF={os.path.basename(sdf_path) if sdf_path else '无'}, "
              f"lucy-ng={os.path.basename(lucy_ng_path) if lucy_ng_path else '无'}",
              file=sys.stderr)

    # 推断化合物名 (已在前面定义)
    output_dir = os.path.join(os.path.abspath(args.output_root), compound_name)
    os.makedirs(output_dir, exist_ok=True)
    out_prefix = os.path.join(output_dir, compound_name)

    # 1. 读 Bruker
    print(f"[1/5] 读取 Bruker 谱图: {folder}", file=sys.stderr)
    specs = read_bruker_folder(folder)
    if not specs:
        print("ERROR: 未在原始数据目录中发现任何 1D 谱图 (peaklist.xml)", file=sys.stderr)
        return 3
    for s in specs:
        print(f"  - exp {s.exp_no}: {s.nucleus} ({s.pulse_program}) "
              f"solvent={s.solvent or '?'} peaks={len(s.peaks)} integrals={len(s.integrals)}",
              file=sys.stderr)

    # 1.5 检测 2D 实验 (--use-2d 开启)
    spectra_2d = {}
    interpretation_2d = None
    if args.use_2d:
        print(f"[1.5] 检测 2D NMR 实验 ...", file=sys.stderr)
        try:
            from .io_bruker_2d import read_bruker_2d_folder
            spectra_2d = read_bruker_2d_folder(folder, quiet=args.quiet)
            if spectra_2d:
                print(f"  -> 检测到 2D 实验: {list(spectra_2d.keys())}", file=sys.stderr)
                for etype, spec in spectra_2d.items():
                    print(f"     {etype}: exp={spec.exp_no} peaks={len(spec.peaks)} "
                          f"({spec.f1_nucleus}-{spec.f2_nucleus})", file=sys.stderr)
            else:
                print(f"  -> 未检测到 2D 数据, 回退到 1D 流程", file=sys.stderr)
        except Exception as ex:
            print(f"WARN: 2D 检测失败, 回退到 1D: {ex}", file=sys.stderr)
            spectra_2d = {}
    else:
        print(f"[1.5] 跳过 2D NMR 检测 (使用 --use-2d 启用)", file=sys.stderr)

    # 2. 构造问题
    problem = _build_problem_from_bruker(specs, keep_solvent=args.keep_solvent)
    print(f"[2/5] 观测峰: 13C={len(problem.carbon_shifts)} 1H={len(problem.proton_shifts)} "
          f"solvent={problem.solvent}", file=sys.stderr)

    # 默认不强制溶剂兼容 (容差内已能吸收常见溶剂差 0.5-1 ppm)
    constraints = Constraints(formula=args.formula)

    # 2.5 DEPT 推断 (用 bridge_nmr2struct 读取并推断碳类型)
    print(f"[2.5/5] DEPT 推断碳类型 ...", file=sys.stderr)
    _tmp_hnmr = out_prefix + '._hnmr_tmp.txt'
    _tmp_cnmr = out_prefix + '._cnmr_tmp.txt'
    bridge_info = convert_bruker_folder_to_nmr2struct(
        folder, _tmp_hnmr, _tmp_cnmr, use_dept=True,
    )
    obs_dept_summary: list[dict] = bridge_info.get('dept_summary', [])
    dept_counts = bridge_info.get('dept_counts', {})
    if obs_dept_summary:
        print(f"  -> 推断 {len(obs_dept_summary)} 个 C 类型: "
              f"CH={dept_counts.get('CH', 0)} "
              f"CH2={dept_counts.get('CH2', 0)} "
              f"CH3={dept_counts.get('CH3', 0)} "
              f"CH2/CH3={dept_counts.get('CH2/CH3', 0)} "
              f"Cq={dept_counts.get('Cq', 0)}", file=sys.stderr)
    else:
        print(f"  -> 无 DEPT 数据, 跳过碳类型推断", file=sys.stderr)
    # 清理临时文件
    for tmp in (_tmp_hnmr, _tmp_cnmr):
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass

    # 2.8 2D NMR 逐步解读 (有 2D 数据时, 在检索前完成)
    hsqc_pairs_for_retrieve = None
    if spectra_2d:
        print(f"[2.8] 2D NMR 逐步解读 (化学家工作流) ...", file=sys.stderr)
        try:
            from .interpret_2d import interpret_2d_sequential, write_step_files
            from .io_bruker import parse_multiplets_from_peaklist
            # 解析 1H 多重峰 (J 耦合). 优先读取 Multiplet 标签/Peak1D+mult,
            # 未做手工多重峰分析时 fallback 到 infer_multiplets_from_peaks.
            multiplets = []
            h1_peaklist_path = None
            h1_field_mhz = 0.0
            for spec in specs:
                if spec.nucleus == '1H' and spec.peaklist_path:
                    multiplets = parse_multiplets_from_peaklist(
                        spec.peaklist_path, field_mhz=spec.field_mhz)
                    h1_peaklist_path = spec.peaklist_path
                    h1_field_mhz = spec.field_mhz
                    break
            interpretation_2d = interpret_2d_sequential(
                spectra_2d=spectra_2d,
                carbon_shifts_1d=problem.carbon_shifts,
                proton_shifts_1d=problem.proton_shifts,
                dept_summary=obs_dept_summary,
                proton_integrals=problem.proton_integrals,
                multiplets=multiplets,
                h1_peaklist_path=h1_peaklist_path,
                h1_field_mhz=h1_field_mhz,
                solvent=problem.solvent or '',
                quiet=args.quiet,
            )
            # 中间结果写入独立 CSV (化学家复核用)
            step_files = write_step_files(interpretation_2d, output_dir, spectra_2d)
            if not args.quiet:
                print(f"  -> 中间结果已写入 {len(step_files)} 个 CSV 文件",
                      file=sys.stderr)
            # 提取 HSQC 对用于配对检索
            if interpretation_2d.ch_pairs:
                hsqc_pairs_for_retrieve = [
                    (p.h_ppm, p.c_ppm) for p in interpretation_2d.ch_pairs
                ]
                print(f"  -> HSQC 配对检索已启用 ({len(hsqc_pairs_for_retrieve)} 个 C-H 对)",
                      file=sys.stderr)
        except Exception as ex:
            print(f"WARN: 2D 解读失败, 回退 1D 检索: {ex}", file=sys.stderr)
            interpretation_2d = None

    # Part 3: --derive 模式 - 从谱图数据直接推导归属表 (不依赖候选检索/打分)
    if args.derive:
        if interpretation_2d is None or not interpretation_2d.ch_pairs:
            # 无 2D 数据时自动降级: 提示用户, 回退到候选检索流程
            print(f"WARN: derive 模式需要 2D 数据 (HSQC), 未检测到. "
                  f"自动回退到候选检索流程.", file=sys.stderr)
            args.derive = False  # 降级

    if args.derive:
        from .assembly import (
            assemble_assignment_table, write_assignment_csv,
            write_assignment_xlsx, write_assignment_meta,
        )
        from .correlation_2d import extract_cosy_pairs, extract_hmbc_pairs
        print(f"[Derive] 从谱图数据直接推导归属表 ...", file=sys.stderr)

        # 提取 COSY/HMBC/ROESY 交叉峰
        cosy_pairs = []
        hmbc_pairs = []
        roesy_pairs = []
        if spectra_2d:
            if 'COSY' in spectra_2d:
                cosy_pairs = extract_cosy_pairs(spectra_2d['COSY'], min_snr=10.0)
            if 'HMBC' in spectra_2d:
                hmbc_result = extract_hmbc_pairs(
                    spectra_2d['HMBC'], min_snr=10.0,
                    hsqc_pairs=[(p.h_ppm, p.c_ppm) for p in interpretation_2d.ch_pairs])
                hmbc_pairs = hmbc_result[0] if isinstance(hmbc_result, tuple) else hmbc_result
            if 'ROESY' in spectra_2d:
                roesy_pairs = extract_cosy_pairs(spectra_2d['ROESY'], min_snr=10.0)
            elif 'NOESY' in spectra_2d:
                roesy_pairs = extract_cosy_pairs(spectra_2d['NOESY'], min_snr=10.0)

        print(f"  COSY: {len(cosy_pairs)} 对, HMBC: {len(hmbc_pairs)} 对, "
              f"ROESY: {len(roesy_pairs)} 对", file=sys.stderr)

        # 组装归属表
        derive_rows = assemble_assignment_table(
            ch_pairs=interpretation_2d.ch_pairs,
            carbon_shifts_1d=problem.carbon_shifts,
            proton_shifts_1d=problem.proton_shifts,
            cosy_pairs=cosy_pairs,
            hmbc_pairs=hmbc_pairs,
            roesy_pairs=roesy_pairs,
        )
        print(f"  归属表: {len(derive_rows)} 行", file=sys.stderr)

        # 输出
        derive_dir = os.path.join(output_dir, 'derive_归属表')
        os.makedirs(derive_dir, exist_ok=True)
        csv_path = os.path.join(derive_dir, 'assignment_table.csv')
        xlsx_path = os.path.join(derive_dir, 'assignment_table.xlsx')
        meta_path = os.path.join(derive_dir, 'meta.txt')
        write_assignment_csv(derive_rows, csv_path)
        # 元信息 (Solvent / Formula / HRESIMS / Field 等) 单独写到 meta.txt
        meta = {}
        if args.formula:
            meta['Formula'] = args.formula
        if problem.solvent:
            meta['Solvent'] = problem.solvent
        if problem.field_mhz:
            meta['Field_MHz'] = f'{problem.field_mhz:.2f}'
        for key in ('hresims', 'optical_rotation', 'ir', 'compound_class'):
            val = getattr(args, key, None) or ''
            if val:
                meta[key.replace('_', ' ').upper()] = val
        write_assignment_xlsx(derive_rows, xlsx_path)
        write_assignment_meta(meta_path, meta)
        print(f"  输出: {csv_path}", file=sys.stderr)
        print(f"  输出: {xlsx_path}", file=sys.stderr)
        if meta:
            print(f"  输出: {meta_path}", file=sys.stderr)

        # 统计
        n_cosy = sum(1 for r in derive_rows if r.cosy_partners)
        n_hmbc = sum(1 for r in derive_rows if r.hmbc_targets)
        n_roesy = sum(1 for r in derive_rows if r.roesy_contacts)
        n_cq = sum(1 for r in derive_rows if r.is_quaternary)
        print(f"  统计: 季碳={n_cq}, COSY非空={n_cosy}, "
              f"HMBC非空={n_hmbc}, ROESY非空={n_roesy}", file=sys.stderr)
        print(f"[Derive] 归属表推导完成.", file=sys.stderr)

    # 已知化合物复核: derive 完成后, 用 SMILES (仅复核用) + 标准归属表做对比
    # 重要: SMILES 不参与推导过程, 仅在结果复核阶段用于化合物结构判等
    if known_table_path and args.derive:
        derive_xlsx_path = xlsx_path  # derive_归属表/assignment_table.xlsx
        derive_meta_path = meta_path  # derive_归属表/meta.txt
        if os.path.isfile(derive_xlsx_path):
            print(f"[复核] 检测到已知化合物, 开始对比标准归属表 ...", file=sys.stderr)
            try:
                std_table = _read_standard_assignment_table(
                    known_table_path, known_smiles=known_smiles,
                )
                comparison = _compare_program_to_standard(
                    derive_xlsx_path, std_table,
                    prog_meta_path=derive_meta_path,
                    c_tol=args.known_c_tol, h_tol=args.known_h_tol,
                    use_smiles=bool(known_smiles),
                    type_placeholder_ok=True,
                )
                # 写出匹配 xlsx + 报告 (到 derive_归属表/ 目录, 与程序 xlsx 同级)
                matched_xlsx = os.path.join(
                    os.path.dirname(derive_xlsx_path),
                    'assignment_table_joint.xlsx',
                )
                comparison_txt = os.path.join(
                    os.path.dirname(derive_xlsx_path),
                    'comparison_report.txt',
                )
                try:
                    # 双 sheet 拼接表: '按标准表为准' + '按程序表为准'
                    _write_joint_xlsx(
                        derive_xlsx_path, std_table, comparison, matched_xlsx,
                    )
                    _write_comparison_report(
                        comparison_txt, derive_xlsx_path, std_table, comparison,
                        c_tol=args.known_c_tol, h_tol=args.known_h_tol,
                        prog_meta_path=derive_meta_path,
                    )
                    print(f"  -> 拼接表 (双sheet): {matched_xlsx}", file=sys.stderr)
                    print(f"  -> 对比报告: {comparison_txt}", file=sys.stderr)
                except Exception as ex:
                    print(f"  WARN: 复核文件写入失败: {ex}", file=sys.stderr)
                print(f"  -> 匹配率: {comparison.match_rate:.1f}% "
                      f"({len(comparison.matches)}/{max(comparison.n_prog, comparison.n_std)})",
                      file=sys.stderr)
                print(f"  -> 判定: {comparison.verdict}", file=sys.stderr)
                if comparison.smiles_equal is not None:
                    eq_str = '一致' if comparison.smiles_equal else '不一致'
                    print(f"  -> SMILES 结构判等: {eq_str}", file=sys.stderr)
            except Exception as ex:
                print(f"WARN: 已知化合物复核失败: {ex}", file=sys.stderr)

    # derive 完成后默认结束 (不进入检索); 仅当 --retrieve 或 --no-derive 时才继续检索
    if args.derive and not args.retrieve:
        print(f"DONE. 归属表已输出. 提示: 加 --retrieve 可继续候选检索+打分.",
              file=sys.stderr)
        return 0

    # 3. 候选检索 (--retrieve 或 --no-derive 时执行)
    if not sdf_path and not lucy_ng_path:
        print(f"[3] 跳过候选检索 (未发现数据库). "
              f"提示: 用 --db 指定数据库目录", file=sys.stderr)
        # 无数据库时直接结束 (derive 已产出归属表)
        if args.derive:
            print(f"DONE. 归属表已输出 (无候选检索).", file=sys.stderr)
            return 0
        return 4
    use_hsqc_mode = hsqc_pairs_for_retrieve is not None
    mode_str = 'HSQC 配对匹配' if use_hsqc_mode else '1D 位移匹配'
    hits = []

    # 3a. nmrshiftdb2 SDF 检索 (优先)
    if sdf_path and os.path.isfile(sdf_path):
        print(f"[3a] nmrshiftdb2 检索 ({mode_str}, {os.path.basename(sdf_path)})", file=sys.stderr)
        progress_every = 0 if args.quiet else 5000
        hits = list(retrieve_candidates(
            sdf_path=sdf_path,
            problem=problem,
            constraints=constraints,
            c_tol=args.c_tol,
            h_tol=args.h_tol,
            max_hits=args.max_hits,
            progress_every=progress_every,
            obs_dept_summary=obs_dept_summary,
            hsqc_pairs=hsqc_pairs_for_retrieve,
        ))
        # HSQC 配对检索无结果时自动回退到 1D 检索
        if not hits and use_hsqc_mode:
            print(f"  -> HSQC 配对检索 0 命中, 回退到 1D 位移匹配 ...", file=sys.stderr)
            use_hsqc_mode = False
            hits = list(retrieve_candidates(
                sdf_path=sdf_path,
                problem=problem,
                constraints=constraints,
                c_tol=args.c_tol,
                h_tol=args.h_tol,
                max_hits=args.max_hits,
                progress_every=progress_every,
                obs_dept_summary=obs_dept_summary,
                hsqc_pairs=None,
            ))
        print(f"  -> nmrshiftdb2 命中 {len(hits)} 条", file=sys.stderr)
    else:
        print(f"[3a] 跳过 nmrshiftdb2 (未找到 SDF 文件)", file=sys.stderr)

    # 3b. lucy-ng 检索 (自动启用: 当发现 lucy-ng 数据库时)
    lucy_ng_hits: List = []
    if lucy_ng_path and os.path.isfile(lucy_ng_path):
        print(f"[3b] lucy-ng 数据库检索 ({os.path.basename(lucy_ng_path)}) ...", file=sys.stderr)
        try:
            from .retrieve_lucy_ng import get_lucy_ng_retriever, LucyNGRetriever
            # Part 3 Phase 7.3 (修复): 实测 13C 峰数 ≠ 候选碳数 (重叠峰会让数字偏多)
            # 用更宽松的碳数估计: 取 min(实测 13C 峰数, lucy_ng_carbon_tol 范围)
            # 如果有 formula 约束则不需要碳数过滤
            if args.formula:
                # 精确分子式约束时, 不限碳数
                est_carbon_count = None
                est_carbon_tol = 0
            else:
                est_carbon_count = len(problem.carbon_shifts)
                est_carbon_tol = max(args.lucy_ng_carbon_tol, 5)  # 至少 ±5 (覆盖重叠峰)
            retriever = get_lucy_ng_retriever(lucy_ng_path, quiet=args.quiet)
            if retriever is not None:
                lucy_ng_hits = retriever.retrieve_by_carbon_count(
                    observed_c13=problem.carbon_shifts,
                    observed_h1=problem.proton_shifts,
                    carbon_count=est_carbon_count,
                    carbon_tol=est_carbon_tol,
                    formula=args.formula,
                    c_tol=args.c_tol,
                    h_tol=args.h_tol,
                    min_c_matches=args.lucy_ng_min_matches,
                    max_results=args.lucy_ng_max_results,
                    progress_every=0 if args.quiet else 2000,
                    quiet=args.quiet,
                )
                # 转成与 retrieve_candidates 兼容的格式 (用 dict 模拟 hit.entry)
                from .retrieve import RetrievalHit
                for lh in lucy_ng_hits:
                    mock_entry = type('Entry', (), {})()
                    mock_entry.canonical_smiles = lh.canonical_smiles or lh.smiles
                    mock_entry.smiles = lh.smiles
                    mock_entry.inchi = lh.inchi
                    mock_entry.inchi_key = lh.inchikey
                    mock_entry.inchikey = lh.inchikey  # evidence.py 也可能用
                    mock_entry.mol_block = ''
                    mock_entry.spectrum_13c = []
                    mock_entry.spectrum_1h = []
                    mock_entry.formula = lh.formula
                    mock_entry.heavy_atoms = lh.carbon_count + sum(
                        1 for s in lh.formula if s.isalpha() and s.upper() not in ('C', 'H')
                    )
                    mock_entry.rec_id = None
                    mock_entry.shiftdb2_id = None
                    mock_entry.source = 'lucy-ng'
                    mock_entry.solvent = None
                    mock_entry.extra_tags = {'source': 'lucy-ng', 'lucy_ng_id': str(lh.compound_id)}
                    # 填充 lib_13c/lib_1h 让 evidence.py 能调用
                    from .io_sdf import SpectrumSignal
                    mock_entry.spectrum_13c = [
                        SpectrumSignal(shift=s, atom_index=i, multiplicity='?')
                        for i, s in enumerate(lh.lib_13c_shifts)
                    ]
                    # lucy-ng 没有 1H, 留空
                    mock_entry.spectrum_1h = []
                    # 提供 mol() 方法 (用 SMILES 解析)
                    from rdkit import Chem
                    def _make_mol(entry_self):
                        try:
                            return Chem.MolFromSmiles(entry_self.smiles)
                        except Exception:
                            return None
                    mock_entry.mol = lambda e=mock_entry: _make_mol(e)
                    mock_hit = RetrievalHit(
                        entry=mock_entry,
                        c_match_count=lh.c_match_count,
                        c_total_obs=lh.c_total_obs,
                        c_total_lib=lh.c_total_lib,
                        c_coverage=lh.c_coverage,
                        c_rms=lh.c_rms,
                        h_match_count=lh.h_match_count,
                        h_total_obs=lh.h_total_obs,
                        h_total_lib=lh.h_total_lib,
                        h_coverage=lh.h_coverage,
                        h_rms=lh.h_rms,
                        formula_ok=True,
                        solvent_ok=True,
                        heavy_atoms=mock_entry.heavy_atoms,
                    )
                    # 在 evidence 阶段需额外信息, 加在 mock_entry 上
                    mock_entry.lucy_score = lh.score
                    mock_entry.lucy_compound_id = lh.compound_id
                    mock_entry.lucy_lib_13c = lh.lib_13c_shifts
                    hits.append(mock_hit)
                print(f"  -> lucy-ng 命中 {len(lucy_ng_hits)} 条 (总计 {len(hits)} 条)",
                      file=sys.stderr)
        except Exception as ex:
            print(f"WARN: lucy-ng 检索失败: {ex}", file=sys.stderr)

    if not hits:
        # 本地库无命中 → 自动扩大到联网数据库
        # 顺序: COCONUT (天然产物专属, 69.5万+ NP) → PubChem (通用兑底, 1亿+)
        from .retrieve import RetrievalHit as _RH

        def _hits_from_online(candidates, source_tag: str):
            """将联网库候选 (COCONUT/PubChem) 转为 RetrievalHit 列表."""
            out = []
            for cand in candidates:
                mock_entry = type('Entry', (), {})()
                mock_entry.canonical_smiles = cand.smiles
                mock_entry.smiles = cand.smiles
                mock_entry.inchi = None
                mock_entry.inchi_key = None
                mock_entry.inchikey = None
                mock_entry.mol_block = ''
                mock_entry.spectrum_13c = []
                mock_entry.spectrum_1h = []
                mock_entry.formula = cand.molecular_formula
                mock_entry.heavy_atoms = cand.heavy_atom_count
                mock_entry.rec_id = None
                mock_entry.shiftdb2_id = None
                mock_entry.source = source_tag
                mock_entry.solvent = None
                # 字段名按 source_tag 区分 (PubChem: cid; COCONUT: coconut_id)
                if source_tag == 'pubchem':
                    mock_entry.pubchem_cid = cand.cid
                    mock_entry.iupac_name = cand.iupac_name
                    mock_entry.extra_tags = {'source': 'pubchem', 'pubchem_cid': str(cand.cid)}
                else:  # coconut
                    mock_entry.coconut_id = cand.coconut_id
                    mock_entry.iupac_name = cand.name
                    mock_entry.extra_tags = {'source': 'coconut', 'coconut_id': str(cand.coconut_id)}
                from rdkit import Chem as _Chem
                def _make_mol(entry_self):
                    try:
                        return _Chem.MolFromSmiles(entry_self.smiles)
                    except Exception:
                        return None
                mock_entry.mol = lambda e=mock_entry: _make_mol(e)
                out.append(_RH(
                    entry=mock_entry,
                    c_match_count=0,
                    c_total_obs=len(problem.carbon_shifts),
                    c_total_lib=0,
                    c_coverage=0.0,
                    c_rms=99.0,
                    h_match_count=0,
                    h_total_obs=len(problem.proton_shifts),
                    h_total_lib=0,
                    h_coverage=0.0,
                    h_rms=99.0,
                    formula_ok=True,
                    solvent_ok=True,
                    heavy_atoms=cand.heavy_atom_count,
                ))
            return out

        # 3c. COCONUT (天然产物专属, 69.5万+ NP 结构, 联网)
        _coconut_auto = (not args.no_coconut)
        if _coconut_auto:
            print(f"[3c] 本地库无命中, 自动扩大到 COCONUT 检索 (天然产物, 联网) ...",
                  file=sys.stderr)
            try:
                from .retrieve_coconut import retrieve_coconut_candidates
                coconut_results = retrieve_coconut_candidates(
                    formula=constraints.formula,
                    c13_shifts=problem.carbon_shifts,
                    max_results=50,
                    quiet=args.quiet,
                )
                if coconut_results:
                    print(f"  -> COCONUT 候选: {len(coconut_results)} 条", file=sys.stderr)
                    hits.extend(_hits_from_online(coconut_results, 'coconut'))
                else:
                    print(f"  -> COCONUT 也无候选", file=sys.stderr)
            except Exception as ex:
                print(f"WARN: COCONUT 检索失败 (网络异常?): {ex}", file=sys.stderr)

        # 3d. PubChem (通用兑底, 1亿+ 结构, 联网)
        _pubchem_auto = (not args.no_pubchem)
        if _pubchem_auto and not hits:
            print(f"[3d] 本地库+COCONUT 仍无命中, 自动扩大到 PubChem 检索 (联网) ...",
                  file=sys.stderr)
            try:
                from .retrieve_pubchem import retrieve_pubchem_candidates
                pubchem_results = retrieve_pubchem_candidates(
                    formula=constraints.formula,
                    c13_shifts=problem.carbon_shifts,
                    max_results=50,
                    quiet=args.quiet,
                )
                if pubchem_results:
                    print(f"  -> PubChem 候选: {len(pubchem_results)} 条", file=sys.stderr)
                    hits.extend(_hits_from_online(pubchem_results, 'pubchem'))
                else:
                    print(f"  -> PubChem 也无候选", file=sys.stderr)
            except Exception as ex:
                print(f"WARN: PubChem 检索失败 (网络异常?): {ex}", file=sys.stderr)

    if not hits:
        print("WARN: 未找到任何候选。可放宽容差或移除分子式约束后重试。", file=sys.stderr)
        # 即使 0 命中也要写报告, 给出未匹配峰
        meta = {
            'problem': {
                'carbon_shifts': problem.carbon_shifts,
                'proton_shifts': problem.proton_shifts,
                'solvent': problem.solvent,
            },
            'constraints': {'formula': constraints.formula},
            'sdf_path': sdf_path,
            'folder': folder,
        }
        # 写空报告（明示无候选）
        write_report_json(out_prefix + '.json', [], meta)
        write_report_markdown(out_prefix + '.md', [], meta, top_n=args.top_n)
        return 4

    # 4. 打分
    print(f"[4/5] 打分 + 证据链 ...", file=sys.stderr)
    scored = score_candidates(hits, use_2d_scoring=use_hsqc_mode)
    evidences = [build_evidence(
        s,
        observed_c=problem.carbon_shifts,
        observed_h=problem.proton_shifts,
        c_tol=args.c_tol,
        h_tol=args.h_tol,
        formula_constraint=constraints.formula,
        obs_dept_summary=obs_dept_summary,
    ) for s in scored]

    # 5. nmrsim 反向校验（可选）
    if args.use_nmrsim:
        print(f"[5/5] nmrsim 反向校验 (cosine 相似度) ...", file=sys.stderr)
        for e in evidences:
            pred_peaks = [s.shift for s in e.smiles and scored[0].hit.entry.spectrum_1h]  # noqa
            # 简单实现：用 entry 的 1H 峰作为预测
            for s in scored:
                if s.hit.entry.smiles == e.smiles:
                    pred_peaks = [sig.shift for sig in s.hit.entry.spectrum_1h]
                    sim = check_candidate_with_nmrsim(
                        observed_peaks_ppm=problem.proton_shifts,
                        predicted_peaks_ppm=pred_peaks,
                        field_mhz=problem.field_mhz or 500.0,
                    )
                    e.notes.append(f"nmrsim 反向校验 cosine = {sim:.3f}")
                    break
    else:
        print(f"[5/5] 跳过 nmrsim 反向校验", file=sys.stderr)

    # 5.5 2D NMR 关联验证 (可选, 依赖 --use-2d + 实际 2D 数据)
    if spectra_2d:
        print(f"[5.5] 2D NMR 关联验证 (HSQC/HMBC/COSY) ...", file=sys.stderr)
        try:
            from .correlation_2d import score_candidate_2d
            for e in evidences[:args.top_n]:
                # 找到对应的 entry (从 scored 中匹配 smiles)
                entry_for_2d = None
                for s in scored:
                    if (s.hit.entry.canonical_smiles or s.hit.entry.smiles) == e.smiles:
                        entry_for_2d = s.hit.entry
                        break
                if entry_for_2d is None:
                    continue
                score_2d = score_candidate_2d(
                    entry_for_2d, spectra_2d,
                    h_tol=args.h_tol, c_tol=args.c_tol,
                )
                e.correlation_2d_score = score_2d
                e.notes.extend(score_2d.notes)
            print(f"  -> 2D 验证完成: {min(args.top_n, len(evidences))} 个候选", file=sys.stderr)
        except Exception as ex:
            print(f"WARN: 2D 关联验证失败: {ex}", file=sys.stderr)
    else:
        if args.use_2d:
            print(f"[5.5] 跳过 2D 验证 (无 2D 数据)", file=sys.stderr)
        else:
            print(f"[5.5] 跳过 2D 验证", file=sys.stderr)

    # 6. HOSE 13C 预测验证（可选）
    hose_verifications = []
    if args.use_hose:
        print(f"[6] HOSE 13C 预测验证 ...", file=sys.stderr)
        # 优先使用 lucy-ng 完整数据库 (如 --lucy-ng-db 指定)
        lucy_provider = None
        try:
            from .lucy_ng_hose import get_lucy_ng_provider, verify_candidate_with_lucy_db
            lucy_provider = get_lucy_ng_provider(lucy_ng_path, quiet=args.quiet)
        except Exception as ex:
            if not args.quiet:
                print(f"  [HOSE] lucy-ng 集成跳过: {ex}", file=sys.stderr)

        if lucy_provider is not None:
            # 用 lucy-ng 数据库预测 (精度高)
            for e in evidences[:args.top_n]:
                vs = verify_candidate_with_lucy_db(
                    lucy_provider, e.smiles, problem.carbon_shifts,
                    c_tol=args.c_tol + 1.0,
                )
                hose_verifications.append(vs)
                e.notes.extend(vs.notes)
            print(f"  -> HOSE 验证完成 (lucy-ng 数据库): {len(hose_verifications)} 个候选",
                  file=sys.stderr)
        else:
            # 回退到本地 SDF 自建表
            try:
                from .nmr_predict import get_predictor, verify_candidate
                hose_cache = args.hose_cache or os.path.join(
                    os.path.dirname(sdf_path), '_hose_lookup_cache.tsv')
                predictor = get_predictor(
                    sdf_path=sdf_path,
                    table_path=hose_cache,
                    max_records=args.hose_max_records,
                    quiet=args.quiet,
                )
                for e in evidences[:args.top_n]:
                    vs = verify_candidate(
                        predictor, e.smiles, problem.carbon_shifts, c_tol=args.c_tol + 1.0,
                    )
                    hose_verifications.append(vs)
                    e.notes.extend(vs.notes)
                print(f"  -> HOSE 验证完成 (本地 SDF 表): {len(hose_verifications)} 个候选",
                      file=sys.stderr)
            except Exception as ex:
                print(f"WARN: HOSE 预测验证失败: {ex}", file=sys.stderr)
    else:
        print(f"[6] 跳过 HOSE 预测验证", file=sys.stderr)

    # 7. COCONUT + PubChem 强制扩展检索 (--use-coconut/--use-pubchem 显式传入时执行)
    # 注: 本地库无命中时的自动触发已在 step 3c/3d 处理
    coconut_candidates = []
    if args.use_coconut:
        print(f"[7a] COCONUT 强制扩展检索 (--use-coconut) ...", file=sys.stderr)
        try:
            from .retrieve_coconut import retrieve_coconut_candidates
            coconut_candidates = retrieve_coconut_candidates(
                formula=constraints.formula,
                c13_shifts=problem.carbon_shifts,
                max_results=50,
                quiet=args.quiet,
            )
            print(f"  -> COCONUT 候选: {len(coconut_candidates)} 条", file=sys.stderr)
        except Exception as ex:
            print(f"WARN: COCONUT 检索失败: {ex}", file=sys.stderr)
    pubchem_candidates = []
    if args.use_pubchem:
        print(f"[7b] PubChem 强制扩展检索 (--use-pubchem) ...", file=sys.stderr)
        try:
            from .retrieve_pubchem import retrieve_pubchem_candidates
            pubchem_candidates = retrieve_pubchem_candidates(
                formula=constraints.formula,
                c13_shifts=problem.carbon_shifts,
                max_results=50,
                quiet=args.quiet,
            )
            print(f"  -> PubChem 候选: {len(pubchem_candidates)} 条", file=sys.stderr)
        except Exception as ex:
            print(f"WARN: PubChem 检索失败: {ex}", file=sys.stderr)
    if not args.use_coconut and not args.use_pubchem:
        print(f"[7] 跳过 COCONUT + PubChem 强制检索", file=sys.stderr)

    # 8. 片段组装（可选）
    fragment_candidates = []
    if args.use_fragment_assembly:
        print(f"[8] 片段组装候选生成 ...", file=sys.stderr)
        try:
            from .fragment_assembler import assemble_candidates, build_constraints
            frag_constraint = build_constraints(
                c13_shifts=problem.carbon_shifts,
                dept_summary=obs_dept_summary,
                formula=constraints.formula,
            )
            fragment_candidates = assemble_candidates(
                frag_constraint, max_candidates=100, quiet=args.quiet,
            )
            print(f"  -> 片段组装候选: {len(fragment_candidates)} 条", file=sys.stderr)
        except Exception as ex:
            print(f"WARN: 片段组装失败: {ex}", file=sys.stderr)
    else:
        print(f"[8] 跳过片段组装", file=sys.stderr)

    # 9. 结构推断置信度评估（可选）
    determination_result = None
    if args.determine:
        print(f"[9] 结构推断置信度评估 ...", file=sys.stderr)
        try:
            from .determination import assess_determination
            # Part 2.F: 提取 NOE / J 命中率
            noe_match_rates = None
            j_match_rates = None
            if noe_verification_results:
                noe_match_rates = []
                for r in noe_verification_results:
                    n_total = r.get('n_total', 0)
                    n_matched = r.get('n_matched', 0)
                    noe_match_rates.append(
                        n_matched / n_total if n_total > 0 else 0.0
                    )
            if j_verification_results:
                j_match_rates = []
                for r in j_verification_results:
                    n_total = r.get('n_total', 0)
                    n_matched = r.get('n_matched', 0)
                    j_match_rates.append(
                        n_matched / n_total if n_total > 0 else 0.0
                    )
            determination_result = assess_determination(
                evidences=evidences[:args.top_n],
                hose_verifications=hose_verifications if hose_verifications else None,
                noe_match_rates=noe_match_rates,
                j_match_rates=j_match_rates,
            )
            print(f"  -> 结论: {determination_result.conclusion_label} "
                  f"(置信度={determination_result.confidence:.3f})", file=sys.stderr)
            if determination_result.top_candidate_smiles:
                print(f"  -> 推断结构: {determination_result.top_candidate_smiles}",
                      file=sys.stderr)
        except Exception as ex:
            print(f"WARN: 置信度评估失败: {ex}", file=sys.stderr)
    else:
        print(f"[9] 跳过结构推断评估 (使用 --determine 启用)", file=sys.stderr)

    # 9.5 为 Top-N 各生成全归属表 (有 2D 数据时)
    full_assignment = None
    full_assignments_topn: List = []  # 存储各 Top-N 候选的归属表 (用于 CSV 输出)
    if interpretation_2d is not None and scored:
        from .correlation_2d import (
            build_full_assignment_table,
            write_assignment_table_csv,
        )
        print(f"[9.5] 为 Top-N 生成完整 NMR 归属表 ...", file=sys.stderr)
        # 新建 完整NMR归属表(论文Table1) 子目录
        assignment_dir = os.path.join(output_dir, '完整NMR归属表(论文Table1)')
        os.makedirs(assignment_dir, exist_ok=True)
        # 为 evidences 中每个候选查找其 entry 并生成归属表
        for ev in evidences[:min(args.top_n, len(evidences))]:
            entry = None
            for s in scored:
                if (s.hit.entry.canonical_smiles or s.hit.entry.smiles) == ev.smiles:
                    entry = s.hit.entry
                    break
            if entry is None:
                continue
            try:
                tbl = build_full_assignment_table(
                    entry, interpretation_2d,
                    observed_c13=problem.carbon_shifts,
                    observed_h1=problem.proton_shifts,
                    h_tol=args.h_tol, c_tol=args.c_tol,
                )
            except Exception as ex:
                print(f"  WARN: Top-{ev.rank} 归属表生成失败: {ex}", file=sys.stderr)
                continue
            # 写 CSV (论文 Table 1 格式)
            safe = re.sub(r'[^A-Za-z0-9_.-]', '_', ev.smiles)[:50] or f'top{ev.rank}'
            csv_path = os.path.join(
                assignment_dir, f'rank{ev.rank}_{safe}.csv')
            try:
                # 若旧文件存在则先删除 (Windows 防病毒可能锁定, 影响覆写)
                if os.path.exists(csv_path):
                    try:
                        os.remove(csv_path)
                    except OSError:
                        pass
                write_assignment_table_csv(tbl, csv_path, rank=ev.rank,
                                           score=getattr(ev, 'unified_score', None))
                full_assignments_topn.append({'rank': ev.rank, 'csv': csv_path})
                if ev.rank == 1:
                    full_assignment = tbl  # Top-1 保留原 meta['full_assignment']
                if not args.quiet:
                    print(f"  -> Top-{ev.rank}: {len(tbl.rows)} 行 -> "
                          f"{os.path.basename(csv_path)}", file=sys.stderr)
            except Exception as ex:
                print(f"  WARN: Top-{ev.rank} CSV 写入失败: {ex}", file=sys.stderr)
        if full_assignment is None and not full_assignments_topn:
            print(f"  -> Top-N 无对应库条目, 跳过归属表", file=sys.stderr)
        if full_assignments_topn and not args.quiet:
            print(f"  -> 归属表目录: {assignment_dir}", file=sys.stderr)

    # 9.7 化学家 Excel 风格归属表 (Part 3 Phase 2-6: 按 H 排列 + 立体氢拆分 + HMBC/COSY/ROESY 回填 + 多重性 + xlsx)
    # 与 9.5 (按 C 排列) 并行, 输出互补.
    assignment_v2_topn: List[dict] = []
    if scored:
        from .assignment_table_v2 import (
            build_assignment_table_v2,
            write_assignment_table_v2_csv,
        )
        from .excel_export import write_assignment_table_v2_xlsx
        print(f"[9.7] 为 Top-N 生成化学家 Excel 风格归属表 (按 H 排列) ...", file=sys.stderr)
        v2_dir = os.path.join(output_dir, '化学家归属表_v2')
        os.makedirs(v2_dir, exist_ok=True)
        multiplets_list = []
        if interpretation_2d is not None:
            multiplets_list = getattr(interpretation_2d, 'multiplet_analyses', []) or []
        for ev in evidences[:min(args.top_n, len(evidences))]:
            entry = None
            for s in scored:
                if (s.hit.entry.canonical_smiles or s.hit.entry.smiles) == ev.smiles:
                    entry = s.hit.entry
                    break
            if entry is None:
                continue
            try:
                tbl_v2 = build_assignment_table_v2(
                    smiles=entry.smiles,
                    mol_block=entry.mol_block,
                    observed_c13=problem.carbon_shifts,
                    observed_h1=problem.proton_shifts,
                    interpretation_2d=interpretation_2d,
                    multiplets=multiplets_list,
                    h_tol=args.h_tol,
                    c_tol=args.c_tol,
                )
                # Part 3 Phase 9.1: 注入化学家元信息
                tbl_v2.hresims = args.hresims
                tbl_v2.alpha = args.optical_rotation
                tbl_v2.ir = args.ir
            except Exception as ex:
                print(f"  WARN: Top-{ev.rank} v2 归属表生成失败: {ex}", file=sys.stderr)
                continue
            safe = re.sub(r'[^A-Za-z0-9_.-]', '_', ev.smiles)[:50] or f'top{ev.rank}'
            csv_path = os.path.join(v2_dir, f'rank{ev.rank}_{safe}_v2.csv')
            xlsx_path = os.path.join(v2_dir, f'rank{ev.rank}_{safe}_v2.xlsx')
            try:
                if os.path.exists(csv_path):
                    try:
                        os.remove(csv_path)
                    except OSError:
                        pass
                write_assignment_table_v2_csv(
                    tbl_v2, csv_path, rank=ev.rank,
                    score=getattr(ev, 'unified_score', None),
                    extra_meta={'Solvent': problem.solvent or '',
                                'Field_MHz': f'{problem.field_mhz:.2f}' if problem.field_mhz else ''},
                )
                write_assignment_table_v2_xlsx(
                    tbl_v2, xlsx_path, rank=ev.rank,
                    score=getattr(ev, 'unified_score', None),
                    extra_meta={'Solvent': problem.solvent or '',
                                'Field_MHz': f'{problem.field_mhz:.2f}' if problem.field_mhz else ''},
                )
                assignment_v2_topn.append({'rank': ev.rank, 'csv': csv_path, 'xlsx': xlsx_path})
                if ev.rank == 1 and not args.quiet:
                    print(f"  -> Top-1: {len(tbl_v2.rows)} 行, 未归属 "
                          f"13C={len(tbl_v2.unassigned_obs_c)} 1H={len(tbl_v2.unassigned_obs_h)}, "
                          f"置信度 {tbl_v2.confidence:.0%}",
                          file=sys.stderr)
                    print(f"     xlsx: {os.path.basename(xlsx_path)}", file=sys.stderr)
            except Exception as ex:
                print(f"  WARN: Top-{ev.rank} v2 表写入失败: {ex}", file=sys.stderr)
        if assignment_v2_topn and not args.quiet:
            print(f"  -> v2 归属表目录: {v2_dir}", file=sys.stderr)
    else:
        print(f"[9.7] 跳过化学家 Excel 风格归属表 (无候选)", file=sys.stderr)

    # 9.6 NOE 3D 构象验证 (有 NOESY 数据 + 归属表时)
    noe_verification_results: List[dict] = []
    if (interpretation_2d is not None
            and interpretation_2d.noe_contacts
            and full_assignments_topn):
        from .correlation_2d import verify_candidate_3d_noe
        print(f"[9.6] NOE 3D 构象验证 ({len(interpretation_2d.noe_contacts)} 个 NOE 对) ...",
              file=sys.stderr)
        for ev in evidences[:min(args.top_n, len(evidences))]:
            # 构建 h_assignments: dH_ppm -> atom_index
            # 从归属表 rows 提取 (h_ppm, atom_label) 映射到 RDKit mol
            h_assignments: dict = {}
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(ev.smiles)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)
                # 简化映射: 用归属表的 h_ppm 作为 key, 用 H 原子 index 作为 value
                # 由于归属表可能不完整, 用顺序匹配
                h_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'H']
                # 从归属表提取已归属的 h_ppm
                for fa_info in full_assignments_topn:
                    if fa_info.get('rank') == ev.rank:
                        # 简化: 用观测 proton_shifts 顺序匹配 H 原子
                        for i, h_ppm in enumerate(problem.proton_shifts):
                            if i < len(h_atoms):
                                h_assignments[round(h_ppm, 3)] = h_atoms[i]
                        break
            except Exception:
                pass

            if not h_assignments:
                continue

            try:
                noe_result = verify_candidate_3d_noe(
                    ev.smiles, interpretation_2d.noe_contacts, h_assignments)
                noe_verification_results.append({
                    'rank': ev.rank,
                    'smiles': ev.smiles,
                    'n_total': noe_result.n_total,
                    'n_matched': noe_result.n_matched,
                    'n_conflict': noe_result.n_conflict,
                    'rms_distance': noe_result.rms_distance,
                    'status': noe_result.verification_status,
                })
            except Exception as ex:
                print(f"  WARN: Top-{ev.rank} NOE 验证失败: {ex}", file=sys.stderr)

        # 写 step6_noe_3d_verification.csv
        if noe_verification_results:
            import csv as csv_mod
            path6 = os.path.join(output_dir, 'step6_noe_3d_verification.csv')
            with open(path6, 'w', newline='', encoding='utf-8') as f:
                w = csv_mod.writer(f)
                w.writerow(['候选rank', 'SMILES', 'NOE总数', '匹配', '冲突',
                            'RMS距离_A', '状态'])
                for r in noe_verification_results:
                    w.writerow([r['rank'], r['smiles'], r['n_total'],
                                r['n_matched'], r['n_conflict'],
                                f"{r['rms_distance']:.3f}", r['status']])
            if not args.quiet:
                print(f"  -> NOE 验证: {len(noe_verification_results)} 个候选 -> "
                      f"step6_noe_3d_verification.csv", file=sys.stderr)

    # 9.6b Karplus 3D J 二面角验证 (Part 2.C)
    j_verification_results: List[dict] = []
    if (interpretation_2d is not None
            and interpretation_2d.multiplet_analyses
            and full_assignments_topn):
        from .correlation_2d import verify_candidate_3d_j
        print(f"[9.6b] Karplus 3D J 验证 (Part 2.C) ...", file=sys.stderr)
        for ev in evidences[:min(args.top_n, len(evidences))]:
            # 构建 h_assignments (与 NOE 验证共用)
            h_assignments_j: dict = {}
            try:
                from rdkit import Chem as _Chem
                mol_j = _Chem.MolFromSmiles(ev.smiles)
                if mol_j is None:
                    continue
                mol_j = _Chem.AddHs(mol_j)
                h_atoms_j = [a.GetIdx() for a in mol_j.GetAtoms() if a.GetSymbol() == 'H']
                # 简化: 用观测 proton_shifts 顺序匹配 H 原子
                for i, h_ppm in enumerate(problem.proton_shifts):
                    if i < len(h_atoms_j):
                        h_assignments_j[round(h_ppm, 3)] = h_atoms_j[i]
            except Exception:
                pass

            if not h_assignments_j:
                continue

            # 构建 j_constraints: 从 multiplet_analyses 提取
            j_constraints: List[dict] = []
            mas = interpretation_2d.multiplet_analyses
            # 对每个 multiplet 找对应的 H label 和 J pattern
            for ma in mas:
                if not ma.j_couplings_Hz or ma.num_subpeaks < 2:
                    continue
                if ma.pattern not in ('d', 't', 'q', 'm'):
                    continue
                # 取第一个 J (从 spin_systems 中找相邻 H 是 Part 1 已做的, 但
                # 当前 Karplus 验证只需要单条 J 与 3D 二面角的预测对比)
                j_obs = ma.j_couplings_Hz[0]
                if j_obs <= 0:
                    continue
                # 简化: partner_h 用该 multiplet 所在 spin_system 中的其他 H
                partner_h = None
                for ss in interpretation_2d.spin_systems:
                    if round(ma.h_ppm, 3) in [round(h, 3) for h in ss.h_ppms]:
                        # 取另一个 H 作为 partner
                        for h in ss.h_ppms:
                            if round(h, 3) != round(ma.h_ppm, 3):
                                partner_h = h
                                break
                        break
                if partner_h is None:
                    # 无 partner: 跳过此 J (单 J 无法做 Karplus 二面角)
                    continue
                j_constraints.append({
                    'h1_label': ma.h_label,
                    'h2_label': ma.h_label + "-neighbor",  # 简化标签
                    'h1_ppm': ma.h_ppm,
                    'h2_ppm': partner_h,
                    'j_obs_hz': j_obs,
                    'multiplicity': ma.pattern,
                })

            if not j_constraints:
                continue

            try:
                j_result = verify_candidate_3d_j(
                    ev.smiles, j_constraints, h_assignments_j)
                j_verification_results.append({
                    'rank': ev.rank,
                    'smiles': ev.smiles,
                    'n_total': j_result.n_total,
                    'n_matched': j_result.n_matched,
                    'n_conflict': j_result.n_conflict,
                    'rms_j_error_hz': j_result.rms_j_error_hz,
                    'best_conformer_idx': j_result.best_conformer_idx,
                    'status': j_result.verification_status,
                    'karplus_constants': j_result.karplus_constants,
                })
            except Exception as ex:
                print(f"  WARN: Top-{ev.rank} Karplus J 验证失败: {ex}",
                      file=sys.stderr)

        # 写 step6b_j_karplus_3d_verification.csv
        if j_verification_results:
            import csv as csv_mod
            path6b = os.path.join(output_dir, 'step6b_j_karplus_3d_verification.csv')
            with open(path6b, 'w', newline='', encoding='utf-8') as f:
                w = csv_mod.writer(f)
                w.writerow(['候选rank', 'SMILES', 'J约束数', '匹配', '冲突',
                            'RMS误差Hz', 'best_conf_idx', '状态'])
                for r in j_verification_results:
                    w.writerow([
                        r['rank'], r['smiles'], r['n_total'],
                        r['n_matched'], r['n_conflict'],
                        f"{r['rms_j_error_hz']:.3f}", r['best_conformer_idx'],
                        r['status'],
                    ])
            if not args.quiet:
                print(f"  -> Karplus J 验证: {len(j_verification_results)} 个候选 -> "
                      f"step6b_j_karplus_3d_verification.csv", file=sys.stderr)

    # 输出 (out_prefix 已在 main() 开头构造, 路径 = output_root / compound_name / compound_name.*)
    meta = {
        'folder': folder,
        'sdf_path': sdf_path,
        'formula_constraint': constraints.formula,
        'solvent_observed': problem.solvent,
        'field_mhz': problem.field_mhz,
        'c_tol': args.c_tol,
        'h_tol': args.h_tol,
        'top_n': args.top_n,
        'used_nmrsim': bool(args.use_nmrsim),
        'used_hose': bool(args.use_hose),
        'used_hose_source': 'lucy-ng-derep.db' if (args.use_hose and lucy_ng_path) else 'local-sdf',
        'used_pubchem': bool(args.use_pubchem) or any(
            getattr(h.entry, 'source', None) == 'pubchem' for h in hits
        ),
        'used_coconut': bool(args.use_coconut) or any(
            getattr(h.entry, 'source', None) == 'coconut' for h in hits
        ),
        'used_fragment_assembly': bool(args.use_fragment_assembly),
        'used_2d': bool(args.use_2d and spectra_2d),
        '2d_experiments': list(spectra_2d.keys()) if spectra_2d else [],
        'dept_counts': dept_counts,
        'observed_carbon_shifts': problem.carbon_shifts,
        'observed_proton_shifts': problem.proton_shifts,
        'obs_dept_summary': obs_dept_summary,
    }
    # 结构推断结论
    if determination_result is not None:
        meta['determination'] = {
            'conclusion': determination_result.conclusion,
            'conclusion_label': determination_result.conclusion_label,
            'confidence': determination_result.confidence,
            'top_candidate_smiles': determination_result.top_candidate_smiles,
            'top_candidate_score': determination_result.top_candidate_score,
            'gap_to_second': determination_result.gap_to_second,
            'evidence_summary': determination_result.evidence_summary,
            'recommendations': determination_result.recommendations,
        }
    # HOSE 验证详情
    if hose_verifications:
        meta['hose_verifications'] = [
            {
                'smiles': hv.smiles,
                'hose_rms': hv.hose_rms,
                'hose_coverage': hv.hose_coverage,
                'hose_match_rate': hv.hose_match_rate,
                'verification_score': hv.verification_score,
                'notes': hv.notes,
            }
            for hv in hose_verifications
        ]
    # PubChem 拓展候选
    if pubchem_candidates:
        meta['pubchem_candidates'] = [
            {
                'cid': pc.cid,
                'smiles': pc.smiles,
                'iupac_name': pc.iupac_name,
                'molecular_formula': pc.molecular_formula,
                'molecular_weight': pc.molecular_weight,
            }
            for pc in pubchem_candidates[:20]  # 限制存储
        ]
    # 片段组装候选
    if fragment_candidates:
        meta['fragment_candidates'] = [
            {'smiles': fc.smiles, 'source': fc.source, 'notes': fc.notes}
            for fc in fragment_candidates[:20]
        ]
    # 2D NMR 关联验证分
    if spectra_2d:
        meta['correlation_2d'] = {}
        for e in evidences[:args.top_n]:
            if e.correlation_2d_score is not None:
                s2d = e.correlation_2d_score
                meta['correlation_2d'][e.smiles] = {
                    'hsqc_match_rate': s2d.hsqc_match_rate,
                    'hmbc_match_rate': s2d.hmbc_match_rate,
                    'cosy_match_rate': s2d.cosy_match_rate,
                    'combined_2d_score': s2d.combined_2d_score,
                    'n_hsqc_obs': s2d.n_hsqc_obs,
                    'n_hmbc_obs': s2d.n_hmbc_obs,
                    'n_cosy_obs': s2d.n_cosy_obs,
                    'n_hmbc_filtered_oneJ': getattr(s2d, 'n_hmbc_filtered_oneJ', 0),  # Part 2.A
                    'notes': s2d.notes,
                }
    # Part 2.B: 溶剂过滤统计 (总是启用, 与有无 2D 无关)
    meta['solvent_filter_stats'] = {
        'observed_solvent': problem.solvent or '',
        'n_h_filtered': problem.n_h_solvent_filtered,
        'n_c_filtered': problem.n_c_solvent_filtered,
    }
    # 2D 逐步解读中间结果 (化学家可审查)
    if interpretation_2d is not None:
        meta['interpretation_2d'] = {
            'ch_pairs': [
                {'h_ppm': p.h_ppm, 'c_ppm': p.c_ppm, 'multiplicity': p.multiplicity,
                 'c_label': p.assigned_c_label, 'h_label': p.assigned_h_label,
                 'h_integral': p.h_integral,
                 'n_h': round(p.h_integral) if p.h_integral is not None else None}
                for p in interpretation_2d.ch_pairs
            ],
            'spin_systems': [
                {'label': ss.label, 'h_ppms': ss.h_ppms,
                 'cosy_links': ss.cosy_links}
                for ss in interpretation_2d.spin_systems
            ],
            'hmbc_connections': [
                {'h_ppm': conn.h_ppm, 'c_ppm': conn.c_ppm,
                 'bond_type': conn.bond_type,
                 'is_oneJ_residue': conn.is_oneJ_residue}  # Part 2.A
                for conn in interpretation_2d.hmbc_connections
            ],
            'n_cq_inferred': interpretation_2d.n_cq_inferred,
            'cq_shifts': interpretation_2d.cq_shifts,
            'cq_evidence': interpretation_2d.cq_evidence,
            'notes': interpretation_2d.notes,
            'overlap_regions': [
                {'h_ppm_center': r.h_ppm_center, 'integral_total': r.integral_total,
                 'n_peaks_in_13c_slice': r.n_peaks_in_13c_slice,
                 'c_peaks': r.c_peaks, 'c_peak_intensities': r.c_peak_intensities,
                 'status': r.status, 'note': r.note}
                for r in interpretation_2d.overlap_regions
            ],
            'tocsy_spin_systems': [
                {'label': ss.label, 'h_ppms': ss.h_ppms}
                for ss in interpretation_2d.tocsy_spin_systems
            ],
            'noe_contacts': [
                {'h1_ppm': nc.h1_ppm, 'h2_ppm': nc.h2_ppm,
                 'intensity': nc.intensity, 'confidence': nc.confidence,
                 'expected_distance': nc.expected_distance}
                for nc in interpretation_2d.noe_contacts
            ],
            # Part 2.C 暴露: TopSpin Manual 多重峰精细化结果
            'multiplet_analyses': [
                {'h_ppm': ma.h_ppm, 'c_ppm': ma.c_ppm,
                 'pattern': ma.pattern, 'num_subpeaks': ma.num_subpeaks,
                 'j_couplings_Hz': ma.j_couplings_Hz,
                 'halfwidth_Hz': ma.halfwidth_Hz,
                 'quality_grade': ma.quality_grade,
                 'c_label': ma.c_label, 'h_label': ma.h_label}
                for ma in interpretation_2d.multiplet_analyses
            ],
            'assignment_validation': {
                'total_observed_c': interpretation_2d.assignment_validation.total_observed_c,
                'assigned_c': interpretation_2d.assignment_validation.assigned_c,
                'total_observed_h': interpretation_2d.assignment_validation.total_observed_h,
                'assigned_h': interpretation_2d.assignment_validation.assigned_h,
                'unassigned_c': interpretation_2d.assignment_validation.unassigned_c,
                'unassigned_h': interpretation_2d.assignment_validation.unassigned_h,
                'warnings': interpretation_2d.assignment_validation.warnings,
                'is_complete': interpretation_2d.assignment_validation.is_complete,
            } if interpretation_2d.assignment_validation else None,
        }
    # NOE 3D 验证结果
    if noe_verification_results:
        meta['noe_verification'] = noe_verification_results

    # Part 2.E: Karplus J 3D 验证结果
    if j_verification_results:
        meta['j_verification'] = j_verification_results

    # Part 2.E: NOE/J 验证结果回喂 score_candidates, 重新打分 (仅 Top-N)
    # 范围限制: 只对 Top-N 个候选 (有 NOE/J 结果) 重打分, 避免无谓的计算
    if (noe_verification_results or j_verification_results) and evidences:
        try:
            noe_dict = {r['smiles']: r for r in noe_verification_results}
            j_dict = {r['smiles']: r for r in j_verification_results}
            # 取有 NOE/J 验证结果的 SMILES 集合
            verified_smiles = set(noe_dict.keys()) | set(j_dict.keys())
            # 仅重打分有 verified 分数的候选
            hits_to_rescore = []
            for ev in evidences:
                if ev.smiles in verified_smiles:
                    hits_to_rescore.append((ev.rank, ev.smiles))
            if hits_to_rescore:
                # 对全部命中候选重打分 (但只有 verified 部分会拿到 NOE/J bonus)
                # 用 scored[i].hit 重建 RetrievalHit 列表 (scored 与 evidences 等长同序)
                final_scored = score_candidates(
                    [s.hit for s in scored],
                    use_2d_scoring=use_hsqc_mode,
                    use_noe_scoring=bool(noe_verification_results),
                    noe_scores=noe_dict if noe_dict else None,
                    use_j_scoring=bool(j_verification_results),
                    j_scores=j_dict if j_dict else None,
                )
                final_scored.sort(key=lambda s: s.score, reverse=True)
                for i, s in enumerate(final_scored, 1):
                    s.rank = i
                # 按 rank 重新排序 evidences (依据 SMILES)
                rank_map = {s.hit.entry.canonical_smiles or s.hit.entry.smiles: s
                            for s in final_scored}
                evidences.sort(key=lambda e: rank_map.get(e.smiles).rank
                               if rank_map.get(e.smiles) else 999)
                if not args.quiet:
                    print(f"  -> NOE/J 验证回喂: 重打分 {len(evidences)} 个候选, "
                          f"最终 Top-1 score={evidences[0].score:.1f}",
                          file=sys.stderr)
                meta['score_breakdown'] = {
                    'base_used_2d': use_hsqc_mode,
                    'noe_bonus_applied': bool(noe_verification_results),
                    'j_bonus_applied': bool(j_verification_results),
                    'n_rescored': len(hits_to_rescore),
                    'final_top1_score': final_scored[0].score if final_scored else None,
                }
        except Exception as ex:
            print(f"  WARN: NOE/J 回喂打分失败: {ex}", file=sys.stderr)
    # 全归属表 (Top-1)
    if full_assignment is not None:
        meta['full_assignment'] = {
            'structure_smiles': full_assignment.structure_smiles,
            'confidence': full_assignment.confidence,
            'rows': [
                {'atom_label': r.atom_label, 'c_ppm': r.c_ppm,
                 'multiplicity': r.multiplicity, 'h_ppm': r.h_ppm,
                 'hsqc_match': r.hsqc_match,
                 'cosy_partners': r.cosy_partners,
                 'hmbc_to_c': r.hmbc_to_c}
                for r in full_assignment.rows
            ],
            'unassigned_obs_c': full_assignment.unassigned_obs_c,
            'unassigned_obs_h': full_assignment.unassigned_obs_h,
        }
    # 结构图生成 (在 HTML 报告之前, 以便 HTML 内嵌 base64 图片)
    print(f"[viz] 生成 Top-N 结构图...", file=sys.stderr)
    try:
        draw_topn_grid(evidences, out_prefix + '_topn.png', top_n=args.top_n)
    except Exception as e:
        print(f"WARN: 结构图生成失败: {e}", file=sys.stderr)
    # 单分子图: 每个 Top-N 候选一张
    for e in evidences[:args.top_n]:
        try:
            draw_single_mol(
                e.smiles,
                f'{out_prefix}_rank{e.rank}.png',
                legend=f'#{e.rank} {e.smiles[:30]}\nscore={e.score:.1f}',
            )
        except Exception as ex:
            print(f"WARN: rank {e.rank} 结构图失败: {ex}", file=sys.stderr)
    # 实测 vs 候选对比图
    try:
        draw_unmatched_overlay(
            observed_c=problem.carbon_shifts,
            observed_h=problem.proton_shifts,
            evidences=evidences,
            out_path=out_prefix + '_overlay.png',
            top_n=min(5, args.top_n),
        )
    except Exception as e:
        print(f"WARN: 对比图生成失败: {e}", file=sys.stderr)

    # 真实谱线形图 (来自 pdata/1/1r; --with-real-spectrum 开启)
    # 价值: peaklist.xml 是竹签, 1r 是完整线形. 叠加展示真实谱 + 峰检测落点是否吻合.
    if args.with_real_spectrum:
        print(f"[real] 读取 pdata/1/1r 真实谱并绘图...", file=sys.stderr)
        try:
            specs_proc = read_processed_folder(folder, specs=specs)
            if specs_proc:
                observed_by_exp = {int(s.exp_no): [p.shift for p in s.peaks]
                                   for s in specs if str(s.exp_no).isdigit()}
                draw_real_spectrum_overlay(specs_proc, observed_by_exp,
                                           out_prefix + '_real_spectrum.png')
                print(f"  -> 真实谱图已生成 ({len(specs_proc)} 个实验)", file=sys.stderr)
            else:
                print(f"  WARN: 未读到任何 1r 真实谱 (检查 pdata/1/1r 是否存在)", file=sys.stderr)
        except ProcessedReadError as e:
            print(f"WARN: 真实谱读取失败 (stage={e.stage}): {e}", file=sys.stderr)

    # 报告输出 (JSON/Markdown/HTML)
    write_report_json(out_prefix + '.json', evidences, meta)
    write_report_markdown(out_prefix + '.md', evidences, meta, top_n=args.top_n)
    # HTML 报告（人类阅读, 含流程/结果/结构图/未匹配峰/DEPT 推断, 图片 base64 内嵌）
    write_report_html(
        out_prefix + '.html',
        evidences=evidences,
        meta=meta,
        bridge_info=bridge_info,
        top_n=args.top_n,
    )

    # HTML 已内嵌所有图片, 删除独立 PNG 文件 (保持输出目录简洁)
    _cleanup_embedded_pngs(out_prefix, args.top_n, evidences)

    # 未匹配峰汇总 CSV
    try:
        write_unmatched_summary_csv(
            evidences, out_prefix + '_unmatched.csv', top_n=args.top_n,
        )
    except Exception as e:
        print(f"WARN: 未匹配峰 CSV 失败: {e}", file=sys.stderr)

    # nmrpy fid 时域二次处理 + 交叉验证 (--validate-with-nmrpy 开启)
    # 价值: 从 fid 第一性原理重出谱, 独立验证 peaklist.xml 的峰是否可复现.
    nmrpy_ok = False
    if args.validate_with_nmrpy:
        print(f"[nmrpy] fid 时域二次处理 + 交叉验证 peaklist.xml...", file=sys.stderr)
        try:
            from .nmrpy_validate import validate_bruker_with_nmrpy
            validate_bruker_with_nmrpy(folder, out_prefix + '_nmrpy')
            nmrpy_ok = True
        except ImportError as e:
            print(f"WARN: nmrpy 未安装, 跳过验证: {e}", file=sys.stderr)
        except Exception as e:
            # 独立可选步骤: 主报告已产出, 不因 nmrpy 内部异常中断整体, 但显式暴露原因
            print(f"WARN: nmrpy 验证失败: {e}", file=sys.stderr)

    print(f"\nDONE. Top-{args.top_n} 已写入:", file=sys.stderr)
    print(f"  {out_prefix}.json  (机器读取)", file=sys.stderr)
    print(f"  {out_prefix}.md    (Markdown 备查)", file=sys.stderr)
    print(f"  {out_prefix}.html  (人类阅读报告, 含全部图片 base64 内嵌)", file=sys.stderr)
    print(f"  {out_prefix}_unmatched.csv (未匹配峰汇总)", file=sys.stderr)
    if nmrpy_ok:
        print(f"  {out_prefix}_nmrpy_validation.{{json,md,png}} (nmrpy fid 交叉验证)", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())