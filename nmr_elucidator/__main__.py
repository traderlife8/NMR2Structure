"""CLI 入口: 端到端运行 NMR 反推结构.

用法:
    python -m nmr_elucidator --folder <原始数据/<sample>> --sdf <SDF path>
                             [--formula <分子式>] [--top-n 10]
                             [--output-prefix <prefix>]
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

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

# 项目 logger (在需要时调用 setup_logging() 初始化)
from .log_setup import setup_logging as _setup_logging
from .log_setup import get_logger as _get_logger


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
    """保留强度 ≥ 第 N 百分位的峰 (去噪声).

    shifts 和 intensities 等长. 如果峰数 ≤ 5 则不过滤.
    """
    if not shifts or not intensities or len(shifts) <= 5:
        return list(shifts)
    paired = sorted(zip(intensities, shifts), reverse=True)
    n = len(paired)
    cut = max(1, int(n * percentile))
    thresh = paired[cut - 1][0]
    return sorted(p for i, p in paired if i >= thresh)


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

    # 强度过滤: 用 20% percentile 阈值去掉噪声峰 (保留真实信号)
    c_shifts = _filter_peaks_by_intensity(c_shifts, [p.intensity for p in c_spec.peaks]) if c_spec else []
    h_shifts = _filter_peaks_by_intensity(h_shifts, [p.intensity for p in h_spec.peaks]) if h_spec else []

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
    parser.add_argument('--sdf', required=True,
                        help='nmrshiftdb2 SDF 文件路径')
    parser.add_argument('--formula', default=None,
                        help='可选：硬约束分子式, e.g. C8H8O2')
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
    parser.add_argument('--use-pubchem', action='store_true',
                        help='启用 PubChem 扩展检索 (按分子式/子结构)')
    parser.add_argument('--use-fragment-assembly', action='store_true',
                        help='启用片段组装候选生成 (库外新骨架)')
    parser.add_argument('--use-2d', action='store_true',
                        help='启用 2D NMR 关联验证 (HSQC/HMBC/COSY); '
                             '无 2D 数据时自动回退到 1D 流程')
    parser.add_argument('--determine', action='store_true',
                        help='输出结构推断结论 (含置信度分级)')
    parser.add_argument('--hose-cache', default=None,
                        help='HOSE 查找表缓存文件路径 (避免每次重建)')
    parser.add_argument('--hose-max-records', type=int, default=50000,
                        help='从 SDF 构建 HOSE 查找表的最大记录数 (默认 50000; '
                             '0=处理全部; 10000~50000 够用, 100000+ 需约 5 分钟)')
    parser.add_argument('--lucy-ng-db', default=None,
                        help='lucy-ng 完整 HOSE 数据库路径 (lucy-ng-derep.db); '
                             '启用后优先用该库查询, 精度更高 (MAE ~3.5 ppm). '
                             '默认查找 E:\\NMR2Structure\\data\\lucy-ng-derep_db\\lucy-ng-derep.db')
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
    args = parser.parse_args(argv)

    # 初始化日志 (默认开启, --no-log 可禁用)
    if not args.no_log:
        log_file = _setup_logging(log_dir=args.log_dir)
        log = _get_logger('main')
        log.info(f'日志文件: {log_file}')

    folder = os.path.abspath(args.folder)
    sdf_path = os.path.abspath(args.sdf)
    if not os.path.isdir(folder):
        print(f"ERROR: 原始数据目录不存在: {folder}", file=sys.stderr)
        return 2
    if not os.path.isfile(sdf_path):
        print(f"ERROR: SDF 文件不存在: {sdf_path}", file=sys.stderr)
        return 2

    # 推断化合物名 (优先 --compound-name, 否则取 folder 路径最后一段)
    compound_name = args.compound_name or os.path.basename(folder.rstrip(r'\/'))
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

    # 3. 检索
    use_hsqc_mode = hsqc_pairs_for_retrieve is not None
    mode_str = 'HSQC 配对匹配' if use_hsqc_mode else '1D 位移匹配'
    print(f"[3] 检索候选 ({mode_str}, SDF={os.path.basename(sdf_path)})", file=sys.stderr)
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
    print(f"  -> 命中 {len(hits)} 条", file=sys.stderr)

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
            lucy_provider = get_lucy_ng_provider(args.lucy_ng_db, quiet=args.quiet)
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

    # 7. PubChem 扩展检索（可选）
    pubchem_candidates = []
    if args.use_pubchem:
        print(f"[7] PubChem 扩展检索 ...", file=sys.stderr)
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
    else:
        print(f"[7] 跳过 PubChem 扩展检索", file=sys.stderr)

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
        'used_hose_source': 'lucy-ng-derep.db' if (args.use_hose and args.lucy_ng_db) else 'local-sdf',
        'used_pubchem': bool(args.use_pubchem),
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