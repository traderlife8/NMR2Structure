"""CLI 入口: 端到端运行 NMR 反推结构.

用法:
    python -m nmr_elucidator --folder <原始数据/<sample>> --sdf <SDF path>
                             [--formula <分子式>] [--top-n 10]
                             [--output-prefix <prefix>]
"""
from __future__ import annotations

import argparse
import os
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


def _build_problem_from_bruker(specs) -> ProblemPeaks:
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

    return ProblemPeaks(
        carbon_shifts=c_shifts,
        proton_shifts=h_shifts,
        proton_integrals=h_integrals,
        solvent=solvent,
        field_mhz=h_spec.field_mhz if h_spec else (c_spec.field_mhz if c_spec else None),
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
    parser.add_argument('--quiet', action='store_true',
                        help='不打印检索进度')
    args = parser.parse_args(argv)

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

    # 2. 构造问题
    problem = _build_problem_from_bruker(specs)
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

    # 3. 检索
    print(f"[3/5] 检索候选 (SDF={sdf_path})", file=sys.stderr)
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
    scored = score_candidates(hits)
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
        'dept_counts': dept_counts,
        'observed_carbon_shifts': problem.carbon_shifts,
        'observed_proton_shifts': problem.proton_shifts,
        'obs_dept_summary': obs_dept_summary,
    }
    write_report_json(out_prefix + '.json', evidences, meta)
    write_report_markdown(out_prefix + '.md', evidences, meta, top_n=args.top_n)
    # HTML 报告（人类阅读, 含流程/结果/结构图/未匹配峰/DEPT 推断）
    write_report_html(
        out_prefix + '.html',
        evidences=evidences,
        meta=meta,
        bridge_info=bridge_info,
        top_n=args.top_n,
    )

    # 结构图与未匹配峰汇总
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
    # 未匹配峰汇总 CSV
    try:
        write_unmatched_summary_csv(
            evidences, out_prefix + '_unmatched.csv', top_n=args.top_n,
        )
    except Exception as e:
        print(f"WARN: 未匹配峰 CSV 失败: {e}", file=sys.stderr)

    # 真实谱线形图 (来自 pdata/1/1r; --with-real-spectrum 开启)
    # 价值: peaklist.xml 是竹签, 1r 是完整线形. 叠加展示真实谱 + 峰检测落点是否吻合.
    real_spectrum_png = None
    if args.with_real_spectrum:
        print(f"[real] 读取 pdata/1/1r 真实谱并绘图...", file=sys.stderr)
        try:
            specs_proc = read_processed_folder(folder, specs=specs)
            if specs_proc:
                # io_bruker 的 exp_no 是 str, 统一转 int 匹配 io_processed 的数字 exp 键
                observed_by_exp = {int(s.exp_no): [p.shift for p in s.peaks]
                                   for s in specs if str(s.exp_no).isdigit()}
                real_spectrum_png = out_prefix + '_real_spectrum.png'
                draw_real_spectrum_overlay(specs_proc, observed_by_exp, real_spectrum_png)
                print(f"  -> {real_spectrum_png} ({len(specs_proc)} 个实验)", file=sys.stderr)
            else:
                real_spectrum_png = None
                print(f"  WARN: 未读到任何 1r 真实谱 (检查 pdata/1/1r 是否存在)", file=sys.stderr)
        except ProcessedReadError as e:
            real_spectrum_png = None
            print(f"WARN: 真实谱读取失败 (stage={e.stage}): {e}", file=sys.stderr)

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
    print(f"  {out_prefix}.html  (人类阅读报告, 含流程+结构图)", file=sys.stderr)
    print(f"  {out_prefix}_topn.png (候选网格)", file=sys.stderr)
    print(f"  {out_prefix}_rank*.png (单分子)", file=sys.stderr)
    print(f"  {out_prefix}_overlay.png (实测 vs 候选对比)", file=sys.stderr)
    print(f"  {out_prefix}_unmatched.csv (未匹配峰汇总)", file=sys.stderr)
    if real_spectrum_png:
        print(f"  {real_spectrum_png} (1r 真实谱线形)", file=sys.stderr)
    if nmrpy_ok:
        print(f"  {out_prefix}_nmrpy_validation.{{json,md,png}} (nmrpy fid 交叉验证)", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())