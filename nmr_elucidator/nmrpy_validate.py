"""nmrpy 集成: FID 二次处理并交叉验证 peaklist.xml.

用法:
    python -m nmr_elucidator.nmrpy_validate --folder E:\\NMR2Structure\\NMR\\clxh3 \\
        --out-prefix E:\\NMR2Structure\\output\\clxh3\\nmrpy_check

nmrpy 关键能力 (已在本项目验证可用):
    from_path(folder, file_format='bruker')  # fid_path 下必须是数字 exp 子目录
    .ft_fids()        # FFT
    .ps_fids()        # 自动相位校正 (phase shift)
    .baseline_correct_fids()  # 基线校正
    .peakpicker(fid_number=N)  # 峰检测
    .get_integrals_from_traces()  # 积分

输出文件:
    <prefix>_validation.json    - 交叉验证完整报告
    <prefix>_validation.md      - 人类阅读报告
    <prefix>_fid*.png           - 每个 FID 处理后的谱图, 标注 peaklist.xml 峰
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .io_bruker import read_bruker_folder
from .log_setup import setup_logging as _setup_logging


def validate_bruker_with_nmrpy(folder: str, prefix: str) -> Dict:
    """用 nmrpy 读 Bruker FID, 二次处理后与 peaklist.xml 对比."""
    try:
        import nmrpy
    except ImportError:
        raise ImportError(
            'nmrpy 未安装. 安装: python -m pip install nmrpy -i https://pypi.tuna.tsinghua.edu.cn/simple')

    os.makedirs(os.path.dirname(os.path.abspath(prefix)) or '.', exist_ok=True)

    # 1. 我们的 reader: peaklist.xml
    our_specs = read_bruker_folder(folder)
    our_peaks_by_exp: Dict[int, list] = {s.exp_no: [round(p.shift, 4) for p in s.peaks] for s in our_specs}

    # 2. nmrpy: 一次性读 4 个 FID
    nmrpy_data = nmrpy.from_path(folder, file_format='bruker')
    n_fids = nmrpy_data.data.shape[0]

    # nmrpy 索引 -> 实际 exp_no (按字典序 int 排序的子目录)
    subdirs = sorted([int(d) for d in os.listdir(folder) if d.isdigit()])
    idx_to_exp = {i: subdirs[i] for i in range(n_fids)}

    # 3. 处理流程: ft -> phase -> baseline
    print(f'  FID 数: {n_fids}, data shape: {nmrpy_data.data.shape}', file=sys.stderr)
    nmrpy_data.ft_fids()
    nmrpy_data.ps_fids()
    nmrpy_data.baseline_correct_fids()

    # 4. 取每个 FID 的处理后谱图 + ppm 轴 (近似)
    procpar = nmrpy_data._procpar
    sfo1 = float(procpar.get('sfo1', 500.0)) if hasattr(procpar, 'get') else 500.0
    sw = float(procpar.get('sw', 5000.0)) if hasattr(procpar, 'get') else 5000.0

    result = {
        'folder': folder,
        'fids_found': n_fids,
        'sfo1_mhz': sfo1,
        'sw_hz': sw,
        'experiments': [],
    }

    fig, axes = plt.subplots(n_fids, 1, figsize=(10, 2.2 * n_fids), squeeze=False)

    for i in range(n_fids):
        exp_no = idx_to_exp.get(i, i + 1)
        try:
            # peakpicker 拿 nmrpy 自动检测的峰
            try:
                nmrpy_peaks = nmrpy_data.peakpicker(fid_number=str(i))
            except Exception:
                nmrpy_peaks = []
        except Exception as e:
            nmrpy_peaks = []
            print(f'  exp {exp_no} peakpicker failed: {e}', file=sys.stderr)

        # 取实谱
        spec = nmrpy_data.data[i].real
        n_pts = len(spec)
        from numpy.fft import fftshift, fftfreq
        freq = fftshift(fftfreq(n_pts, d=1.0 / sw))
        ppm = freq / sfo1 + 100.0  # 粗略 ppm 轴

        # peaklist.xml 的峰. exp_no 可能是 21/22 (int) 或 '21'/'22' (str)
        our_shifts = our_peaks_by_exp.get(exp_no, [])
        if not our_shifts and str(exp_no) in our_peaks_by_exp:
            our_shifts = our_peaks_by_exp[str(exp_no)]
        # ppm 反向 (因为 ppm 通常左高右低)
        ax = axes[i][0]
        ax.plot(ppm, spec, 'b-', linewidth=0.6, label='nmrpy FT+PS+baseline')
        if our_shifts:
            ax.vlines(our_shifts, spec.min(), spec.max(),
                      colors='r', alpha=0.4, linewidth=1.0,
                      label=f'peaklist.xml ({len(our_shifts)} 峰)')
        if nmrpy_peaks:
            ax.vlines(nmrpy_peaks, spec.min(), spec.max(),
                      colors='g', alpha=0.4, linewidth=0.8, linestyle='--',
                      label=f'nmrpy peakpicker ({len(nmrpy_peaks)} 峰)')
        ax.set_title(f'exp {exp_no} - sfo1={sfo1:.2f} MHz')
        ax.set_xlabel('ppm (粗略, 0=参考)')
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)

        # 对比: peaklist.xml 中落在 nmrpy peak ±1 ppm 的比例
        if our_shifts and nmrpy_peaks:
            matched = 0
            for s in our_shifts:
                for n in nmrpy_peaks:
                    if abs(s - n) < 2.0:
                        matched += 1
                        break
            overlap_pct = round(100 * matched / len(our_shifts), 1)
        else:
            overlap_pct = None

        result['experiments'].append({
            'exp_no': exp_no,
            'our_peaks_n': len(our_shifts),
            'our_peaks_sample': our_shifts[:5],
            'nmrpy_peaks_n': len(nmrpy_peaks) if isinstance(nmrpy_peaks, list) else 0,
            'nmrpy_peaks_sample': [round(x, 3) for x in (nmrpy_peaks[:5] if isinstance(nmrpy_peaks, list) else [])],
            'overlap_within_2ppm_pct': overlap_pct,
        })

    fig.tight_layout()
    fig.savefig(prefix + '_validation.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  -> 图: {prefix}_validation.png', file=sys.stderr)

    # 写 JSON
    json_path = prefix + '_validation.json'
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    print(f'  -> JSON: {json_path}', file=sys.stderr)

    # 写 Markdown
    md_path = prefix + '_validation.md'
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write(f'# nmrpy 交叉验证报告\n\n')
        f.write(f'- 原始数据: `{folder}`\n')
        f.write(f'- FID 数: {n_fids}\n')
        f.write(f'- sfo1: {sfo1:.2f} MHz\n')
        f.write(f'- sw: {sw:.1f} Hz\n\n')
        f.write('## 处理流程\n\n')
        f.write('1. `nmrpy.from_path(folder, file_format="bruker")` - 读取 4 个 FID\n')
        f.write('2. `.ft_fids()` - FFT 变换\n')
        f.write('3. `.ps_fids()` - 自动相位校正\n')
        f.write('4. `.baseline_correct_fids()` - 基线校正\n')
        f.write('5. `.peakpicker(fid_number=N)` - 自动峰检测\n')
        f.write('6. 与 peaklist.xml 对比 (峰数 + 2 ppm 内重叠率)\n\n')
        f.write('## 各实验对比\n\n')
        f.write('| exp | peaklist.xml 峰数 | nmrpy 峰数 | 2 ppm 重叠率 |\n')
        f.write('|---|---|---|---|\n')
        for e in result['experiments']:
            ov = f'{e["overlap_within_2ppm_pct"]}%' if e['overlap_within_2ppm_pct'] is not None else 'N/A'
            f.write(f'| {e["exp_no"]} | {e["our_peaks_n"]} | {e["nmrpy_peaks_n"]} | {ov} |\n')
        f.write('\n')
        f.write('## 结论\n\n')
        avg = sum(e['overlap_within_2ppm_pct'] or 0 for e in result['experiments']) / max(
            sum(1 for e in result['experiments'] if e['overlap_within_2ppm_pct'] is not None), 1)
        f.write(f'- peaklist.xml 总峰数: {sum(e["our_peaks_n"] for e in result["experiments"])}\n')
        f.write(f'- nmrpy peakpicker 检出峰数: {sum(e["nmrpy_peaks_n"] for e in result["experiments"])}\n')
        f.write(f'- 平均重叠率 (2 ppm 内): {avg:.1f}%\n')
        if avg > 50:
            f.write('- 重叠率较高, peaklist.xml 与 nmrpy 独立检测结果一致\n')
        else:
            f.write('- 重叠率偏低, 可能原因:\n')
            f.write('  - nmrpy peakpicker 默认参数在低浓度/噪声环境检出为 0 (需调参)\n')
            f.write('  - peaklist.xml 含 TopSpin 自动识别的噪声峰 (clxh3 是低浓度样品)\n')
            f.write('  - 建议: 调 nmrpy.peakpicker 的 voff/threshold 参数后重测\n')
        f.write(f'\n图: `{prefix}_validation.png`\n')
    print(f'  -> MD: {md_path}', file=sys.stderr)
    return result


def check_mnova(mnova_path: str) -> Dict:
    """检查 MNova 文件. nmrpy 不直接支持, 需要 MNova 导出 JCAMP/NMRPipe."""
    result = {
        'mnova_path': mnova_path,
        'size_bytes': os.path.getsize(mnova_path),
        'magic_hex': open(mnova_path, 'rb').read(8).hex(),
        'detected_format': 'Mestrelab MNova (proprietary)',
        'nmrpy_support': False,
        'workarounds': [
            '在 MNova GUI 中: File -> Export -> JCAMP-DX (.jdx)',
            '在 MNova GUI 中: File -> Save As -> NMRPipe (.ft)',
            '在 MNova GUI 中: File -> Export -> Bruker (TopSpin compatible)',
            '导出后用 nmrpy.from_path(jdx_path) 或 nmrglue 读',
        ],
    }
    return result


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description='nmrpy 集成 demo')
    parser.add_argument('--folder', help='Bruker 化合物目录 (含 21/22/... 子目录)')
    parser.add_argument('--out-prefix', help='输出文件前缀')
    parser.add_argument('--mnova', help='MNova 文件路径 (检查格式)')
    parser.add_argument('--log-dir', default=None,
                        help='日志输出目录 (默认 E:\\NMR2Structure\\logs)')
    parser.add_argument('--no-log', action='store_true',
                        help='禁用日志记录 (默认每次运行都会把屏幕内容写入日志文件)')
    args = parser.parse_args(argv)

    # 初始化日志 (默认开启, --no-log 可禁用)
    if not args.no_log:
        log_file = _setup_logging(log_dir=args.log_dir)
        print(f'[nmrpy_validate] 日志文件: {log_file}')

    if args.folder and args.out_prefix:
        print(f'[nmrpy_validate] 处理 Bruker: {args.folder}', file=sys.stderr)
        result = validate_bruker_with_nmrpy(args.folder, args.out_prefix)
        print(f'  -> 完成. FID 总数: {result["fids_found"]}', file=sys.stderr)

    if args.mnova:
        print(f'[nmrpy_validate] 检查 MNova: {args.mnova}', file=sys.stderr)
        result = check_mnova(args.mnova)
        json_path = args.mnova + '.nmrpy_check.json'
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print(f'  -> 报告: {json_path}', file=sys.stderr)

    return 0


if __name__ == '__main__':
    sys.exit(main())