"""一键把 Bruker 实测谱图通过 NMR2Struct 反推结构.

用法:
    # 阶段 1: 仅生成输入文件 (不需要 torch)
    python -m nmr_elucidator.run_bridge --folder <bruker_dir> --out-prefix <output_dir> --stage 1

    # 阶段 2: 调用 NMR2Struct 推理 (需要 torch)
    python -m nmr_elucidator.run_bridge --folder <bruker_dir> --out-prefix <output_dir> --stage 2

    # 阶段 3: 解析 h5 输出为 Markdown (不需要 torch)
    python -m nmr_elucidator.run_bridge --folder <bruker_dir> --out-prefix <output_dir> --stage 3

    # 阶段 4: 阶段 1+2+3 全流程
    python -m nmr_elucidator.run_bridge --folder <bruker_dir> --out-prefix <output_dir> --stage all
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Optional


def stage1(folder: str, prefix: str) -> dict:
    """把 Bruker 转成 NMR2Struct 期望的 hnmr/cnmr .txt."""
    from .bridge_nmr2struct import convert_bruker_folder_to_nmr2struct

    h_out = prefix + '_hnmr.txt'
    c_out = prefix + '_cnmr.txt'
    info = convert_bruker_folder_to_nmr2struct(folder, h_out, c_out)
    print(f'[stage1] H 谱写出 {info["n_1h_written"]} 峰 -> {info["hnmr_file"]}')
    print(f'[stage1] C 谱写出 {info["n_13c_written"]} 峰 -> {info["cnmr_file"]}')
    print(f'[stage1] 溶剂: {info["solvent"]}  频率: {info["field_mhz"]:.4f} MHz')
    return info


def stage2(hnmr_file: str, cnmr_file: str, savedir: str, num_pred: int = 15, ckpt: str = None) -> str:
    """调 NMR2Struct 推理. 返回 h5 路径.

    num_pred: 候选 SMILES 数量 (默认 15, 调试时设 3-5 加速)
    """
    try:
        import torch  # noqa: F401
    except ImportError:
        print('[stage2] torch 未安装, 跳过推理. 安装命令:')
        print('         python -m pip install torch -i https://pypi.tuna.tsinghua.edu.cn/simple')
        raise
    from .run_nmr2struct import run_nmr2struct_inference
    out_h5 = run_nmr2struct_inference(hnmr_file, cnmr_file, savedir, ckpt=ckpt, num_pred_per_tgt=num_pred)
    print(f'[stage2] 推理 h5 输出: {out_h5}')
    return out_h5


def stage3(h5_path: str, md_path: str) -> list:
    """把 h5 解析为 Markdown 报告."""
    from .parse_nmr2struct_output import parse_predictions_h5, write_candidates_markdown
    cands = parse_predictions_h5(h5_path)
    write_candidates_markdown(cands, md_path)
    print(f'[stage3] 共 {len(cands)} 个 example, 写入 {md_path}')
    return cands


def main() -> int:
    ap = argparse.ArgumentParser(description='通过 NMR2Struct 反推分子结构')
    ap.add_argument('--folder', required=True, help='Bruker 原始数据目录')
    ap.add_argument('--out-prefix', required=True, help='输出文件前缀 (不含扩展名), '
                    '例如 E:\\NMR2Structure\\output\\clxh3\\clxh3')
    ap.add_argument('--stage', default='all', choices=['1', '2', '3', 'all'],
                    help='运行阶段: 1=生成输入, 2=推理, 3=解析, all=全部')
    ap.add_argument('--h5', default=None, help='阶段 3 单独用: h5 文件路径')
    ap.add_argument('--num-pred', type=int, default=15, help='候选 SMILES 数量 (默认 15, 调试设 3-5)')
    ap.add_argument('--ckpt', default=None, help='自定义权重路径, 默认 multitask_checkpoint.pt')
    args = ap.parse_args()

    folder = os.path.abspath(args.folder)
    prefix = os.path.abspath(args.out_prefix)
    os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
    savedir = prefix + '_nmr2struct_savedir'

    info: Optional[dict] = None
    h5_path: Optional[str] = None

    if args.stage in ('1', 'all'):
        info = stage1(folder, prefix)

    if args.stage in ('2', 'all'):
        # 单跑 stage 2 时复用之前 stage 1 写出的文件
        if info is None:
            h_path = prefix + '_hnmr.txt'
            c_path = prefix + '_cnmr.txt'
            if not (os.path.exists(h_path) and os.path.exists(c_path)):
                raise FileNotFoundError(
                    f'stage 2 需要先有输入文件: {h_path} 和 {c_path} 不存在, 请先跑 stage 1 或 all')
            info = {
                'hnmr_file': h_path,
                'cnmr_file': c_path,
                'n_1h_written': -1,
                'n_13c_written': -1,
            }
        h5_path = stage2(info['hnmr_file'], info['cnmr_file'], savedir,
                         num_pred=args.num_pred, ckpt=args.ckpt)

    if args.stage in ('3', 'all'):
        h5_path = h5_path or args.h5 or (prefix + '_nmr2struct_savedir' + '/predictions_dataset_0_0.h5')
        if not os.path.exists(h5_path):
            # 在 savedir 找
            cand_dir = os.path.dirname(h5_path)
            if os.path.isdir(cand_dir):
                hits = [f for f in os.listdir(cand_dir) if f.startswith('predictions_') and f.endswith('.h5')]
                if hits:
                    h5_path = os.path.join(cand_dir, hits[0])
                    print(f'[stage3] 找到 h5: {h5_path}')
        md_path = prefix + '_nmr2struct.md'
        stage3(h5_path, md_path)

    print('\nDONE.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
