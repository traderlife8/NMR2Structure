"""NMR2Struct 单谱图推理封装 (调用模板).

完整流程 (用户在本机装好 torch 后可直接跑):
    1. python -m nmr_elucidator.bridge_nmr2struct  (无需 torch)
    2. python -m nmr_elucidator.run_nmr2struct      (需要 torch + 已下载 NMR2Struct)

NMR2Struct 命令行 (单谱图推理):
    nmr_infer_single_spectrum \\
        --config <yaml> \\
        --hnmr_file <hnmr.txt> \\
        --cnmr_file <cnmr.txt> \\
        --hnmr_shifts <HNMR_shifts.p> \\
        --cnmr_shifts <CNMR_shifts.p> \\
        --ckpt <multitask_checkpoint.pt> \\
        --normalize

输入格式: 详 bridge_nmr2struct.py
输出: h5 文件含 15 个 SMILES 候选 + 分数, 解析见 parse_nmr2struct_output.py
"""
from __future__ import annotations

import os
import shutil
import sys
from typing import Optional


NMR2STRUCT_ROOT = r'E:\NMR2Structure\external\NMR2Struct'
EXAMPLE_CONFIGS = os.path.join(NMR2STRUCT_ROOT, 'example_configs')
CHECKPOINTS = os.path.join(NMR2STRUCT_ROOT, 'checkpoints')
MULTITASK_CKPT = os.path.join(CHECKPOINTS, 'multitask_checkpoint.pt')
HNMR_SHIFTS = os.path.join(EXAMPLE_CONFIGS, 'HNMR_shifts.p')
CNMR_SHIFTS = os.path.join(EXAMPLE_CONFIGS, 'CNMR_shifts.p')
ALPHABET = os.path.join(EXAMPLE_CONFIGS, 'alphabet.npy')

# token 约定: alphabet size = 21 -> tgt_pad=21, tgt_start=22, tgt_stop=23
# 28045 = H+C 特征维度 (28000 + 80 + offset)
ALPHABET_SIZE = 21
TGT_START = ALPHABET_SIZE + 1  # 22
TGT_STOP = ALPHABET_SIZE + 2   # 23


def build_inference_config(savedir: str, num_pred_per_tgt: Optional[int] = None) -> str:
    """生成一份新的 yaml 配置文件, 适配单谱图推理.

    复制 example_configs/training_spectrum_to_struct.yaml,
    替换关键字段: 路径、token、savedir 等
    返回 yaml 文件路径
    num_pred_per_tgt: 覆盖默认 15 候选, 调试时设 3-5 可加速
    """
    src = os.path.join(EXAMPLE_CONFIGS, 'training_spectrum_to_struct.yaml')
    if not os.path.exists(src):
        raise FileNotFoundError(f'NMR2Struct example config not found: {src}')

    os.makedirs(savedir, exist_ok=True)
    dst = os.path.join(savedir, 'nmr2struct_inference_config.yaml')

    with open(src, 'r', encoding='utf-8') as f:
        text = f.read()

    # 替换关键字段
    replacements = {
        "savedir: \"checkpoints\"": f"savedir: \"{savedir.replace(chr(92), '/')}\"",
        "tgt_start_token: null": f"tgt_start_token: {TGT_START}",
        "tgt_stop_token: null": f"tgt_stop_token: {TGT_STOP}",
        "alphabet: metadatas/alphabet.npy": f"alphabet: {ALPHABET.replace(chr(92), '/')}",
        "substructure_model_ckpt: null": f"substructure_model_ckpt: \"{os.path.join(CHECKPOINTS, 'cnmr_only_checkpoint.pt').replace(chr(92), '/')}\"",
        "structure_model_ckpt: #Populate": f"structure_model_ckpt: \"{os.path.join(CHECKPOINTS, 'transformer_checkpoint.ckpt').replace(chr(92), '/')}\"",
        "n_hnmr_features: 28000": "n_hnmr_features: 28000",  # 保持
        "n_cnmr_features: 80": "n_cnmr_features: 80",        # 保持
        "source_size: null #Populate": f"source_size: 28045",  # H+C+offset
        "target_size: null #Populate": f"target_size: {ALPHABET_SIZE + 3}",
        "src_pad_token: null #Populate": "src_pad_token: 0",
        "tgt_pad_token: null #Populate": f"tgt_pad_token: {ALPHABET_SIZE}",
    }
    if num_pred_per_tgt is not None:
        replacements["num_pred_per_tgt: 15"] = f"num_pred_per_tgt: {num_pred_per_tgt}"
    for old, new in replacements.items():
        text = text.replace(old, new)

    with open(dst, 'w', encoding='utf-8') as f:
        f.write(text)
    return dst


def patch_model_args_for_strict_load(model_args: dict) -> None:
    """在 create_model 之前把 structure_model_ckpt / substructure_model_ckpt 置 None,
    避免 partial load 阶段打印大量 key 不匹配的 warning (因为 multitask ckpt 才是真权重).

    修改 dict 原地.
    """
    ma = model_args.get('model_args', {})
    ma['structure_model_ckpt'] = None
    ma['substructure_model_ckpt'] = None


def run_nmr2struct_inference(
    hnmr_file: str,
    cnmr_file: str,
    savedir: str,
    ckpt: Optional[str] = None,
    num_pred_per_tgt: Optional[int] = None,
    use_nmrsim: bool = False,  # 占位, NMR2Struct 不用 nmrsim
) -> str:
    """调用 NMR2Struct 库做单谱图推理.

    hnmr_file: 由 bridge_nmr2struct.convert_bruker_folder_to_nmr2struct 产生
    cnmr_file: 同上
    savedir: NMR2Struct 会把 h5 输出写到 global_args.savedir
    ckpt: 默认 multitask_checkpoint.pt
    num_pred_per_tgt: 覆盖 yaml 中的候选数, 默认 15; 调试时可设 3-5 加速

    返回 h5 输出文件路径
    """
    if not os.path.exists(NMR2STRUCT_ROOT):
        raise FileNotFoundError(f'NMR2Struct 未找到: {NMR2STRUCT_ROOT}')
    if ckpt is None:
        ckpt = MULTITASK_CKPT
    if not os.path.exists(ckpt):
        raise FileNotFoundError(f'权重未找到: {ckpt}')

    # 1. 生成 yaml 配置
    config_path = build_inference_config(savedir)
    print(f'生成配置: {config_path}')

    # 2. 把 NMR2Struct 加入 sys.path
    if NMR2STRUCT_ROOT not in sys.path:
        sys.path.insert(0, NMR2STRUCT_ROOT)

    # 3. 调用 single_spectrum 入口
    from nmr.scripts import run_nmr_inference_single_spectrum as rni
    # patch sys.argv 让 argparse 读到正确参数
    sys.argv = [
        'run_nmr_inference_single_spectrum.py',
        '--config', config_path,
        '--hnmr_file', hnmr_file,
        '--cnmr_file', cnmr_file,
        '--hnmr_shifts', HNMR_SHIFTS,
        '--cnmr_shifts', CNMR_SHIFTS,
        '--ckpt', ckpt,
        '--normalize',
    ]
    print('开始推理...')
    rni.main()
    print('推理完成')

    # 4. 返回 h5 输出路径
    out_h5 = os.path.join(savedir, 'predictions_dataset_0_0.h5')
    if not os.path.exists(out_h5):
        # 有些版本会写其他名字, 列目录找
        cands = [f for f in os.listdir(savedir) if f.startswith('predictions_') and f.endswith('.h5')]
        if cands:
            out_h5 = os.path.join(savedir, cands[0])
        else:
            raise FileNotFoundError(f'推理完成但未找到 h5 输出 in {savedir}')
    return out_h5
