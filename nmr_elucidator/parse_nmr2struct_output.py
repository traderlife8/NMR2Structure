"""解析 NMR2Struct 单谱图推理输出 (h5) -> 候选 SMILES 列表 + 分数.

NMR2Struct 输出的 h5 结构 (从 inference_fxns.infer_transformer_model 推断):
  'test' (group, 来自 saved_inference_predictions 的 savetag='dataset_0')
    - smiles_or_targets: (num_examples,) 字符串数组
    - num_pred_per_tgt: 数字 15
    - '0' (sub-group, 每个 example 一个)
        - predictions: (15,) 字符串数组, 候选 SMILES
        - scores: (15,) 数组, log 概率 (越高越好)

实际 h5 结构还需打开看; 此模块兼容两种布局.
"""
from __future__ import annotations

import os
from typing import List, Optional, Tuple

import h5py
import numpy as np


def parse_predictions_h5(h5_path: str) -> List[dict]:
    """读 h5 推理结果, 返回每个 example 的候选列表.

    返回格式:
    [
        {
            'index': 0,
            'target': '...' (如果有),
            'candidates': [
                {'rank': 1, 'smiles': 'CCO', 'score': -3.45},
                ...
            ]
        }
    ]
    """
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f'h5 输出未找到: {h5_path}')
    out: list[dict] = []
    with h5py.File(h5_path, 'r') as f:
        print(f'h5 keys (top-level): {list(f.keys())}')
        # 找 test / val / train / 数字索引 group
        top_keys = list(f.keys())
        # 找包含 'predictions' 或 'scores' 的 group
        for k in top_keys:
            g = f[k]
            print(f'  [{k}] type={type(g).__name__}, sub-keys={list(g.keys())[:10]}')
            if isinstance(g, h5py.Group):
                # 单个 example
                preds, scores = _extract_predictions_from_group(g)
                if preds is not None:
                    out.append({
                        'index': _try_int(k),
                        'target': None,
                        'candidates': _rank_candidates(preds, scores),
                    })
    return out


def _extract_predictions_from_group(g) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """从一个 group 提取 predictions + scores 数组."""
    keys = set(g.keys())
    # 直接在 group
    if 'predictions' in keys or 'prediction_strings' in keys:
        preds = g.get('predictions') or g.get('prediction_strings')
        scores = g.get('scores') or g.get('prediction_scores')
        if isinstance(preds, h5py.Dataset):
            preds = preds[()]
        if isinstance(scores, h5py.Dataset):
            scores = scores[()]
        return preds, scores
    # 在子 group '0' '1' ...
    for sub in sorted(g.keys()):
        try:
            sub_g = g[sub]
        except Exception:
            continue
        if isinstance(sub_g, h5py.Group):
            sub_keys = set(sub_g.keys())
            if 'predictions' in sub_keys or 'prediction_strings' in sub_keys:
                preds = sub_g.get('predictions') or sub_g.get('prediction_strings')
                scores = sub_g.get('scores') or sub_g.get('prediction_scores')
                if isinstance(preds, h5py.Dataset):
                    preds = preds[()]
                if isinstance(scores, h5py.Dataset):
                    scores = scores[()]
                return preds, scores
    return None, None


def _flatten_2d(arr):
    """如果 arr 是 (1, N) 2D, 摊成 (N,) 1D; 否则保持."""
    import numpy as np
    if hasattr(arr, 'shape') and len(arr.shape) == 2 and arr.shape[0] == 1:
        return arr[0]
    return arr


def _try_int(s) -> int:
    try:
        return int(s)
    except (ValueError, TypeError):
        return -1


def _rank_candidates(preds, scores) -> list[dict]:
    """把 (smiles_array, score_array) 排序成 [{rank, smiles, score}, ...]."""
    if preds is None:
        return []
    # 摊平 (1, N) -> (N,)
    preds = _flatten_2d(preds)
    if scores is not None:
        scores = _flatten_2d(scores)
    # decode bytes -> str
    smiles_list = []
    for p in preds:
        if isinstance(p, bytes):
            smiles_list.append(p.decode('utf-8', errors='replace'))
        else:
            smiles_list.append(str(p))
    if scores is None:
        scores = np.zeros(len(smiles_list))
    else:
        scores = np.asarray(scores)
    # 按分数降序 (log prob 越大越好)
    order = np.argsort(-scores)
    out = []
    for i, idx in enumerate(order):
        out.append({
            'rank': i + 1,
            'smiles': smiles_list[idx],
            'score': float(scores[idx]),
        })
    return out


def write_candidates_markdown(candidates_list: list[dict], out_path: str) -> None:
    """把候选 SMILES 列表写为 Markdown 报告, 含结构图 (RDKit)."""
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    # 中文字体 (与 drawing.py 同步)
    _ZH = [('SimHei', r'C:\Windows\Fonts\simhei.ttf'),
           ('Microsoft YaHei', r'C:\Windows\Fonts\msyh.ttc')]
    for _name, _path in _ZH:
        if os.path.exists(_path):
            try:
                font_manager.fontManager.addfont(_path)
            except Exception:
                pass
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or '.', exist_ok=True)
    lines = ['# NMR2Struct 推理候选\n']
    img_dir = os.path.splitext(out_path)[0] + '_imgs'
    os.makedirs(img_dir, exist_ok=True)

    for ex in candidates_list:
        idx = ex.get('index', '?')
        lines.append(f'## Example {idx}\n')
        cands = ex.get('candidates', [])
        for c in cands:
            smi = c['smiles']
            mol = Chem.MolFromSmiles(smi)
            valid = mol is not None
            c['valid'] = valid
            c['formula'] = rdMolDescriptors.CalcMolFormula(mol) if mol else None
            c['mw'] = Descriptors.MolWt(mol) if mol else None
            img_path = os.path.join(img_dir, f'ex{idx}_rank{c["rank"]}.png')
            if mol is not None:
                try:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    img.save(img_path)
                except Exception:
                    img_path = None
            lines.append(f'### Rank {c["rank"]}')
            lines.append(f'- SMILES: `{smi}`')
            if valid:
                lines.append(f'- 分子式: {c["formula"]}  MW: {c["mw"]:.2f}')
                if img_path and os.path.exists(img_path):
                    rel = os.path.relpath(img_path, os.path.dirname(out_path))
                    lines.append(f'- 结构图: ![]({rel})')
            else:
                lines.append(f'- (解析失败, 跳过)')
            lines.append(f'- Score: {c["score"]:.4f}\n')

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'报告: {out_path}')
    print(f'图片: {img_dir}/')
