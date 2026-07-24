"""COCONUT 候选检索模块.

COCONUT (COlleCtion of Open Natural prodUcTs) 是最大的开放天然产物数据库之一,
收录超过 69.5 万 NP 结构 (整合自 53 个公开数据源).

通过 REST API 查询, 无需 API Key, 支持:
  - 按分子式搜索
  - 按 SMILES 子结构/相似度搜索
  - 按名字/InChI/InChIKey 搜索

API 文档: https://coconut.naturalproducts.net/api-documentation

设计要点:
  - 仅结构信息, 无 NMR 位移数据
  - 检索后用 HOSE 预测验证
  - 网络失败时静默降级, 不阻断主流程
  - 限速: 避免高频请求 (单次最多 100 条)
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem


COCONUT_API_BASE = 'https://coconut.naturalproducts.net'


@dataclass
class CoconutCandidate:
    """COCONUT 检索到的候选."""

    coconut_id: str
    smiles: str
    name: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: float = 0.0
    heavy_atom_count: int = 0
    # 验证分数 (由 nmr_predict 填充)
    verification_score: float = 0.0
    source: str = 'coconut'


@dataclass
class CoconutSearchResult:
    """COCONUT 检索结果."""

    candidates: List[CoconutCandidate] = field(default_factory=list)
    query_formula: Optional[str] = None
    query_smiles: Optional[str] = None
    total_found: int = 0
    error: Optional[str] = None
    notes: List[str] = field(default_factory=list)


def _coconut_api_get(endpoint: str, params: dict, timeout: float = 30.0) -> Optional[dict]:
    """调用 COCONUT REST API GET 端点.

    Args:
        endpoint: API 路径 (如 '/api/search/formula')
        params: 查询参数
        timeout: 超时秒数

    Returns:
        解析后的 JSON dict, 失败返回 None
    """
    url = f"{COCONUT_API_BASE}{endpoint}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={
        'Accept': 'application/json',
        'User-Agent': 'NMR2Structure/1.0',
    })
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
            return json.loads(data.decode('utf-8'))
    except urllib.error.HTTPError as e:
        return {'_error': f'HTTP {e.code}: {e.reason}'}
    except urllib.error.URLError as e:
        return {'_error': f'URL 错误: {e.reason}'}
    except json.JSONDecodeError as e:
        return {'_error': f'JSON 解析失败: {e}'}
    except Exception as e:
        return {'_error': f'请求失败: {e}'}


def search_coconut_by_formula(
    formula: str,
    max_results: int = 100,
    min_heavy_atoms: int = 3,
    max_heavy_atoms: int = 50,
    quiet: bool = False,
) -> CoconutSearchResult:
    """按分子式从 COCONUT 检索候选.

    Args:
        formula: 分子式 (如 'C8H8O2')
        max_results: 最大返回数
        min_heavy_atoms: 最小重原子数
        max_heavy_atoms: 最大重原子数
        quiet: 静默模式

    Returns:
        CoconutSearchResult
    """
    result = CoconutSearchResult(query_formula=formula)

    if not quiet:
        print(f"  [COCONUT] 按分子式检索: {formula} ...", file=sys.stderr)

    # COCONUT 高级搜索端点 (分子式)
    data = _coconut_api_get('/api/advanced_search', {
        'formula': formula,
        'size': min(max_results, 200),  # API 上限
    })

    if data is None:
        result.error = 'COCONUT API 无响应'
        return result

    if '_error' in data:
        result.error = data['_error']
        return result

    # 解析结果 (COCONUT 返回 hits + total 结构)
    hits = data.get('hits', data.get('results', data.get('data', [])))
    if isinstance(data, list):
        hits = data

    if not isinstance(hits, list):
        result.notes.append(f'COCONUT 返回格式异常: {type(hits).__name__}')
        return result

    result.total_found = len(hits)

    # 预筛: 重原子数范围
    for item in hits:
        if len(result.candidates) >= max_results:
            break
        try:
            # COCONUT 字段名标准化
            smiles = item.get('smiles') or item.get('SMILES') or item.get('canonical_smiles')
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            heavy = mol.GetNumHeavyAtoms()
            if heavy < min_heavy_atoms or heavy > max_heavy_atoms:
                continue
            cid = (item.get('coconut_id') or item.get('id')
                   or item.get('_id') or item.get('uri') or '')
            result.candidates.append(CoconutCandidate(
                coconut_id=str(cid),
                smiles=smiles,
                name=item.get('name') or item.get('trivial_name'),
                molecular_formula=item.get('molecular_formula') or item.get('formula'),
                molecular_weight=float(item.get('molecular_weight') or item.get('mw') or 0.0),
                heavy_atom_count=heavy,
            ))
        except Exception:
            continue

    if not quiet:
        print(f"  [COCONUT] 找到 {result.total_found} 条, "
              f"预筛后保留 {len(result.candidates)} 条", file=sys.stderr)

    return result


def search_coconut_by_substructure(
    smiles: str,
    max_results: int = 50,
    min_heavy_atoms: int = 3,
    max_heavy_atoms: int = 50,
    quiet: bool = False,
) -> CoconutSearchResult:
    """按子结构 (SMILES) 从 COCONUT 检索候选.

    Args:
        smiles: 子结构 SMILES (如 'c1ccccc1' 苯环)
        max_results: 最大返回数
        quiet: 静默模式

    Returns:
        CoconutSearchResult
    """
    result = CoconutSearchResult(query_smiles=smiles)

    if not quiet:
        print(f"  [COCONUT] 按子结构检索: {smiles} ...", file=sys.stderr)

    data = _coconut_api_get('/api/search/substructure', {
        'smiles': smiles,
        'size': min(max_results, 200),
    })

    if data is None or '_error' in data:
        result.error = (data or {}).get('_error', 'COCONUT API 无响应')
        return result

    hits = data.get('hits', data.get('results', data.get('data', [])))
    if isinstance(data, list):
        hits = data

    if not isinstance(hits, list):
        result.notes.append('COCONUT 返回格式异常')
        return result

    result.total_found = len(hits)

    for item in hits:
        if len(result.candidates) >= max_results:
            break
        try:
            smi = item.get('smiles') or item.get('SMILES') or item.get('canonical_smiles')
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            heavy = mol.GetNumHeavyAtoms()
            if heavy < min_heavy_atoms or heavy > max_heavy_atoms:
                continue
            cid = (item.get('coconut_id') or item.get('id')
                   or item.get('_id') or item.get('uri') or '')
            result.candidates.append(CoconutCandidate(
                coconut_id=str(cid),
                smiles=smi,
                name=item.get('name'),
                molecular_formula=item.get('molecular_formula'),
                molecular_weight=float(item.get('molecular_weight') or 0.0),
                heavy_atom_count=heavy,
            ))
        except Exception:
            continue

    if not quiet:
        print(f"  [COCONUT] 子结构检索: 找到 {result.total_found} 条, "
              f"预筛后保留 {len(result.candidates)} 条", file=sys.stderr)

    return result


def retrieve_coconut_candidates(
    formula: Optional[str] = None,
    c13_shifts: Optional[List[float]] = None,
    max_results: int = 50,
    quiet: bool = False,
) -> List[CoconutCandidate]:
    """综合 COCONUT 检索入口.

    优先按分子式检索; 若无分子式或结果不足, 按 13C 推断的子结构检索.

    Args:
        formula: 分子式 (可选)
        c13_shifts: 实测 13C 化学位移 (用于子结构推断)
        max_results: 最大返回数
        quiet: 静默模式

    Returns:
        CoconutCandidate 列表
    """
    all_candidates: List[CoconutCandidate] = []
    seen_smiles: set = set()

    # 1. 按分子式检索
    if formula:
        result = search_coconut_by_formula(
            formula, max_results=max_results, quiet=quiet,
        )
        if result.error and not quiet:
            print(f"  [COCONUT] WARN: {result.error}", file=sys.stderr)
        for c in result.candidates:
            if c.smiles not in seen_smiles:
                seen_smiles.add(c.smiles)
                all_candidates.append(c)

    # 2. 按子结构检索 (补充)
    if c13_shifts and len(all_candidates) < max_results // 2:
        # 复用 PubChem 的子结构推断逻辑
        from .retrieve_pubchem import infer_substructures_from_c13
        smarts_list = infer_substructures_from_c13(c13_shifts)
        for smarts in smarts_list[:2]:
            if len(all_candidates) >= max_results:
                break
            result = search_coconut_by_substructure(
                smarts, max_results=max_results // 4, quiet=quiet,
            )
            if result.error and not quiet:
                print(f"  [COCONUT] WARN: {result.error}", file=sys.stderr)
            for c in result.candidates:
                if c.smiles not in seen_smiles:
                    seen_smiles.add(c.smiles)
                    all_candidates.append(c)

    return all_candidates[:max_results]