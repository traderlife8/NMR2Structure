"""PubChem 候选检索模块.

当 nmrshiftdb2 无命中或命中不足时, 从 PubChem 按分子式/子结构检索候选.
使用已安装的 pubchempy 库, 处理网络失败、超时、空结果等异常.

设计要点:
  - 按分子式检索 (get_compounds(formula, 'formula'))
  - 按 13C 特征位移推断的子结构做子结构检索 (可选)
  - 对候选用 HOSE 预测验证 (依赖 nmr_predict 模块)
  - 限制: 单次最多返回 100 条, 按分子量/重原子数预筛
  - 网络失败时静默降级, 不阻断主流程
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


@dataclass
class PubChemCandidate:
    """PubChem 检索到的候选."""

    cid: int
    smiles: str
    iupac_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: float = 0.0
    heavy_atom_count: int = 0
    # 验证分数 (由 nmr_predict 填充)
    verification_score: float = 0.0
    source: str = 'pubchem'


@dataclass
class PubChemSearchResult:
    """PubChem 检索结果."""

    candidates: List[PubChemCandidate] = field(default_factory=list)
    query_formula: Optional[str] = None
    query_substructure: Optional[str] = None
    total_found: int = 0
    error: Optional[str] = None
    notes: List[str] = field(default_factory=list)


def search_pubchem_by_formula(
    formula: str,
    max_results: int = 100,
    min_heavy_atoms: int = 3,
    max_heavy_atoms: int = 30,
    quiet: bool = False,
) -> PubChemSearchResult:
    """按分子式从 PubChem 检索候选.

    Args:
        formula: 分子式 (如 'C8H8O2')
        max_results: 最大返回数
        min_heavy_atoms: 最小重原子数
        max_heavy_atoms: 最大重原子数
        quiet: 静默模式

    Returns:
        PubChemSearchResult
    """
    result = PubChemSearchResult(query_formula=formula)

    try:
        import pubchempy as pcp
    except ImportError:
        result.error = 'pubchempy 未安装'
        result.notes.append('pip install pubchempy 后可启用 PubChem 检索')
        return result

    if not quiet:
        print(f"  [PubChem] 按分子式检索: {formula} ...", file=sys.stderr)

    try:
        compounds = pcp.get_compounds(formula, 'formula')
    except pcp.ServerBusyError:
        result.error = 'PubChem 服务器繁忙, 请稍后重试'
        return result
    except pcp.TimeoutError:
        result.error = 'PubChem 请求超时'
        return result
    except pcp.NotFoundError:
        result.notes.append(f'PubChem 中未找到分子式 {formula} 的化合物')
        return result
    except Exception as e:
        result.error = f'PubChem 检索失败: {e}'
        return result

    if not compounds:
        result.notes.append(f'PubChem 中未找到分子式 {formula} 的化合物')
        return result

    result.total_found = len(compounds)

    # 预筛: 重原子数范围
    for c in compounds:
        if len(result.candidates) >= max_results:
            break
        try:
            mol = Chem.MolFromSmiles(c.smiles)
            if mol is None:
                continue
            heavy = mol.GetNumHeavyAtoms()
            if heavy < min_heavy_atoms or heavy > max_heavy_atoms:
                continue
            result.candidates.append(PubChemCandidate(
                cid=c.cid,
                smiles=c.smiles,
                iupac_name=c.iupac_name,
                molecular_formula=c.molecular_formula,
                molecular_weight=c.molecular_weight or 0.0,
                heavy_atom_count=heavy,
            ))
        except Exception:
            continue

    if not quiet:
        print(f"  [PubChem] 找到 {result.total_found} 条, "
              f"预筛后保留 {len(result.candidates)} 条", file=sys.stderr)

    return result


def search_pubchem_by_substructure(
    smarts: str,
    max_results: int = 50,
    min_heavy_atoms: int = 3,
    max_heavy_atoms: int = 30,
    quiet: bool = False,
) -> PubChemSearchResult:
    """按子结构 (SMARTS) 从 PubChem 检索候选.

    Args:
        smarts: SMARTS 查询 (如 'c1ccccc1' 苯环)
        max_results: 最大返回数
        min_heavy_atoms: 最小重原子数
        max_heavy_atoms: 最大重原子数
        quiet: 静默模式

    Returns:
        PubChemSearchResult
    """
    result = PubChemSearchResult(query_substructure=smarts)

    try:
        import pubchempy as pcp
    except ImportError:
        result.error = 'pubchempy 未安装'
        return result

    if not quiet:
        print(f"  [PubChem] 按子结构检索: {smarts} ...", file=sys.stderr)

    try:
        compounds = pcp.get_compounds(smarts, 'smiles', searchtype='substructure')
    except pcp.ServerBusyError:
        result.error = 'PubChem 服务器繁忙'
        return result
    except pcp.TimeoutError:
        result.error = 'PubChem 请求超时'
        return result
    except Exception as e:
        result.error = f'PubChem 子结构检索失败: {e}'
        return result

    if not compounds:
        result.notes.append(f'PubChem 中未找到含子结构 {smarts} 的化合物')
        return result

    result.total_found = len(compounds)

    for c in compounds:
        if len(result.candidates) >= max_results:
            break
        try:
            mol = Chem.MolFromSmiles(c.smiles)
            if mol is None:
                continue
            heavy = mol.GetNumHeavyAtoms()
            if heavy < min_heavy_atoms or heavy > max_heavy_atoms:
                continue
            result.candidates.append(PubChemCandidate(
                cid=c.cid,
                smiles=c.smiles,
                iupac_name=c.iupac_name,
                molecular_formula=c.molecular_formula,
                molecular_weight=c.molecular_weight or 0.0,
                heavy_atom_count=heavy,
            ))
        except Exception:
            continue

    if not quiet:
        print(f"  [PubChem] 子结构检索: 找到 {result.total_found} 条, "
              f"预筛后保留 {len(result.candidates)} 条", file=sys.stderr)

    return result


def infer_substructures_from_c13(c13_shifts: List[float]) -> List[str]:
    """从 13C 化学位移推断可能存在的子结构 (SMARTS).

    基于经验化学位移范围:
      - 160-220 ppm: 羰基 (C=O)
      - 110-160 ppm: 芳香/烯烃碳
      - 50-90 ppm: 连杂原子脂肪碳 (C-O, C-N)
      - 0-50 ppm: 饱和脂肪碳

    Returns:
        SMARTS 列表
    """
    smarts_list: List[str] = []

    has_carbonyl = any(160 <= s <= 220 for s in c13_shifts)
    has_aromatic = any(110 <= s <= 160 for s in c13_shifts)
    has_hetero_aliphatic = any(50 <= s <= 90 for s in c13_shifts)

    if has_carbonyl:
        smarts_list.append('[CX3]=[OX1]')  # 羰基
    if has_aromatic:
        smarts_list.append('c1ccccc1')  # 苯环 (最常见的芳香体系)
    if has_hetero_aliphatic:
        smarts_list.append('[CX4]-[OX2]')  # C-O 单键

    return smarts_list


def retrieve_pubchem_candidates(
    formula: Optional[str] = None,
    c13_shifts: Optional[List[float]] = None,
    max_results: int = 100,
    quiet: bool = False,
) -> List[PubChemCandidate]:
    """综合 PubChem 检索入口.

    优先按分子式检索; 若无分子式, 按 13C 推断的子结构检索.

    Args:
        formula: 分子式 (可选)
        c13_shifts: 实测 13C 化学位移 (用于子结构推断)
        max_results: 最大返回数
        quiet: 静默模式

    Returns:
        PubChemCandidate 列表
    """
    all_candidates: List[PubChemCandidate] = []
    seen_smiles: set = set()

    # 1. 按分子式检索
    if formula:
        result = search_pubchem_by_formula(
            formula, max_results=max_results, quiet=quiet,
        )
        if result.error and not quiet:
            print(f"  [PubChem] WARN: {result.error}", file=sys.stderr)
        for c in result.candidates:
            if c.smiles not in seen_smiles:
                seen_smiles.add(c.smiles)
                all_candidates.append(c)

    # 2. 按子结构检索 (补充)
    if c13_shifts and len(all_candidates) < max_results // 2:
        smarts_list = infer_substructures_from_c13(c13_shifts)
        for smarts in smarts_list[:2]:  # 最多 2 个子结构查询
            if len(all_candidates) >= max_results:
                break
            result = search_pubchem_by_substructure(
                smarts, max_results=max_results // 4, quiet=quiet,
            )
            if result.error and not quiet:
                print(f"  [PubChem] WARN: {result.error}", file=sys.stderr)
            for c in result.candidates:
                if c.smiles not in seen_smiles:
                    seen_smiles.add(c.smiles)
                    all_candidates.append(c)

    return all_candidates[:max_results]
