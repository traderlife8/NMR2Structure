"""lucy-ng 数据库候选检索器 (Part 3 Phase 7).

与 `nmrshiftdb2 SDF` 检索的区别:
  - 数据规模: 928K 化合物 + 24M 化学位移 (vs nmrshiftdb2 58K + 159MB)
  - 检索粒度: SQLite 直接 SQL 查询 (vs SDF 流式文本扫描)
  - 谱库来源: 文献汇总 (Coconut 等) + 计算 13C 谱 (vs 实测实验谱)

策略:
  1. 硬筛: 碳数范围 (carbon_count ± 2) + 分子式约束 (formula, 可选)
  2. 软筛: 13C 化学位移贪心匹配 + 1H 化学位移贪心匹配 + RMS 计算
  3. 输出: 与 retrieve.py RetrievalHit 接口兼容的候选列表

依赖:
  - lucy-ng-derep.db (E:\\NMR2Structure\\data\\lucy-ng-derep_db\\lucy-ng-derep.db)
  - sqlite3 (标准库)
"""
from __future__ import annotations

import math
import os
import sqlite3
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


@dataclass
class LucyNGHit:
    """lucy-ng 检索的单条命中.

    与 retrieve.RetrievalHit 字段兼容, 但简化版(无 DEPT/HSQC 证据).
    """

    compound_id: int
    name: str
    smiles: str
    canonical_smiles: str = ''
    inchi: str = ''
    inchikey: str = ''
    formula: str = ''
    carbon_count: int = 0
    lib_13c_shifts: List[float] = field(default_factory=list)
    lib_1h_shifts: List[float] = field(default_factory=list)
    # 匹配打分
    c_match_count: int = 0
    c_total_obs: int = 0
    c_total_lib: int = 0
    c_rms: float = 0.0
    c_coverage: float = 0.0
    h_match_count: int = 0
    h_total_obs: int = 0
    h_total_lib: int = 0
    h_rms: float = 0.0
    h_coverage: float = 0.0
    # 综合分
    score: float = 0.0


class LucyNGRetriever:
    """lucy-ng 数据库检索器."""

    DEFAULT_DB_PATH = r"E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db"

    def __init__(self, db_path: Optional[str] = None):
        if db_path is None:
            db_path = self.DEFAULT_DB_PATH
        if not os.path.isfile(db_path):
            raise FileNotFoundError(
                f"lucy-ng 数据库不存在: {db_path}\n"
                f"请用 util/download_lucy_hose_db.py 下载并解压"
            )
        self._db_path = db_path
        # read-only URI 模式
        try:
            uri = f"file:{db_path}?mode=ro"
            self._conn = sqlite3.connect(uri, uri=True)
        except sqlite3.OperationalError:
            self._conn = sqlite3.connect(db_path)
        self._conn.row_factory = sqlite3.Row

    def close(self):
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False

    def count_candidates(self, carbon_count: int, tol: int = 2) -> int:
        """统计碳数 ±tol 范围内的候选数."""
        cur = self._conn.execute(
            "SELECT COUNT(*) FROM compounds WHERE carbon_count BETWEEN ? AND ?",
            (carbon_count - tol, carbon_count + tol),
        )
        return cur.fetchone()[0]

    def retrieve_by_carbon_count(
        self,
        observed_c13: List[float],
        observed_h1: List[float],
        carbon_count: Optional[int] = None,
        carbon_tol: int = 2,
        formula: Optional[str] = None,
        c_tol: float = 2.0,
        h_tol: float = 0.3,
        min_c_matches: int = 10,
        max_results: Optional[int] = None,
        progress_every: int = 2000,
        quiet: bool = True,
    ) -> List[LucyNGHit]:
        """按碳数范围 + 谱图匹配检索.

        Args:
            observed_c13: 实测 13C 化学位移
            observed_h1: 实测 1H 化学位移
            carbon_count: 目标碳数 (None = 不指定)
            carbon_tol: 碳数容差 (±)
            formula: 精确分子式约束 (None = 不指定)
            c_tol: 13C 匹配容差 (ppm)
            h_tol: 1H 匹配容差 (ppm)
            min_c_matches: 最少 13C 匹配数 (过滤无意义候选)
            max_results: 限制返回数 (None = 不限制)
            progress_every: 进度打印间隔
            quiet: 静默模式
        """
        if carbon_count is not None:
            sql = "SELECT id, name, smiles, inchi, inchi_key, formula, carbon_count " \
                  "FROM compounds WHERE carbon_count BETWEEN ? AND ?"
            params: tuple = (carbon_count - carbon_tol, carbon_count + carbon_tol)
            if formula:
                sql += " AND formula = ?"
                params = params + (formula,)
        elif formula:
            sql = "SELECT id, name, smiles, inchi, inchi_key, formula, carbon_count " \
                  "FROM compounds WHERE formula = ?"
            params = (formula,)
        else:
            sql = "SELECT id, name, smiles, inchi, inchi_key, formula, carbon_count FROM compounds"
            params = ()

        cur = self._conn.execute(sql, params)
        candidates = cur.fetchall()

        if not quiet:
            print(f"  lucy-ng: 按碳数 {carbon_count}±{carbon_tol if carbon_count else '无'} "
                  f"+ 分子式 {formula or '无'} → {len(candidates)} 个候选", flush=True)

        obs_c = sorted(observed_c13, reverse=True)
        obs_h = sorted(observed_h1, reverse=True)
        hits: List[LucyNGHit] = []

        for i, row in enumerate(candidates):
            if progress_every and not quiet and i % progress_every == 0:
                print(f"  lucy-ng: 扫描 {i}/{len(candidates)}...", flush=True)
            hit = self._match_one(
                row=row,
                obs_c=obs_c,
                obs_h=obs_h,
                c_tol=c_tol,
                h_tol=h_tol,
            )
            if hit and hit.c_match_count >= min_c_matches:
                hits.append(hit)

        # 综合打分: 13C 匹配数 + 覆盖率 + (-RMS)
        # 设计: 主分主要看 13C (lucy-ng 主谱库是 13C), 1H 缺失不扣分
        for hit in hits:
            h_pen = 0.0
            if hit.h_total_lib > 0:
                # 只有 1H 库存在时才用 1H RMS 扣分 (避免 lucy-ng 候选被 999 拖毙)
                h_pen = 5.0 * hit.h_rms
            hit.score = (
                hit.c_match_count * 4.0
                + hit.c_coverage * 30.0
                - hit.c_rms * 3.0
                + hit.h_match_count * 1.0
                - h_pen
            )
        hits.sort(key=lambda h: h.score, reverse=True)

        if max_results is not None:
            hits = hits[:max_results]

        if not quiet:
            print(f"  lucy-ng: 命中 {len(hits)} 条 "
                  f"(min_c_matches={min_c_matches})", flush=True)

        return hits

    def _match_one(
        self,
        row: sqlite3.Row,
        obs_c: List[float],
        obs_h: List[float],
        c_tol: float,
        h_tol: float,
    ) -> Optional[LucyNGHit]:
        """对单个候选做 13C + 1H 贪心匹配."""
        cid = row['id']
        # 拉 13C shifts
        cur = self._conn.execute(
            "SELECT shift_ppm FROM shifts WHERE compound_id=? ORDER BY shift_ppm DESC",
            (cid,),
        )
        lib_c = [s['shift_ppm'] for s in cur.fetchall()]
        # 1H shifts 在 lucy-ng 中不一定有 (主要存 13C); 留空
        lib_h: List[float] = []

        if not lib_c:
            return None

        # 13C 一对一贪心匹配
        used_c = set()
        c_pairs = []
        for o in obs_c:
            best_d = c_tol + 1
            best_j = -1
            for j, c in enumerate(lib_c):
                if j in used_c:
                    continue
                d = abs(c - o)
                if d < best_d:
                    best_d = d
                    best_j = j
            if best_d <= c_tol and best_j >= 0:
                used_c.add(best_j)
                c_pairs.append(best_d)

        n_c_matched = len(c_pairs)
        c_rms = math.sqrt(sum(d * d for d in c_pairs) / n_c_matched) if c_pairs else 999.0
        c_coverage = n_c_matched / max(1, len(lib_c))

        # 1H 匹配 (若无 lib_h, 跳过)
        h_pairs: List[float] = []
        if lib_h and obs_h:
            used_h = set()
            for o in obs_h:
                best_d = h_tol + 1
                best_j = -1
                for j, h in enumerate(lib_h):
                    if j in used_h:
                        continue
                    d = abs(h - o)
                    if d < best_d:
                        best_d = d
                        best_j = j
                if best_d <= h_tol and best_j >= 0:
                    used_h.add(best_j)
                    h_pairs.append(best_d)
        h_rms = math.sqrt(sum(d * d for d in h_pairs) / len(h_pairs)) if h_pairs else 999.0
        h_coverage = len(h_pairs) / max(1, len(lib_h)) if lib_h else 0.0

        # 规范化 SMILES
        smi = row['smiles'] or ''
        canon = ''
        try:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                canon = Chem.MolToSmiles(m)
        except Exception:
            pass

        return LucyNGHit(
            compound_id=cid,
            name=row['name'] or '',
            smiles=smi,
            canonical_smiles=canon,
            inchi=row['inchi'] or '',
            inchikey=row['inchi_key'] or '',
            formula=row['formula'] or '',
            carbon_count=row['carbon_count'] or 0,
            lib_13c_shifts=lib_c,
            lib_1h_shifts=lib_h,
            c_match_count=n_c_matched,
            c_total_obs=len(obs_c),
            c_total_lib=len(lib_c),
            c_rms=c_rms,
            c_coverage=c_coverage,
            h_match_count=len(h_pairs),
            h_total_obs=len(obs_h),
            h_total_lib=len(lib_h),
            h_rms=h_rms,
            h_coverage=h_coverage,
        )


# ---- 便捷接口 ----

_global_retriever: Optional[LucyNGRetriever] = None


def get_lucy_ng_retriever(
    db_path: Optional[str] = None,
    quiet: bool = False,
) -> Optional[LucyNGRetriever]:
    """获取全局 LucyNGRetriever (懒加载)."""
    global _global_retriever
    if _global_retriever is not None:
        return _global_retriever
    if db_path is None:
        db_path = LucyNGRetriever.DEFAULT_DB_PATH
    if not os.path.isfile(db_path):
        if not quiet:
            print(f"  [lucy-ng] WARN: 数据库不存在: {db_path}", flush=True)
        return None
    try:
        _global_retriever = LucyNGRetriever(db_path)
        if not quiet:
            print(f"  [lucy-ng] 已连接: {_global_retriever._db_path}", flush=True)
        return _global_retriever
    except Exception as e:
        if not quiet:
            print(f"  [lucy-ng] WARN: 连接失败: {e}", flush=True)
        return None