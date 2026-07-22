"""lucy-ng 完整 HOSE 数据库的查询接口.

封装 `external/lucy-ng/data/reference/lucy-ng-derep.db` (~2.8 GB, 928K 化合物, 7.9M HOSE 统计).
本模块实现:

  1. 用 `hosegen` (Ratsemaat/HOSE_code_generator) 生成与 lucy-ng 数据库兼容的 HOSE 编码
  2. 直接查 SQLite `hose_stats` 表, O(1) 命中, 支持 radius 1~6 回退
  3. 输出与本项目 `nmr_predict.py` 兼容的 `PredictionResult` / `PredictedShift`

依赖:
  pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps
"""
from __future__ import annotations

import math
import os
import sqlite3
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Mol


# ---------- 可选依赖 ----------

try:
    from hosegen import HoseGenerator
    HOSEGEN_AVAILABLE = True
except ImportError:
    HOSEGEN_AVAILABLE = False
    HoseGenerator = None  # type: ignore[misc, assignment]


# ---------- 数据结构 ----------

@dataclass
class LNHOSEStats:
    """lucy-ng `hose_stats` 表的一行."""

    hose_code: str
    radius: int
    mean: float
    std: float
    count: int


@dataclass
class LNPredictedShift:
    """单个碳原子的预测结果 (lucy-ng 风格)."""

    atom_index: int
    predicted_ppm: float
    confidence: float
    hose_code: str
    radius_used: int
    match_count: int
    std_dev: float


@dataclass
class LNPredictionResult:
    """整个分子的 13C 预测结果 (lucy-ng 风格)."""

    smiles: str
    predictions: List[LNPredictedShift]
    carbon_count: int = 0
    success_count: int = 0

    @property
    def success_rate(self) -> float:
        if self.carbon_count == 0:
            return 0.0
        return self.success_count / self.carbon_count

    def get_shifts_sorted(self) -> List[float]:
        """返回预测位移列表 (降序)."""
        return sorted([p.predicted_ppm for p in self.predictions], reverse=True)


# ---------- SQLite 查询器 ----------

class LucyNGHOSEProvider:
    """lucy-ng-derep.db 的 HOSE 统计查询接口.

    用法::

        provider = LucyNGHOSEProvider(r'E:/path/lucy-ng-derep.db')
        result = provider.predict_smiles('CC(=O)Oc1ccccc1C(=O)O')
        for ps in result.predictions:
            print(ps.atom_index, ps.predicted_ppm, ps.radius_used, ps.match_count)

    内部用 `hosegen.HoseGenerator` 生成 HOSE 编码 (格式与 lucy-ng 库一致),
    然后按 radius 从大到小查 `hose_stats` 表.
    """

    DEFAULT_DB_PATH = r'E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db'

    def __init__(
        self,
        db_path: Optional[str] = None,
        max_radius: int = 6,
        min_radius: int = 1,
    ):
        if not HOSEGEN_AVAILABLE:
            raise ImportError(
                "hosegen 包未安装. 请运行: "
                "pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps"
            )
        if db_path is None:
            db_path = self.DEFAULT_DB_PATH
        if not os.path.isfile(db_path):
            raise FileNotFoundError(
                f"lucy-ng 数据库不存在: {db_path}\n"
                f"请用 util/download_lucy_hose_db.py 下载并解压"
            )
        self._db_path = db_path
        self._conn: Optional[sqlite3.Connection] = None
        self._hose_gen = HoseGenerator()
        self._max_radius = max_radius
        self._min_radius = min_radius
        self._open_db()

    def _open_db(self):
        """打开 SQLite 连接 (read-only 模式)."""
        # URI mode: file:...?mode=ro
        uri = f"file:{self._db_path}?mode=ro"
        try:
            self._conn = sqlite3.connect(uri, uri=True)
        except sqlite3.OperationalError:
            # 退回普通模式 (部分老版 sqlite3 不支持 URI)
            self._conn = sqlite3.connect(self._db_path)
        self._conn.row_factory = sqlite3.Row
        # 预编译查询 (sqlite3 >= 3.21; 低版本退化为 execute)
        try:
            self._stmt = self._conn.prepare(
                "SELECT mean, std, count FROM hose_stats WHERE hose_code=? AND radius=?"
            )
            self._use_prepared = True
        except AttributeError:
            self._stmt = None
            self._use_prepared = False

    def close(self):
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __repr__(self) -> str:
        try:
            cnt = self._conn.execute("SELECT COUNT(*) FROM hose_stats").fetchone()[0]
            return f"LucyNGHOSEProvider(stats={cnt:,}, db={self._db_path})"
        except Exception:
            return f"LucyNGHOSEProvider(db={self._db_path})"

    # --- 核心查询 ---

    def lookup_stats(self, hose_code: str, radius: int) -> Optional[LNHOSEStats]:
        """查 hose_stats 表. 返回 None 表示该 (code, radius) 无命中."""
        if not hose_code:
            return None
        if self._use_prepared:
            row = self._stmt.execute((hose_code, radius)).fetchone()
        else:
            row = self._conn.execute(
                "SELECT mean, std, count FROM hose_stats WHERE hose_code=? AND radius=?",
                (hose_code, radius),
            ).fetchone()
        if row is None:
            return None
        return LNHOSEStats(
            hose_code=hose_code,
            radius=radius,
            mean=row['mean'],
            std=row['std'],
            count=row['count'],
        )

    def predict_smiles(self, smiles: str) -> LNPredictionResult:
        """预测 SMILES 分子所有碳的 13C 化学位移."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return LNPredictionResult(smiles=smiles, predictions=[])
        # lucy-ng 内部 RemoveHs 后生成
        mol = Chem.RemoveHs(mol)
        predictions = self.predict_mol(mol)
        carbon_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
        return LNPredictionResult(
            smiles=smiles,
            predictions=predictions,
            carbon_count=carbon_count,
            success_count=len(predictions),
        )

    def predict_mol(self, mol: Mol) -> List[LNPredictedShift]:
        """预测 RDKit Mol 对象所有碳的 13C 化学位移."""
        predictions: List[LNPredictedShift] = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'C':
                continue
            atom_idx = atom.GetIdx()
            ps = self._predict_one(mol, atom_idx)
            if ps is not None:
                predictions.append(ps)
        return predictions

    def _predict_one(self, mol: Mol, atom_idx: int) -> Optional[LNPredictedShift]:
        """预测单个碳: radius 从大到小回退查询."""
        for radius in range(self._max_radius, self._min_radius - 1, -1):
            try:
                hose_code = self._hose_gen.get_Hose_codes(
                    mol, atom_idx, max_radius=radius,
                )
            except Exception:
                continue
            stats = self.lookup_stats(hose_code, radius)
            if stats is None:
                continue
            confidence = self._calc_confidence(radius, stats.count, stats.std)
            return LNPredictedShift(
                atom_index=atom_idx,
                predicted_ppm=stats.mean,
                confidence=confidence,
                hose_code=hose_code,
                radius_used=radius,
                match_count=stats.count,
                std_dev=stats.std,
            )
        return None

    @staticmethod
    def _calc_confidence(radius: int, match_count: int, std_dev: float) -> float:
        """计算置信度 (0-1, 与 lucy-ng C13Predictor._calculate_confidence 一致)."""
        radius_score = radius / 6.0
        match_score = min(1.0, 0.5 + 0.5 * math.log10(max(1, match_count)) / 2)
        std_score = max(0.5, 1.0 - std_dev / 20)
        return round(0.5 * radius_score + 0.3 * match_score + 0.2 * std_score, 3)


# ---------- 验证分数 (兼容 nmr_predict.verify_candidate 接口) ----------

def verify_candidate_with_lucy_db(
    provider: LucyNGHOSEProvider,
    smiles: str,
    observed_c13: List[float],
    c_tol: float = 3.0,
) -> 'LucyNGVerificationScore':
    """用 lucy-ng 数据库预测 SMILES 的 13C, 与实测对比.

    返回值结构与 `nmr_predict.VerificationScore` 类似, 但加注来源标识.

    Args:
        provider: LucyNGHOSEProvider 实例
        smiles: 候选 SMILES
        observed_c13: 实测 13C 位移列表 (ppm)
        c_tol: 匹配容差 (ppm), 默认 3.0
    """
    from .nmr_predict import VerificationScore

    result = provider.predict_smiles(smiles)
    if not result.predictions:
        return VerificationScore(
            smiles=smiles,
            hose_rms=999.0,
            hose_coverage=0.0,
            hose_match_rate=0.0,
            verification_score=0.0,
            predicted_shifts=[],
            matched_pairs=[],
            notes=['lucy-ng 预测: 无碳原子或全部查询失败'],
        )

    predicted = result.get_shifts_sorted()
    obs_sorted = sorted(observed_c13, reverse=True)

    # 一对一贪心匹配
    pairs: List[Tuple[float, float, float]] = []  # (pred, obs, diff)
    used_obs = set()
    matched_count = 0
    for p in predicted:
        best_i = -1
        best_diff = c_tol + 1
        for i, o in enumerate(obs_sorted):
            if i in used_obs:
                continue
            d = abs(p - o)
            if d < best_diff:
                best_diff = d
                best_i = i
        if best_i >= 0 and best_diff <= c_tol:
            used_obs.add(best_i)
            pairs.append((p, obs_sorted[best_i], best_diff))
            matched_count += 1

    # RMS
    if pairs:
        rms = math.sqrt(sum(d * d for _, _, d in pairs) / len(pairs))
    else:
        rms = 999.0

    # 覆盖率 (成功预测碳数 / 总碳数)
    coverage = result.success_rate
    # 匹配率 (匹配对数 / min(预测数, 实测数))
    match_rate = matched_count / max(1, min(len(predicted), len(obs_sorted)))

    # 综合验证分 (0-1)
    if matched_count > 0:
        rms_score = max(0.0, 1.0 - rms / 10.0)  # RMS=0→1, RMS>=10→0
        coverage_score = coverage
        match_score = match_rate
        verification_score = round(0.4 * rms_score + 0.3 * coverage_score + 0.3 * match_score, 3)
    else:
        verification_score = 0.0

    return VerificationScore(
        smiles=smiles,
        hose_rms=round(rms, 3),
        hose_coverage=round(coverage, 3),
        hose_match_rate=round(match_rate, 3),
        verification_score=verification_score,
        predicted_shifts=predicted,
        matched_pairs=pairs,
        notes=[
            f"lucy-ng 数据库预测: {result.success_count}/{result.carbon_count} 碳, "
            f"匹配 {matched_count}/{min(len(predicted), len(obs_sorted))} 峰, "
            f"RMS={rms:.2f} ppm, 综合分={verification_score}",
        ],
    )


# ---------- 便捷入口 ----------

_global_provider: Optional[LucyNGHOSEProvider] = None


def get_lucy_ng_provider(
    db_path: Optional[str] = None,
    quiet: bool = False,
) -> Optional[LucyNGHOSEProvider]:
    """获取全局 LucyNGHOSEProvider (懒加载, 单例)."""
    global _global_provider
    if _global_provider is not None:
        return _global_provider
    if not HOSEGEN_AVAILABLE:
        if not quiet:
            print(
                '  [lucy-ng] WARN: hosegen 未安装, 跳过 lucy-ng 数据库集成\n'
                '         安装: pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps',
                file=sys.stderr,
            )
        return None
    if db_path is None:
        db_path = LucyNGHOSEProvider.DEFAULT_DB_PATH
    if not os.path.isfile(db_path):
        if not quiet:
            print(
                f'  [lucy-ng] WARN: 数据库不存在: {db_path}\n'
                f'         下载: python util/download_lucy_hose_db.py --strategy colab',
                file=sys.stderr,
            )
        return None
    try:
        _global_provider = LucyNGHOSEProvider(db_path)
        if not quiet:
            print(f'  [lucy-ng] 已连接数据库: {_global_provider}', file=sys.stderr)
        return _global_provider
    except Exception as e:
        if not quiet:
            print(f'  [lucy-ng] WARN: 连接失败: {e}', file=sys.stderr)
        return None