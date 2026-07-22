"""13C NMR 化学位移预测与候选验证.

核心思路:
  1. 从 nmrshiftdb2 SDF 构建 HOSE 编码 → 13C 化学位移统计查找表
  2. 对候选结构生成 HOSE 编码, 查表预测每个碳的化学位移
  3. 将预测位移与实测对比, 计算 RMS 和覆盖率作为验证分数

不依赖外部 hosegen 包或 2.8 GB 数据库, 完全基于本地 nmrshiftdb2 SDF 自建.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

from .hose_generator import generate_hose_code, prepare_mol


# ---------- 数据结构 ----------

@dataclass
class HOSEStats:
    """单个 HOSE 编码对应的 13C 化学位移统计."""

    mean: float
    std: float
    count: int
    min_shift: float
    max_shift: float


@dataclass
class PredictedShift:
    """单个碳原子的预测结果."""

    atom_index: int
    predicted_ppm: float
    confidence: float  # 0-1
    radius_used: int
    match_count: int


@dataclass
class PredictionResult:
    """整个分子的 13C 预测结果."""

    smiles: str
    predictions: List[PredictedShift] = field(default_factory=list)
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


@dataclass
class VerificationScore:
    """候选结构验证分数."""

    smiles: str
    # HOSE 预测 vs 实测
    hose_rms: float = 999.0  # 预测-实测 RMS (ppm), 越小越好
    hose_coverage: float = 0.0  # 成功预测的碳占比
    hose_match_rate: float = 0.0  # 预测位移在容差内匹配实测的比例
    # 综合验证分 (0-1, 越高越可信)
    verification_score: float = 0.0
    # 详细
    predicted_shifts: List[float] = field(default_factory=list)
    matched_pairs: List[Tuple[float, float, float]] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


# ---------- HOSE 查找表 ----------

class HOSELookupTable:
    """HOSE 编码 → 13C 化学位移统计查找表.

    支持多级 radius 回退: 先查 radius=6, 无命中则降到 5, 4, ...
    """

    def __init__(self):
        # {hose_code: [shift1, shift2, ...]}
        self._raw: Dict[str, List[float]] = {}
        # 缓存统计
        self._stats_cache: Dict[str, HOSEStats] = {}
        self._built = False

    @property
    def size(self) -> int:
        return len(self._raw)

    def add_observation(self, hose_code: str, shift_ppm: float):
        """添加一条观测 (HOSE 编码, 化学位移)."""
        if not hose_code:
            return
        if hose_code not in self._raw:
            self._raw[hose_code] = []
        self._raw[hose_code].append(shift_ppm)
        # 清除缓存
        self._stats_cache.pop(hose_code, None)

    def build_from_sdf(
        self,
        sdf_path: str,
        max_records: int = 0,
        progress_every: int = 5000,
        quiet: bool = False,
    ):
        """从 nmrshiftdb2 SDF 构建查找表.

        对每条记录的每个碳原子:
          1. 解析分子
          2. 生成 HOSE 编码 (radius 1~6)
          3. 将对应 13C 化学位移加入查找表

        Args:
            sdf_path: nmrshiftdb2 SDF 文件路径
            max_records: 最大处理记录数 (0=全部)
            progress_every: 每 N 条打印进度
            quiet: 不打印进度
        """
        from .io_sdf import stream_nmrshiftdb2_sdf

        count = 0
        n_carbons = 0
        for entry in stream_nmrshiftdb2_sdf(sdf_path):
            count += 1
            if max_records > 0 and count > max_records:
                break
            if not quiet and progress_every > 0 and count % progress_every == 0:
                print(f"  [HOSE 建表] 已处理 {count} 条, {n_carbons} 个碳 ...",
                      file=sys.stderr)

            # 需要有 13C 数据
            if not entry.spectrum_13c:
                continue

            # 解析分子
            mol = entry.mol()
            if mol is None:
                continue
            mol = Chem.RemoveHs(mol)
            if mol is None:
                continue

            # 获取碳原子列表
            carbon_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == 'C']
            if not carbon_atoms:
                continue

            # 13C 位移列表
            c_shifts = [sig.shift for sig in entry.spectrum_13c]
            if not c_shifts:
                continue

            # 策略: 将 13C 位移按顺序分配给碳原子
            # nmrshiftdb2 中 spectrum_13c 的 atom_index 字段可能不可靠,
            # 用贪心匹配: 对每个碳生成 HOSE 编码, 将最近的位移分配给它
            # 简化: 如果碳数 == 位移数, 按索引对应; 否则跳过
            n_c_atoms = len(carbon_atoms)
            n_c_shifts = len(c_shifts)

            if n_c_atoms == n_c_shifts:
                # 直接对应
                for i, atom in enumerate(carbon_atoms):
                    shift = c_shifts[i]
                    for radius in range(1, 7):
                        code = generate_hose_code(mol, atom.GetIdx(), max_radius=radius)
                        if code:
                            self.add_observation(code, shift)
                    n_carbons += 1
            elif n_c_shifts > 0:
                # 碳数 != 位移数: 用贪心匹配
                # 为每个碳生成 radius=3 的编码, 按顺序分配位移
                # (保守策略: 只用 radius<=3 避免过度拟合)
                for i, atom in enumerate(carbon_atoms):
                    if i < n_c_shifts:
                        shift = c_shifts[i]
                        for radius in range(1, 4):
                            code = generate_hose_code(mol, atom.GetIdx(), max_radius=radius)
                            if code:
                                self.add_observation(code, shift)
                        n_carbons += 1

        self._built = True
        if not quiet:
            print(f"  [HOSE 建表] 完成: {count} 条记录, {n_carbons} 个碳, "
                  f"{self.size} 个唯一编码", file=sys.stderr)

    def lookup(self, hose_code: str) -> Optional[HOSEStats]:
        """查找 HOSE 编码对应的统计."""
        if not hose_code:
            return None
        if hose_code in self._stats_cache:
            return self._stats_cache[hose_code]
        shifts = self._raw.get(hose_code)
        if not shifts:
            return None
        n = len(shifts)
        mean = sum(shifts) / n
        if n > 1:
            var = sum((s - mean) ** 2 for s in shifts) / (n - 1)
            std = math.sqrt(var)
        else:
            std = 0.0
        stats = HOSEStats(
            mean=mean, std=std, count=n,
            min_shift=min(shifts), max_shift=max(shifts),
        )
        self._stats_cache[hose_code] = stats
        return stats

    def save(self, path: str):
        """保存查找表到文件 (简单文本格式: code\\tshift1,shift2,...)."""
        with open(path, 'w', encoding='utf-8') as f:
            for code, shifts in self._raw.items():
                shifts_str = ','.join(f'{s:.2f}' for s in shifts)
                f.write(f'{code}\t{shifts_str}\n')

    def load(self, path: str):
        """从文件加载查找表."""
        self._raw.clear()
        self._stats_cache.clear()
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or '\t' not in line:
                    continue
                code, shifts_str = line.split('\t', 1)
                shifts = [float(x) for x in shifts_str.split(',') if x]
                if shifts:
                    self._raw[code] = shifts
        self._built = True


# ---------- 预测器 ----------

class C13Predictor:
    """13C 化学位移预测器.

    使用 HOSE 查找表, 对候选结构的每个碳原子预测化学位移.
    支持多级 radius 回退 (6→5→4→...→1).
    """

    def __init__(self, lookup: HOSELookupTable, max_radius: int = 6, min_radius: int = 1):
        self._lookup = lookup
        self._max_radius = max_radius
        self._min_radius = min_radius

    def predict_smiles(self, smiles: str) -> PredictionResult:
        """预测 SMILES 对应分子的 13C 化学位移."""
        mol = prepare_mol(smiles)
        if mol is None:
            return PredictionResult(smiles=smiles, carbon_count=0, success_count=0)

        predictions: List[PredictedShift] = []
        carbon_count = 0

        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'C':
                continue
            carbon_count += 1
            pred = self._predict_atom(mol, atom.GetIdx())
            if pred is not None:
                predictions.append(pred)

        return PredictionResult(
            smiles=smiles,
            predictions=predictions,
            carbon_count=carbon_count,
            success_count=len(predictions),
        )

    def _predict_atom(self, mol, atom_idx: int) -> Optional[PredictedShift]:
        """预测单个碳原子的化学位移 (多级回退)."""
        for radius in range(self._max_radius, self._min_radius - 1, -1):
            code = generate_hose_code(mol, atom_idx, max_radius=radius)
            if not code:
                continue
            stats = self._lookup.lookup(code)
            if stats is not None and stats.count >= 1:
                confidence = self._calc_confidence(radius, stats)
                return PredictedShift(
                    atom_index=atom_idx,
                    predicted_ppm=round(stats.mean, 2),
                    confidence=confidence,
                    radius_used=radius,
                    match_count=stats.count,
                )
        return None

    def _calc_confidence(self, radius: int, stats: HOSEStats) -> float:
        """计算预测置信度 (0-1)."""
        # radius 越大越特异, 置信度越高
        radius_score = radius / self._max_radius
        # 匹配数越多越可靠
        match_score = min(1.0, 0.5 + 0.5 * math.log10(max(1, stats.count)) / 2)
        # 标准差越小越可靠
        std_score = max(0.5, 1.0 - stats.std / 20.0)
        return round(0.5 * radius_score + 0.3 * match_score + 0.2 * std_score, 3)


# ---------- 验证函数 ----------

def verify_candidate(
    predictor: C13Predictor,
    smiles: str,
    observed_c13: List[float],
    c_tol: float = 3.0,
) -> VerificationScore:
    """验证候选结构: 预测 13C 位移并与实测对比.

    Args:
        predictor: C13Predictor 实例
        smiles: 候选 SMILES
        observed_c13: 实测 13C 化学位移列表
        c_tol: 匹配容差 (ppm)

    Returns:
        VerificationScore 验证结果
    """
    result = predictor.predict_smiles(smiles)
    score = VerificationScore(smiles=smiles)

    if not result.predictions or not observed_c13:
        score.notes.append('预测或实测数据为空, 无法验证')
        return score

    predicted = result.get_shifts_sorted()
    observed = sorted(observed_c13, reverse=True)
    score.predicted_shifts = predicted
    score.hose_coverage = result.success_rate

    # 贪心匹配: 对每个实测峰找最近的预测峰
    used_pred = [False] * len(predicted)
    matched_pairs: List[Tuple[float, float, float]] = []  # (obs, pred, delta)
    deltas: List[float] = []

    for obs in observed:
        best_i = -1
        best_d = float('inf')
        for i, pred in enumerate(predicted):
            if used_pred[i]:
                continue
            d = abs(obs - pred)
            if d < best_d:
                best_d = d
                best_i = i
        if best_i >= 0 and best_d <= c_tol:
            used_pred[best_i] = True
            matched_pairs.append((obs, predicted[best_i], best_d))
            deltas.append(best_d)

    # 计算指标
    n_matched = len(matched_pairs)
    n_obs = len(observed)
    score.matched_pairs = matched_pairs
    score.hose_match_rate = n_matched / n_obs if n_obs > 0 else 0.0

    if deltas:
        rms = math.sqrt(sum(d ** 2 for d in deltas) / len(deltas))
        score.hose_rms = round(rms, 3)
    else:
        score.hose_rms = 999.0

    # 综合验证分 (0-1)
    # 权重: match_rate 0.5 + (1 - rms/10) 0.3 + coverage 0.2
    rms_score = max(0.0, 1.0 - score.hose_rms / 10.0)
    score.verification_score = round(
        0.5 * score.hose_match_rate + 0.3 * rms_score + 0.2 * score.hose_coverage,
        4,
    )

    score.notes.append(
        f'HOSE 验证: 匹配 {n_matched}/{n_obs} 峰, '
        f'RMS={score.hose_rms:.2f} ppm, '
        f'覆盖率={score.hose_coverage:.0%}, '
        f'综合分={score.verification_score:.3f}'
    )

    return score


# ---------- 便捷入口 ----------

_global_predictor: Optional[C13Predictor] = None
_global_table: Optional[HOSELookupTable] = None


def get_predictor(
    sdf_path: Optional[str] = None,
    table_path: Optional[str] = None,
    max_records: int = 0,
    quiet: bool = False,
) -> C13Predictor:
    """获取全局 C13Predictor (懒加载, 单例).

    优先从 table_path 加载已缓存的查找表;
    否则从 sdf_path 构建 (耗时较长, 建议缓存).

    Args:
        sdf_path: nmrshiftdb2 SDF 路径
        table_path: 已缓存的查找表文件路径
        max_records: 建表时最大记录数
        quiet: 静默模式
    """
    global _global_predictor, _global_table
    if _global_predictor is not None:
        return _global_predictor

    table = HOSELookupTable()

    if table_path and os.path.isfile(table_path):
        if not quiet:
            print(f"  [HOSE] 加载缓存查找表: {table_path}", file=sys.stderr)
        table.load(table_path)
    elif sdf_path and os.path.isfile(sdf_path):
        if not quiet:
            print(f"  [HOSE] 从 SDF 构建查找表 (首次较慢) ...", file=sys.stderr)
        table.build_from_sdf(sdf_path, max_records=max_records, quiet=quiet)
        # 缓存
        if table_path:
            os.makedirs(os.path.dirname(table_path) or '.', exist_ok=True)
            table.save(table_path)
            if not quiet:
                print(f"  [HOSE] 查找表已缓存: {table_path}", file=sys.stderr)
    else:
        if not quiet:
            print("  [HOSE] WARN: 无 SDF 或缓存文件, 预测不可用", file=sys.stderr)

    _global_table = table
    _global_predictor = C13Predictor(table)
    return _global_predictor
