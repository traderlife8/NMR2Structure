# 完整实现化学家解谱工作流缺失步骤 (修正版)

## Summary

按 `doc/MestReNova_解谱指南.md` 5.1/5.4 节化学家实操, 当前程序有 8 个缺失步骤, 本计划全部实现:

| 步骤 | 缺失项 | 优先级 |
|---|---|---|
| 1.2 | 1H 多重峰分析 (J 耦合) | 高 |
| 2.3 | HSQC-DEPT 颜色区分 | 中 |
| 3.2 | TOCSY 长程自旋系统 | 中 |
| 4.3 | HMBC 2J/3J 键数判断 | 高 |
| 5.1 | NOESY 谱读取 | 高 |
| 5.2 | NOESY 空间近邻 H-H 提取 | 高 |
| 5.3 | 候选 3D 构象 NOE 验证 | 高 |
| 6.2 | 归属完整性自动核对 | 高 |

附加: README 对照表 + AGENTS.md + MestReNova 指南同步.

---

## Part 1: README.md 文档补充

### 修改文件: `NMR2Structure/README.md`

新增 "## 解谱工作流对标" 章节, 含 20 行对照表 (5.4 六步 + 5.1 四步 + 0.3 五步).

**插入位置**: "## 主流程 CLI" 之后, "### 必填" 之前.

---

## Part 2: 模块 A -- 1H 多重峰分析 (J 耦合)

### 修改文件: `nmr_elucidator/io_bruker.py`

新增 `Multiplet` 数据类 + `parse_multiplets_from_peaklist()`:

```python
@dataclass
class Multiplet:
    center_ppm: float
    integral: float
    j_couplings: List[float] = field(default_factory=list)
    multiplicity: str = ''   # 's'/'d'/'t'/'q'/'m'
    pattern_type: str = ''
```

读取策略: Bruker peaklist.xml 的 `<multiplet>` 标签, regex + 容错 try-except 兼容 TopSpin 2/3/4.

### 修改文件: `nmr_elucidator/interpret_2d.py`

**CHPair 新增字段**: `j_couplings: List[float]`, `h_neighbors: List[float]`.

**`_step2_cosy_spin_systems` 增强**:
- 对每个 COSY 交叉峰 (dHa, dHb), 检查 Ha 的 multiplet 中 J 值是否指向 Hb
- 仅当 |J| > 4 Hz 才视为真实邻位耦合 (COSY 3J 通常 4-12 Hz)
- J 值来自 1D peaklist multiplet, COSY 交叉峰来自 2D 谱, 两者通过 dH 匹配关联

**输出**: `step1_hsqc_ch_pairs.csv` 新增 `J_couplings_Hz` / `h_neighbors_ppm` 列.

---

## Part 3: 模块 B -- HSQC-DEPT 颜色区分

### 修改文件: `nmr_elucidator/interpret_2d.py`

`_draw_hsqc_contour_all_pairs` 函数:
- 对每个 C-H 对, 根据其 multiplicity 单独设 contour 颜色:
  - CH/CH3 -> 红色 (传统约定, 正峰)
  - CH2 -> 蓝色 (反相, 负峰)
  - Cq/unknown -> 灰色
- 实现方式: 不用 cmap (cmap 是连续映射), 改为对每个 C-H 对的 contour 调用设 `colors=` 参数
- 图例标注: 右上角添加 "CH/CH3=红, CH2=蓝" 图例

**输出**: `step1_hsqc_ch_pairs/{overlap,normal}/*.png` 颜色升级.

---

## Part 4: 模块 C -- TOCSY 长程自旋系统

### 修改文件: `nmr_elucidator/interpret_2d.py`

复用 `_step2_cosy_spin_systems` 框架, 新增 `_step2_tocsy_spin_systems`:

**算法**: 与 COSY 同 (Union-Find), 但容差放大:
- COSY h_tol = 0.3 ppm (仅 3J 直接邻位)
- TOCSY h_tol = 0.6 ppm (整个自旋系统, 通常 4-7 个 H)

**自动识别**: 实验类型为 'TOCSY' 时调用, COSY 时不调用.

**与 COSY 合并策略**:
- 若同一自旋系统 COSY 和 TOCSY 都检出, 合并为一条, 加 `source` 列 ('COSY'/'TOCSY'/'COSY+TOCSY')
- TOCSY 检出的额外 H (COSY 未覆盖) 追加到自旋系统

**输出**: `step2_tocsy_spin_systems.csv` (与 step2_cosy 格式一致, 含 source 列).

---

## Part 5: 模块 D -- HMBC 2J/3J 键数判断

### 修改文件: `nmr_elucidator/interpret_2d.py`

**HMBCConnection 新增字段**:
- `bond_type: str` -- '2J' / '3J' / 'long_range' / ''
- `inferred_distance: int` -- 推断键数 (2/3)
- `inference_basis: str` -- 'topology' / 'unresolved'

**新增函数 `_annotate_hmbc_bond_type`**:
```python
def _annotate_hmbc_bond_type(
    hmbc_connections: List[HMBCConnection],
    ch_pairs: List[CHPair],
    spin_systems: List[SpinSystem],
) -> List[HMBCConnection]:
    """基于自旋系统拓扑推断 2J/3J.

    规则:
    - 若 H 的 HSQC C 在同一自旋系统 -> 距离 1-2 -> likely 2J
    - 若 H 与 C 在不同自旋系统 -> 跨系统 -> likely 3J
    - 自旋系统不确定 -> long_range
    """
```

**集成**: `_step3_hmbc_connections` 末尾调用.

**输出**: `step3_hmbc_connections.csv` 新增 `bond_type` 列.

**打分**: `evidence.score_candidates` 中 `bond_type=2J` 权重略高于 `3J` (仅记录, 不直接改分).

---

## Part 6: 模块 E -- NOESY 立体化学

### E.1 NOESY 谱读取

**修改文件**: `nmr_elucidator/io_bruker_2d.py`

实验类型识别加入 'NOESY' / 'ROESY' (轴均为 1H x 1H):
- `experiment_type` 推断: PULPROG 含 'noesy'/'roesy' 关键字
- 数据结构沿用 `Spectrum2D` (f1_ppm / f2_ppm 均为 1H 维)

**新增函数** `_classify_noesy_2d()`:
- 从 PULPROG 文本匹配 'noesygpph' / 'roesyph' 等关键字
- 区分 NOESY vs ROESY

### E.2 NOE 空间近邻提取

**修改文件**: `nmr_elucidator/interpret_2d.py`

**新增数据结构**:
```python
@dataclass
class NOEContact:
    h1_ppm: float
    h2_ppm: float
    intensity: float
    confidence: str        # 'strong' / 'medium' / 'weak'
    expected_distance: float  # 期望距离 (A), 基于强度-距离映射
    is_symmetric: bool
```

**NOE 强度-距离映射** (化学家约定):
| 置信度 | 强度范围 | 期望距离 (A) |
|---|---|---|
| strong | top 30% | < 2.5 |
| medium | 30-70% | 2.5 - 3.5 |
| weak | bottom 30% | 3.5 - 5.0 |

**新增函数 `_step5_extract_noe_contacts`**:
- 对角线附近裁剪 (排除主对角线 +/- 0.1 ppm)
- 检测交叉峰, 每个交叉峰生成 NOEContact
- 对称化: (h1, h2) 和 (h2, h1) 合并为一条

**输出**: `step5_noe_spatial_contacts.csv`

### E.3 候选 3D 构象 NOE 验证

**修改文件**: `nmr_elucidator/correlation_2d.py`

**新增数据结构**:
```python
@dataclass
class NOEVerificationResult:
    n_total: int = 0
    n_matched: int = 0
    n_conflict: int = 0
    rms_distance: float = 0.0
    conflicts: List[dict] = field(default_factory=list)
    best_conformer_idx: int = 0
    verification_status: str = ''  # 'good' / 'partial' / 'conflict' / 'fallback_2d'
```

**Correlation2DScore 新增字段**:
- `noe_match_rate: Optional[float] = None`
- `noe_conflict_count: int = 0`
- `noe_rms_distance: Optional[float] = None`

**新增函数 `verify_candidate_3d_noe`**:
```python
def verify_candidate_3d_noe(
    candidate_smiles: str,
    noe_contacts: List[NOEContact],
    h_assignments: dict,  # {dH_ppm: atom_index} 从归属表构建
    max_distance: float = 5.0,
    num_confs: int = 10,
) -> NOEVerificationResult:
    """对候选结构用 RDKit 嵌入多构象, 验证 NOE 约束.

    算法:
    1. RDKit Chem.AddHs -> AllChem.EmbedMultipleConfs(mol, numConfs=10)
    2. 对每个构象: MMFFOptimizeMolecule -> 力场优化
    3. 遍历每个 NOEContact:
       - 根据 h_assignments 找到两个 H 原子 index
       - 计算 3D 距离
       - 与 expected_distance 比较 (容差 +/- 1.0 A)
       - 距离 > max_distance -> 标为冲突
    4. 取 NOE 命中率最高的构象为 best_conformer
    5. 输出: 命中率, RMS 距离偏差, 冲突列表

    h_assignments 构建方式:
    - 从 FullAssignmentTable.rows 提取 (h_ppm, atom_label)
    - atom_label (如 'H-3') 映射到 RDKit mol 的 atom index
    - 映射方法: mol.GetAtomWithIdx(i).GetSymbol() == 'H' 且
      其重原子邻居的 atom_label 匹配
    """
```

**输出**: `step6_noe_3d_verification.csv` (Top-N 各一行)

### E.4 集成到打分公式

**修改文件**: `nmr_elucidator/evidence.py`

`score_candidates` 新增 `use_noe_scoring: bool = False` 参数:

```python
if use_noe_scoring and s2d.noe_match_rate is not None:
    noe_bonus = s2d.noe_match_rate * 10 - 5 * s2d.noe_conflict_count
else:
    noe_bonus = 0
final_score = base_score + noe_bonus
```

### E.5 集成到管线 (修正时序)

**修改文件**: `nmr_elucidator/__main__.py`

**正确时序** (关键修正):
```
1. interpret_2d_sequential -> 含 _step5_extract_noe_contacts
2. retrieve_candidates
3. build_full_assignment_table (Top-N)  <-- 归属表在此生成
4. 构建 h_assignments 映射 (从归属表 rows)
5. verify_candidate_3d_noe (Top-N)     <-- NOE 验证在归属表之后
6. score_candidates(use_noe_scoring=True)
7. meta 新增 noe_contacts / noe_verification
```

**h_assignments 构建**:
```python
# 从 FullAssignmentTable 构建 dH -> atom_index 映射
h_assignments = {}
for row in table.rows:
    if row.h_ppm is not None and row.atom_label:
        # atom_label 如 'H-3', 找到 RDKit mol 中对应 H 的 index
        h_assignments[row.h_ppm] = _find_h_atom_index(mol, row.atom_label)
```

---

## Part 7: 模块 F -- 归属完整性自动核对

### 修改文件: `nmr_elucidator/interpret_2d.py`

**新增 `AssignmentValidation` 数据类** + `validate_assignment` 函数.

核对规则:
| 规则 | 触发警告 |
|---|---|
| 13C 中 CH/CH2/CH3 无 HSQC 对应 | "可能的季碳未被识别" |
| 1H 中有峰但无 HSQC 对应 | "活泼氢: NH/OH" |
| HSQC 中 C 无 1H 积分对应 | "HSQC 假阳峰" |
| 总积分 != 分子式 H 数 | "积分与分子式不符" |

**输出文件**: `step7_assignment_validation.csv`

---

## step 编号规划 (修正)

| 文件 | 内容 |
|---|---|
| step1_hsqc_ch_pairs.csv | HSQC C-H 对 (含 J 耦合列) |
| step1_hsqc_ch_pairs/{overlap,normal}/ | 等高线图 (含 DEPT 颜色) |
| step1b_overlap_analysis.csv | 重叠峰分解 |
| step2_cosy_spin_systems.csv | COSY 自旋系统 |
| step2_tocsy_spin_systems.csv | TOCSY 自旋系统 (若有) |
| step3_hmbc_connections.csv | HMBC 远程连接 (含 bond_type) |
| step4_quaternary_carbons.csv | 季碳推断 |
| step5_noe_spatial_contacts.csv | NOE 空间近邻 (若有 NOESY) |
| step6_noe_3d_verification.csv | NOE 3D 构象验证 (Top-N) |
| step7_assignment_validation.csv | 归属完整性核对 |
| 完整NMR归属表(论文Table1)/ | Top-N 归属表 CSV |

---

## 修改文件清单

| 文件 | 改动 |
|---|---|
| `NMR2Structure/README.md` | 新增"解谱工作流对标"章节 |
| `nmr_elucidator/io_bruker.py` | `Multiplet` + `parse_multiplets_from_peaklist` |
| `nmr_elucidator/io_bruker_2d.py` | 实验类型识别加入 NOESY/ROESY |
| `nmr_elucidator/interpret_2d.py` | 模块 A/C/D/E.2/F: CHPair 加 J 字段, TOCSY, HMBC bond_type, NOE 提取, validate_assignment |
| `nmr_elucidator/correlation_2d.py` | 模块 E.3: `NOEVerificationResult` + `verify_candidate_3d_noe` + Correlation2DScore 新增 NOE 字段 |
| `nmr_elucidator/evidence.py` | `score_candidates` 新增 NOE 加分 |
| `nmr_elucidator/__main__.py` | 接线所有新模块 + 修正 NOE 验证时序 + meta 字段扩展 |
| `nmr_elucidator/report.py` | 2.5 节新增 step5/6/7 路径索引 |
| `nmr_elucidator/__init__.py` | 导出新 API |
| `AGENTS.md` | 数据流 + 模块职责表同步 |
| `doc/MestReNova_解谱指南.md` | 5.5 节补充 NOESY/TOCSY 程序化实现说明 |

---

## 实施顺序

1. README 文档补充
2. 模块 A: 1H 多重峰分析 (io_bruker + interpret_2d)
3. 模块 B: HSQC-DEPT 颜色 (interpret_2d)
4. 模块 C: TOCSY (interpret_2d)
5. 模块 D: HMBC 2J/3J (interpret_2d)
6. 模块 E.1: NOESY 谱读取 (io_bruker_2d)
7. 模块 E.2: NOE 提取 (interpret_2d)
8. 模块 E.3: 3D 构象 NOE 验证 (correlation_2d)
9. 模块 E.4-E.5: 集成 + 打分 + 时序修正 (__main__ + evidence)
10. 模块 F: 归属完整性核对 (interpret_2d)
11. report.py: 2.5 节路径索引
12. __init__.py: 导出
13. AGENTS.md + MestReNova 指南同步
14. 验证 clxh100 (有 2D) + clxh3 (无 2D 回归)

---

## 关键设计决策

| 决策 | 选择 | 理由 |
|---|---|---|
| peaklist.xml multiplet 解析 | regex 容错 | 跨 TopSpin 2/3/4 版本兼容 |
| J 与 COSY 关联 | 对每个 COSY 交叉峰 (dHa, dHb), 检查 Ha multiplet 中 J 是否指向 Hb | 1D J 值验证 2D COSY 交叉峰真实性 |
| TOCSY 容差 | 0.6 ppm | 比 COSY 大, 覆盖整个自旋系统 |
| TOCSY/COSY 合并 | 同一自旋系统合并, 加 source 列 | 避免重复计数 |
| HMBC 2J/3J 判定 | 自旋系统拓扑 | 候选未定时唯一可用信息 |
| HSQC-DEPT 颜色 | 对每个 C-H 对单独设 contour colors 参数 | cmap 是连续映射, 不适用于分类着色 |
| NOE 距离映射 | strong < 2.5A, medium 2.5-3.5A, weak 3.5-5.0A | 化学家 r^-6 约定 |
| 3D 构象生成 | EmbedMultipleConfs(numConfs=10) + 取最佳 | 柔性分子需多构象 |
| NOE 验证时序 | retrieve -> build_assignment -> h_assignments -> verify_noe -> score | 归属表是 NOE 验证的前提 |
| NOE 加分公式 | 命中率*10 - 冲突*5 | 与现有打分同尺度 |
| 完整性核对分子式 | 不强制 | 用户可能未知分子式 |
| step 编号 | step5=NOE, step6=3D验证, step7=完整性核对 | 避免 step5/5b/5_noe_3d 混乱 |

---

## 测试验证

### clxh100 (有 2D, --use-2d --top-n 3)
- step1 CSV: 含 J_couplings_Hz 列 (若有 multiplet 数据)
- step2_cosy: 49 个系统 + J 验证
- step2_tocsy (若有 TOCSY 谱)
- step3: 含 bond_type 列 (2J/3J/long_range)
- step5: NOE 空间近邻表 (若有 NOESY 谱)
- step6: NOE 3D 验证 (Top-N)
- step7: 归属完整性核对
- 完整NMR归属表/rank*.csv: 新增 NOE 验证列
- HTML 2.5 节: 10 个文件路径索引

### clxh3 (无 2D)
- 行为完全不变 (新增字段为空)

### 单元测试
- Multiplet 解析: 构造 mock peaklist.xml
- NOE 提取: 构造 mock NOESY 谱
- 3D 验证: 已知小分子 (甲醇) -> RDKit 嵌入 -> 距离矩阵正确
- h_assignments 构建: 从 FullAssignmentTable 构建映射

---

## 风险与缓解

| 风险 | 缓解 |
|---|---|
| peaklist.xml 跨版本差异 | regex + try-except 静默失败, 不阻断主流程 |
| NOESY 谱可能不含 | 检测到缺失时跳过 E.1-E.5 |
| 3D 嵌入失败 (RDKit) | fallback: 用 2D 拓扑距离估算, 标 verification_status='fallback_2d' |
| TOCSY 与 COSY 冲突 | 合并到 Step 2 单一表, 加 source 列区分 |
| 归属核对与候选检索时序 | validation 在 interpret_2d 末尾, retrieval 前 |
| NOE 验证计算开销 | EmbedMultipleConfs(10) ~1s/分子; 仅 Top-N (10 个) 验证 |
| 1H 多重峰版本不兼容 | 找不到 multiplet 标签时留空, 不报错 |
| NOE 验证依赖归属表 | 时序修正: 先 build_assignment 再 verify_noe |
| h_assignments 映射不完整 | 未归属 H 跳过 NOE 验证, 标 'unassigned_h_skipped' |