# Part 3: 从谱图数据逐步推导归属表

## 核心思路转变

**旧流程**: 检索候选 → 打分排名 → 归属表（拟合）
**新流程**: 谱图信号提取 → 逐步推导每行归属 → 数据库仅提供原子编号

化学家 Excel 每行数据来源:
| 列 | 数据来源 | 推导步骤 |
|---|---|---|
| No. (编号) | 分子结构 + 编号体系 | 数据库匹配或用户指定 |
| 13C NMR | 13C peaklist (exp22) | Step 1 |
| Types | DEPT-90/135 (exp23/24) | Step 2 |
| 1H (mult., J) | 1H peaklist (exp21) 多重峰分析 | Step 3 |
| COSY | COSY 交叉峰 (exp25) | Step 4 |
| HMBC | HMBC 交叉峰 (exp27) | Step 5 |
| ROESY | ROESY 交叉峰 (exp28) | Step 6 |

## 诊断发现的关键问题

1. **HSQC 峰质量差**: 274 个原始峰 → 52 个过滤后, 仅 17/36 匹配化学家
   - 假阳性: 溶剂区 (dH 7-9 ppm) 和噪声峰未过滤
   - 假阴性: 拥挤区 (dH 1-2 ppm) 的真实 C-H 对未检出
   - 根因: `io_bruker_2d.py` 的 SNR 阈值和峰聚类参数不适配

2. **DEPT 分类失效**: 几乎所有 C-H 对标记为 "Cq"
   - 根因: DEPT peaklist 的 sign 推断逻辑有 bug (peaklist.xml 无 sign 字段)
   - 需要: 用 DEPT-90 有无 + DEPT-135 相对强度推断

3. **1H 多重峰不足**: 仅 19 个 multiplet (应 ~36)
   - 根因: `parse_multiplets_from_peaklist` 依赖 peaklist.xml 的 Multiplet 标签, 化学家未做手工多重峰分析时 fallback 不足

4. **COSY/HMBC/ROESY 列未回填**: 有原始交叉峰但未映射到具体原子标签

---

## Phase 1: 修复 HSQC C-H 对提取 (P0)

**目标**: 从 HSQC 谱可靠地提取 ~36 个真实 C-H 对

### 修改: `nmr_elucidator/io_bruker_2d.py`

- 峰检测改用**自适应阈值**: `threshold = median + k * MAD` (k=5, MAD-based robust noise)
- 新增**溶剂区过滤**: 排除已知溶剂 HSQC 交叉峰 (CD3OD: dH 3.31/dC 49.0, DMSO: 2.50/39.5, CDCl3: 7.26/77.16)
- 新增**13C 峰表约束**: HSQC 的 F1 坐标必须与 1D 13C peaklist 中的某个峰匹配 (±2 ppm), 否则视为噪声
- 新增**1H 峰表约束**: HSQC 的 F2 坐标必须与 1D 1H peaklist 中的某个峰匹配 (±0.1 ppm)

### 修改: `nmr_elucidator/interpret_2d.py` `_step1_hsqc_to_ch_pairs`

- 接收 `carbon_shifts_1d` 和 `proton_shifts_1d` 作为约束
- 对每个 HSQC 交叉峰: 仅当 F1 在 13C 峰表 ±2ppm 内 AND F2 在 1H 峰表 ±0.1ppm 内时保留
- 去重: 同一 (dC, dH) 对只保留最强峰

### 验收
clxh100: HSQC C-H 对数量 30-40 个, 与化学家 36 对匹配率 >= 80%

---

## Phase 2: 修复 DEPT 碳类型推断 (P0)

**目标**: 正确分类每个 13C 为 CH/CH2/CH3/Cq

### 修改: `nmr_elucidator/bridge_nmr2struct.py` `assign_carbon_multiplicity`

当前 bug: peaklist.xml 无 sign 字段, 导致所有碳标记为 Cq。

修复策略 (优先级从高到低):
1. **DEPT-90 有峰 → CH** (DEPT-90 只出 CH)
2. **DEPT-135 有峰 + DEPT-90 无峰 → CH2 或 CH3**
   - 用 1H 积分区分: 积分 ~2H → CH2, 积分 ~3H → CH3
   - 或用 DEPT-135 相对强度: 负峰 → CH2, 正峰 → CH/CH3
3. **DEPT-90 和 DEPT-135 均无峰 → Cq** (季碳)

### 关键: DEPT-135 sign 推断

peaklist.xml 不记录 sign。改用:
- 读取 DEPT-135 的 `pdata/1/1r` 实谱数据 (nmrglue)
- 在对应 ppm 位置取强度值: 正 → CH/CH3, 负 → CH2
- 或: 对比 DEPT-90 和 DEPT-135 峰表, DEPT-90 有而 DEPT-135 也有 → CH; DEPT-90 无而 DEPT-135 有 → CH2/CH3

### 验收
clxh100: 41 个 13C 中, 正确识别 5 个 Cq (C-5/10/13/20/22), 其余为 CH/CH2/CH3

---

## Phase 3: 增强 1H 多重峰分析 (P1)

**目标**: 为每个 C-H 对中的 H 确定 multiplicity (s/d/t/q/m) 和 J 值

### 修改: `nmr_elucidator/io_bruker.py`

当前 `parse_multiplets_from_peaklist` 仅读取 peaklist.xml 中的 `<Multiplet>` 标签。
化学家通常不做手工多重峰分析, 所以 peaklist 中无此标签。

新增 `infer_multiplicity_from_peaklist`:
- 对每个 1H 信号, 检查 peaklist.xml 中该峰的 `Peak1D` 子峰列表
- 若有多个子峰 (sub-peaks): 计算间距 → J 耦合 → 推断 d/t/q/m
- 若只有单峰: 标记为 s 或 m (根据半高宽)
- 结合 HSQC C-H 对的碳级数: CH3 → 预期 q/d, CH2 → 预期 t/m, CH → 预期 d/m

### 输出
每个 C-H 对增加: `h_multiplicity: str`, `h_j_values: List[float]`

### 验收
clxh100: H-6 (d, J=5.0), Glc-1' (d, J=7.7), H3-18 (s), H3-21 (s) 正确识别

---

## Phase 4: COSY 回填 (P1)

**目标**: 对每个 C-H 对, 列出其 COSY 相关的 H 标签

### 修改: `nmr_elucidator/interpret_2d.py`

现有 `_step2_cosy_spin_systems` 已建立自旋系统。需要新增**回填步骤**:

```python
def backfill_cosy_column(
    ch_pairs: List[CHPair],
    cosy_pairs: List[Tuple[float, float]],  # 从 COSY 谱提取的 (dH1, dH2) 对
    h_tol: float = 0.05,
) -> Dict[int, List[str]]:
    """对每个 C-H 对 (by index), 找 COSY 中与其 H 相关的其他 H, 返回标签列表."""
```

逻辑:
- 对 C-H 对 i (dH_i), 在 cosy_pairs 中找所有 (dH_i, dH_j) 或 (dH_j, dH_i)
- 将 dH_j 映射回对应的 C-H 对 j 的标签 (如 "H-2")
- 排除对角线峰 (|dH_i - dH_j| < 0.02)
- 排除溶剂峰

### 输出
归属表 COSY 列: 如 "H-1/H-2" 表示 H-1 与 H-2 有 COSY 相关

### 验收
clxh100: H-6 的 COSY 列 = "H-7" (化学家: H-6/H-7)

---

## Phase 5: HMBC 回填 (P1)

**目标**: 对每个 C-H 对, 列出其 HMBC 远程相关的 C 标签

### 修改: `nmr_elucidator/interpret_2d.py`

```python
def backfill_hmbc_column(
    ch_pairs: List[CHPair],
    hmbc_pairs: List[Tuple[float, float]],  # (dH, dC) 远程对
    carbon_shifts_1d: List[float],
    h_tol: float = 0.05,
    c_tol: float = 2.0,
) -> Dict[int, List[str]]:
    """对每个 C-H 对, 找 HMBC 中从该 H 出发的远程 C, 返回 C 标签列表."""
```

逻辑:
- 对 C-H 对 i (dH_i), 在 hmbc_pairs 中找所有 (dH_i, dC_remote)
- 排除 1J 残留: 若 dC_remote 与该 C-H 对自身的 dC 匹配 (±2ppm), 跳过
- 将 dC_remote 映射到 13C 峰表中最接近的 C 标签
- 格式: "H-4/C-2/3/5/6/10" (化学家格式: 从 H-4 看到 C-2,3,5,6,10)

### 验收
clxh100: H3-18 的 HMBC 列 = "C-12/13/14/17" (化学家: H3-18/C-12/13/14/17)

---

## Phase 6: ROESY 回填 (P2)

**目标**: 对每个 C-H 对, 列出其 ROESY 空间相关的 H 标签

### 修改: `nmr_elucidator/interpret_2d.py`

```python
def backfill_roesy_column(
    ch_pairs: List[CHPair],
    roesy_pairs: List[Tuple[float, float]],  # (dH1, dH2) 空间对
    h_tol: float = 0.05,
) -> Dict[int, List[str]]:
    """对每个 C-H 对, 找 ROESY 中与其 H 空间接近的其他 H."""
```

### 验收
clxh100: H-8 的 ROESY 列 = "H3-18/H3-19" (化学家: H-8/H3-18/H3-19)

---

## Phase 7: 归属表组装 (P0)

**目标**: 将 Phase 1-6 的结果组装成化学家格式的归属表

### 新建: `nmr_elucidator/assembly.py`

```python
@dataclass
class AssignmentRow:
    """归属表一行 (对标化学家 Excel)."""
    label: str           # 原子标签 (如 "1a", "Glc-1'") 或匿名 (如 "C-H-01")
    c_ppm: float         # 13C 化学位移
    carbon_type: str     # CH/CH2/CH3/C
    h_ppm: Optional[float]  # 1H 化学位移
    h_mult: str          # 多重性 (s/d/t/q/m/dd/...)
    h_j: List[float]     # J 耦合值 (Hz)
    cosy_partners: List[str]   # COSY 相关 H 标签
    hmbc_targets: List[str]    # HMBC 远程 C 标签
    roesy_contacts: List[str]  # ROESY 空间 H 标签
    is_ch2_stereo: bool = False  # 是否为 CH2 的立体氢 (a/b)

def assemble_assignment_table(
    ch_pairs: List[CHPair],
    dept_types: Dict[float, str],
    h_multiplets: Dict[float, Tuple[str, List[float]]],
    cosy_map: Dict[int, List[str]],
    hmbc_map: Dict[int, List[str]],
    roesy_map: Dict[int, List[str]],
    cq_shifts: List[float],  # 季碳列表
    atom_labels: Optional[Dict[int, str]] = None,  # 原子编号 (来自结构匹配)
) -> List[AssignmentRow]:
    """组装完整归属表."""
```

组装规则:
1. 每个 C-H 对生成一行 (CH2 拆为 a/b 两行)
2. 季碳单独成行 (无 1H 列)
3. 按 dC 降序排列 (或按原子编号排序, 若有)
4. 糖基识别: dC 60-105 ppm 区域的连续 CH 信号 → 标记为 Glc-1'/2'/...

### 输出
- CSV: `assignment_table.csv` (通用中间产物)
- XLSX: `assignment_table.xlsx` (化学家格式, 用 excel_export.py)

### 验收
clxh100: 输出 41 行归属 (含季碳), 与化学家 Excel 行数一致

---

## Phase 8: 原子编号 (P2, 可选)

**目标**: 为归属表中的匿名标签赋予化学家编号 (如 "1a", "Glc-1'")

### 策略

编号需要分子结构。来源:
1. **用户指定**: `--smiles <SMILES>` 或 `--mol <molfile>`
2. **数据库匹配**: lucy-ng 中 C33H52O11 有 29 个候选, 用 13C RMS 选最佳
3. **匿名模式**: 无结构时用 "C1-H1", "C2-H2a/b" 等匿名标签

### 修改: `nmr_elucidator/assembly.py`

```python
def assign_atom_labels(
    ch_pairs: List[CHPair],
    mol: Mol,  # RDKit Mol (含原子编号)
    c_tol: float = 2.0,
    h_tol: float = 0.1,
) -> Dict[int, str]:
    """将 C-H 对映射到分子结构中的原子编号."""
```

逻辑:
- 对候选结构做 13C 预测 (HOSE 或数据库值)
- 将预测 dC 与实测 dC 一对一匹配 (Hungarian 算法)
- 匹配成功后, 用候选结构的原子编号标记每个 C-H 对

---

## 实施顺序

1. **Phase 1**: 修复 HSQC C-H 对提取 (最关键, 决定后续所有步骤)
2. **Phase 2**: 修复 DEPT 碳类型推断
3. **Phase 7**: 归属表组装框架 (先搭骨架, 用匿名标签)
4. **Phase 3**: 1H 多重峰分析
5. **Phase 4-6**: COSY/HMBC/ROESY 回填
6. **Phase 8**: 原子编号 (可选)

---

## 修改文件清单

| 文件 | Phase | 改动 |
|---|---|---|
| `nmr_elucidator/io_bruker_2d.py` | 1 | 自适应阈值 + 1D 峰表约束 + 溶剂过滤 |
| `nmr_elucidator/interpret_2d.py` | 1,4,5,6 | HSQC 约束过滤 + COSY/HMBC/ROESY 回填函数 |
| `nmr_elucidator/bridge_nmr2struct.py` | 2 | DEPT sign 推断修复 |
| `nmr_elucidator/io_bruker.py` | 3 | 1H 多重峰 fallback 推断 |
| `nmr_elucidator/assembly.py` (新) | 7,8 | 归属表组装 + 原子编号映射 |
| `nmr_elucidator/__main__.py` | 全部 | 集成新流程, 新增 `--derive` 模式 |
| `nmr_elucidator/excel_export.py` | 7 | 适配新 AssignmentRow 格式 |

---

## 关键设计决策

| 决策 | 选择 | 理由 |
|---|---|---|
| HSQC 峰过滤策略 | 1D 峰表双向约束 | 比纯 SNR 阈值更可靠, 利用已有 1D 数据 |
| DEPT sign 来源 | DEPT-90 有无 + 1H 积分 | peaklist.xml 无 sign 字段, 必须间接推断 |
| 归属表排序 | 默认按 dC 降序; 有编号时按编号 | 化学家习惯按结构编号, 无结构时按 dC |
| CH2 立体氢拆分 | 仅当 1H 峰表中有两个对应 H 时拆分 | 避免过度拆分 |
| 数据库角色 | 仅提供原子编号, 不参与归属推导 | 用户明确要求"不通过打分拟合" |
| 泛化性 | 所有步骤不依赖特定化合物 | 输入为 Bruker 文件夹, 输出为通用归属表 |

---

## 验收标准 (clxh100)

| 检查项 | 期望 |
|---|---|
| HSQC C-H 对数 | 30-40 (化学家 36) |
| 与化学家匹配率 | >= 30/36 (83%) |
| DEPT 分类正确率 | >= 35/41 (85%) |
| 季碳识别 | 5 个 (C-5/10/13/20/22) |
| COSY 列非空行 | >= 15 行 |
| HMBC 列非空行 | >= 8 行 |
| ROESY 列非空行 | >= 5 行 |
| 归属表总行数 | 41 (含季碳行) |
| 输出格式 | CSV + XLSX, 8 列布局与化学家一致 |
