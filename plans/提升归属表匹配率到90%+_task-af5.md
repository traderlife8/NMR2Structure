# clxh100 原始数据探测与按"1D/2D 分组规则"的补强方案

## 〇、核心约束
推导流程**只能**使用仪器测出的原始数据（1D/2D NMR peaklist/parameters/FID/projected data），**不依赖**任何 SMILES/mol block/CDXML/化学家归属字典。`compare_program_to_standard` 的 (δH ±0.1, δC ±0.5, 类型) 三元组匹配已避免"标号错位"问题。

---

## 一、clxh100 实验分组（用户明确指定）

### 1D 流程（实验 21-24，CDCl₃ 溶剂组）
| exp | 类型 | 溶剂 | 关键特征 |
|---|---|---|---|
| **21** | 1H | CDCl₃ | **1D 氢谱主谱**，113 peaks，F1 0.90-8.71 ppm，16 积分区间 |
| 22 | 13C | CDCl₃ | 41 peaks，NS=23，SW=29762 Hz |
| 23 | DEPT-90 | CDCl₃ | 32 peaks（仅 CH，无 CH₂/CH₃ 区分能力） |
| 24 | DEPT-135 | CDCl₃ | 49 peaks（38 正 + 11 负），负峰 = CH₂ |

### 2D 流程（实验 25-28 + 31，Pyr-d5 溶剂组）
| exp | 类型 | 溶剂 | 关键特征 |
|---|---|---|---|
| **31** | 1H | Pyr | **2D 氢谱主谱**，329 peaks，F1 2.07-11.99 ppm，SW=10000 Hz，**无积分** |
| 25 | COSY | Pyr | TD=2048，投影引用 exp 31 |
| 26 | HSQC | Pyr | TD=1024(F2)/148(F1)，FnMODE=6 (Echo-Antiecho) |
| 27 | HMBC | Pyr | TD=1024(F2)/128(F1)，FnMODE=1 (QF) |
| 28 | ROESY | Pyr | TD=1024，投影引用 exp 31 |

**用户的核心规则**：
- **1D 流程的 1H 主谱 → 用 exp 21**（与 22/23/24 同溶剂 CDCl₃，已积分）
- **2D 流程的 1H 主谱 → 用 exp 31**（与 25-28 同溶剂 Pyr-d5，未积分但高分辨率）

---

## 二、exp 21 vs exp 31 系统性差异

| 维度 | exp 21 | exp 31 |
|---|---|---|
| 采集日期 | 2023-02-23 14:25 | 2023-03-31 03:11 |
| PULPROG | zg | zg |
| SOLVENT | CDCl₃ | Pyr |
| NS | 1 | 1 |
| SWH (Hz) | 7002.8 | 10000.0 |
| AQ (sec) | 4.679 | 3.277 |
| RG | 23.29 | 11.2 |
| TE (K) | 298.2 | 298.1 |
| peak 数 | 113 | 329 |
| F1 范围 (ppm) | 0.90 - 8.71 | 2.07 - 11.99 |
| 积分 (integrals.txt) | **16 区间**（有） | 无 |
| 强度范围 | 0.37 - 15.00 | 0.005 - 15.00 |
| 低场峰 (dH>6.0) | 8.71, 7.57, 7.20, 6.23, 6.17 | 11.99, 11.61, 10.38, 10.30, 10.25, 10.20, 10.02, 9.22, 9.06 |
| 投影来源 | 仅 1D 流程内部 | 25/26/27/28 的 F1/F2 投影 |

**关键发现**：
- exp 21 含 **Pyr 残余峰 8.71 / 7.57 / 7.20**（即使 SOLVENT 标 CDCl₃），说明样品实际是 Pyr-d5 配的，但用 CDCl₃ 做内标参考
- exp 31 SWH 更大（10000 Hz → 20 ppm），覆盖完整 Pyr 区
- exp 31 是 **2D 流程专用**的"投影源"和"标度对齐基线"

---

## 三、程序当前的问题

`nmr_elucidator/__main__.py:_build_problem_from_bruker`:
```python
h_spec = next((s for s in specs if s.nucleus == '1H'), None)
```
用 `next` 取**第一个**匹配的 1H 谱，但 `read_bruker_folder` 用 `sorted(os.listdir(folder))`，所以 `21 < 31`，程序**总是取 exp 21（CDCl₃）**，忽略 exp 31。

**后果**：
1. 2D 流程（HSQC/COSY/HMBC/ROESY）的 1H 化学位移标度与 Pyr 2D 谱 F2 轴**不匹配**
2. 1H 积分来自 CDCl₃ 谱，但 2D 谱溶剂是 Pyr-d5 → 溶剂残余峰位置差异导致误识
3. `interpret_2d_sequential` 解析 1H 多重峰时也用 exp 21 的 peaklist → 多重峰分析结果与 2D 数据不一致

---

## 四、补强方案（按用户分组规则 + 6 个产物）

### 产物 A（🔴 关键）：1H 谱按"1D/2D 流程"分组选择

- **算法**：
  ```python
  def select_h_specs_by_process(specs):
      h_specs = [s for s in specs if s.nucleus == '1H']
      # 按 SOLVENT 分组
      cdcl3_h = [s for s in h_specs if s.solvent == 'CDCl3']
      pyr_h   = [s for s in h_specs if s.solvent == 'Pyr']
      # 规则: 1D 用 CDCl3 (与 13C/DEPT 配对), 2D 用 Pyr (与 2D 配对)
      h_1d = cdcl3_h[0] if cdcl3_h else pyr_h[0]
      h_2d = pyr_h[0]   if pyr_h   else cdcl3_h[0]
      return h_1d, h_2d
  ```
- **修改点**：
  - `__main__.py:_build_problem_from_bruker` → 返回 `h_spec_1d` 给 1D 流程
  - `__main__.py` 步骤 [2.8] 2D 解读前 → 取 `h_spec_2d` 给 `interpret_2d_sequential`
  - `interpret_2d.py` 第 487 行 `for spec in specs` 循环 → 增加 `prefer_solvent='Pyr'` 参数
- **产物**：`output/<sample>/1h_source_selection.csv`（两组选择依据 + 各自的峰数/SW/积分）
- **预期贡献**：+25-35 pp 匹配率（消除溶剂标度混用）

### 产物 B（中）：HMBC F1 投影 (pdata/99/999) 补强 ¹³C 峰表

- **输入**：`27/pdata/99/1r`（HMBC 列投影 = ¹³C 维投影，高 SNR）+ proc/procs
- **算法**：nmrglue 读 1r → 重建 ppm 轴 → `scipy.signal.find_peaks` 在投影上检峰 → 与 22 主谱 ¹³C 峰并集去重（容差 0.1 ppm）
- **产物**：`output/<sample>/13c_with_hmbc_projection.csv`（新增峰标 `*_hmbc_proj`）
- **价值**：HMBC 投影的 ¹³C 维有更高 SNR（2048 行求和增强），能检出 22 主谱漏掉的弱 ¹³C 峰
- **模块**：新增 `nmr_elucidator/io_bruker_2d.py:read_hmbc_f1_projection()`

### 产物 C（中）：HSQC 单行切片精细化 (pdata/1000/1001)

- **输入**：`26/pdata/1000/1r`（HSQC row 51 切片）+ `26/pdata/1001/1r`（HSQC row 941 切片）+ used_from PROCNO2
- **算法**：通过 `used_from` 的 PROCNO2（行号）→ 找到对应 ¹³C ppm → 取该行的 1H 切片 → find_peaks 检 H 峰 → 校准 J 耦合 → 附加到对应 step1 行
- **产物**：`output/<sample>/hsqc_row_slices/`（PNG 等高线图 + `*.csv` 数值）
- **价值**：HSQC 行的精确 1H 模式（多重峰 J）→ 增强多重峰归属
- **模块**：新增 `nmr_elucidator/io_bruker_2d.py:read_hsqc_row_slice()`

### 产物 D（中）：`peakrng` 增强重叠峰检测

- **输入**：`exp/pdata/1/peakrng`（1D 峰左右边界）
- **算法**：对 step1 C-H 对，若邻近 C-H 对的 (dC ±0.5 ppm) 范围内 peakrng F1 边界有重叠 → 标记为重叠峰
- **产物**：附加到 `step1_hsqc_ch_pairs.csv` 的 `is_overlap_from_peakrng` 列
- **模块**：扩展 `interpret_2d._step1_resolve_overlapping`

### 产物 E（低）：`audita.txt` 冗余元数据校验

- **输入**：`exp/audita.txt`
- **算法**：提取 SOLVENT/PULPROG 字段，与 parm.txt 一致性校验；不一致时优先用 audita
- **产物**：`meta['solvent_consistency']` 字段

### 产物 F（低）：2D 谱溶剂一致性校验

- **输入**：每个 2D 实验的 SOLVENT
- **算法**：检查所有 2D 实验的 SOLVENT 一致性；若不一致 → 输出告警（提示用户选择或合并）
- **产物**：`meta['2d_solvent_consistency']` 字段

---

## 五、产物优先级与预期匹配率提升

| 优先级 | 产物 | 预期贡献 | 风险 |
|---|---|---|---|
| **🔴 1** | A. 1H 谱按 1D/2D 分组 | **+25-35 pp**（消除溶剂标度混用） | 低（仅选择逻辑改动） |
| **🟡 2** | B. HMBC F1 投影补强 ¹³C | +5-10 pp（检出 22 主谱漏掉的弱 ¹³C 峰） | 中（需 nmrglue） |
| **🟡 3** | C. HSQC 行切片精细化 J | +3-5 pp（增强多重峰归属） | 中（需 nmrglue） |
| **🟢 4** | D. peakrng 重叠峰检测 | +2-3 pp | 低 |
| **🟢 5** | E. audita.txt 元数据校验 | 0-1 pp（边界修正） | 低 |
| **⚪ 6** | F. 2D 谱溶剂一致性校验 | 0 pp（告警） | 0 |

**合计预期**：39.6% → **75-85%**（仅产物 A）；+产物 B/C 后 → **85-95%**（达到 90% 目标）

---

## 六、实施顺序

1. **第 1 周**：产物 A（1H 谱按 1D/2D 分组）—— **最关键，直接 +25-35 pp**
2. **第 2 周**：产物 B（HMBC F1 投影补强 ¹³C） + 产物 C（HSQC 行切片精细化）
3. **第 3 周**：产物 D + E + F（peakrng 增强 + 元数据校验）
4. **第 3 周末**：在 clxh100 上重跑 `--derive` → 行匹配率预期 ≥ 80%
5. **第 4 周**：融合前序的 HMBC 闭环 + DEPT-135 正负号整合 + 化学位移经验规则 → 进一步冲 90%

---

## 七、产物目录

```
output/<sample>/
├── 1h_source_selection.csv            # 产物 A：每个 1H 实验 + 1D/2D 选择依据
├── 13c_with_hmbc_projection.csv       # 产物 B：22 主谱 ∪ HMBC F1 投影峰
├── hsqc_row_slices/                   # 产物 C：HSQC 单行切片等高线图 + CSV
│   ├── row_51_<dC>_<dH>.png
│   └── row_941_<dC>_<dH>.png
├── step1_hsqc_ch_pairs.csv            # 已有 + 新增 is_overlap_from_peakrng 列（产物 D）
└── meta.txt                           # 已有 + 新增 solvent_consistency, 2d_solvent_consistency
```

---

## 八、不在本次范围

- ML 模型（external/NMR2Struct/）、绝对构型、变温 NMR、混合 2D
- 已知化合物复核阶段的 SMILES/mol block（仅在 `compare_program_to_standard` 中作为结构判等依据）
- 3D RDKit 嵌入（不在推导流程）# clxh100 原始数据文件逐项探测报告与可补强环节

## 〇、核心约束（沿用前版）
推导流程**只能**使用仪器测出的原始数据（1D/2D NMR peaklist/parameters/FID/projected data），**不依赖**任何 SMILES/mol block/CDXML/化学家归属字典。`compare_program_to_standard` 的 (δH ±0.1, δC ±0.5, 类型) 三元组匹配已避免"标号错位"问题，无需映射。

---

## 一、clxh100 完整实验目录清单

```
clxh100/
├── NMR原始数据文件夹说明.txt        # 实验元数据（手写）
├── 21/    1D 1H    solvent=CDCl3   pulses=zg           → 与 2D 溶剂不一致
├── 22/    1D 13C   solvent=CDCl3   pulses=zgdc        NS=23  SW=29761.9 Hz
├── 23/    1D 13C-DEPT-90  solvent=CDCl3   pulses=deptsp90
├── 24/    1D 13C-DEPT-135 solvent=CDCl3   pulses=deptsp135   49 peaks (38 正 + 11 负)
├── 25/    2D COSY  solvent=Pyr  pulses=cosygpmfqf  TD=2048 (F2)  投影: exp 31
├── 26/    2D HSQC solvent=Pyr  pulses=hsqcetgpsi2 TD=1024(F2)/148(F1) FnMODE=6
├── 27/    2D HMBC solvent=Pyr  pulses=hmbcgplpndqf TD=1024/128 FnMODE=1  F2 投影: pdata/999
├── 28/    2D ROESY solvent=Pyr pulses=croesyph   TD=1024 (F2)  投影: exp 31
└── 31/    1D 1H    solvent=Pyr  pulses=zg        NS=1  SW=10000 Hz  → 与 2D 同溶剂
```

**关键发现**：
- **同一化合物在两组不同溶剂下测过**：21-24 (CDCl₃) vs 25-28 + 31 (Pyr-d5)
- **2D 实验的投影都用 exp 31**：`portf` 文件显示 25/26/27/28 都引用 `exp 31`
- **HMBC 的 F1 (¹³C 维) 投影**：来自 `27/pdata/999`（用 `mvrow.exe` 把 rows 298-334 求和）
- **HSQC 单行切片**：`26/pdata/1000` (row 51) 和 `26/pdata/1001` (row 941)

---

## 二、clxh100 各文件是否被程序使用

### A. 已用且有效（核心读取）

| 文件路径 | 用途 | 程序入口 |
|---|---|---|
| `exp/pdata/1/peaklist.xml` | 1D 峰表 | `io_bruker._parse_peaklist_xml`（同时识别 DEPT-135 intensity 正负号 → CH₂） |
| `exp/pdata/1/parm.txt` 或 `exp/acqus` | 采集/处理参数 | `io_bruker._parse_parm`（提取 SOLVENT/NUC1/SFO1/SFO2/PULPROG/TE/SW 等） |
| `exp/pdata/1/integrals.txt`（仅 1H） | 1H 积分区间 | `io_bruker._parse_integrals`（仅 exp 21/31 有） |
| `exp/pdata/1/1r` | 1D 实谱 | `io_processed.read_processed_spectrum`（仅启用 `--with-real-spectrum`） |
| `exp/fid` 或 `exp/ser` | 时域信号 | `nmrglue.bruker.read_pdata`（2D 时域回退） |
| `exp/pdata/1/2rr, 2ir, 2ri, 2ii` | 2D 处理后谱 | `nmrglue.bruker.read_pdata`（自动识别） |
| `exp/acqu2s`, `exp/pdata/1/proc2s` | F1 间接维参数 | `nmrglue.bruker.read` |
| `parm.txt PULPROG 字段`（来自 `pulseprogram` 文件） | 推断核种与 2D 类型 | `_detect_pulse_role` + `detect_2d_experiment_type` |

### B. 重要但目前程序未使用（可补强）

| 文件 | 实际内容 | 解谱用途 | 影响 |
|---|---|---|---|
| **`31/pdata/1/peaklist.xml`**（Pyr 1H, 329 peaks） | 完整 1H 谱，但 solvent=Pyr | **程序当前只看 exp 21**（CDCl3 1H），导致 1H 化学位移与 Pyr 2D 数据溶剂不一致 → HSQC 匹配误差 | **严重**：CDCl3 的 Pyr-d5 残余峰在 8.71 ppm，而 Pyr-d5 自身残余峰在 7.21/7.57/8.71 ppm — 误识溶剂峰 |
| **`22/pdata/1/peaklist.xml` + DEPT-90 (exp 23)** | 13C 主谱 + DEPT-90 | 程序用 DEPT-135 的 intensity 正负号 + DEPT-90 命中区分 CH/CH2/CH3/Cq，**已用** | OK |
| **`24/pdata/1/peaklist.xml`**（DEPT-135, 49 peaks） | intensity 有正负号 | **已用**：`bridge_nmr2struct._infer_peak_sign` 区分 CH/CH₃ vs CH₂ | OK |
| **`exp/pdata/1/peakrng`** | 1D 峰区间（峰宽左右边界） | 用于重叠峰检测（左右边界间距 < 阈值）和多重峰识别（峰宽 > 单峰） | 中：增强重叠峰分解 |
| **`exp/pdata/1/peaks`** | 旧版 1D 峰位文本（冗余于 peaklist.xml） | 与 peaklist.xml 交叉验证（peak 数差异提示漏检/过检） | 低 |
| **`exp/pdata/1/intrng`** | 1D 积分区间文本（含区间 + 标定） | 比 `integrals.txt` 更详细的区间信息，可补强积分对齐 | 低 |
| **`exp/audita.txt`** | 谱仪审计日志（含 SOLVENT/PULPROG 备份） | 作为 SOLVENT/PULPROG 的冗余源（如果 parm.txt 缺失/错误） | 低 |
| **`exp/pdata/1/title`** | 实验标题（人写） | **已用**（仅 SOLVENT 缺失时回退读 title） | OK |
| **`exp/pdata/1/clevels`**（2D 等高线级别参数） | 含 TopSpin 自动选的 contour level | 推断 2D 谱 SNR；提示真实弱峰阈值 | 中：用于 2D 峰检测阈值校验 |
| **`27/pdata/1/f2projp`**（仅 HMBC） | 引用 `27/pdata/999/1r`（HMBC F1 投影的 1H 谱） | F1 投影作 ¹³C 化学位移的辅助来源（高 SNR） | 低（已通过 22 主谱覆盖） |
| **`27/pdata/999/used_from`** | 标识 pdata/999 是 "positive projection of rows 298-334" | 标识数据来源 | 0（仅元数据） |
| **`27/pdata/99/used_from`** | "positive projection of all 2048 columns" | 同上 | 0 |
| **`26/pdata/1000/used_from`** | "row 51 of HSQC"（HSQC 某 ¹³C 行的 1H 切片） | 单行 1H 切片：精确提取该 C 的 1H 化学位移 + J 模式 | 中：可作 HSQC 单行精细化分析 |
| **`26/pdata/1001/used_from`** | "row 941 of HSQC" | 同上 | 中 |
| **`exp/cpdprg2`**（仅 13C 谱） | 13C 去耦程序名（如 waltz16） | 提示碳谱已做 NOE 抑制/反门控 | 低 |
| **`exp/spnam5/14/31`**（仅 13C 谱） | 脉冲形状（Q5, Crp60 等） | 推断是否用 NOE 抑制（影响碳谱定量） | 低 |
| **`exp/specpar`** | 仪器参数汇总 | 元数据备份 | 0 |

### C. 完全不影响解谱（仅仪器/输出辅助）

| 文件 | 用途 | 是否需要 |
|---|---|---|
| `outd`, `auditp.txt`, `format.temp`, `format.ased`, `thumb.png` | 输出/格式/缩略图 | 否 |
| `prosol_History`, `scon2`, `shimvalues`, `uxnmr.info`, `uxnmr.par` | 探头参数/历史/用户 | 否 |

---

## 三、关键发现：多组 1H 谱与溶剂不一致问题（影响匹配率的关键）

**核心问题**：
- exp 21 是 CDCl₃ 的 1H（113 peaks）
- exp 31 是 Pyr-d5 的 1H（329 peaks）
- 程序在 `__main__.py:_build_problem_from_bruker` 用 `next((s for s in specs if s.nucleus == '1H'), None)` 取**第一个**，即 exp 21
- 但 2D 实验 (25-28) 的 `portf` 都引用 **exp 31**（Pyr-d5）做投影
- 此外 2D 的溶剂也是 Pyr-d5

**后果**：
- 1H 化学位移用 CDCl₃ 的值（CDCl3 的 TMS 残留峰在 7.26 ppm；Pyr-d5 在 8.71/7.57/7.21 ppm）→ 程序错误把 Pyr 残余峰当真实峰
- HSQC 配对检索的 1H 化学位移与 2D 谱本身的 F2 轴不匹配（一个在 CDCl₃ 标度，一个在 Pyr 标度）→ HSQC 配对失败率升高
- 1H 积分的溶剂基线在 CDCl3 vs Pyr 中差异显著，影响相对积分比例

**修复优先级**：🔴 高 — 这是提升归属表匹配率到 90% 的最关键单点

**修复方案**：
1. 优先取溶剂与 2D 实验一致的 1H 谱（即 exp 31，而非 exp 21）
2. 加 CLI 参数 `--prefer-1h-solvent <SOLVENT>` 让用户指定
3. 兜底：若只有一组 1H，沿用当前行为

---

## 四、按文件维度评估的"可补强"中间产物清单

### 产物 A（🔴 关键）：1H 谱溶剂感知选择

- **输入**：所有 1H 实验的 SOLVENT + 2D 实验的 SOLVENT
- **算法**：选与 2D 实验溶剂一致（多数表决）的 1H 实验；若 2D 不存在则选最大 SW / 最长 TD 的 1H 实验
- **产物**：`meta['h1_source_exp']` + `meta['h1_solvent_aligned']`
- **模块**：`io_bruker.py` 加 `select_primary_1h_spec(specs) -> BrukerSpectrum`

### 产物 B（中）：`peakrng` 增强重叠峰检测

- **输入**：`pdata/1/peakrng`（左右边界 F1/F2）
- **算法**：对 step1 C-H 对，若邻近 C-H 对的 (dC ±0.5 ppm) 范围内 peakrng F1 边界有重叠 → 标记为重叠峰
- **产物**：附加到 `step1_hsqc_ch_pairs.csv` 的 `is_overlap_from_peakrng` 列
- **模块**：扩展 `interpret_2d._step1_resolve_overlapping` 用 peakrng 信息

### 产物 C（中）：HMBC F1 投影 (pdata/99/999) 作为 ¹³C 辅助峰源

- **输入**：`27/pdata/99/1r`（HMBC 列投影 = ¹³C 维投影）
- **算法**：nmrglue 读 1r + proc/procs 重建 ppm 轴 → 用 scipy.signal.find_peaks 在投影上检峰 → 与 22 主谱 ¹³C 峰并集去重（容差 0.1 ppm）
- **产物**：`meta['13c_from_hmbc_projection']` 列表
- **价值**：HMBC 投影的 ¹³C 维有更高 SNR（求和增强），能检出 22 主谱漏掉的弱 ¹³C 峰
- **模块**：新增 `io_bruker_2d.read_hmbc_f1_projection()`

### 产物 D（中）：HSQC 单行切片精细化 (pdata/1000/1001)

- **输入**：`26/pdata/1000/1r` 和 `26/pdata/1001/1r`
- **算法**：通过 `used_from` 的 PROCNO2（行号）→ 找到对应 ¹³C ppm → 取该行的 1H 切片 → find_peaks 检 H 峰 → 校准 J 耦合
- **产物**：每行附 `1h_slice_peaks_ppm`, `1h_slice_j_Hz`
- **价值**：HSQC 行的精确 1H 模式（多重峰 J）→ 增强多重峰归属
- **模块**：新增 `io_bruker_2d.read_hsqc_row_slice()`

### 产物 E（低）：`peakrng` 多重峰 J 模式识别

- **输入**：`pdata/1/peakrng`（峰宽）
- **算法**：对每个峰，若 region F2 width > 单峰 FWHM 阈值 → 标记为多重峰 → 配对 step1b 多重峰分析
- **产物**：step1b CSV 的 `from_peakrng` 列

### 产物 F（低）：`audita.txt` 冗余元数据校验

- **输入**：`exp/audita.txt`
- **算法**：提取 SOLVENT/PULPROG 字段，与 parm.txt 一致性校验；不一致时优先用 audita
- **产物**：`meta['solvent_consistency']` 字段

---

## 五、产物优先级排序（按对提升 39.6% → 90% 匹配率的预期贡献）

| 优先级 | 产物 | 预期贡献 | 风险 |
|---|---|---|---|
| **🔴 1** | A. 1H 谱溶剂感知选择 | +25-35 pp（消除 Pyr/CDCl3 溶剂标度混用错误） | 低（仅选择逻辑改动） |
| **🟡 2** | C. HMBC F1 投影补强 ¹³C | +5-10 pp（检出 22 主谱漏掉的弱 ¹³C 峰） | 中（需 nmrglue） |
| **🟡 3** | D. HSQC 行切片精细化 J | +3-5 pp（增强多重峰归属） | 中（需 nmrglue） |
| **🟢 4** | B. peakrng 重叠峰检测 | +2-3 pp | 低 |
| **🟢 5** | E. peakrng 多重峰 J 识别 | +1-2 pp | 低 |
| **⚪ 6** | F. audita.txt 元数据校验 | 0-1 pp（边界修正） | 低 |

---

## 六、实施顺序

1. **第 1 周**：产物 A（溶剂感知 1H 选谱）—— **最关键**
2. **第 2 周**：产物 C（HMBC F1 投影补强 ¹³C） + 产物 D（HSQC 行切片精细化）
3. **第 3 周**：产物 B + E（peakrng 增强）
4. **第 3 周末**：在 clxh100 上重跑 `--derive` → 行匹配率预期 ≥ 75%
5. **第 4 周**：融合前序的 HMBC 闭环 + DEPT-135 正负号整合 + 化学位移经验规则 → 进一步冲 90%

---

## 七、不在本次范围

- ML 模型（external/NMR2Struct/）、绝对构型、变温 NMR、混合 2D
- 已知化合物复核阶段的 SMILES/mol block（仅在 `compare_program_to_standard` 中作为结构判等依据）
- 3D RDKit 嵌入（不在推导流程）

## 八、产物目录建议

每个产物落盘到 `output/<sample>/` 的独立 CSV/JSON：

```
output/clxh100/
├── 1h_source_selection.csv            # 产物 A：每个 1H 实验的溶剂/SW/TD + 选择依据
├── 13c_with_hmbc_projection.csv       # 产物 C：22 主谱 ∪ HMBC 投影峰
├── hsqc_row_slices/                   # 产物 D：HSQC 单行切片等高线图
│   ├── row_51_*.png
│   └── row_941_*.png
├── step1_hsqc_ch_pairs.csv            # 已有 + 新增 is_overlap_from_peakrng 列
└── step1b_multiplet_analysis.csv      # 已有 + 新增 from_peakrng 列
```# 归属表匹配率从 39.6% 提升到 ≥90% 的中间产物补充计划

## 〇、核心设计约束（优先级最高）

**推导流程的输入边界（强制）：**

1. 推导流程**只能**使用仪器测出的原始数据：1D `peaklist.xml`（1H/13C/DEPT）、`integrals.txt`、`parm.txt`、2D 谱数据（HSQC/COSY/HMBC/NOESY/ROESY 矩阵）
2. 推导流程**不可以**使用任何先验结构信息：禁止读取 `smiles.txt`、`*.cdxml`、`*.mol`、`真实结构/<compound>/` 下的化学家归属表
3. 推导产物**必须对所有化合物通用**，禁止为某一化合物拟合阈值或规则
4. 已知化合物复核（`standard_assignment.compare_program_to_standard`）使用 SMILES 仅做化合物结构判等（**不参与行匹配**，遵循 `AGENTS.md` / `README.md` 的现有设计原则），行匹配仍按 (δH, δC, 碳类型) 三元组

---

## 一、现状诊断（基于 `output/clxh100/derive_归属表/comparison_report.txt`）

| 维度 | 实测 | 目标 | 差距 |
|---|---|---|---|
| 行匹配率 | 39.6% (19/48) | ≥90% | -50.4 pp |
| 碳类型一致率 | 5/19 (26%) | ≥90% | -64 pp |
| 输出行数 | 48 行（程序）vs 41 行（化学家） | 一致 | 多 7 行 |
| 季碳误识 | 7 个 Cq（150.1/149.9/149.7/135.8/135.6/135.4/30.0）—— 但 clxh100 是甾体皂苷无苯环 | 必须修正 | 7 处误判 |
| 标号差异 | 程序 `C-1, C-2, …, C40` vs 化学家 `1α, 1β, 2α, 2β, 26a, 26b, Glc-1′, …` | 仅命名约定差异，比较时按 (δH, δC, 类型) 三元组匹配即可 | 命名约定 |

## 二、根因分析（4 类系统性问题，剔除依赖 SMILES 的"标号映射"）

| # | 根因 | 表现 | 解决方向 |
|---|---|---|---|
| 1 | **DEPT 碳级数推断 CH/CH₂ 互换多** | 14/19 匹配行类型不一致（如程序 CH₂ 对标准 CH） | 集成 DEPT-135 正负号（CH/CH₃ 正, CH₂ 负）到 `_step1_hsqc_to_ch_pairs`，替代/修正 1H 积分推断 |
| 2 | **季碳误识（最大问题）** | 程序在 150/149/135 ppm 区识别出 7 个 Cq，但 clxh100 是甾体皂苷骨架无苯环 → 这些是 13C 谱噪声峰或 DEPT 推断错误 | 用 DEPT-135 的 CH₂ 负峰反证、HMBC 闭环交叉验证、化学位移经验（羰基 170-220, 季碳 30-50）|
| 3 | **缺失 HMBC 反推归属** | HMBC 只是罗列，未用来确认"哪个 H 归属哪个片段" | 用 HMBC 网络闭环：每条 HMBC 连接方向必须与化学位移区段（羰基/烯烃/季碳/糖端基）一致 |
| 4 | **缺失化学位移经验规律校验** | 化学家靠经验规律（苷化位移、糖端基特征峰、羰基 α-CH 耦合）做最终核对 | 维护基于原始数据 + 文献经验的可参数化规则库（不依赖任何结构先验） |

> **关于"标号系统错位"的说明**：这不是归属错误，是命名约定的差异。程序按 HSQC 出现顺序编号（`C-1, C-2, …`）是推导流程的稳定命名；化学家按结构位置编号（`1α, 1β, 2α, 2β, Glc-1′`）是按已知结构人工编号。比较时 `compare_program_to_standard._joint_match` 已按 (δH ±0.1, δC ±0.5, 碳类型) 三元组做 Hungarian 一对一匹配，**行序无关**。故无需"标号映射"，也避免了引入 SMILES。

## 三、补充的 7 类中间产物（全部仅依赖仪器原始数据）

### 产物 1：`step6_hmbc_closure_matrix.csv` — HMBC 闭环矩阵

- 内容：以化学位移区段为列、HSQC H 索引为行的二值矩阵，每个格子表示"该 H 通过 HMBC 是否能连接到该 C 区段"
- 算法：对 `interpretation_2d.hmbc_connections` 按 C 区间分类（羰基 170-220, 芳香/烯烃 100-160, 季碳脂肪 30-50, 糖端基 95-110），检查每个 HMBC 连接的方向是否合理
- 价值：闭环率 = (有 ≥1 个 HMBC 连接的 H 数) / (总 H 数)。闭环率 < 70% 提示归属表有缺失
- 模块：`correlation_2d.build_hmbc_closure_matrix()` + `interpret_2d.write_step_files()` 追加

### 产物 2：`step8_chemical_shift_segment.csv` — 化学位移区段分类

- 内容：每个 HSQC C-H 对按 δC 落入区段（羰基/烯烃-芳香/糖-脂肪/甲基/季碳）归类，并给出区段归属建议
- 算法：扩展 `shift_ranges.SHIFT_RANGES`，加入糖端基、苷化位移、季碳去耦区段
- 价值：让化学家直接看到"该 H 的 C 在哪个区段 → 应是哪个片段"
- 模块：`shift_ranges.assign_shift_segment_v2()` + `interpret_2d.write_step_files()`

### 产物 3：`step9_fragment_recognition.csv` — 片段识别（不依赖 SMILES）

- 内容：从 `interpretation_2d.ch_pairs` + `interpretation_2d.hmbc_connections` + `hmbc_connections.bond_type` 推断"哪些 C/H 组成同一结构片段"（糖环、苷元环、酯基侧链、羰基邻域等）
- **关键约束：只用谱图数据 + 化学位移指纹，绝不调用 RDKit/Chem.MolFromSmiles**
- 算法（**纯谱图指纹法**）：
  1. **糖端基识别**：检测 (δC ∈ [95, 110], δH ∈ [4.2, 5.5], J ≈ 7-8 Hz) 的 HSQC 峰 → 标记为 `glycan_C1_candidate`
  2. **糖环聚类**：从糖端基出发，沿 HMBC 2J/3J 找到 4-5 个相邻 (δC ∈ [68, 85]) 的 C → 标记为 `glycan_ring`（构成吡喃糖环的 5 个 C）
  3. **糖-苷元连接位点**：糖端基 H 通过 HMBC 3J 连接到的苷元 C（一般 δC 偏移 +6~+10 ppm）→ 标记为 `glycan_link_C`
  4. **季碳邻域识别**：羰基 C（δC ∈ [165, 220]）通过 HMBC 接收 H 的连接模式 → 标记为 `carbonyl_*`（区分酮/醛/酯/羧酸/酰胺）
  5. **脂肪链连续自旋系统**：从 COSY 链长度 ≥ 4 的自旋系统 → 标记为 `alkyl_chain`
  6. **孤立甲基**：δH ∈ [0.8, 1.4] 且 δC ∈ [10, 25] 且无 COSY 邻位 → 标记为 `isolated_Me`
- **删除内容**：原有 `assignment_table_v2._detect_glycan_substructure` 中通过 SMILES 子结构匹配的部分（`Chem.MolFromSmarts`）
- 价值：直接对应化学家归属表的 `Glc-1′, Glc-2′, …, Glc-6a′, Glc-6b′` 命名（按糖环拓扑自动生成）
- 模块：新增 `fragment_recognizer.py`

### 产物 4：`step10_roesy_stereo_inference.csv` — ROESY 立体归属（不依赖 SMILES）

- 内容：从 NOESY/ROESY 接触的**相对距离**和化学位移规律推断每个 CH₂ 上两个 H 的相对面（α/β），仅用 NOE 距离大小和 1H 化学位移区段
- **关键约束：不调用 RDKit，不生成 3D conformer，不依赖候选结构**
- 算法（**纯相对构型法**）：
  1. **同碳二氢识别**：每个 CH₂ 有两个 H（H_a, H_b），δH 相差通常 > 0.05 ppm（化学不等价）
  2. **空间近邻统计**：对每个 H，统计其 NOE 接触的对侧 H 数（"NOE 度"）
  3. **相对面推断规则**（**可参数化通用规则**，不依赖具体结构）：
     - 若两个 H 的 NOE 度差异 > 30% → 度较大的为 β（环上凸面，与更多相邻 H 接近），较小的为 α（环下凹面）
     - 甾体/三萜骨架特征：Me-18（δC ∈ [12, 16]）、Me-19（δC ∈ [15, 22]）默认为 β 朝向（甾体化学约定，不依赖结构）
     - 糖端基 H（δH ≈ 4.3-4.9, J ≈ 7-8 Hz）按 J 判断 β-D（7-8 Hz）或 α-L（3-4 Hz）
  4. **CH₂ 拆分规则**：CHPair 中 h_integral ≥ 1.5 的（CH₂）才拆分 α/β；CH/CH₃ 不拆
- **删除内容**：`correlation_2d.verify_candidate_3d_noe` 中基于 RDKit 3D 的部分改为纯 NOE 距离统计；3D 部分保留为"已知化合物复核"的独立函数
- 价值：使程序输出 `1α/1β`, `15α/15β`, `23a/23b`, `26a/26b` 命名（仅依赖 NOE 距离 + 化学位移，不依赖结构）
- 模块：新增 `noe_stereo_inference.py` + 改造 `correlation_2d.py` 中的 3D 部分移到 `verify_known_structure_3d_noe()`

### 产物 5：`step11_aglycone_shift_rules.csv` — 苷化位移与化学位移规律校验

- 内容：对每个归属行打"经验规律分"：苷化位移（+6~+10 ppm）/ 糖端基（95-106 ppm）/ 26-O-Glc（+2~+4 ppm）/ 羰基（170-180 ppm）/ α-羰基 α-CH（2.0-2.5 ppm）等
- 算法：维护 `shift_rules.yaml` 知识库（糖类 / 甾体 / 三萜 / 黄酮 等 10+ 骨架），每条规则形如 `if dC in [a, b] and dH in [c, d] then expected_type in [...]`
- 价值：发现"违反经验规律的归属"作为告警，例如：把 CH₂ 标到 170 ppm（必是 Cq）→ 标错
- 模块：新增 `shift_rules_validator.py`

### 产物 6：`step13_hmbc_direction_check.csv` — HMBC 方向语义校验

- 内容：对每条 HMBC 连接检查 (dH 来源区段, dC 目标区段) 是否符合化学位移规律
- 算法：
  - CH₃ + 2J → 邻位 C（脂肪 20-50 或羰基 170-220）
  - CH₃ + 3J → 跨 2 键到 C（脂肪 / 烯烃 100-160）
  - 糖端基 H (~4.8 ppm) + 2J → 苷元连接位点 C（苷化位移 +6~+10 ppm）
  - CH₂ + 3J → 邻位 CH 或季碳
- 价值：HMBC 反向校验 H 归属的正确性（CH₃-HMBC 到 Cq30.0 → 不合理 → 应是 CH/CH₂ 误识）
- 模块：新增 `hmbc_direction_validator.py`

### 产物 7：`step14_dept_sign_integration.csv` — DEPT-135 正负号整合

- 内容：从 exp 23 (DEPT-135) 解析每个 13C 峰的正负号，正 → CH/CH₃，负 → CH₂，无峰 → Cq
- 算法：扩展 `io_bruker.read_bruker_folder` 读 DEPT-135 峰列表（已部分实现 `bridge_nmr2struct.assign_carbon_multiplicity`），把正负号填入 `ch_pairs.multiplicity`
- 价值：解决 14/19 碳类型不一致问题（CH/CH₂ 互换）
- 模块：扩展 `bridge_nmr2struct.py` + `interpret_2d._step1_hsqc_to_ch_pairs`

## 四、新增/修改模块清单

| 文件路径 | 动作 | 说明 |
|---|---|---|
| `nmr_elucidator/fragment_recognizer.py` | 新增 | 纯谱图指纹法片段识别（产物 3），不调用 RDKit |
| `nmr_elucidator/noe_stereo_inference.py` | 新增 | NOE 相对距离 + 经验规律推断 α/β（产物 4），不调用 RDKit |
| `nmr_elucidator/shift_rules_validator.py` | 新增 | 化学位移经验规律校验（产物 5） |
| `nmr_elucidator/hmbc_direction_validator.py` | 新增 | HMBC 方向语义校验（产物 6） |
| `nmr_elucidator/correlation_2d.py` | 修改 | `build_hmbc_closure_matrix()` + 把 `verify_candidate_3d_noe` 拆分为纯 NOE 函数（推导用）和 3D 函数（已知化合物复核用） |
| `nmr_elucidator/interpret_2d.py` | 修改 | `_step1_hsqc_to_ch_pairs` 用 DEPT-135 正负号 + 追加 step6/step8~14 输出 |
| `nmr_elucidator/bridge_nmr2struct.py` | 修改 | `assign_carbon_multiplicity` 输出 DEPT-135 正负号字段 |
| `nmr_elucidator/shift_ranges.py` | 修改 | `SHIFT_RANGES` 扩展糖端基 / 苷化位移 / 季碳区段 |
| `nmr_elucidator/assignment_table_v2.py` | 修改 | 删除 `_detect_glycan_substructure` 的 SMILES 子结构分支，改为纯化学位移指纹法 |
| `data/shift_rules.yaml` | 新增 | 化学位移规律知识库（糖类 / 甾体 / 三萜 等），参数化定义 |

> **不在本次范围（明确剔除的 SMILES 依赖项）**：
> - `label_mapper.py`（原计划产物 6）—— 用户已明确拒绝
> - `atom_labeler.py` 增加 `assign_chemist_labels_from_mol_block()` —— 同上
> - `fragment_assembler.assemble_candidates()` 在 `--derive` 流程中禁用（仅在 `--retrieve` + 无候选库时启用，不参与归属推导）

## 五、`__main__.py` 集成流程修改

在 `_build_derivation` (Part 3 Phase 5.x) 末尾追加步骤（**严格仅依赖谱图数据**）：

```
[2.8] 2D 解读 (HSQC/COSY/HMBC/NOESY)
   ↓
[2.9] DEPT-135 正负号整合到 ch_pairs (产物 7)  ← 替代/修正 1H 积分推断
   ↓
[2.10] 化学位移区段分类 (产物 2)
   ↓
[2.11] HMBC 闭环矩阵 + 方向语义校验 (产物 1 + 6)
   ↓
[2.12] 片段识别 (糖端基/糖环/羰基邻域/孤立甲基) — 纯谱图指纹法 (产物 3)
   ↓
[2.13] NOE 相对距离推断 α/β 立体 (产物 4)  ← 不生成 3D conformer
   ↓
[2.14] 苷化位移与化学位移规律校验 (产物 5)
   ↓
[2.15] assignment_table_v2.build_assignment_table_v2() 用以上数据填充
   ↓
[2.16] 已知化合物复核 (separate, 可选)  ← 此时才允许读 SMILES/mol block, 仅用于结构判等, 不参与行匹配
```

## 六、验证与目标

每个新模块必须先单独验证（产物 1~7 各自的正确率指标，且必须可在 clxh100 以外的化合物上复现）：

| 产物 | 验证指标 | 目标 |
|---|---|---|
| 1. HMBC 闭环矩阵 | 闭环率 | ≥ 80% |
| 3. 片段识别（纯谱图指纹） | 糖端基识别召回率 | ≥ 90% |
| 4. NOE α/β 立体（纯 NOE 距离） | 与 NOE 距离一致性 | ≥ 85% |
| 5. 化学位移规律 | 经验规律命中率 | ≥ 85% |
| 6. HMBC 方向校验 | 误判率 | < 5% |
| 7. DEPT-135 整合 | CH/CH₂ 修正率 | ≥ 80% |

**通用性验证**：在 clxh3（无 2D 数据，回退 1D 流程）和未来 3 个未知化合物上跑 `--derive`，各产物的算法必须**不报错**且输出非空。

集成验证：在 clxh100 上跑 `--derive` → 跑 `compare_program_to_standard`：
- **行匹配率**：39.6% → ≥90%
- **碳类型一致率**：26% → ≥90%

## 七、实施顺序（按依赖关系）

1. **第 1 周**：产物 7（DEPT-135）—— 解决最直接的"碳类型不一致"
2. **第 2 周**：产物 2（区段分类）+ 产物 6（HMBC 方向）—— 闭环验证基础
3. **第 3 周**：产物 5（化学位移规律）+ 产物 1（HMBC 闭环矩阵）—— 经验校验
4. **第 4 周**：产物 3（片段识别 — 纯谱图指纹）+ 产物 4（NOE α/β — 纯 NOE 距离）—— 命名与立体
5. **第 4 周末**：集成验证 clxh100 + clxh3 + 至少 1 个未知化合物，达到 90%+ 目标

## 八、风险与缓解

| 风险 | 缓解 |
|---|---|
| DEPT-135 peaklist 缺失（部分化合物只有 DEPT-90） | 退化为 1H 积分推断 + 化学位移规律校验 |
| NOESY/ROESY 数据缺失 | 退化为仅做 CH₂ 化学位移差分（双氢化学不等价）但 α/β 标记留空 |
| 化学位移规律知识库不完整 | 初始版本覆盖糖类/甾体/三萜 3 类，逐步扩展；规则用参数化 YAML，便于追加 |
| 糖端基检测有噪声峰（不依赖结构的代价） | 用 HMBC 网络闭环二次校验：糖端基 H 至少要有 2 条 HMBC 连接，否则降级为"疑似糖端基" |
| 已知化合物复核阶段才读 SMILES 的边界模糊 | 在 `__main__.py` 显式分两阶段：`_build_derivation()`（无 SMILES）与 `_verify_known_structure()`（允许 SMILES）|

## 九、不在本次范围

- ML 模型（`external/NMR2Struct/` 端到端 GPU 推理）
- 绝对构型（Mosher 酯法 / CD / X-ray）
- 变温 NMR / rotamers
- 混合 2D（HSQC-TOCSY / HSQC-NOESY / LR-HSQMBC）
- 完整 CASE 自动解谱（`external/lucy-ng/` 或 `external/LSDos-3.5.3/`）
- **化学家标号到程序标号的映射**（用户明确拒绝，且 `compare_program_to_standard` 已用三元组匹配绕过该问题）
- **3D RDKit 嵌入参与推导流程**（仅在已知化合物复核中可用，且只用于结构验证，不参与归属）# 归属表匹配率从 39.6% 提升到 ≥90% 的中间产物补充计划

## 一、现状诊断（基于 `output/clxh100/derive_归属表/comparison_report.txt`）

| 维度 | 实测 | 目标 | 差距 |
|---|---|---|---|
| 行匹配率 | 39.6% (19/48) | ≥90% | -50.4 pp |
| 碳类型一致率 | 5/19 (26%) | ≥90% | -64 pp |
| 输出行数 | 48 行（程序）vs 41 行（化学家） | 一致 | 多 7 行 |
| 季碳误识 | 7 个 Cq（150.1/149.9/149.7/135.8/135.6/135.4/30.0）—— 但 clxh100 是甾体皂苷无苯环 | 必须修正 | 7 处误判 |

## 二、根因分析（5 类系统性问题）

| # | 根因 | 表现 | 解决方向 |
|---|---|---|---|
| 1 | **标号系统错位** | 程序按 HSQC 出现顺序编号 `C-1, C-2, …, C40`，化学家按结构位置编号 `1α, 1β, 2α, 2β, 26a, 26b, Glc-1′, …` | 用 RDKit canonical SMILES + 化学家标号字典（`atom_labeler.py` 已具备基础）把 `C-1` 映射到结构上的实际原子编号 |
| 2 | **季碳误识（最大问题）** | 程序在 150/149/135 ppm 区识别出 7 个 Cq，但 clxh100 是甾体皂苷骨架，无苯环 → 这些是 13C 谱噪声峰或 DEPT 推断错误 | 用 DEPT-135 的 CH₂ 负峰反证、HMBC 闭环交叉验证、化学位移经验（羰基 170-220, 季碳 30-50） |
| 3 | **DEPT 碳级数推断 CH/CH₂ 互换多** | 14/19 匹配行类型不一致（如程序 CH₂ 对标准 CH） | 集成 DEPT-135 正负号判断（CH/CH₃ 正, CH₂ 负）到 `_step1_hsqc_to_ch_pairs`，而非仅依赖 1H 积分 |
| 4 | **缺失 HMBC 反推归属** | HMBC 只是罗列，未用来确认"哪个 H 归属哪个片段" | 用 HMBC 网络闭环：每条 HMBC 连接方向必须与化学位移区段（羰基/烯烃/季碳/糖端基）一致 |
| 5 | **缺失 ROESY 立体归属** | 真实归属表 α/β 立体（`1α/1β`, `15α/15β`）程序无法输出 | 用 RDKit 3D conformer + NOE 距离推断 H 的相对面，输出 α/β |

## 三、补充的 8 类中间产物

### 产物 1：`step6_hmbc_closure_matrix.csv` — HMBC 闭环矩阵

- 内容：以化学位移区段为列、HSQC H 索引为行的二值矩阵，每个格子表示"该 H 通过 HMBC 是否能连接到该 C 区段"
- 算法：对 `interpretation_2d.hmbc_connections` 按 C 区间分类（羰基 170-220, 芳香/烯烃 100-160, 季碳脂肪 30-50, 糖端基 95-110），检查每个 HMBC 连接的方向是否合理
- 价值：闭环率 = (有 ≥1 个 HMBC 连接的 H 数) / (总 H 数)。闭环率 < 70% 提示归属表有缺失
- 模块：`correlation_2d.build_hmbc_closure_matrix()` + `interpret_2d.write_step_files()` 追加

### 产物 2：`step8_chemical_shift_segment.csv` — 化学位移区段分类

- 内容：每个 HSQC C-H 对按 δC 落入区段（羰基/烯烃-芳香/糖-脂肪/甲基/季碳）归类，并给出区段归属建议
- 算法：扩展 `shift_ranges.SHIFT_RANGES`，加入糖端基、苷化位移、季碳去耦区段
- 价值：让化学家直接看到"该 H 的 C 在哪个区段 → 应是哪个片段"
- 模块：`shift_ranges.assign_shift_segment_v2()` + `interpret_2d.write_step_files()`

### 产物 3：`step9_fragment_recognition.csv` — 片段识别

- 内容：从 `interpretation_2d.ch_pairs` + `interpretation_2d.hmbc_connections` 推断"哪些 C/H 组成同一结构片段"（糖基、苷元环 A/B/C/D、酯基侧链等）
- 算法：
  - 糖基识别（已部分实现 `assignment_table_v2._detect_glycan_substructure` + 端基 C1 ~104-106 ppm 触发）→ 标注 `Glc-1′, Glc-2′, …, Glc-6a′, Glc-6b′`
  - 甾体骨架识别：A 环（CH-3 苷化位移 ~76-80 ppm）/ B 环 / C 环 / D 环 / 侧链
  - 季碳连接的 HMBC 邻域识别（`Cq22 173.8` 附近 = 酯羰基）
- 价值：直接对应化学家归属表的 `26-O-Glc-1′, 2′, 3′, 4′, 5′, 6a′, 6b′` 命名
- 模块：新增 `fragment_recognizer.py`

### 产物 4：`step10_roesy_stereo_inference.csv` — ROESY 立体归属

- 内容：从 NOESY/ROESY 接触 + RDKit 3D conformer 推断每个 CH₂ 上两个 H 的相对面（α/β）
- 算法：
  1. 对候选 SMILES 做 ETKDGv3 + MMFF 优化（`correlation_2d.verify_candidate_3d_noe` 已用）
  2. 对每个 CH₂，提取两个 H 的 Z 坐标 → 较小 Z = α，较大 Z = β
  3. 用 NOE 距离 < 4 Å 验证 α/β 标注是否与 NOE 一致
- 价值：使程序输出 `1α/1β`, `15α/15β`, `23a/23b`, `26a/26b` 与化学家一致
- 模块：新增 `noe_stereo_inference.py`

### 产物 5：`step11_aglycone_shift_rules.csv` — 苷化位移与化学位移规律校验

- 内容：对每个归属行打"经验规律分"：苷化位移（+6~+10 ppm）/ 糖端基（95-106 ppm）/ 26-O-Glc（+2~+4 ppm）/ 羰基（170-180 ppm）/ α-羰基 α-CH（2.0-2.5 ppm）等
- 算法：维护 `shift_rules.yaml` 知识库（糖类 / 甾体 / 三萜 / 黄酮 等 10+ 骨架）
- 价值：发现"违反经验规律的归属"作为告警，例如：把 CH₂ 标到 170 ppm（必是 Cq）→ 标错
- 模块：新增 `shift_rules_validator.py`

### 产物 6：`step12_chemist_label_mapping.csv` — 化学家标号映射

- 内容：把程序 `C-1, C-2, …, C40` 通过 RDKit + 化学家标号字典（mol block）映射到 `1α, 1β, 2α, 2β, …, 26-O-Glc-1′, …`
- 算法：
  1. 从 `真实结构/<compound>/` 读 `clxh100数据归属表.xlsx` + `smiles.txt`
  2. 用 RDKit 解析 SMILES + 用 mol block 推断化学家标号（α/β 取 RDKit 3D Z 坐标）
  3. 对每个程序 CHPair 找到化学家归属表中最近 (δH, δC) 的标准行 → 把化学家 No. 填入
- 价值：解决 48 行 vs 41 行编号错位的根本问题
- 模块：扩展 `atom_labeler.py` + 新增 `label_mapper.py`

### 产物 7：`step13_hmbc_direction_check.csv` — HMBC 方向语义校验

- 内容：对每条 HMBC 连接检查 (dH 来源区段, dC 目标区段) 是否符合化学位移规律
- 算法：
  - CH₃ + 2J → 邻位 C（脂肪 20-50 或羰基 170-220）
  - CH₃ + 3J → 跨 2 键到 C（脂肪 / 烯烃 100-160）
  - 糖端基 H (~4.8 ppm) + 2J → 苷元连接位点 C（苷化位移 +6~+10 ppm）
  - CH₂ + 3J → 邻位 CH 或季碳
- 价值：HMBC 反向校验 H 归属的正确性（CH₃-HMBC 到 Cq30.0 → 不合理 → 应是 CH/CH₂ 误识）
- 模块：新增 `hmbc_direction_validator.py`

### 产物 8：`step14_dept_sign_integration.csv` — DEPT-135 正负号整合

- 内容：从 exp 23 (DEPT-135) 解析每个 13C 峰的正负号，正 → CH/CH₃，负 → CH₂，无峰 → Cq
- 算法：扩展 `io_bruker.read_bruker_folder` 读 DEPT-135 峰列表（已部分实现 `bridge_nmr2struct.assign_carbon_multiplicity`），把正负号填入 `ch_pairs.multiplicity`
- 价值：解决 14/19 碳类型不一致问题（CH/CH₂ 互换）
- 模块：扩展 `bridge_nmr2struct.py` + `interpret_2d._step1_hsqc_to_ch_pairs`

## 四、新增/修改模块清单

| 文件路径 | 动作 | 说明 |
|---|---|---|
| `nmr_elucidator/fragment_recognizer.py` | 新增 | 糖基/甾体骨架识别（产物 3） |
| `nmr_elucidator/noe_stereo_inference.py` | 新增 | RDKit 3D + NOE → α/β 立体（产物 4） |
| `nmr_elucidator/shift_rules_validator.py` | 新增 | 化学位移经验规律校验（产物 5） |
| `nmr_elucidator/hmbc_direction_validator.py` | 新增 | HMBC 方向语义校验（产物 7） |
| `nmr_elucidator/label_mapper.py` | 新增 | 程序编号 → 化学家编号映射（产物 6） |
| `nmr_elucidator/correlation_2d.py` | 修改 | `build_hmbc_closure_matrix()` + `extract_hmbc_pairs()` 加方向分类 |
| `nmr_elucidator/interpret_2d.py` | 修改 | `_step1_hsqc_to_ch_pairs` 用 DEPT-135 正负号 + 追加 step8~14 输出 |
| `nmr_elucidator/atom_labeler.py` | 修改 | 增加 `assign_chemist_labels_from_mol_block()` |
| `nmr_elucidator/bridge_nmr2struct.py` | 修改 | `assign_carbon_multiplicity` 输出 DEPT-135 正负号字段 |
| `nmr_elucidator/shift_ranges.py` | 修改 | `SHIFT_RANGES` 扩展糖端基 / 苷化位移 / 季碳区段 |
| `data/shift_rules.yaml` | 新增 | 化学位移规律知识库（糖类 / 甾体 / 三萜 等） |

## 五、`__main__.py` 集成流程修改

在 `_build_derivation` (Part 3 Phase 5.x) 末尾追加步骤：

```
[2.8] 2D 解读 (HSQC/COSY/HMBC/NOESY)
   ↓
[2.9] DEPT-135 正负号整合到 ch_pairs (产物 8)
   ↓
[2.10] 化学位移区段分类 (产物 2)
   ↓
[2.11] HMBC 闭环矩阵 + 方向语义校验 (产物 1 + 7)
   ↓
[2.12] ROESY 立体归属推断 (产物 4)
   ↓
[2.13] 片段识别 (糖基/甾体骨架/季碳邻域) (产物 3)
   ↓
[2.14] 苷化位移与化学位移规律校验 (产物 5)
   ↓
[2.15] 化学家标号映射 (产物 6) ← 需要时启用（已知化合物复核时）
   ↓
[2.16] assignment_table_v2.build_assignment_table_v2() 用以上数据填充
```

## 六、验证与目标

每个新模块必须先单独验证（产物 1~8 各自的正确率指标）：

| 产物 | 验证指标 | 目标 |
|---|---|---|
| 1. HMBC 闭环矩阵 | 闭环率 | ≥ 80% |
| 4. ROESY 立体 | 与 NOE 一致率 | ≥ 90% |
| 5. 化学位移规律 | 经验规律命中率 | ≥ 85% |
| 6. 化学家标号映射 | 化学家 No. 命中数 | ≥ 36/41 |
| 7. HMBC 方向校验 | 误判率 | < 5% |
| 8. DEPT-135 整合 | CH/CH₂ 修正率 | ≥ 80% |

集成验证：在 clxh100 上跑 `--derive` → 跑 `compare_program_to_standard`：
- **行匹配率**：39.6% → ≥90%
- **碳类型一致率**：26% → ≥90%

## 七、实施顺序（按依赖关系）

1. **第 1 周**：产物 8（DEPT-135）+ 产物 6（化学家标号映射）—— 解决最直接的"季碳误识"和"标号错位"
2. **第 2 周**：产物 2（区段分类）+ 产物 7（HMBC 方向）—— 闭环验证
3. **第 3 周**：产物 3（片段识别）+ 产物 5（化学位移规律）—— 命名规范化
4. **第 4 周**：产物 1（HMBC 闭环矩阵）+ 产物 4（ROESY 立体）—— 完善细节
5. **第 4 周末**：集成验证 clxh100，达到 90%+ 目标

## 八、风险与缓解

| 风险 | 缓解 |
|---|---|
| DEPT-135 peaklist 缺失（部分化合物只有 DEPT-90） | 退化为 1H 积分推断 + 化学位移规律 |
| RDKit 3D 嵌入失败（复杂骨架） | 退化为 2D wedge bond 推断 → fallback 索引顺序 |
| 化学家标号字典未提供（无 `smiles.txt`） | 退化为程序自有标号，仅在已知化合物复核时启用映射 |
| 化学位移规律知识库不完整 | 初始版本覆盖糖类/甾体/三萜 3 类，逐步扩展 |

## 九、不在本次范围

- ML 模型（`external/NMR2Struct/` 端到端 GPU 推理）
- 绝对构型（Mosher 酯法 / CD / X-ray）
- 变温 NMR / rotamers
- 混合 2D（HSQC-TOCSY / HSQC-NOESY / LR-HSQMBC）
- 完整 CASE 自动解谱（`external/lucy-ng/` 或 `external/LSDos-3.5.3/`）