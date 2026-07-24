# NMR 反推结构工具 (`nmr_elucidator`)

基于 **nmrshiftdb2 (58K) + lucy-ng (928K) 双数据源** + **nmrglue Bruker 解析** + **RDKit 化学信息学** 的 NMR 反推结构工具。
读取 1H/13C 谱图 → 在双数据源检索 (优先 lucy-ng 928K 化合物 + lucy-ng-derep.db, 回退 nmrshiftdb2 SDF) → 输出 Top-N 候选 + 证据链 + 未匹配峰标注。
可选附加：**fid 交叉验证**（独立算法路径）、**1r 真实谱可视化**、**nmrsim 反向校验**、**NMR2Struct 端到端 ML**、**HOSE 13C 预测验证**（本地表 + lucy-ng 928K 双数据库）、**PubChem 扩展检索**、**片段组装候选生成**、**2D NMR 关联验证 (HSQC/HMBC/COSY)**、**统一日志记录**。

**核心目标**：从"输出候选列表"升级为"**推断真实结构并报告置信度**"——通过 HOSE code 预测验证 + 多源证据融合打分，输出 determined / probable / ambiguous / no_match 四级置信度结论。

## 工作流

```
Bruker 原始谱 (NMR/<sample>/)
   ├─ peaklist.xml / parm.txt / integrals.txt
   ├─ pdata/1/1r           (--with-real-spectrum 启用)
   └─ fid                  (--validate-with-nmrpy 启用)
        ↓
[1] read_bruker_folder → BrukerSpectrum[]
        ↓
[2] ProblemPeaks (1H + 13C, 20% 强度阈值去噪) + Constraints
[2.5] DEPT 推断 (CH/CH2/CH3/Cq, 用 exp 23+24)
        ↓
[2.7] **已知化合物检测**: 如 E:\NMR2Structure\真实结构/<sample>/ 下存在 *数据归属*.xlsx
        → 加载标准归属表 + 已知 SMILES (仅复核用, 不参与推导)
[2.8] 2D NMR 逐步解读 (化学家工作流 HSQC→COSY→HMBC→季碳)
        ↓
[3] stream_nmrshiftdb2_sdf → retrieve_candidates (窗口匹配 + DEPT 软约束)
        ↓
[4] score_candidates (coverage + RMS penalty) + build_evidence
        ↓
[5] (可选) nmrsim 反向校验 cosine
        ↓
[5.5] (可选) 2D NMR 关联验证 (HSQC/HMBC/COSY, --use-2d 启用)
        ↓
[6] (可选) HOSE code 13C 预测验证 (**双数据源**: 优先 lucy-ng 928K SQLite，回退 nmrshiftdb2 本地表, MAE ~3.5 ppm / ~5 ppm)
        ↓
[3c] (自动) 本地库无命中时, 自动扩大到 PubChem 联网检索
        ↓
[7] (可选) PubChem 强制检索 (--use-pubchem, 即使本地有命中也检索)
        ↓
[8] (可选) 片段组装候选生成 (库外新骨架, DEPT + 化学位移范围驱动)
        ↓
[9] (可选) 结构推断置信度评估 (assess_determination: determined/probable/ambiguous/no_match)
        ↓
[9.8] **已知化合物复核**: 如存在标准归属表, 自动对 derive 输出做行级联合匹配
        (13C ±0.5 ppm / 1H ±0.1 ppm, 兼容 CH2/CH3 占位符).
        → 写出 assignment_table_joint.xlsx + comparison_report.txt
        (SMILES 仅用于化合物结构判等, 不参与行匹配)
        ↓
output/<sample>/<sample>.{json,md,html}  (HTML 内嵌所有图片 base64)
            <sample>_unmatched.csv
            step1_hsqc_ch_pairs.csv      (2D 解读中间结果, 化学家复核用)
            step1_hsqc_ch_pairs/overlap/ (HSQC 重叠峰等高线图, Mnova 风格)
            step1_hsqc_ch_pairs/normal/  (HSQC 非重叠峰等高线图)
            step1b_overlap_analysis.csv  (5.9 节重叠峰分解详情)
            step2_cosy_spin_systems.csv
            step3_hmbc_connections.csv
            step4_quaternary_carbons.csv (含 DEPT 证据来源)
            完整NMR归属表(论文Table1)/   (--use-2d 启用，生成 Top-N 个 CSV, 论文 Table 1 格式)
            <sample>_nmrpy_validation.{json,md,png}  (--validate-with-nmrpy)
            derive_归属表/assignment_table.xlsx       (derive 主输出, 第 1 行表头, 后续数据行, 贴近化学家 Excel)
            derive_归属表/meta.txt                    (元信息: Solvent / Formula / Field_MHz / HRESIMS 等, 独立存放)
            derive_归属表/assignment_table_joint.xlsx (已知化合物复核时生成, 双 sheet 拼接: '按标准表为准' + '按程序表为准')
            derive_归属表/comparison_report.txt        (已知化合物复核时生成, 行级联合匹配明细 + 结论)
```

详见 `从NMR谱图反推分子结构代码库.md`（调研 + 升级路径），本文档只讲**如何使用**。

## 安装

**一键安装** (推荐):

```powershell
# 在项目根目录执行, 自动安装所有依赖
powershell -ExecutionPolicy Bypass -File E:\NMR2Structure\util\install_deps.ps1
```

**手动安装**:

```bash
# 核心依赖
pip install nmrglue nmrsim rdkit pubchempy numpy scipy matplotlib pandas openpyxl -i https://pypi.tuna.tsinghua.edu.cn/simple
# 项目增强
pip install sympy seaborn -i https://pypi.tuna.tsinghua.edu.cn/simple
# 可选 (FID 交叉验证)
pip install nmrpy -i https://pypi.tuna.tsinghua.edu.cn/simple
# 可选 (lucy-ng HOSE 集成)
pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps
```

依赖分层:
- **必需**: rdkit / pubchempy / numpy / scipy / matplotlib / pandas / nmrglue / openpyxl
- **增强**: sympy (公式校验) / seaborn (统计可视化) / nmrsim (反向校验)
- **可选**: nmrpy (FID 交叉验证) / HOSE_code_generator (lucy-ng 集成)

## 数据准备

下载 **nmrshiftdb2** 数据文件（仅 1 个，约 159 MB）：

| 文件 | 大小 | 链接 | 用途 |
|---|---|---|---|
| `nmrshiftdb2withsignals.sd` | 158.6 MB | [SourceForge](https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd/download) | 结构 + 13C/1H 化学位移（反推检索必需） |

放置：`E:\NMR2Structure\data\nmrshiftdb2\nmrshiftdb2withsignals.sd`
不需要下载：无信号的 `nmrshiftdb2.sd`、1.1 GB 的 SQL dump、`predictor*.jar`（Java 工具用不上）。

## 未来升级: 天然产物专属数据库 (留待下载)

以下两个数据库包含天然产物专属 NMR 位移数据，命中精度高于通用库，但需手工下载后集成:

| 数据库 | 规模 | 位移类型 | 下载来源 | 预计体积 |
|---|---|---|---|---|
| **LOTUS+13C** | 21.8万 NP | ACD 预测 13C | Zenodo 搜 `LOTUS predicted 13C NMR Nuzillard` (文件名 `acd_lotusv9.sdf.zip`) | ~200-400 MB |
| **NP-MRD** | 28万 NP, 550万谱 | 实验 13C/1H/J/2D | https://np-mrd.org/downloads → 2D Structures SDF (7卷) + Peaklist TXT (7卷) + Assignments CSV (453 KB) | ~510 MB |

目标放置位置 (暂未实现自动发现):
```
data/lotus_13c/acd_lotusv9.sdf
data/np-mrd/structures/
data/np-mrd/peaklists/
data/np-mrd/assignments.csv
```

计划升级: 下载后实现 `retrieve_lotus_13c.py` (复用现有 SDF 解析器) 和 `retrieve_npmrd.py` (解析 SDF + TXT + CSV), 接入检索级联 [3a-3b] 之间。

## 快速开始

```bash
cd E:/NMR2Structure/NMR2Structure
# 默认模式 (derive + 2D): 从谱图直接推导归属表, 不进入检索
python -m nmr_elucidator --folder "E:/NMR2Structure/NMR/clxh100"

# 推导 + 检索: derive 完成后继续候选检索+打分
python -m nmr_elucidator --folder "E:/NMR2Structure/NMR/clxh100" --retrieve

# 纯检索模式: 跳过 derive, 直接走数据库检索+打分
python -m nmr_elucidator --folder "E:/NMR2Structure/NMR/clxh3" --no-derive
```

程序内部已配置 UTF-8 编码, 无需 `python -X utf8`。
输出落到 `E:\NMR2Structure\output\clxh3\`（默认输出根目录 + 从 folder 路径末段推断化合物名）。
化合物名显式指定用 `--compound-name my_compound`，输出根目录改 `--output-root D:/my_output`。

## 主流程 CLI

```
python -m nmr_elucidator --folder <原始数据/<sample>> [选项...]
```

默认行为:
- `--derive` 开启: 从谱图直接推导归属表 (HSQC/COSY/HMBC/ROESY), 默认仅推导不检索
- `--use-2d` 开启: 自动检测 2D 数据, 无则回退 1D
- `--formula` 不传: 适用于未知化合物
- `--db` 默认: NMR2Structure 外层平级 `data/` 目录 (自动发现 nmrshiftdb2 + lucy-ng)
- UTF-8 编码: 程序内部配置, 无需 `-X utf8`

关闭默认: `--no-derive` (跳过推导, 直接检索) / `--no-2d` (强制 1D)
开启检索: `--retrieve` (derive 后继续检索) / `--no-derive` (纯检索模式)

检索级联: nmrshiftdb2 SDF → lucy-ng 928K → COCONUT 69.5万 (本地无命中自动联网, `--no-coconut` 禁用) → PubChem 1亿+ (`--no-pubchem` 禁用)

### 辅助工具脚本（util/）

| 脚本 | 用途 |
|---|---|
| `util/install_deps.ps1` | 一键安装所有项目依赖 (PowerShell) |
| `util/download_lucy_hose_db.py` | 下载 lucy-ng HOSE 完整数据库 (~830 MB / ~2.8 GB)；支持 Colab / aria2 / curl / proxy 4 种策略 |

### 测试

无 pytest; `tests/run_all.py` 串行 `import` 各 `test_*.py`.

```powershell
cd E:\NMR2Structure\NMR2Structure
python tests\run_all.py                # 全部
python tests\test_bridge_dept.py       # 单文件 (有 __main__ 块)
python tests\test_io_processed.py
python tests\test_new_features.py
```

`test_io_bruker.py` / `test_io_sdf.py` / `test_retrieve.py` / `test_score.py` 只有函数定义 (无 `__main__`), 仅通过 `run_all` 触发.

### 清理

```powershell
Get-ChildItem -Path . -Recurse -Filter __pycache__ | Remove-Item -Recurse -Force
```

## 解谱工作流对标

按 `doc/MestReNova_解谱指南.md` 5.4 节化学家实操 6 步 + 5.1 节拼图策略, 当前程序能力对照:

| 步骤 | 化学家动作 | 程序实现 | 状态 |
|---|---|---|---|
| 1.1 | 1H 谱积分 + 标峰 | `io_bruker.py` 读 peaklist.xml + integrals.txt | 已做 |
| 1.2 | 1H 多重峰分析 (J 耦合) | `io_bruker.parse_multiplets_from_peaklist` (含 `infer_multiplets_from_peaks` 自动反推 fallback) | 已做 |
| 1.2b | **TopSpin Manual 多重峰精细化分析** | `io_bruker.analyze_multiplet_at_ppm` 对每个 CHPair 单独做局部窗口多重峰拟合 (pattern/J/halfwidth/质量等级) → `step1b_multiplet_analysis.csv` | 已做 |
| 1.3 | 13C 标峰 + DEPT 推断级数 | `bridge_nmr2struct.assign_carbon_multiplicity` | 已做 |
| 1.4 | 化学位移经验区段初判 | `shift_ranges` | 已做 |
| 2.1 | HSQC 建立 (dH, dC) 对 | `interpret_2d._step1_hsqc_to_ch_pairs` | 已做 |
| 2.2 | HSQC 投影叠加 | 等高线图 1H/13C 双投影 | 已做 |
| 2.3 | HSQC-DEPT 颜色区分 | 等高线图按 multiplicity 着色 (CH/CH3=红, CH2=蓝) | 已做 |
| 2.4 | 季碳从 13C 扣除 HSQC 对应 | `_step4_infer_quaternary_carbons` DEPT 优先 | 已做 |
| 2.5 | 重叠峰 13C 切片分解 | `_step1b_resolve_overlapping` (5.9 节) → `step1b_overlap_regions.csv` (fallback) | 已做 |
| 2.5b | **精细化 1H 多重峰分析** | `_step1b_multiplet_analysis_for_ch_pairs` 对每个非季碳 C-H 对单独取窗口内 Peak1D, 拟合 pattern/J/半高宽/强度均匀性/质量等级 → `step1b_multiplet_analysis.csv` (主输出) | 已做 |
| 3.1 | COSY 拼自旋系统 | `_step2_cosy_spin_systems` Union-Find + J 验证 | 已做 |
| 3.2 | TOCSY 长程自旋系统 | `_step2_tocsy_spin_systems` (h_tol=0.6) | 已做 |
| 4.1 | HMBC 找 2-3 键外 C | `_step3_hmbc_connections` | 已做 |
| 4.2 | HMBC 排除 1J 干扰 | c_tol 过滤 | 已做 |
| 4.3 | HMBC 区分 2J vs 3J | `_annotate_hmbc_bond_type` 自旋系统拓扑 | 已做 |
| 5.1 | NOESY/ROESY 谱读取 | `io_bruker_2d.detect_2d_experiment_type` | 已做 |
| 5.2 | NOESY 空间近邻 H-H 提取 | `_step5_extract_noe_contacts` | 已做 |
| 5.3 | 候选 3D 构象 NOE 验证 | `correlation_2d.verify_candidate_3d_noe` | 已做 |
| 6.1 | 完整归属表 (论文 Table 1) | `build_full_assignment_table` Top-N CSV | 已做 |
| 6.2 | 归属完整性自动核对 | `validate_assignment` | 已做 |
| 6.3 | 候选结构检索 | `retrieve.retrieve_candidates` nmrshiftdb2 | 已做 |
| 6.4 | 综合打分排序 | `evidence.score_candidates` (含 NOE 加分) | 已做 |
| 6.5 | 置信度四级评估 | `determination.assess_determination` | 已做 |

### step1b 多重峰数据在后续流程的应用规划

`step1b_multiplet_analysis.csv` 提供每个非季碳 C-H 对的多重峰精细化拟合
(pattern/J/半高宽/强度均匀性/质量等级), 在项目后续环节具有以下增强潜力:

| 环节 | 当前位置 | step1b 增强方式 | 状态 |
|---|---|---|---|
| **HSQC 配对检索** | `__main__.py` 提取 `[(h_ppm, c_ppm), ...]` | 扩展为 `[(h_ppm, c_ppm, multiplicity, j_pattern), ...]`, 候选中对应 H 的 J 模式应匹配实测 | 未做 |
| **候选打分 (J-RMS 维度)** | `evidence.score_candidates` (13C/1H RMS + DEPT coverage + NOE) | 加 `J_RMS` 维度: 实测 J 与候选结构 RDKit 3D 构象二面角预测的 J 对比 (Karplus 方程) | 未做 |
| **COSY 交叉峰验证** | `_step2_cosy_spin_systems` 用 `\|J\|>4 Hz` 过滤 | 进一步: 如果实测 H 是 s (无 J) 但 COSY 显示交叉峰, 视为假阳性并降权 | 未做 |
| **3D 构象 J 验证** | `correlation_2d.verify_candidate_3d_noe` 仅 NOE 距离 | 加 `verify_candidate_3d_j`: ETKDG 嵌入 + MMFF 优化 + 二面角扫描 + Karplus 预测 J, 与实测 J_RMS 比对 | 未做 |
| **归属置信度** | `determination.assess_determination` (基于 1D/2D 一致性) | 加 J 一致性维度: 实测 pattern (s/d/t/q/m) 应与候选 CH/CH2/CH3 结构相符 | 未做 |
| **邻位 H 拓扑反推** | `_fill_h_neighbors_from_spin_systems` 仅用 COSY | 用 J 反推: `\|J\|>5 Hz` → 3J 直接邻位, `1-5 Hz` → 邻位或 2J 长程, `<1 Hz` → 4J+ 长程 | 未做 |
| **化学家人工复核** | 看 step1 CSV 数值 | 重点复核 **Very Poor** 行 (强度剖面严重偏离理论), 这类 C-H 对归属需手动验证 | 已支持 (按 quality_grade 排序) |
| **Karplus 方程常数表** | 无 | 需维护: ³J(H-C-C-H) Karplus A/B/C 系数 (典型 7.0/-1.0/0.0; 烯烃顺/反 7.0/-1.5/2.0; 芳烃 8.5/-1.0/0.5) | 未做 |

**优先级建议**:
1. **Karplus 3D J 验证** (优先级高): 与 NOE 3D 验证互补, 二面角 + 距离双维度交叉, 显著提高候选结构可靠性
2. **HSQC 配对检索加入 J 模式**: 减少误匹配 (CH2 vs CH3 不能仅靠 ppm 区分)
3. **归属置信度 J 一致性**: 简单加权即可, 让 determination 输出更稳定

**实施示例** (Karplus 3D J 验证):

```python
def verify_candidate_3d_j(candidate_smiles, j_observations, h_assignments,
                          field_mhz, num_confs=10):
    """对候选 SMILES 做 3D 构象采样 + Karplus J 预测 vs 实测.

    Args:
        candidate_smiles: 候选 SMILES
        j_observations: [MultipletAnalysis, ...] 实测 J (Hz)
        h_assignments: {h_ppm: atom_idx} H 归属 (从 build_assignment 推断)
        field_mhz: 实测频率 (Hz/ppm 转换)

    Returns:
        J_RMS_Hz, J_match_count, J_conflict_count
    """
    mol = Chem.MolFromSmiles(candidate_smiles)
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs)
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    # 对每个 conf 提取 H-C-C-H 二面角, 用 Karplus 计算 ³J
    # 取与实测 J_RMS 最小的 conf 为最佳
```

### 必填

| 参数 | 说明 |
|---|---|
| `--folder` | Bruker 化合物根目录（含数字子目录 21/22/...，每个下有 `pdata/1/peaklist.xml`）|

### 数据库

| 参数 | 默认 | 说明 |
|---|---|---|
| `--db` | `E:\NMR2Structure\data` | 数据库路径: 传入总目录自动发现所有数据库, 或传入单个文件路径 |
| `--lucy-ng-db` | 自动发现 | lucy-ng 数据库路径覆盖 (优先于 --db 自动发现) |

`--db` 支持两种模式:
- **目录**: 自动扫描 `nmrshiftdb2/nmrshiftdb2withsignals.sd` + `lucy-ng-derep_db/lucy-ng-derep.db`
- **文件**: 按扩展名判断 (.sd → nmrshiftdb2, .db → lucy-ng)

检索级联: nmrshiftdb2 SDF → lucy-ng 928K → COCONUT 69.5万 (天然产物) → PubChem 1亿+ (本地无命中自动联网)

### 检索控制

| 参数 | 默认 | 说明 |
|---|---|---|
| `--formula` | 无 | 硬约束分子式（如 `C8H8O2`），命中候选必须含该分子式 |
| `--top-n` | 10 | 输出 Top-N 候选数 |
| `--max-hits` | 不限 | 限制候选命中数（快速验证用，如 30 表示只取前 30 个） |
| `--c-tol` | 2.0 | 13C 化学位移窗口容差 (ppm) |
| `--h-tol` | 0.3 | 1H 化学位移窗口容差 (ppm) |

### 输出路径

| 参数 | 默认 | 说明 |
|---|---|---|
| `--output-root` | `E:\NMR2Structure\output` | 输出根目录，每个化合物一个子目录 |
| `--compound-name` | 从 `--folder` 路径末段推断 | 子目录与文件名前缀 |

例：`--folder .../clxh3` 等价于 `--compound-name clxh3`，产物全部落到 `E:\NMR2Structure\output\clxh3\`。

### 可选附加步骤

| 参数 | 说明 | 产物 |
|---|---|---|
| `--use-nmrsim` | nmrsim 反向校验每个候选 1H 谱的 cosine 相似度 | `.md` 报告里追加 `nmrsim 反向校验 cosine = X.XXX` |
| `--validate-with-nmrpy` | 用 nmrpy 从 `fid` 时域信号二次处理（FT/相位/基线/峰检测），与 peaklist.xml 交叉验证 | `<prefix>_nmrpy_validation.{json,md,png}` |
| `--with-real-spectrum` | 读取 `pdata/1/1r` 真实处理后谱，绘制真实谱线形图（叠加 peaklist 峰位） | 图片 base64 内嵌到 HTML 报告 |
| `--use-hose` | HOSE code 13C 预测验证；**双数据源策略**：默认自动查找 `data/lucy-ng-derep_db/lucy-ng-derep.db` (928K 化合物, MAE ~3.5 ppm)；不存在时自动从 nmrshiftdb2 SDF 自建轻量本地表 (默认 50K, MAE ~5 ppm) | JSON/MD 报告里追加 `hose_verifications[]` + `meta.used_hose_source` (`lucy-ng-derep.db` 或 `local-sdf`) |
| `--use-pubchem` | PubChem **强制**扩展检索：即使本地库有命中也联网检索 PubChem (1亿+ 化合物)。默认行为: 本地库无命中时自动触发, 无需此参数 | JSON/MD 报告里追加 `pubchem_candidates[]` 列表 |
| `--no-pubchem` | 禁用 PubChem 检索 (即使本地库无命中也不联网) | — |
| `--use-coconut` | COCONUT **强制**扩展检索：即使本地库有命中也联网检索 COCONUT (69.5万+ 天然产物结构)。默认行为: 本地库无命中时自动触发 (天然产物类样本优先) | JSON/MD 报告里追加 `coconut_candidates[]` 列表 |
| `--no-coconut` | 禁用 COCONUT 检索 (即使本地库无命中也不联网) | — |
| `--use-fragment-assembly` | 片段组装候选生成：当数据库检索均无命中时，从 13C 位移范围 + DEPT 约束组装结构（仅 1D 谱可走） | JSON 报告里追加 `fragment_candidates[]` |
| `--use-2d` | 2D NMR 关联验证（**默认开启**, `--no-2d` 关闭）：自动检测 COSY/HSQC/HMBC/ROESY 实验；按化学家工作流顺序解读（HSQC → COSY → HMBC → 季碳）；含重叠峰分解；Top-1 生成完整 NMR 归属表；无 2D 数据时自动回退到 1D 流程 | JSON/MD/HTML 报告里追加 `interpretation_2d` + `full_assignment` 章节；中间结果独立 CSV 落盘 |
| `--determine` | 结构推断置信度评估：多源证据融合打分，输出 determined / probable / ambiguous / no_match 四级结论 | JSON/MD/HTML 报告顶部追加 `determination` 章节 |
| `--hose-cache <path>` | HOSE 查找表缓存路径（默认 `<sdf 同目录>/_hose_lookup_cache.tsv`，避免重复建表） | — |
| `--hose-max-records <N>` | 从 SDF 构建 HOSE 查找表的最大记录数（默认 **50000**；0=全部；10000 轻度，50000 够用，100000+ 需约 5 分钟） | — |
| `--lucy-ng-db <path>` | lucy-ng 完整 HOSE 数据库路径（`lucy-ng-derep.db`, ~2.8 GB SQLite）; **默认已自动查找** `data/lucy-ng-derep_db/lucy-ng-derep.db`, 无需手动指定; 用此选项可指向备份位置 | — |
| `--log-dir <path>` | 日志输出目录 (默认 `E:\NMR2Structure\logs`); 每次运行自动生成 `NMR2Structure_YYYY-MM-DD_HH-MM-SS.log` 文件, 包含屏幕上的所有内容 (print + logging) | — |
| `--no-log` | 禁用日志记录 (默认每次运行都写日志) | — |
| `--no-derive` | 禁用 derive 模式, 直接走候选检索+打分流程 | — |
| `--retrieve` | derive 完成后继续执行候选检索+打分 (默认仅推导不检索) | — |
| `--no-2d` | 禁用 2D NMR 检测, 强制使用 1D 流程 (默认 2D 开启) | — |
| `--real-structure-dir <path>` | 已知化合物根目录 (默认 `E:\NMR2Structure\真实结构`). 当该目录下存在 `<sample>/数据归属*.xlsx` 时自动切换到已知化合物流程 (derive 不使用 SMILES, 复核阶段才用于结构判等); 根目录/子目录/归属表缺失时自动跳过复核 | `derive_归属表/assignment_table_joint.xlsx` + `comparison_report.txt` |
| `--known-c-tol <ppm>` | 已知化合物复核对比 13C 匹配容差 (默认 0.5 ppm, 覆盖实际推导偏差) | — |
| `--known-h-tol <ppm>` | 已知化合物复核对比 1H 匹配容差 (默认 0.1 ppm) | — |

### 已知化合物复核流程

如果 `E:\NMR2Structure\真实结构/<compound_name>/` 目录下存在化学家 Excel 标准归属表
(`*数据归属*.xlsx`, 同时推荐同目录放 `smiles.txt`), 运行主流程时会自动:

1. **推导阶段不依赖 SMILES**: 仍然从 1D + 2D 谱图直接推导归属表 (与未知化合物流程一致).
   SMILES **不参与推导**, 只用于复核阶段的化合物结构判等 (RDKit canonical SMILES 一致性).
2. **复核对比**: derive 完成后, 自动以 (13C ±0.5 ppm, 1H ±0.1 ppm, 碳类型含 CH2/CH3 占位符兼容) 做行级联合匹配.
3. **产出文件**:
   - `derive_归属表/assignment_table_joint.xlsx` 双 sheet 左右拼接表:
     - Sheet 1 **按标准表为准**: 行为标准表所有行 (行序与标准表一致);匹配上的程序输出行拼接在右侧并填 ΔdC / ΔdH / 类型一致;未匹配行的程序侧列填 `—`.
     - Sheet 2 **按程序表为准**: 行为程序输出所有行 (行序与程序 xlsx 一致);匹配上的标准表行拼接在右侧并填 ΔdC / ΔdH / 类型一致;未匹配行的标准侧列填 `—`.
     - 列结构 (左右两例): [表 A 全部列] + '拼接' + [表 B 全部列] + ΔdC (ppm) + ΔdH (ppm) + 类型一致
   - `derive_归属表/comparison_report.txt` (含化合物结构信息 + 逐行匹配明细 + 综合判定)
4. **判定等级**: 匹配率 ≥80% & 碳类型全一致 → "一致"; ≥60% → "部分一致"; <60% → "不一致".

设计原则: **化学家标号 = ground truth**. SMILES 仅作为最终结构判等依据 (InChIKey / canonical SMILES),
不参与行匹配. 推导过程的化合物结构验证完全依赖谱图推断.

**Python API**:

```python
from nmr_elucidator.standard_assignment import (
    discover_known_compound, read_standard_assignment_table,
    compare_program_to_standard, write_matched_xlsx, write_comparison_report,
)

# 1. 发现
result = discover_known_compound('clxh100', r'E:\NMR2Structure\真实结构')
table_path, smiles = result  # table_path, smiles 或 None

# 2. 读取标准归属表
std_table = read_standard_assignment_table(table_path, known_smiles=smiles)

# 3. 对比程序输出
comparison = compare_program_to_standard(
    'output/clxh100/derive_归属表/assignment_table.xlsx', std_table,
    c_tol=0.5, h_tol=0.1, use_smiles=True,
)

# 4. 写匹配 xlsx + 报告
write_matched_xlsx('output/.../assignment_table.xlsx',
                   'output/.../assignment_table_matched.xlsx', comparison)
write_comparison_report('output/.../comparison_report.txt',
                        'output/.../assignment_table.xlsx', std_table, comparison)
```

### 运行日志

每次运行 CLI / 子命令时, 项目会在 `E:\NMR2Structure\logs\` 下生成一个独立的日志文件:

```
E:\NMR2Structure\logs\
├── NMR2Structure_2026-07-22_22-09-13.log   # 此次运行的完整记录
├── NMR2Structure_2026-07-22_22-12-59.log
└── ...
```

日志文件格式 (UTF-8):

```
==============================================================================
日志文件: E:\NMR2Structure\logs\NMR2Structure_2026-07-22_22-09-13.log
启动时间: 2026-07-22T22:09:13
命令行参数: '...'
Python: 3.13.5  平台: win32
工作目录: E:\NMR2Structure\NMR2Structure
==============================================================================
[1/5] 读取 Bruker 谱图: E:\NMR2Structure\NMR\clxh3
  - exp 21: 1H (zg) solvent=Pyr peaks=103 integrals=13
  ...
DONE.
```

- **TEI 机制**: `_TeeTextIO` 包装 sys.stdout / sys.stderr, 所有 print() 原文与 logging 输出同步进入文件, 无需修改现有 25+ 处 print 调用
- **DEBUG 级日志**: 仅写入文件, 不进屏幕 (控制台为 INFO 级)
- **未捕获异常**: traceback 自动写入日志
- **覆盖子命令**: `__main__` / `run_bridge` / `nmrpy_validate` 三个入口都启用
- **环境变量覆盖**: `NMR_ELUCIDATOR_LOG_DIR` / `NMR_ELUCIDATOR_LOG_NAME` 可改默认路径和文件名前缀

**Python API**:

```python
from nmr_elucidator.log_setup import setup_logging, get_logger, get_current_log_file

# 自定义目录 + 文件名前缀
log_file = setup_logging(log_dir='D:/mylogs', name_prefix='myexperiment')

# 获取项目 logger (子模块化)
log = get_logger('evidence')  # -> nmr_elucidator.evidence
log.info('开始计算证据链')
log.debug('matched 24/30 peaks')
```

### 其他

| 参数 | 说明 |
|---|---|
| `--quiet` | 不打印检索进度（每 5000 条打印一行） |

## 子模块（独立运行）

主流程会自动调用，但也可单独用：

### `nmrpy_validate` — FID 二次处理 + 交叉验证

```bash
python -m nmr_elucidator.nmrpy_validate \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmrpy_check"
```

输出 `<prefix>_validation.{json,md,png}` 三件套。**当前局限**：nmrpy peakpicker 默认参数在低浓度样品（如 clxh3）上检出 0 峰，需手工调 voff/threshold；ppm 轴是粗略近似（`freq/sfo1 + 100`）。

### `run_bridge` — NMR2Struct 端到端 ML

需要 `external/NMR2Struct` 目录已下载权重（见主仓库 `从NMR谱图反推分子结构代码库.md`）。

```bash
# stage 1: 数据转换 (无需 torch)
python -m nmr_elucidator.run_bridge \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmr2struct" --stage 1

# stage 2: NMR2Struct 推理 (需 torch + 已下载权重)
python -m nmr_elucidator.run_bridge \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmr2struct" --stage 2

# stage 3: h5 -> Markdown (无需 torch)
python -m nmr_elucidator.run_bridge \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmr2struct" --stage 3
```

## Python API

包内符号可直接 import（见 `nmr_elucidator/__init__.py` 的 `__all__`）。

### 读 Bruker + 检索 + 打分（pipeline 复用）

```python
from nmr_elucidator import (
    read_bruker_folder,
    Constraints,
    retrieve_candidates,
    score_candidates,
    build_evidence,
)

specs = read_bruker_folder('E:/NMR2Structure/NMR/clxh3')
h_spec = next(s for s in specs if s.nucleus == '1H')
c_spec = next(s for s in specs if s.nucleus == '13C' and 'DEPT' not in s.nucleus)

# 简化版 ProblemPeaks 构造（真实程序用 _build_problem_from_bruker）
problem = type('PP', (), {
    'carbon_shifts': [p.shift for p in c_spec.peaks],
    'proton_shifts': [p.shift for p in h_spec.peaks],
    'proton_integrals': [],
    'solvent': h_spec.solvent,
    'field_mhz': h_spec.field_mhz,
})()

constraints = Constraints(formula=None)
hits = list(retrieve_candidates(
    sdf_path='E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd',
    problem=problem, constraints=constraints,
    c_tol=2.0, h_tol=0.3, max_hits=200,
))
scored = score_candidates(hits)
evidences = [build_evidence(s, problem.carbon_shifts, problem.proton_shifts) for s in scored]
```

### 读 1r 真实谱 + 重采样到 ML 网格

```python
from nmr_elucidator import read_bruker_folder
from nmr_elucidator.io_processed import read_processed_folder, resample_to_grid

specs = read_bruker_folder('.../clxh3')
specs_proc = read_processed_folder('.../clxh3', specs=specs)

# 1H 谱重采样到 ChefNMR h_10k 网格 (-2 ~ 10 ppm, 10000 维)
h_grid = resample_to_grid(specs_proc[21].ppm, specs_proc[21].intensity,
                          grid_min=-2.0, grid_max=10.0, n_points=10000)
# 13C 谱重采样到 ChefNMR c_80 网格 (3.42 ~ 231.3 ppm, 80 维)
c_grid = resample_to_grid(specs_proc[22].ppm, specs_proc[22].intensity,
                          grid_min=3.42, grid_max=231.3, n_points=80)
# 同样适用于 NMR2Struct CNMR 80 维 [3.4, 231.3]
```

### 写报告 + 绘图

```python
from nmr_elucidator import write_report_json, write_report_markdown
from nmr_elucidator.drawing import draw_topn_grid, draw_unmatched_overlay

write_report_json('out.json', evidences, meta={'top_n': 10})
write_report_markdown('out.md', evidences, meta={'top_n': 10}, top_n=10)
draw_topn_grid(evidences, 'out_topn.png', top_n=10)
draw_unmatched_overlay(observed_c=problem.carbon_shifts,
                       observed_h=problem.proton_shifts,
                       evidences=evidences, out_path='out_overlay.png', top_n=5)
```

### HOSE 预测

```python
from nmr_elucidator.hose_generator import prepare_mol, generate_hose_codes_for_carbons
from nmr_elucidator.nmr_predict import get_predictor, verify_candidate

# 从 nmrshiftdb2 SDF 自建 HOSE 查找表 (只需一次, 后序读缓存)
predictor = get_predictor(
    sdf_path='E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd',
    table_path='E:/NMR2Structure/data/nmrshiftdb2/_hose_lookup_cache.tsv',
    max_records=10000,
)

# 对任意 SMILES 生成每个碳的 HOSE 编码并预测化学位移
mol = prepare_mol('CC(=O)Oc1ccccc1C(=O)O')  # 阿司匹林
codes = generate_hose_codes_for_carbons(mol, max_radius=4)
predictions = predictor.predict_mol(mol, max_radius=4)  # list[PredictedShift]

# 与实测 13C 验证 (返回 RMS / coverage / matched_carbon_ratio)
observed_c13 = [170.5, 150.2, 133.8, 129.5, 126.3, 123.7, 21.0]
vs = verify_candidate(predictor, 'CC(=O)Oc1ccccc1C(=O)O', observed_c13, c_tol=3.0)
print(f'RMS={vs.hose_rms:.2f} coverage={vs.hose_coverage:.2f} matched={vs.hose_match_rate:.2f}')
```

### lucy-ng 完整数据库查询（默认已启用，需本地 2.8 GB SQLite）

```python
from nmr_elucidator.lucy_ng_hose import (
    LucyNGHOSEProvider, verify_candidate_with_lucy_db,
)

# 打开 lucy-ng SQLite 数据库 (默认查找 data/lucy-ng-derep_db/lucy-ng-derep.db)
provider = LucyNGHOSEProvider()  # 也可传 db_path=...
result = provider.predict_smiles('CC(=O)Oc1ccccc1C(=O)O')
for ps in result.predictions:
    print(f'C{ps.atom_index}: {ps.predicted_ppm:.2f} ppm  radius={ps.radius_used}  matches={ps.match_count}')
# 与本地 50K 表对比: lucy-ng RMS 约 1.4-1.7 ppm (本地约 1.75-2.0 ppm)

vs = verify_candidate_with_lucy_db(
    provider, 'CC(=O)Oc1ccccc1C(=O)O',
    observed_c13=[170.5, 150.2, 133.8, 129.5, 126.3, 123.7, 21.0],
    c_tol=3.0,
)
print(f'lucy-ng RMS={vs.hose_rms} ppm, verification_score={vs.verification_score}')
provider.close()
```

### 结构推断置信度评估

```python
from nmr_elucidator.determination import assess_determination, compute_unified_score

# 为每个候选计算融合分 (HOSE RMS + coverage + DEPT + 分子式一致性)
unified_scores = [compute_unified_score(e, hose_verification=vs) for e, vs in zip(evidences, hose_verifications)]

# 融合打分 + 输出四级结论
result = assess_determination(
    evidences=evidences,
    hose_verifications=hose_verifications,
    unified_scores=unified_scores,
)
print(f'结论: {result.confidence_label} ({result.conclusion})')
print(f'Top-1 SMILES: {result.top_candidate_smiles}')
print(f'证据链: {result.evidence_summary}')
```

### PubChem 扩展检索 + 片段组装

```python
from nmr_elucidator.retrieve_pubchem import retrieve_pubchem_candidates
from nmr_elucidator.fragment_assembler import build_constraints, assemble_candidates
from nmr_elucidator.shift_ranges import infer_functional_groups

# 从 13C 推断官能团
observed_c13 = [170.5, 150.2, 133.8, 129.5, 126.3, 123.7, 21.0]
groups = infer_functional_groups(observed_c13)
# → ['酯羰基', '芳香 C-O', '芳香 C', '脂肪 CH3', ...]

# PubChem 检索 (按分子式, 可选按子结构)
pubchem_cands = retrieve_pubchem_candidates(
    formula='C9H8O4',  # 阿司匹林分子式
    max_results=50,
)

# 片段组装 (从 DEPT + 化学位移范围 + 分子式生成库外结构)
constraints = build_constraints(
    dept_summary={'CH': 5, 'CH2': 0, 'CH3': 1, 'Cq': 3},
    proton_total=8,
    formula='C9H8O4',
)
assembled = assemble_candidates(constraints, max_heavy_atoms=12, max_candidates=200)
```

## 输出文件清单

每个化合物目录（如 `E:\NMR2Structure\output\clxh3\`）下：

### 始终产出（主流程）

| 文件 | 说明 |
|---|---|
| `<sample>.json` | 完整证据链（matched/unmatched/conflict），程序可读 |
| `<sample>.md` | Markdown 报告，含未匹配峰专项章节 + pandas 候选表（InChIKey 去重）|
| `<sample>.html` | HTML 报告（人类阅读，含流程 + 结构图 + DEPT 推断 + seaborn 统计可视化，所有图片 base64 内嵌）|
| `<sample>_unmatched.csv` | 未匹配峰汇总（Excel 友好，UTF-8 BOM） |

> 注：所有 PNG 图片（topn / rank / overlay / heatmap / scatter / matrix / overview / real_spectrum）
> 已 base64 内嵌到 HTML 中，输出目录不再保留独立 PNG 文件。

### 四种输出格式内容对比

| 内容 | JSON | MD | HTML | CSV |
|---|:---:|:---:|:---:|:---:|
| 完整 meta（含 observed 原始位移数组） | ✅ | 部分 | 部分（图表展示） | — |
| **全部**候选数据（不限 top_n） | ✅ | ❌ | ❌ | ❌ |
| DEPT 推断表（碳类型 + 方向） | ✅ | ✅ | ✅ | — |
| 2D NMR 关联验证（HSQC/HMBC/COSY） | ✅ | ✅ | ✅ | — |
| **2D 解读路径索引**（步骤 1/1b/2/3/4 + Top-1 指向各独立文件） | ❌ | ❌ | ✅ | — |
| **HSQC 等高线图**（Mnova 风格, overlap/normal 分文件夹） | ❌ | ❌ | ❌ (不嵌入) | ✅ (PNG 子目录) |
| **Top-N 完整 NMR 归属表**（论文 Table 1 格式，HTML 中显示几个 top 就生成几个 CSV） | ❌ | ❌ | ❌ (仅索引目录) | ✅ (专属子目录) |
| 每候选 13C/1H 逐峰匹配明细 | ✅ | ✅ | ✅ | — |
| 每候选未匹配峰（逐条） | ✅ | ✅ | ✅ | ✅ |
| InChIKey 去重候选表（pandas） | ❌ | ✅ | ❌ | — |
| 结构图 / 谱图可视化 | ❌ | ❌ | ✅（base64 内嵌） | — |
| seaborn 统计图（热力图/散点/矩阵） | ❌ | ❌ | ✅（base64 内嵌） | — |
| 未匹配峰汇总（Excel 可筛选排序） | — | — | — | ✅ |
| 程序可读结构化证据链 | ✅ | ❌ | ❌ | — |

**定位说明**：
- **JSON** = 程序接口（完整数据，含全部候选 + 原始位移数组，供下游自动化）
- **MD** = 纯文本备查（含 InChIKey 去重表，适合 diff / 版本管理）
- **HTML** = 人类阅读（可视化 + 单文件可移植，图片全部内嵌）
- **CSV** = Excel 工作流（化学家按候选筛选/排序未匹配峰）

### `--with-real-spectrum` 附加

不新增独立文件，真实谱图（4 子图：每个实验的真实谱曲线 + peaklist 峰位）已 base64 内嵌到 HTML 报告中。

### `--validate-with-nmrpy` 附加

| 文件 | 说明 |
|---|---|
| `<sample>_nmrpy_validation.json` | FID 交叉验证数据（各实验峰数 + 2ppm 重叠率） |
| `<sample>_nmrpy_validation.md` | 人类阅读报告 |
| `<sample>_nmrpy_validation.png` | 4 子图：每个 FID 的 nmrpy 处理结果 |

### `--use-hose` 附加

不新增独立文件，**丰富 JSON/MD/HTML 报告**：
- `meta.hose_verifications[]`：每个候选的 `{smiles, hose_rms, hose_coverage, hose_match_rate, verification_score, notes}`
- `meta.used_hose: true` 标记启用
- Top-1 候选的 HOSE RMS 嵌入证据链

### HOSE 查找表进阶：lucy-ng 928K (默认) vs nmrshiftdb2 本地表 (回退)

**默认路径**：自动使用 lucy-ng 928K SQLite (~2.8 GB, MAE ~3.5 ppm):

```bash
python -m nmr_elucidator --folder ... --use-hose
# 默认自动查找 data/lucy-ng-derep_db/lucy-ng-derep.db
```

**回退路径**：如果 lucy-ng 数据库不可用, 从本项目 `nmrshiftdb2 SDF` 自建轻量本地表 (默认 50K 条, 约 1-2 分钟, MAE ~5 ppm):

```bash
python -m nmr_elucidator --folder ... --use-hose --hose-max-records 50000
```

**升级路径**：下载 lucy-ng 完整 HOSE 库（928K 化合物，~830 MB 下载，~2.8 GB 解压，MAE ~3.5 ppm）：

```bash
# 1. 下载（4 种策略可选）
python util/download_lucy_hose_db.py --strategy colab    # 生成 Google Colab Notebook
python util/download_lucy_hose_db.py --strategy aria2    # 生成 aria2 多线程脚本

# 2. 解压到 E:\NMR2Structure\data\lucy-ng-derep_db\
#    pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps

# 3. 使用完整库（精度提升明显）
python -m nmr_elucidator --folder ... --use-hose
# 默认自动查找 data/lucy-ng-derep_db/lucy-ng-derep.db
# 也可指定：--lucy-ng-db D:\my\lucy-ng-derep.db
```

**优先级**：`--lucy-ng-db` 指定的库存在时，优先用它；否则回退本地 SDF 表。

**实测效果**（clxh3，Top-1 候选 `CC1CC(C)c2ccccc21`）：

| 库 | RMS | 覆盖率 | 验证分 |
|---|---|---|---|
| 本地 10K 表 | 1.85 ppm | 100% | 0.611 |
| 本地 50K 表 | ~1.75 ppm | 100% | ~0.68 |
| **lucy-ng 928K** | **1.70 ppm** | 100% | **0.741** |

figshare 在国内无法直接访问，Google Colab 是官方推荐的曲线下载方案。

### `--use-pubchem` 附加

- `meta.pubchem_candidates[]`：PubChem 检索到的候选列表（含 `cid, smiles, iupac_name, molecular_formula, molecular_weight`）

### `--use-fragment-assembly` 附加

- `meta.fragment_candidates[]`：片段组装的库外结构候选（含 `smiles, fragment_set, source`）

### `--determine` 附加

在 JSON/MD/HTML 报告**顶部追加** `determination` 章节：
- `conclusion`：推断结构 SMILES
- `confidence`：四级（determined / probable / ambiguous / no_match）
- `confidence_label`：中文标签（确定 / 很可能 / 不确定 / 无匹配）
- `top_candidate_smiles`：排名第一的 SMILES
- `evidence_summary`：证据链文本摘要
- `unified_score`：Top-1 融合分
- `gap_to_runner_up`：与第二名差距

### `--use-2d` 附加（化学家工作流详解）

启用后，额外产出 **6 个独立中间文件**（化学家人工复核用，UTF-8 CSV/PNG），报告新增内容：

| 文件 | 内容 |
|---|---|
| `step1_hsqc_ch_pairs.csv` | HSQC C-H 直接连接对（52 行示例）。列：编号/C_label/H_label/δH/δC/1H积分/H个数/碳级数/J_couplings_Hz/h_neighbors_ppm（精细化 1H 多重峰分析覆盖的 J 与 COSY 反推的邻位 H）。 |
| `step1_hsqc_ch_pairs/overlap/*.png` | 积分 ≥ 4 的 C-H 对 HSQC 等高线图（Mnova 风格：2D等高线+1H投影+13C投影+13C切片检峰），可目视检查重叠峰。 |
| `step1_hsqc_ch_pairs/normal/*.png` | 非重叠峰等高线图（同上布局）。 |
| `step1b_multiplet_analysis.csv` | **精细化 1H 多重峰分析**（主输出，对标 TopSpin Manual Multiple Peak Analysis）。每个非季碳 C-H 对一行：编号/C_label/H_label/δH/δC/碳级数/multiplet_pattern (s/d/t/q/m)/子峰数/J1_Hz/J2_Hz/J3_Hz/多重峰中心_ppm/半高宽_Hz/强度均匀性 (0-1)/质量等级 (Excellent/Good/Poor/Very Poor/N/A)/说明。质量等级基于子峰强度剖面与理论 1:1/1:2:1/1:3:3:1 的吻合度，用于化学家重点复核 "Very Poor" 行。 |
| `step2_cosy_spin_systems.csv` | COSY 沿耦合链连接的自旋系统（Union-Find 聚类）。 |
| `step3_hmbc_connections.csv` | HMBC 远程 H→C 连接（排除 1J，键类型标记 2J/3J）。完整不截断。 |
| `step4_quaternary_carbons.csv` | 推断的季碳。**推断依据**列标明 `DEPT判定为Cq`（优先级 1）或 `13C有峰但无HSQC对应`（优先级 2）。 |

**HTML 报告 2.5 节**：仅提供一张步骤路径索引表（步骤 1 / 1b / 2 / 3 / 4 + Top-1），
列出每步的具体文件路径（如 `<code>step1_hsqc_ch_pairs.csv</code>` 与目录
`<code>step1_hsqc_ch_pairs/overlap/</code>`）。**不嵌入详细表格与 base64 等高线图**，
详细原始数据在同名独立文件中供化学家单独复核。Top-N 完整归属表
位于 `<code>完整NMR归属表(论文Table1)/</code>` 子目录，每个 Top 一个 CSV。
这将 HTML 总大小从 ~8 MB 减到 ~0.8 MB（以 clxh100 为例）。

**关键设计决策**：
- DEPT 推断优先级 > HSQC 匹配：避免 13C 与 HSQC 微小位移偏差 (< 0.2 ppm) 导致季碳被误判
- 重叠峰聚类容差 0.3 ppm（与 h_tol 一致）
- 13C 切片宽度 ±0.15 ppm（约为典型 1H 线宽 3 倍）
- 13C 峰最小间距 2.0 ppm（避免重复检峰）
- 分解后 H 个数 = 总积分 / 13C 切片峰数（假设各峰贡献相等）

**JSON/MD 报告新增章节**：

## 典型工作流示例

### 示例 1：基础检索（最少命令）

```bash
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3"
# 默认 top-n=10, c-tol=2.0, h-tol=0.3; 数据库从 --db 默认目录自动发现
```

### 示例 2：分子式硬约束 + 1H 谱可视化

```bash
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --formula "C8H8O2" --top-n 5 --with-real-spectrum
```

### 示例 3：完整流水线（检索 + 1r 可视化 + FID 交叉验证 + 反向校验）

```bash
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --top-n 10 --max-hits 200 \
    --with-real-spectrum --validate-with-nmrpy --use-nmrsim
```

### 示例 4：快速验证（小样本 + 跳过校验）

```bash
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --max-hits 30 --top-n 5 --quiet
# 仅 30 秒内可跑完
```

### 示例 5：ML 升级路径（先用主流程出未匹配峰，再用 NMR2Struct）

```bash
# 步骤 1: 检索 + 未匹配峰标注
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3"
# 看 output/clxh3/clxh3_unmatched.csv

# 步骤 2: NMR2Struct 端到端反推（需 GPU + 已下载权重）
python -m nmr_elucidator.run_bridge --folder ... --out-prefix ... --stage 1
python -m nmr_elucidator.run_bridge --folder ... --out-prefix ... --stage 2
python -m nmr_elucidator.run_bridge --folder ... --out-prefix ... --stage 3
```

### 示例 6：CPU 路径结构推断（HOSE 验证 + 置信度评估，无 GPU）

```bash
# 检索 + HOSE 预测验证 + 结构推断置信度（无需 GPU）
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --top-n 10 \
    --use-hose \
    --determine
# 查看 HTML 报告顶部 "结构推断结论" 章节
# 例: "结论: determined (确定) | Top-1: C9H8O4 (阿司匹林) | 融合分 0.899"
```

### 示例 7：库外新骨架检索（PubChem + 片段组装）

```bash
# 当 nmrshiftdb2 命中不足时 (如仅 2-3 个命中)，启用 PubChem 扩展 + 片段组装
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --formula "C9H8O4" \
    --top-n 10 \
    --use-pubchem \
    --use-fragment-assembly \
    --determine
# JSON 报告里 meta.pubchem_candidates[] + meta.fragment_candidates[] 含补充候选
```

### 示例 8：完整流水线（全部验证 + 推断结论）

```bash
# 最完整的验证 + 推断路径（检索 + 1r 可视化 + FID 交叉 + 反向校验 + HOSE + PubChem + 推断结论）
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --top-n 10 --max-hits 200 \
    --with-real-spectrum --validate-with-nmrpy --use-nmrsim \
    --use-hose --use-pubchem --determine
```

### 示例 9：2D NMR 关联验证（有 COSY/HSQC/HMBC 数据时）

```bash
# 启用 2D NMR 关联验证 (自动检测 COSY/HSQC/HMBC/ROESY)
python -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh100" \
    --top-n 10 \
    --determine
# 输出报告中包含 HSQC/HMBC/COSY 匹配率
# 无 2D 数据时 (clxh3) 自动回退到 1D 流程
```

## clxh3 实测样本（参考产物）

`output/clxh3/` 下 4 个产物文件（`clxh3.{json,md,html}` + `clxh3_unmatched.csv`）。关键观察：

- **DEPT 推断**：6 CH + 2 CH2 + 1 CH3 + 2 CH2/CH3 + **1 Cq** + 0 ?（详见 `clxh3.html` 报告）
- **1r 真实谱**：exp 24 (DEPT-135) 的 CH2 负峰在真实谱里清晰可见（向下），peaklist.xml 竹签完全看不出
- **peaklist.xml 过度识别**：clxh3 是低浓度样品，peaklist 含 268 个峰，但真实谱只显 3 个真实大峰（exp 22: 13C），其余都在基线 — 印证 20% 强度过滤的必要性

## 已实现要点

- ✅ Bruker 解析（peaklist.xml + parm.txt + integrals.txt + DEPT 谱识别）
- ✅ 自定义 nmrshiftdb2 SDF 解析器（绕过 RDKit 中文路径 bug）
- ✅ 流式检索（nmrshiftdb2 58K + lucy-ng 928K, 30 秒内完成）
- ✅ DEPT 推断（CH/CH2/CH3/Cq 四态 + 占位 CH2/CH3）
- ✅ 候选打分（coverage + RMS penalty + 峰数差异 penalty）
- ✅ 证据链（matched/unmatched/conflict 三态 + DEPT 软约束）
- ✅ **未匹配峰专项高亮**（重点关注 → 留待 ML 模型 + 化学家人工复核）
- ✅ 结构图（Top-N 网格 + 单分子 + 实测/候选对比 + 真实谱线形）
- ✅ nmrpy FID 交叉验证（独立算法路径）
- ✅ nmrsim 反向校验（cosine 相似度，可选）
- ✅ JSON + Markdown + HTML 三种报告
- ✅ CSV 导出（Excel 兼容）
- ✅ 1r 真实谱重采样到 ML 网格（ChefNMR/NMR2Struct 共用对接点）
- ✅ **scipy.cKDTree 化学位移索引**（`match_index.py`，单条候选内部 O(m log m) 加速，失败回退纯 Python 暴力）
- ✅ **scipy.signal.find_peaks 二次找峰**（`io_bruker.enhance_peaks_with_scipy`，从 1r 实谱补漏 peaklist.xml，依赖未装时静默退化）
- ✅ **pandas 候选表**（`candidate_table.py`，DataFrame 表示 + InChIKey 去重 + Markdown/CSV 导出）
- ✅ **SymPy 分子式校验**（`formula_solver.py`，不饱和度/H-C 比/分子量，辅助候选结构判定）
- ✅ **seaborn 统计可视化**（`viz.py`，评分热力图 + ppm 对比散点 + 约束匹配矩阵，自动嵌入 HTML 第 8 节）
- ✅ **HOSE code 13C 预测验证** (`hose_generator.py` + `nmr_predict.py`); **双数据源**: lucy-ng 928K SQLite (默认自动查找) + 本地 nmrshiftdb2 50K 表 (回退)
- ✅ **lucy-ng 完整 HOSE 数据库集成** (`lucy_ng_hose.py`, ~2.8 GB SQLite, 928K 化合物 + 7.9M HOSE 统计, MAE ~3.5 ppm; 自动查找 `E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db`, 不存在时回退本地表)
- ✅ **PubChem 扩展检索**（`retrieve_pubchem.py`，按分子式 + 13C 推断子结构检索，处理网络失败/速率限制）
- ✅ **片段组装候选生成**（`shift_ranges.py` + `fragment_assembler.py`，从 13C 位移范围 + DEPT 约束 + 分子式生成库外结构，控制组合爆炸）
- ✅ **结构推断置信度评估**（`determination.py`，多源证据融合打分，四级结论：确定/很可能/不确定/无匹配）
- ✅ **2D NMR 关联验证**（`io_bruker_2d.py` + `correlation_2d.py`，自动检测 COSY/HSQC/HMBC/ROESY，MAD-SNR 自适应检峰，候选结构 2D 匹配率打分，无 2D 数据时自动回退）
- ✅ **CPU 路径完整闭环**：从 1D 谱 + 已知分子式可输出推断结论 + 证据链，无需 GPU / ML 权重
- ✅ **与 MestReNova 功能对齐**（`doc/MestReNova_解谱指南.md`，本项目增强后覆盖 Database / NMRPredict / Verify / 部分 CASE 能力）

## 已知限制与边界

| 限制 | 说明 | 缓解 |
|---|---|---|
| 候选库规模 | nmrshiftdb2 58K + lucy-ng 928K 化合物（默认 lucy-ng 优先, HOSE-l4-l6 覆盖宽) | 需同时安装两个数据库; `util/download_lucy_hose_db.py` 下载辅助 |
| 分子式无推断 | 当前需用户提供（可选） | TODO: 接 HRMS 峰自动推分子式 |
| 不处理耦合常数 | 仅用化学位移窗口匹配，未用 J 值 | `nmrsim_check` 模块已预留接口 |
| clxh3 实测峰含大量噪声 | peaklist 268 个峰 → 20% 强度过滤去噪 | `--c-tol`/`--h-tol` 控制匹配容差 |
| nmrpy peakpicker 默认参数 | 在低浓度样品上检出 0 峰 | 需手工调 voff/threshold |
| nmrpy ppm 轴 | 粗略近似（`freq/sfo1+100`） | 精确 ppm 轴用 `--with-real-spectrum` 路径 |
| SDF 中文路径 | RDKit `SDMolSupplier` C++ 层 OSError | 自定义 bytes 解析器绕开（见 `io_sdf.py`） |
| ChefNMR / NMR2Struct ML 模型 | 需 GPU + Zenodo checkpoints | 见主仓库 `从NMR谱图反推分子结构代码库.md` |
| HOSE 查找表精度 | **lucy-ng 928K (MAE ~3.5 ppm, 默认)** / 本地 nmrshiftdb2 50K (MAE ~5 ppm, 回退) | lucy-ng 已默认集成; 如未下载则自动回退本地表; 最大 radius = 6 |
| 片段组装仅 1D 谱 | 仅有 13C + DEPT + 1H 积分时，结构可能不唯一（异构体区分需要 2D） | 报告里明确标注 "需 2D NMR 进一步确认" |
| PubChem 检索速率限制 | API 5 req/s，超额会报错 | 检索函数内部含 try/except 静默降级，主流程不阻断 |

## 结构推断能力边界（与 MestReNova CASE 对比）

本项目增强后的能力对应到 MestReNova 的高级模块（详见 `doc/MestReNova_解谱指南.md`）：

| 能力维度 | 本项目 | MestReNova |
|---|---|---|
| 数据库检索 | nmrshiftdb2 (58K) + PubChem (需网络) | DB (多数据库，含商用) |
| 结构→谱图预测 | HOSE code 13C; lucy-ng 928K (MAE ~3.5 ppm, 默认) / 本地 nmrshiftdb2 50K (MAE ~5 ppm, 回退) + nmrsim 1H | NMRPredict（HOSE + ML，MAE ~3 ppm） |
| 谱图→结构生成 | 片段组装（仅 1D,最多 500 个候选,涵盖几乎所有常见骨架） | Structure Elucidation CASE（含 2D,完整 CASE 算法） |
| 候选打分 | 多源融合（6 维度 + HOSE + lucy-ng 928K） | Verify（独立打分引擎） |
| HOSE 查找表覆盖 | lucy-ng 928K (7.9M HOSE 统计) | 内置闭源库 (估计百万级) |
| 置信度评估 | determined / probable / ambiguous / no_match | 内置评分报告 |
| 需 2D NMR 才能唯一确定 | ✅ 明确报告该情况 | 严重依赖 2D NMR |
| 无 GPU 可用 | ✅ 完整 CPU 路径 | ✅ |

### 什么情况下结论是 determined？

需同时满足：
1. Top-1 综合分 > 0.85
2. Top-1 与 Top-2 差距 > 0.15
3. 至少 2 个验证维度高分（13C coverage ≥ 0.80，HOSE RMS ≤ 5 ppm，DEPT 完全匹配）

### 什么情况下应该报告 ambiguous？

- 仅 1D 谱（无 HSQC / HMBC）且分子式存在多个异构体
- nmrshiftdb2 + PubChem 均无命中，片段组装生成 > 50 个候选
- HOSE 预测对 Top-N 候选分数接近（差距 < 0.08）

### 下一步提升路径（按优先级）

1. **优先级 1**：获取 2D NMR（HSQC + HMBC），从 1H-13C 直接关联 + 长程偶联，能区分 90% 异构体
2. **优先级 2**：接 HRMS（高分辨质谱）精确分子式，减少候选搜索空间 100 倍以上
3. **优先级 3**：接 J 值耦合常数提取（HSQC/HMBC/COSY），接入 `nmrsim_check` 现有接口提高打分维度
4. **优先级 4**：容器化部署 (Docker + MCP 服务公开接口，让外部 CLM/Drug-Agent 调用)

## 模块索引

| 模块 | 职责 | 关键 API |
|---|---|---|
| `io_bruker.py` | 读 Bruker 原始谱（peaklist.xml / parm.txt / integrals.txt） | `read_bruker_folder`, `BrukerSpectrum`, `enhance_peaks_with_scipy` |
| `io_processed.py` | 读 Bruker 已处理谱（pdata/1/1r）+ 重采样到 ML 网格 | `read_processed_folder`, `resample_to_grid` |
| `io_sdf.py` | 自定义 nmrshiftdb2 SDF 解析器 | `stream_nmrshiftdb2_sdf`, `NmrShiftDb2Entry` |
| `constraints.py` | 硬约束（分子式、溶剂、重原子数） | `Constraints`, `formula_match` |
| `retrieve.py` | 流式窗口匹配检索（贪心 + DEPT 软约束，cKDTree fast path） | `retrieve_candidates`, `RetrievalHit` |
| `evidence.py` | 打分 + 证据链生成（matched/unmatched/conflict，含 DoU/H-C 校验） | `score_candidates`, `build_evidence` |
| `report.py` | JSON + Markdown + HTML 报告（含 pandas 候选表 + seaborn 统计图） | `write_report_json`, `write_report_markdown`, `write_report_html` |
| `drawing.py` | 结构图 + 谱图对比 + 真实谱可视化 | `draw_topn_grid`, `draw_unmatched_overlay`, `draw_real_spectrum_overlay` |
| `match_index.py` | scipy.cKDTree 化学位移索引（retrieve 内部加速） | `build_shift_index`, `ShiftIndex` |
| `candidate_table.py` | pandas 候选表 + InChIKey 去重 + Markdown/CSV 导出 | `evidences_to_dataframe`, `dedupe_by_inchikey`, `to_markdown`, `to_csv` |
| `formula_solver.py` | SymPy 分子式解析 + 元素平衡 + DoU + 分子量 | `parse_formula`, `formula_balance`, `formula_dou`, `formula_molecular_weight`, `validate_formula_match` |
| `viz.py` | seaborn 候选评分热力图 + ppm 对比散点 + 约束匹配矩阵 | `plot_candidate_score_heatmap`, `plot_ppm_comparison_scatter`, `plot_constraint_match_matrix` |
| `nmrsim_check.py` | nmrsim 反向校验（cosine 相似度） | `check_candidate_with_nmrsim` |
| `nmrpy_validate.py` | nmrpy FID 二次处理 + 交叉验证 | `validate_bruker_with_nmrpy` |
| `bridge_nmr2struct.py` | Bruker → NMR2Struct 输入格式转换 + DEPT 推断 | `convert_bruker_folder_to_nmr2struct` |
| `run_bridge.py` | NMR2Struct 数据转换 + 推理编排 | `main` (stage 1/2/3) |
| `run_nmr2struct.py` | NMR2Struct 模型推理 | `run_nmr2struct_inference` |
| `parse_nmr2struct_output.py` | NMR2Struct h5 输出解析 + 写 Markdown | `parse_predictions_h5`, `write_candidates_markdown` |
| `hose_generator.py` | HOSE-like 原子环境编码生成器（纯 RDKit，1~6 层） | `generate_hose_code`, `generate_hose_codes_for_carbons`, `prepare_mol` |
| `nmr_predict.py` | 13C 化学位移预测 + 候选验证（从 nmrshiftdb2 SDF 自建 HOSE 查找表） | `get_predictor`, `verify_candidate`, `HOSELookupTable`, `C13Predictor` |
| `lucy_ng_hose.py` | lucy-ng 完整 HOSE 数据库集成（928K 化合物，~2.8 GB SQLite，MAE ~3.5 ppm） | `LucyNGHOSEProvider`, `predict_smiles`, `verify_candidate_with_lucy_db`, `get_lucy_ng_provider` |
| `determination.py` | 结构推断置信度评估（多源证据融合 + 四级结论） | `assess_determination`, `compute_unified_score`, `DeterminationResult` |
| `io_bruker_2d.py` | 2D Bruker 谱图读取与峰检测 (COSY/HSQC/HMBC/ROESY) | `read_bruker_2d_folder`, `read_bruker_2d_experiment`, `pick_peaks_2d`, `Spectrum2D`, `Peak2D` |
| `interpret_2d.py` | 2D NMR 逐步解谱引擎 (化学家工作流: HSQC→COSY→HMBC→季碳推断+5.9 节重叠峰分解，每步产出独立 CSV + 等高线图供化学家复核) | `interpret_2d_sequential`, `write_step_files`, `CHPair`, `SpinSystem`, `HMBCConnection`, `OverlapRegion`, `Interpretation2D` |
| `correlation_2d.py` | 2D NMR 关联提取、候选验证打分、全归属表生成 | `score_candidate_2d`, `build_full_assignment_table`, `extract_hsqc_pairs`, `Correlation2DScore`, `FullAssignmentTable` |
| `retrieve_pubchem.py` | PubChem 扩展检索（分子式 + 子结构，错误静默降级） | `retrieve_pubchem_candidates`, `PubChemCandidate` |
| `shift_ranges.py` | 13C 化学位移经验范围表 + 官能团推断 | `infer_functional_groups`, `assign_shift_ranges`, `SHIFT_RANGES` |
| `fragment_assembler.py` | 片段组装候选生成器（库外新骨架，从 1D 约束出发） | `build_constraints`, `assemble_candidates`, `FragmentConstraint` |
| `log_setup.py` | 项目统一日志初始化（Tee stdout/stderr + logging，多入口共用） | `setup_logging`, `get_logger`, `get_current_log_file`, `auto_setup` |
| `standard_assignment.py` | **已知化合物流程**：读取化学家 Excel 标准归属表 + 自动发现 + 行级联合匹配 + SMILES 结构判等（SMILES 仅复核用） | `read_standard_assignment_table`, `discover_known_compound`, `compare_program_to_standard`, `write_matched_xlsx`, `write_comparison_report`, `StandardAssignmentTable` |
| `atom_labeler.py` | SMILES → RDKit mol → 化学家标号映射（仅在人工复核阶段用于结构验证，不参与推导） | `assign_chemist_labels`, `parse_smiles_to_mol`, `format_h_cell` |
| `__main__.py` | CLI 入口（端到端主流程） | `main` |

## 未匹配峰处理（重点关注）

按设计目标：**与已知结构匹配不上的氢谱和碳谱峰要重点记录，可能是新结构**。

报告中的处理：
1. **JSON 报告**：`CandidateEvidence.unmatched_carbon_shifts` / `unmatched_proton_shifts` 字段
2. **Markdown 报告**：开篇“重点关注”章节，列出 Top-N 候选都未能解释的实测化学位移
3. **HTML 报告**：同 Markdown，含交互 + overlay 对比图（base64 内嵌）
4. **CSV 报告**：`<sample>_unmatched.csv` 每行一个候选，列 `rank/smiles/unmatched_*.shifts + note`

后续可对接 ML 模型（chefnmr/NMR2Struct/lucy-ng）输入：`SMILES` + `unmatched_*.shifts`。

## 不在项目范围 (Part 2.G)

补全 8 项核心解谱步骤后, 以下能力**明确不在 CPU-only 本项目范围内**, 需要外部工具补充:

| 不支持项 | 需要补足 | 外部路径 |
|---|---|---|
| **绝对构型** | Mosher 酯法 / CD 圆二色谱 / X-ray 单晶 | external 工具 (Mnova CD) |
| **变温 NMR / rotamers** | 变温实验 + 动态核磁分析 | 温度序列 (180-340K) + 谱图叠合 |
| **混合 2D (HSQC-TOCSY/HSQC-NOESY/LR-HSQMBC)** | 高级 2D 实验 | Bruker 脉冲程序 + 后续解析 |
| **²J(CH) / ³J(CH) 直接读取** | J-resolved / J-HMBC 实验 | peaklist.xml 不含此字段 |
| **qNMR 定量** | ERETIC / 校准曲线 | nmrpy 已支持 ERETIC, 需单独接入 |
| **完整 CASE 自动解谱 (无分子式已知)** | 2D 约束驱动的结构生成 | `external/lucy-ng/` 或 `external/LSDos-3.5.3/` |
| **ML 谱图预测** | HSQC + SMILES 端到端 | `external/NMR2Struct/` (GPU) / `external/chefnmr/` (GPU) |
| **蛋白-配体相互作用 / HSQC 滴定** | 配体浓度梯度 HSQC | 自由衰减 / 化学位移扰动分析 |

### 已实现: Part 2 改进矩阵

| Part | 能力 | 模块 | 文档 |
|---|---|---|---|
| 2.A | HMBC ¹J 残留过滤 (与 HSQC 比对, `c_tol=1.0 ppm` / `h_tol=0.05 ppm`) | `extract_hmbc_pairs(hsqc_pairs=...)` | `meta['hmbc_stats']` |
| 2.B | 溶剂残余峰过滤 (CDCl3/DMSO/MeOD/H2O/C6D6/CD2Cl2/Pyr/TFA/THF/MeCN, 默认启用, `--keep-solvent` 关闭) | `solvent_filter.py` | `step1*.csv` 新列 `is_solvent_residual` |
| 2.C | **Karplus 3D J 二面角验证** (距离 + 二面角双维度, 文献经典常数) | `verify_candidate_3d_j` | `step6b_j_karplus_3d_verification.csv` |
| 2.E | NOE + J 验证回喂 `score_candidates` (最终打分时序修正) | `__main__.py` 步骤 [9.6b] | `meta['score_breakdown']` |
| 2.F | J 一致性降级维度 (Top-1 `j_match_rate<0.5` 时 determined → probable) | `assess_determination(j_match_rates=...)` | `evidence_summary[]` |
| 2.G | 本节文档归档 | `AGENTS.md` / 本文档 / `MestReNova_解谱指南.md` 5.10 节 | 同左 |
