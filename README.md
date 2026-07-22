# NMR 反推结构工具 (`nmr_elucidator`)

基于 **nmrshiftdb2** 数据库 + **nmrglue Bruker 解析** + **RDKit 化学信息学** 的 NMR 反推结构工具。
读取 1H/13C 谱图 → 在 58187 个库分子中检索 → 输出 Top-N 候选 + 证据链 + 未匹配峰标注。
可选附加：**fid 交叉验证**（独立算法路径）、**1r 真实谱可视化**、**nmrsim 反向校验**、**NMR2Struct 端到端 ML**。

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
[3] stream_nmrshiftdb2_sdf → retrieve_candidates (窗口匹配 + DEPT 软约束)
        ↓
[4] score_candidates (coverage + RMS penalty) + build_evidence
        ↓
[5] (可选) nmrsim 反向校验 cosine
        ↓
output/<sample>/<sample>.{json,md,html}
            <sample>_topn.png / _rank*.png / _overlay.png / _unmatched.csv
            <sample>_real_spectrum.png        (--with-real-spectrum)
            <sample>_nmrpy_validation.{json,md,png}  (--validate-with-nmrpy)
```

详见 `从NMR谱图反推分子结构代码库.md`（调研 + 升级路径），本文档只讲**如何使用**。

## 安装

```bash
pip install nmrglue nmrsim -i https://pypi.tuna.tsinghua.edu.cn/simple
# nmrpy 可选（FID 交叉验证才需要）:
pip install nmrpy -i https://pypi.tuna.tsinghua.edu.cn/simple
# 预装依赖: rdkit / pubchempy / numpy / scipy / matplotlib / pandas
# 以下为该项目增强: sympy (公式校验) / seaborn (统计可视化)
```

## 数据准备

下载 **nmrshiftdb2** 数据文件（仅 1 个，约 159 MB）：

| 文件 | 大小 | 链接 | 用途 |
|---|---|---|---|
| `nmrshiftdb2withsignals.sd` | 158.6 MB | [SourceForge](https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd/download) | 结构 + 13C/1H 化学位移（反推检索必需） |

放置：`E:\NMR2Structure\data\nmrshiftdb2\nmrshiftdb2withsignals.sd`
不需要下载：无信号的 `nmrshiftdb2.sd`、1.1 GB 的 SQL dump、`predictor*.jar`（Java 工具用不上）。

## 快速开始

```bash
cd E:/NMR2Structure/NMR2Structure
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --top-n 10
```

输出落到 `E:\NMR2Structure\output\clxh3\`（默认输出根目录 + 从 folder 路径末段推断化合物名）。
化合物名显式指定用 `--compound-name my_compound`，输出根目录改 `--output-root D:/my_output`。

## 主流程 CLI

```
python -X utf8 -m nmr_elucidator --folder <原始数据/<sample>> --sdf <SDF>
                              [选项...]
```

### 必填

| 参数 | 说明 |
|---|---|
| `--folder` | Bruker 化合物根目录（含数字子目录 21/22/...，每个下有 `pdata/1/peaklist.xml`）|
| `--sdf` | nmrshiftdb2 SDF 文件路径 |

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
| `--with-real-spectrum` | 读取 `pdata/1/1r` 真实处理后谱，绘制真实谱线形图（叠加 peaklist 峰位） | `<prefix>_real_spectrum.png` |

### 其他

| 参数 | 说明 |
|---|---|
| `--quiet` | 不打印检索进度（每 5000 条打印一行） |

## 子模块（独立运行）

主流程会自动调用，但也可单独用：

### `nmrpy_validate` — FID 二次处理 + 交叉验证

```bash
python -X utf8 -m nmr_elucidator.nmrpy_validate \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmrpy_check"
```

输出 `<prefix>_validation.{json,md,png}` 三件套。**当前局限**：nmrpy peakpicker 默认参数在低浓度样品（如 clxh3）上检出 0 峰，需手工调 voff/threshold；ppm 轴是粗略近似（`freq/sfo1 + 100`）。

### `run_bridge` — NMR2Struct 端到端 ML

需要 `external/NMR2Struct` 目录已下载权重（见主仓库 `从NMR谱图反推分子结构代码库.md`）。

```bash
# stage 1: 数据转换 (无需 torch)
python -X utf8 -m nmr_elucidator.run_bridge \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmr2struct" --stage 1

# stage 2: NMR2Struct 推理 (需 torch + 已下载权重)
python -X utf8 -m nmr_elucidator.run_bridge \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --out-prefix "E:/NMR2Structure/output/clxh3/nmr2struct" --stage 2

# stage 3: h5 -> Markdown (无需 torch)
python -X utf8 -m nmr_elucidator.run_bridge \
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

## 输出文件清单

每个化合物目录（如 `E:\NMR2Structure\output\clxh3\`）下：

### 始终产出（主流程）

| 文件 | 说明 |
|---|---|
| `<sample>.json` | 完整证据链（matched/unmatched/conflict），程序可读 |
| `<sample>.md` | Markdown 报告，含未匹配峰专项章节 + pandas 候选表（InChIKey 去重）|
| `<sample>.html` | HTML 报告（人类阅读，含流程 + 结构图 + DEPT 推断 + seaborn 统计可视化）|
| `<sample>_topn.png` | Top-N 候选 2D 结构网格图 |
| `<sample>_rank<N>.png` | 每个候选的单分子图（N=1..top-n） |
| `<sample>_overlay.png` | 实测谱 vs Top-5 候选预测谱对比（红色=实测，散点=候选） |
| `<sample>_unmatched.csv` | 未匹配峰汇总（Excel 友好，UTF-8 BOM） |
| `<sample>_heatmap.png` | 候选评分热力图（seaborn，行=候选，列=维度，归一化 [0, 1]）|
| `<sample>_scatter.png` | 实测 vs 候选 ppm 散点对比（seaborn，13C + 1H 双子图）|
| `<sample>_matrix.png` | 约束匹配矩阵（seaborn，6 维度满足度热力图）|

### `--with-real-spectrum` 附加

| 文件 | 说明 |
|---|---|
| `<sample>_real_spectrum.png` | 4 子图：每个实验的真实谱曲线（蓝）+ peaklist 峰位（红） |

### `--validate-with-nmrpy` 附加

| 文件 | 说明 |
|---|---|
| `<sample>_nmrpy_validation.json` | FID 交叉验证数据（各实验峰数 + 2ppm 重叠率） |
| `<sample>_nmrpy_validation.md` | 人类阅读报告 |
| `<sample>_nmrpy_validation.png` | 4 子图：每个 FID 的 nmrpy 处理结果 |

## 典型工作流示例

### 示例 1：基础检索（最少命令）

```bash
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd"
# 默认 top-n=10, c-tol=2.0, h-tol=0.3
```

### 示例 2：分子式硬约束 + 1H 谱可视化

```bash
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --formula "C8H8O2" --top-n 5 --with-real-spectrum
```

### 示例 3：完整流水线（检索 + 1r 可视化 + FID 交叉验证 + 反向校验）

```bash
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --top-n 10 --max-hits 200 \
    --with-real-spectrum --validate-with-nmrpy --use-nmrsim
```

### 示例 4：快速验证（小样本 + 跳过校验）

```bash
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --max-hits 30 --top-n 5 --quiet
# 仅 30 秒内可跑完
```

### 示例 5：ML 升级路径（先用主流程出未匹配峰，再用 NMR2Struct）

```bash
# 步骤 1: 检索 + 未匹配峰标注
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd"
# 看 output/clxh3/clxh3_unmatched.csv

# 步骤 2: NMR2Struct 端到端反推（需 GPU + 已下载权重）
python -X utf8 -m nmr_elucidator.run_bridge --folder ... --out-prefix ... --stage 1
python -X utf8 -m nmr_elucidator.run_bridge --folder ... --out-prefix ... --stage 2
python -X utf8 -m nmr_elucidator.run_bridge --folder ... --out-prefix ... --stage 3
```

## clxh3 实测样本（参考产物）

`output/clxh3/` 下 27 个产物文件（~30 MB）。关键观察：

- **DEPT 推断**：6 CH + 2 CH2 + 1 CH3 + 2 CH2/CH3 + **1 Cq** + 0 ?（详见 `clxh3.html` 报告）
- **1r 真实谱**：exp 24 (DEPT-135) 的 CH2 负峰在真实谱里清晰可见（向下），peaklist.xml 竹签完全看不出
- **peaklist.xml 过度识别**：clxh3 是低浓度样品，peaklist 含 268 个峰，但真实谱只显 3 个真实大峰（exp 22: 13C），其余都在基线 — 印证 20% 强度过滤的必要性

## 已实现要点

- ✅ Bruker 解析（peaklist.xml + parm.txt + integrals.txt + DEPT 谱识别）
- ✅ 自定义 nmrshiftdb2 SDF 解析器（绕过 RDKit 中文路径 bug）
- ✅ 流式检索（58187 条记录，30 秒内完成）
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

## 已知限制与边界

| 限制 | 说明 | 缓解 |
|---|---|---|
| 候选库规模 | nmrshiftdb2 58187 条（不含 COCONUT 等） | `lucy-ng` 可扩到 928K（待集成） |
| 分子式无推断 | 当前需用户提供（可选） | TODO: 接 HRMS 峰自动推分子式 |
| 不处理耦合常数 | 仅用化学位移窗口匹配，未用 J 值 | `nmrsim_check` 模块已预留接口 |
| clxh3 实测峰含大量噪声 | peaklist 268 个峰 → 20% 强度过滤去噪 | `--c-tol`/`--h-tol` 控制匹配容差 |
| nmrpy peakpicker 默认参数 | 在低浓度样品上检出 0 峰 | 需手工调 voff/threshold |
| nmrpy ppm 轴 | 粗略近似（`freq/sfo1+100`） | 精确 ppm 轴用 `--with-real-spectrum` 路径 |
| SDF 中文路径 | RDKit `SDMolSupplier` C++ 层 OSError | 自定义 bytes 解析器绕开（见 `io_sdf.py`） |
| ChefNMR / NMR2Struct ML 模型 | 需 GPU + Zenodo checkpoints | 见主仓库 `从NMR谱图反推分子结构代码库.md` |

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
| `__main__.py` | CLI 入口（端到端主流程） | `main` |

## 未匹配峰处理（重点关注）

按设计目标：**与已知结构匹配不上的氢谱和碳谱峰要重点记录，可能是新结构**。

报告中的处理：
1. **JSON 报告**：`CandidateEvidence.unmatched_carbon_shifts` / `unmatched_proton_shifts` 字段
2. **Markdown 报告**：开篇"重点关注"章节，列出 Top-N 候选都未能解释的实测化学位移
3. **HTML 报告**：同 Markdown，含交互
4. **CSV 报告**：`<sample>_unmatched.csv` 每行一个候选，列 `rank/smiles/unmatched_*.shifts + note`
5. **结构图**：`<sample>_overlay.png` 红色竖线标记实测峰，散点显示候选预测峰

后续可对接 ML 模型（chefnmr/NMR2Struct/lucy-ng）输入：`SMILES` + `unmatched_*.shifts`。