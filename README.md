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
[7] (可选) PubChem 扩展检索 (按分子式 / 13C 推断子结构)
        ↓
[8] (可选) 片段组装候选生成 (库外新骨架, DEPT + 化学位移范围驱动)
        ↓
[9] (可选) 结构推断置信度评估 (assess_determination: determined/probable/ambiguous/no_match)
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

### 辅助工具脚本（util/）

| 脚本 | 用途 |
|---|---|
| `util/download_lucy_hose_db.py` | 下载 lucy-ng HOSE 完整数据库 (~830 MB / ~2.8 GB)；支持 Colab / aria2 / curl / proxy 4 种策略 |

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
| `--use-hose` | HOSE code 13C 预测验证；**双数据源策略**：默认自动查找 `E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db` (928K 化合物, MAE ~3.5 ppm)；不存在时自动从 nmrshiftdb2 SDF 自建轻量本地表 (默认 50K, MAE ~5 ppm) | JSON/MD 报告里追加 `hose_verifications[]` + `meta.used_hose_source` (`lucy-ng-derep.db` 或 `local-sdf`) |
| `--use-pubchem` | PubChem 扩展检索：当 nmrshiftdb2 命中不足时，按分子式 / 13C 推断子结构去 PubChem 查候选 | JSON/MD 报告里追加 `pubchem_candidates[]` 列表 |
| `--use-fragment-assembly` | 片段组装候选生成：当数据库检索均无命中时，从 13C 位移范围 + DEPT 约束组装结构（仅 1D 谱可走） | JSON 报告里追加 `fragment_candidates[]` |
| `--use-2d` | 2D NMR 关联验证：自动检测 COSY/HSQC/HMBC/ROESY 实验，提取交叉峰关联，验证候选结构；无 2D 数据时自动回退到 1D 流程 | JSON/MD/HTML 报告里追加 `correlation_2d` 章节 |
| `--determine` | 结构推断置信度评估：多源证据融合打分，输出 determined / probable / ambiguous / no_match 四级结论 | JSON/MD/HTML 报告顶部追加 `determination` 章节 |
| `--hose-cache <path>` | HOSE 查找表缓存路径（默认 `<sdf 同目录>/_hose_lookup_cache.tsv`，避免重复建表） | — |
| `--hose-max-records <N>` | 从 SDF 构建 HOSE 查找表的最大记录数（默认 **50000**；0=全部；10000 轻度，50000 够用，100000+ 需约 5 分钟） | — |
| `--lucy-ng-db <path>` | lucy-ng 完整 HOSE 数据库路径（`lucy-ng-derep.db`, ~2.8 GB SQLite）; **默认已自动查找** `E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db`, 无需手动指定; 用此选项可指向备份位置 | — |
| `--log-dir <path>` | 日志输出目录 (默认 `E:\NMR2Structure\logs`); 每次运行自动生成 `NMR2Structure_YYYY-MM-DD_HH-MM-SS.log` 文件, 包含屏幕上的所有内容 (print + logging) | — |
| `--no-log` | 禁用日志记录 (默认每次运行都写日志) | — |

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

# 打开 lucy-ng SQLite 数据库 (默认查找 E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db)
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

### `--use-hose` 附加

不新增独立文件，**丰富 JSON/MD/HTML 报告**：
- `meta.hose_verifications[]`：每个候选的 `{smiles, hose_rms, hose_coverage, hose_match_rate, verification_score, notes}`
- `meta.used_hose: true` 标记启用
- Top-1 候选的 HOSE RMS 嵌入证据链

### HOSE 查找表进阶：lucy-ng 928K (默认) vs nmrshiftdb2 本地表 (回退)

**默认路径**：自动使用 lucy-ng 928K SQLite (~2.8 GB, MAE ~3.5 ppm):

```bash
python -X utf8 -m nmr_elucidator --folder ... --sdf ... --use-hose
# 默认自动查找 E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db
```

**回退路径**：如果 lucy-ng 数据库不可用, 从本项目 `nmrshiftdb2 SDF` 自建轻量本地表 (默认 50K 条, 约 1-2 分钟, MAE ~5 ppm):

```bash
python -X utf8 -m nmr_elucidator --folder ... --sdf ... --use-hose --hose-max-records 50000
```

**升级路径**：下载 lucy-ng 完整 HOSE 库（928K 化合物，~830 MB 下载，~2.8 GB 解压，MAE ~3.5 ppm）：

```bash
# 1. 下载（4 种策略可选）
python util/download_lucy_hose_db.py --strategy colab    # 生成 Google Colab Notebook
python util/download_lucy_hose_db.py --strategy aria2    # 生成 aria2 多线程脚本

# 2. 解压到 E:\NMR2Structure\data\lucy-ng-derep_db\
#    pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps

# 3. 使用完整库（精度提升明显）
python -X utf8 -m nmr_elucidator --folder ... --sdf ... --use-hose
# 默认自动查找 E:\NMR2Structure\data\lucy-ng-derep_db\lucy-ng-derep.db
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

### 示例 6：CPU 路径结构推断（HOSE 验证 + 置信度评估，无 GPU）

```bash
# 检索 + HOSE 预测验证 + 结构推断置信度（无需 GPU）
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --top-n 10 \
    --use-hose \
    --determine
# 查看 HTML 报告顶部 "结构推断结论" 章节
# 例: "结论: determined (确定) | Top-1: C9H8O4 (阿司匹林) | 融合分 0.899"
```

### 示例 7：库外新骨架检索（PubChem + 片段组装）

```bash
# 当 nmrshiftdb2 命中不足时 (如仅 2-3 个命中)，启用 PubChem 扩展 + 片段组装
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
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
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh3" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --top-n 10 --max-hits 200 \
    --with-real-spectrum --validate-with-nmrpy --use-nmrsim \
    --use-hose --use-pubchem --determine
```

### 示例 9：2D NMR 关联验证（有 COSY/HSQC/HMBC 数据时）

```bash
# 启用 2D NMR 关联验证 (自动检测 COSY/HSQC/HMBC/ROESY)
python -X utf8 -m nmr_elucidator \
    --folder "E:/NMR2Structure/NMR/clxh100" \
    --sdf "E:/NMR2Structure/data/nmrshiftdb2/nmrshiftdb2withsignals.sd" \
    --top-n 10 \
    --use-2d \
    --determine
# 输出报告中包含 HSQC/HMBC/COSY 匹配率
# 无 2D 数据时 (clxh3) 自动回退到 1D 流程
```

## clxh3 实测样本（参考产物）

`output/clxh3/` 下 27 个产物文件（~30 MB）。关键观察：

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
| `correlation_2d.py` | 2D NMR 关联提取与候选验证打分 | `score_candidate_2d`, `extract_hsqc_pairs`, `extract_hmbc_pairs`, `extract_cosy_pairs`, `Correlation2DScore` |
| `retrieve_pubchem.py` | PubChem 扩展检索（分子式 + 子结构，错误静默降级） | `retrieve_pubchem_candidates`, `PubChemCandidate` |
| `shift_ranges.py` | 13C 化学位移经验范围表 + 官能团推断 | `infer_functional_groups`, `assign_shift_ranges`, `SHIFT_RANGES` |
| `fragment_assembler.py` | 片段组装候选生成器（库外新骨架，从 1D 约束出发） | `build_constraints`, `assemble_candidates`, `FragmentConstraint` |
| `log_setup.py` | 项目统一日志初始化（Tee stdout/stderr + logging，多入口共用） | `setup_logging`, `get_logger`, `get_current_log_file`, `auto_setup` |
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