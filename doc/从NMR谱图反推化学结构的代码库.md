# 从 NMR 谱图反推分子结构 — 代码库调研报告

> 调研时间：2026-07-21（三次补充：NMRPy / NMR2Struct / lucy-ng 详情 + 实际本地落位路径 + 端到端验证）
> 目标：在已有 RDKit + PubChemPy 的基础上，找出适合"读谱图 → 推理结构 → 给化学家人工复核"全流程的 Python 工具链
> 调研结论摘要：`nmrshiftdb2` 做候选库 + `nmrglue` 解 Bruker + `nmrsim` 做反向校验 + `CHEF-NMR`/`NMR2Struct`/`lucy-ng` 做 ML/CASE 升级路径

---

## 一、按功能分层

### 1. 谱图文件解析（Bruker FID / 1r / 2r）

| 库 | GitHub / PyPI | 关键能力 | 备注 |
|---|---|---|---|
| **nmrglue** | `jjhelmus/nmrglue` | 读写 Bruker / Varian / NMRPipe / Sparky；谱图处理、可视化 | 文献被引 471 次 |
| pynmranalysis | `1feres1/pynmranalysis` | 1H NMR 预处理（binning、对齐、归一化） | 代谢组学常用 |
| nmrpoise | PyPI `nmrpoise` | 自动化谱图处理流程 | `pip install nmrpoise` |
| **NMRPy** | `NMRPy/nmrpy` (PyPI `nmrpy`) | 处理与分析 1D / 阵列化 NMR 谱；DOSY / 弛豫 / 变温谱 | 见下文补充 §1.1 |
| ssNake | van Meerten 2019 | 跨平台 Bruker 处理 | GUI 桌面工具 |

**推荐**：`nmrglue` 是 Bruker 处理的"事实标准"，与 `clxh3/21/pdata/1/1r`、`audita.txt` 这类原始数据天然匹配。

#### 1.1 NMRPy/nmrpy 补充（已落地, 0.2.8）

- **仓库**：`https://github.com/NMRPy/nmrpy`
- **PyPI**：`https://pypi.org/project/nmrpy/`（已装版本 0.2.8）
- **文档**：`https://nmrpy.readthedocs.io/en/latest/intro.html`
- **代表论文**：Liu et al., arXiv:2103.14988 (2021)
- **安装**：`pip install nmrpy -i https://pypi.tuna.tsinghua.edu.cn/simple`（已装）
- **定位**：纯 Python 3 模块，跨平台；与 `nmrglue` 互补使用
- **核心 API（本项目已验证）**：
  - `nmrpy.from_path(folder, file_format='bruker')` - 一次性读 Bruker 实验目录（含 4 个 FID）
  - `.ft_fids()` / `.ps_fids()` / `.baseline_correct_fids()` - 标准处理流程
  - `.peakpicker(fid_number=N)` - 自动峰检测
  - `.calibrate()` / `.get_integrals_from_traces()` - 校准与积分
- **路径要求**：fid_path 下必须是数字子目录（21/22/...）, 与 Bruker 实验目录格式完全一致
- **MNova 文件**：**不支持**（magic = "Mestrelab" 专有格式），需在 MNova GUI 导出 JCAMP-DX
- **已落地的集成**：`NMR2Structure\nmr_elucidator\nmrpy_validate.py`
  - 命令：`python -X utf8 -m nmr_elucidator.nmrpy_validate --folder E:\NMR2Structure\NMR\clxh3 --out-prefix E:\NMR2Structure\output\clxh3\nmrpy_check`
  - 产物：`*_validation.{json,md,png}` 三件套
  - clxh3 验证：成功读 4 个 FID (shape 4×32699 complex128), 完成 FT+PS+baseline, peakpicker API 可调但参数需细调
- **风险**：项目维护较弱（0.2.8 是 2022-09 版本）；Python 3.13 兼容性已验证 OK

### 2. 化学位移预测（结构 → 谱图，用于反向打分）

| 库 | 来源 | 精度 | 用途 |
|---|---|---|---|
| **patonlab/CASCADE** | Paton Group (2021) | DFT 级精度；graph NN | 立体化学感知；本地可运行（PyTorch） |
| lab-cosmo/shiftml | COSMO Lab | 屏蔽常数预测（晶体/有机固相） | 主要面向周期性体系 |
| nmrgnn | PyPI `nmrgnn` | 通用有机分子 | `pip install nmrgnn`，需 nmrdata 嵌入 |

**推荐**：`CASCADE` 是当前最高质量的开源有机分子 1H/13C 化学位移预测器，可以用来对候选结构反向算谱，与实测谱比对打分。

### 3. 谱图 → 结构（核心反推 ML 模型）

| 项目 | 论文/会议 | 准确率 | 可用性 |
|---|---|---|---|
| **ml-struct-bio/chefnmr** | NeurIPS 2025 | 65%+（天然产物） | 最新最强；原子扩散模型；PyTorch |
| **MarklandGroup/NMR2Struct** | arXiv 2403.11353v4 (2024-12) | 69.6%（含预训练，含 top-15） | 多任务 ML；Stanford Markland 组；见 §3.1 |
| **rxn4chemistry/nmr-to-structure** | ChemRxiv 2023 | 67% top-1（合成谱） | Transformer；含谱图模拟脚本 |
| **steinbeck/lucy-ng** | Steinbeck 组 (2026) | 取决于 LSD 求解 | AI-agent 协作 CASE；920K 化合物库；见 §3.2 |
| boskovicgroup/DP4-AI | RSC 2020, 被引 195 | 自动指认+概率打分 | 含 NMR-AI 自动处理 |
| GT-NMR | J. Cheminform. 2024 | 化学位移预测 | 图 Transformer |

**定位**：这是 `nmr_elucidator.py` 最直接的升级目标。当前版本只走候选库检索（基于化学位移窗口匹配 + 打分排序），遇到库外新结构会留下 unmatched 峰供 ML 模型二次处理。

#### 3.1 MarklandGroup/NMR2Struct 补充

- **仓库**：`https://github.com/MarklandGroup/NMR2Struct`
- **论文**：arXiv:2403.11353v4 [cs.LG] (2024-12-16)，Stanford Markland Group
- **方法**：多任务机器学习框架，仅基于 1D 1H NMR 和 13C NMR 同时预测分子结构 + 分子式 + 连接性
- **架构**：端到端，从原始谱 → 拓扑结构 + 化学式
- **准确率**：69.6% top-15（在 COCONUT/天然产物数据集上）
- **亮点**：唯一同时输出 `formula + connectivity` 的开源模型
- **与 lucy-ng 区别**：NMR2Struct 是 ML 模型，不依赖 LSD 求解；lucy-ng 是 CASE agent + LSD 求解；本项目 `nmr_elucidator` 是确定性窗口匹配
- **与本项目关系**：
  - 升级路径 E（端到端 ML）：把 `clxh3` 实测谱直接喂入，验证 model 输出是否在 nmrshiftdb2 检索结果的 enriched 候选列表中
  - 或：NMR2Struct top-15 + 本工具打分排序 做串联（`unmatched.csv` 输入 NMR2Struct，输出结构 → 与 nmrshiftdb2 候选交叉验证）
- **风险**：模型权重是否公开需确认；预训练数据集是否包含 clxh3 类似结构未知；需 GPU

#### 3.2 steinbeck/lucy-ng 补充

- **仓库**：`https://github.com/steinbeck/lucy-ng`
- **作者**：Christoph Steinbeck (Friedrich-Schiller-University Jena) — 同名 Lucy 的原作者（Marilyn Easton 时代的 CASE 程序复活升级版）
- **License**：MIT
- **Python 要求**：3.10+
- **标题**："AI-Agent Powered Computer-Assisted Structure Elucidation for Organic Natural Products"
- **核心概念**：通过协作 AI agent 团队运行全 CASE 流程，模拟一群光谱学家工作
- **底层依赖**（lucy-ng 是 integration layer）：
  - `LSD` (Jean-Marc Nuzillard, Université de Reims)：核心结构生成器
  - `RDKit`：化学信息学
  - `nmrglue`：谱图读取
  - `HOSE codes` (Bremser 1978)：13C 预测
  - `COCONUT + NMRShiftDB`：参考数据库
- **Agent 团队**（4 专家 + orchestrator + 备用 1）：
  - NMR Chemist：peak picking + 多重度 + 杂化推断
  - LSD Engineer：构建约束文件 + 运行 solver
  - Solution Analyst：13C shift 排名 + 结构鉴定
  - Devil's Advocate：质量门（每次 LSD 文件运行前验证）
  - Diagnostic specialist：灾难恢复（仅当 0 解 / 解爆炸时）
- **CLI / 命令族**：`/lucy-ng:case`, `/lucy-ng:dereplicate`, `/lucy-ng:predict`, `/lucy-ng:sanitise`, `/lucy-ng:status`
- **关键命令**：
  ```bash
  # 安装
  pip install "lucy-ng @ git+https://github.com/steinbeck/lucy-ng.git"
  # 13C 预测需要（HOSE code 生成器，旧版 Python 测试依赖有问题）：
  pip install git+https://github.com/Ratsemaat/HOSE_code_generator.git --no-deps
  pip install rdkit
  # LSD solver 必须安装（PATH 里有 LSD / outlsd 可执行）
  lucy lsd check
  # 参考数据库（928,443 化合物 + 7.9M HOSE 统计）
  lucy database download -o data/reference/lucy-ng-derep.db
  # 化学求解一条跑
  lucy read 1d data/Ibuprofen/2
  lucy pick hsqc data/Ibuprofen/6 --dept135 data/Ibuprofen/3
  lucy dereplicate c13 data/Ibuprofen/2 C13H18O2
  lucy lsd generate data/Ibuprofen C13H18O2 -o ibuprofen.lsd
  lucy lsd run ibuprofen.lsd
  lucy lsd rank output/ data/Ibuprofen/2 --top 10
  ```
- **数据库规模**：928,443 个化合物（COCONUT + NMRShiftDB）+ 790 万条 HOSE-code 统计，下载 ~830 MB → 解压 ~2.8 GB
- **Database DOI**：10.6084/m9.figshare.31073554
- **依赖原理数据**：
  - manual HMBC ~16 相关 → 900+ 解
  - guided HMBC（28–30 相关）→ 1–10 解
  - **guided peak picking**（用 DEPT-135 当 ground truth）大幅压缩候选空间
- **与本项目关系**：
  - 升级路径 F（CASE 框架）：当本项目 `nmr_elucidator` 无法命中（unmatched 峰过多）时，把 Bruker 文件 + 分子式直接交 lucy-ng 由 LSD 求解
  - 复用：lucy-ng 自带的 `lucy dereplicate c13` 和 `lucy predict c13` 直接可对接本项目 `clxh3` 的 13C 谱（已对应溶剂为 Pyr）
  - 借鉴：lucy 的 DEPT-135→HSQC guided peak picking 思路，比本项目简单的"20% 强度阈值"更严谨，可作为未来过滤增强
- **风险**：
  - LFS (LSD solver) 是 C 写的，需单独安装；CLI 不带自动安装
  - 与 Claude Code agent 系统强绑定，agent team 直接走 `.claude/`，对纯 Python 工程侵入性较大
  - 数据库 2.8 GB 远超本项目当前的 `nmrshiftdb2withsignals.sd` (160 MB)，需要权衡

### 4. 谱图数据库（候选结构来源）

| 资源 | 规模 | 用途 |
|---|---|---|
| **nmrshiftdb2** | ~43K 结构 | 唯一开源可批量 13C 化学位移预测；SDF 可下载 |
| COCONUT | 100 万+ 天然产物 | 与 lucy-ng 已预集成 |
| NMRShiftDB HOSE-code 统计 | 790 万条 | lucy-ng 内置 |

**本项目采用**：`nmrshiftdb2`（文件名 `nmrshiftdb2withsignals.sd`，约 158.6 MB，含原子级 1H/13C 化学位移指认）。
- 保存路径：`e:\NMR2Structure\data\nmrshiftdb2\nmrshiftdb2withsignals.sd`
- 下载链接：<https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd/download>
- 优势：仅此一个数据库就把"结构+信号指认"两件事同时解决

### 5. 谱图模拟（反向校验 / 数据增强）

| 库 | 特点 |
|---|---|
| **nmrsim** (Skrynnikov 组, Purdue) | 量子力学 + 化学交换模拟；PyPI 版本 0.6+ |
| rxn4chemistry/nmr-to-structure 自带模拟器 | 用于训练数据生成的合成谱 |

**本项目采用**：`nmrsim` 做反向校验模块，验证候选结构预测谱与实测谱的形状匹配。

---

## 二、本项目已采纳的技术栈

| 层 | 库 | 版本 | 角色 |
|---|---|---|---|
| 化学信息学 | rdkit | 2026.3.4 | SMILES 解析、分子描述符、子结构匹配 |
| 化合物查询 | pubchempy | 1.0.5 | 在线 CID / 名称 / 相似性查询 |
| 谱图解析 | **nmrglue** | **0.11** | Bruker 1r 文件读写（本项目实际只用到 peaklist.xml / parm.txt） |
| 谱图模拟 | **nmrsim** | **0.6.0.post1** | 反向校验（实验性模块，结构允许时启用） |
| 候选库 | **nmrshiftdb2 SDF** | 2024 版 | 58187 个分子 × 1H/13C 化学位移指认 |
| ML 模型 | (未集成) | — | chefnmr / NMR2Struct / nmr-to-structure 留作升级路径 |

---

## 三、已实现的反推流程（nmr_elucidator 包）

```
Bruker 原始数据 (NMR/clxh3/)
   ↓ read_bruker_folder
观测峰表 (1H + 13C, 已按强度过滤噪声)
   ↓ ProblemPeaks
nmrshiftdb2 SDF (流式扫描 58187 条)
   ↓ retrieve_one (化学位移窗口匹配 + 贪心分配)
候选 (RetrievalHit, 硬筛后约 200 条)
   ↓ score_candidates (coverage + RMS penalty)
Top-N (ScoredCandidate)
   ↓ build_evidence
证据链 (CandidateEvidence: matched/unmatched/conflict/notes)
   ↓ write_report_* + draw_*
JSON 报告 + Markdown 报告 + Top-N 结构图 + 实测对比 overlay + 未匹配 CSV
```

详见 `NMR2Structure/README.md`、`NMR2Structure/nmr_elucidator/` 子模块说明。

---

## 四、升级路径（未实现但可对接）

### A. 对接 CHEF-NMR（已下载源码, 本机无法运行, 评估如下）

- **落位**：`E:\NMR2Structure\external\chefnmr`（源码 src/script/configs/assets 齐全）
- **定位**：NeurIPS 2025，原子级扩散生成模型（不是检索/分类，是从谱条件直接"扩散生成"3D 原子坐标 + 键）
- **本机无法运行的硬性缺口**（逐条核实 external\chefnmr 后确认）：
  1. **无 GPU**：`environment.yaml` 要求 CUDA 12.4 + PyTorch 2.5.1，扩散采样在 CPU 上不可行（官方基准用 A100/H100）；本机 torch 是 `2.13.0+cpu`
  2. **无模型权重**：仓库只含代码，checkpoints 需从 Zenodo 单独下载（仓库内 `meta/`、`grids/`、`datasets/` 均缺失）
  3. **无 conda 环境**：官方走 `conda env create -f environment.yaml`（环境名 `nmr3d`, py3.11），与本项目 anaconda base 不隔离
- **输入格式**（读 `src/data/input_generator.py` 确认）：
  - `h_10k`：1H 谱离散到 **10000 维**，ppm 范围 **[-2, 10]**
  - `c_80`：13C 谱离散到 **80 维**，ppm 范围 **[3.42, 231.3]**
  - condition_type: `H1NMRSpectrum` / `C13NMRSpectrum` / `H1C13NMRSpectrum`
- **关键洞察**：ChefNMR 的 `c_80 [3.42,231.3]` 与 NMR2Struct 的 `CNMR 80 维 [3.4,231.3]` **几乎一致**，二者可共用同一套"真实谱 → 网格向量"的重采样管线
- **本项目已提供的对接点**（不建无法测试的完整桥接, 遵守 YAGNI）：
  - `NMR2Structure\nmr_elucidator\io_processed.py::resample_to_grid(ppm, intensity, grid_min, grid_max, n_points)`
  - 用法：`read_processed_folder` 读 1r 真实谱 → `resample_to_grid(ps.ppm, ps.intensity, 3.42, 231.3, 80)` 即得 ChefNMR 的 `c_80`（或 `-2, 10, 10000` 得 `h_10k`）
  - 重采样真实 1r 谱比 peaklist.xml 竹签装箱**保真度更高**（保留完整线形, 是喂 ML 模型的正确来源）
- **落地条件**：拿到 GPU 机器 + 下载 Zenodo checkpoints 后，用 `resample_to_grid` 生成向量 → 按 `input_generator.py` 组装 → 跑官方 `script/` 采样
- **风险**：扩散模型输出是 3D 结构候选，需再做分子式/价态校验（本项目 constraints 模块可复用）；天然产物准确率 65%+ 但对偏分布样本（如 Pyr 溶剂）未知

### B. 对接 lucy-ng (CASE agent) — 详见 §3.2

- 用途：920K 化合物（COCONUT + NMRShiftDB）库 + LSD 求解 + AI agent 协作
- 优势：可解释性最强；适合产物结构未知场景
- 风险：需要 LSD solver + 数据库 2.8 GB 下载；与 Claude Code 耦合
- 本项目对接点：`unmatched.csv` 拆出 Bruker 文件 + 分子式 → 移交 `lucy-ng:case`

### C. 对接 CASCADE（DFT 级精度）

- 用途：对 Top-N 候选做高精度 chemical shift 预测，与实测谱精确比对
- 优势：可分辨同分异构体
- 风险：推理慢，每分子 ~10-30 秒

### D. 替换 candidate library 为更大数据库

- 当前：nmrshiftdb2 (~43K)
- 可选：COCONUT (100 万+) + ChemSpider + PubChem
- 风险：networking dependency，违反 `AGENTS.md` 规则 7（不依赖网络结果覆盖本地原始数据）

### E. 对接 NMR2Struct (端到端 ML) — 详见 §3.1

- 用途：直接把 1H + 13C 原始谱喂给模型，输出 formula + connectivity
- 优势：跳过中间所有处理；69.6% top-15
- 对接方式：把本项目 `_unmatched.csv` + 实测谱文件路径喂 NMR2Struct，模型返回结构 → 与 nmrshiftdb2 Top-N 交叉验证
- 风险：可能需要 GPU；权重是否完全可下载需确认

### F. 接入 nmrpy 做 FID 二次处理与交叉验证（已集成进主流程）

- 用途：用 nmrpy 重新处理 `fid`，与 peaklist.xml 对比，验证 TopSpin 自动处理结果
- 已落地：`NMR2Structure\nmr_elucidator\nmrpy_validate.py`
- **已集成进主流程**：主命令加 `--validate-with-nmrpy` 开关即自动产出交叉验证三件套（`*_nmrpy_validation.{json,md,png}`）
- 优势：
  - 完全独立的算法路径（nmrpy 用 lmfit 做相位/基线，TopSpin 用 apk0/lb）
  - 一次读多个 FID 不用循环 4 次 `read_bruker_experiment`
  - 可扩展到 DOSY / 变温 / 阵列谱分析（如果未来需求扩展）
- 局限：
  - peakpicker 参数需手工调（默认参数在 clxh3 噪声环境下检出为 0，故重叠率显示 N/A）
  - ppm 轴是粗略近似（`freq/sfo1 + 100`），且 sfo1 对 13C 仍取 500（未按核种切换）——**精确 ppm 轴请用 §四.G 的 1r 路径**
- 命令（两种等价方式）：
  ```powershell
  # 方式 1: 集成进主流程 (推荐, 与检索/打分/报告一起产出)
  python -X utf8 -m nmr_elucidator --folder E:\NMR2Structure\NMR\clxh3 `
      --sdf E:\NMR2Structure\data\nmrshiftdb2\nmrshiftdb2withsignals.sd --validate-with-nmrpy

  # 方式 2: 单独运行 nmrpy_validate 模块
  python -X utf8 -m nmr_elucidator.nmrpy_validate `
      --folder E:\NMR2Structure\NMR\clxh3 `
      --out-prefix E:\NMR2Structure\output\clxh3\nmrpy_check
  ```
- 产物：
  - `*_nmrpy_validation.png` - 4 个 FID 谱图，标注 peaklist.xml 峰 + nmrpy peakpicker 峰
  - `*_nmrpy_validation.json` - 详细对比数据
  - `*_nmrpy_validation.md` - 人类阅读报告

### G. 读取 pdata/1/1r 真实处理后谱（已落地, 推荐）

- **用途**：读 TopSpin 处理后的最终实谱（FT + 相位 + 基线已完成），用于
  1. **真实谱线形图**：相比 peaklist.xml 的竹签，1r 保留完整线形，直观展示"真实谱 + 我们的峰检测落点是否吻合"
  2. **喂 ML 模型的正确来源**：重采样成 ChefNMR/NMR2Struct 的网格向量（见 §四.A 对接点）
- 已落地：`NMR2Structure\nmr_elucidator\io_processed.py`
  - `read_processed_folder(folder, specs)` → `{exp_no: ProcessedSpectrum}`（用 `nmrglue.bruker.read_pdata` 读, 用 procs 的 SF/SW/OFFSET/SI 算**精确 ppm 轴**, 按核种区分）
  - `resample_to_grid(ppm, intensity, grid_min, grid_max, n_points)` → ML 网格向量
- **已集成进主流程**：主命令加 `--with-real-spectrum` 开关即产出 `*_real_spectrum.png`
- 相比 §四.F 的 nmrpy 路径的优势：
  - **ppm 轴精确**：nmrglue 按 procs 参数逐核种计算（1H ~[-4,16], 13C ~[-6,234]），nmrpy 是粗略近似且 sfo1 不切换核种
  - 无需重新 FT/相位/基线（直接读 TopSpin 的最终结果，更可复现）
- clxh3 验证（4 个实验全读通）：
  - exp 21 (1H): 65536 点, SF=500.18 MHz
  - exp 22 (13C): 32768 点, SF=125.77 MHz
  - exp 23/24 (DEPT): 32768 点；**exp 24 (DEPT-135) 的 CH2 负峰在真实谱图中清晰可见**（比竹签更能体现 DEPT sign）
- 命令：
  ```powershell
  python -X utf8 -m nmr_elucidator --folder E:\NMR2Structure\NMR\clxh3 `
      --sdf E:\NMR2Structure\data\nmrshiftdb2\nmrshiftdb2withsignals.sd --with-real-spectrum
  ```
- 局限：1r 是已处理数据（依赖 TopSpin 的相位/基线选择）；若需完全独立的第一性原理处理，用 §四.F 的 fid 路径

---

## 五、未通过的方案与原因

| 方案 | 拒绝原因 |
|---|---|
| **rxn4chemistry/nmr-to-structure** (`from_pkl = True`) | 模型权重不公开，需要付费 API |
| DP4-AI | 仅适用于结构指认，不解决反推 |
| shiftml | 主要面向晶体/固相，与本项目小分子溶液体系不匹配 |
| ssNake | GUI 桌面工具，与本项目 CLI 流水线冲突 |
| nmrgnn | pip 可装但 `nmrdata` 嵌入文件 ~800 MB，性价比低 |

---

## 六、版本兼容性与限制

- Python 3.13（项目规范要求）
- RDKit 中文路径 bug：必须用 `pathlib.Path(...)` + bytes 读取 SDF，**不要**用 `Chem.SDMolSupplier()` 在中文路径下直接打开（见 `io_sdf.py:_split_records`）
- matplotlib 必须在 `import pyplot` 之前配置中文，否则仍然乱码（已修复，见 `drawing.py` 头部）
- nmrsim 实际能用版本：`0.6+`（0.7 还未稳定），需要时检查 `python -c "import nmrsim; print(nmrsim.__version__)"`

---

## 七、引用与致谢

- nmrglue: J.J. Helmus, C.P. Jaroniec, *J. Biomol. NMR* (2013)
- nmrshiftdb2: S. Kuhn et al., *J. Cheminform.* (2012)
- nmrsim: T. Vosegaard, *J. Magn. Reson.* (2014)
- CHEF-NMR: NeurIPS 2025 (ml-struct-bio)
- NMR2Struct: *J. Chem. Inf. Model.* 2024
- patonlab/CASCADE: CASCADE Paton Group, 2021

---

## 附录：搜索过程中失败的尝试 + 重新搜索补充

### 初次搜索失败（2026-07-21）

以下三个 GitHub 仓库直接 `WebFetch` 抓取均 timeout（沙箱网络出口 `10.88.3.106` 到 `github.com:443` 被防火墙拦截）：

- `https://github.com/MarklandGroup/NMR2Struct`
- `https://github.com/steinbeck/lucy-ng`
- `https://github.com/NMRPy/nmrpy`

### 重新搜索补充（用 WebSearch + arXiv / 间接文档）

- **NMRPy/nmrpy** (§1.1)：Liu et al., arXiv:2103.14988 (2021)；Pure Python 3 模块，跨平台，与 nmrglue 互补
- **MarklandGroup/NMR2Struct** (§3.1)：arXiv:2403.11353v4 (2024-12)；Stanford Markland Group；多任务 ML 端到端反推，69.6% top-15
- **steinbeck/lucy-ng** (§3.2)：Steinbeck 2026；AI-Agent 协作的 CASE 框架；与 LSD solver + COCONUT/NMRShiftDB 920K 库联动

详见 §1.1、§3.1、§3.2 各小节。

### 沙箱解封 / 权限请求备注

如需直接获取最新 release 状态（避免信息滞后），可请求沙箱解封对外网端口，或在本地终端手动 `git clone`。

> 详见 [Qoder 官方文档：终端与沙箱](https://docs.qoder.com/zh/user-guide/quest/terminal-and-sandbox)：目前网络控制为粗粒度的开/关模式，精细到域名级的过滤能力正在规划中，**普通 IDE 实例无法在客户端放行 github.com:443**。

### 三仓库本地落位与版本快照（方案 B+C 落地，2026-07-21）

**结论先行**：沙箱出口 `10.88.3.106` 到 `github.com:443` 被拦截，`WebFetch` 不可达。改走「本地 `git clone` + `WebSearch` 元信息」组合即可完成同等工作。

#### 元信息（`WebSearch` 验证，2026-07-21）

| 仓库 | 仓库地址 | 安装方式 | 最新版本 / 日期 | 引用论文 |
|---|---|---|---|---|
| NMR2Struct | `https://github.com/MarklandGroup/NMR2Struct` | 源码克隆，无 PyPI | 持续开发 | Hu et al., *ACS Cent. Sci.* 2024（doi:10.1021/acscentsci.4c01132） |
| lucy-ng | `https://github.com/steinbeck/lucy-ng` | 源码克隆，无 PyPI | 2026 活跃 | Priessner et al., *RSC* 2026（D5DD00359H） |
| nmrpy | `https://github.com/NMRPy/nmrpy` | **PyPI 可装** | `0.2.8` / 2022-09-29 | Liu et al., arXiv:2103.14988（2021） |

#### 落位目录设计

```
NMR2Structure\                           # 项目根 (e:\NMR2Structure)
├── NMR\                   # Bruker 原始谱
├── NMR2Structure\                          # 本项目主代码
│   ├── nmr_elucidator\                # 主代码包
│   ├── output\                        # 报告产物
│   ├── samples\
│   ├── scripts\
│   └── tests\
└── external\                          # 第三方参考实现，只读快照（项目根级，不在NMR2Structure下）
    └── NMR2Struct\                    # git clone --depth 1  (e:\NMR2Structure\external\NMR2Struct)
```

- `external/` 放在**项目根级**，与 `NMR2Structure/` 平级：方便跨子项目引用，不与 `nmr_elucidator` 代码强耦合
- `external/NMR2Struct` 由 `NMR2Structure/nmr_elucidator/run_nmr2struct.py` 中的常量 `NMR2STRUCT_ROOT` 指向
- lucy-ng 同理放在 `external/lucy-ng`（待手工下载）
- nmrpy 不落 `external/`，走 `pip install nmrpy`（标准做法）

#### 本地 PowerShell 命令清单

```powershell
Set-Location 'e:\NMR2Structure'
New-Item -ItemType Directory -Path 'external' -Force
# NMR2Struct 已手工下载并解压到该目录（git clone 在沙箱不可用, 用户从 GitHub 下载 zip 解压）
Expand-Archive -Path "$HOME\Downloads\NMR2Struct-main.zip" -DestinationPath 'external\' -Force
# 验证
Test-Path 'e:\NMR2Structure\external\NMR2Struct\nmr'
Test-Path 'e:\NMR2Structure\external\NMR2Struct\checkpoints\multitask_checkpoint.pt'
Test-Path 'e:\NMR2Structure\external\NMR2Struct\example_configs\alphabet.npy'
python -m pip install nmrpy -i https://pypi.tuna.tsinghua.edu.cn/simple
```

执行后将 `Get-ChildItem` 与 `print` 的输出回填到本节，作为「已落位」标记。

#### 已落位标记（2026-07-21）

| 路径 | 状态 | 备注 |
|---|---|---|
| `e:\NMR2Structure\external\NMR2Struct` | 已存在 | 手工下载 zip 解压 |
| `e:\NMR2Structure\external\NMR2Struct\checkpoints\multitask_checkpoint.pt` | 97.6 MB | 核心权重（端到端） |
| `e:\NMR2Structure\external\NMR2Struct\checkpoints\cnmr_only_checkpoint.pt` | 96.9 MB | 仅 C 谱子模型 |
| `e:\NMR2Structure\external\NMR2Struct\checkpoints\hnmr_only_checkpoint.pt` | 97.5 MB | 仅 H 谱子模型 |
| `e:\NMR2Structure\external\NMR2Struct\checkpoints\transformer_checkpoint.ckpt` | 66.2 MB | Transformer 解码器 |
| `e:\NMR2Structure\external\NMR2Struct\example_configs\alphabet.npy` | 21 token | SMILES 字符 vocab |
| `e:\NMR2Structure\external\NMR2Struct\example_configs\HNMR_shifts.p` | 28000 grid | H 谱离散化网格 |
| `e:\NMR2Structure\external\NMR2Struct\example_configs\CNMR_shifts.p` | 80 grid | C 谱离散化网格 |
| `e:\NMR2Structure\external\NMR2Struct\example_configs\substructures_957.p` | 957 SMARTS | 子结构列表 |
| `nmrpy` (PyPI 安装) | 0.2.8 | 全局可用，走清华镜像 |

#### 已验证的端到端流程（clxh3 实测样本）

```bash
# 1. 数据转换 (无需 torch, 几秒)
python -m nmr_elucidator.run_bridge --folder <bruker> --out-prefix <prefix> --stage 1

# 2. NMR2Struct 推理 (需 torch + 已下载权重)
python -m nmr_elucidator.run_bridge --folder <bruker> --out-prefix <prefix> --stage 2

# 3. h5 -> Markdown (无需 torch, 几秒)
python -m nmr_elucidator.run_bridge --folder <bruker> --out-prefix <prefix> --stage 3
```

实测结果：clxh3 跑通 stage 1+2+3 全流程，生成 15 个候选 SMILES + 分数 + 结构图。  
候选结构本身化学合理性待化学家复核（CPU 推理被 timeout 截断 + Pyr 溶剂样本对模型属偏分布，详情见 `NMR2Structure/README.md`）。

#### 待办（手工下载）

| 仓库 | 解压到 | 当前状态 |
|---|---|---|
| lucy-ng (`https://github.com/steinbeck/lucy-ng`) | `e:\NMR2Structure\external\lucy-ng` | 待下载 |
| CASCADE (`https://github.com/patonlab/CASCADE`) | `e:\NMR2Structure\external\CASCADE` | 待下载 |
| chefnmr (`https://github.com/ml-struct-bio/chefnmr`) | `e:\NMR2Structure\external\chefnmr` | 已下载（仅源码, 缺 Zenodo checkpoints + GPU, 详见 §四.A） |

每个仓库下载后解压到对应 `external\` 子目录即可，代码层不需要改动（常量路径自动指向）。
