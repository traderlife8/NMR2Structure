# MestReNova 解谱指南

> 本文档面向本项目（`NMR2Structure`）解谱与结构推断工程背景，整理 MestReNova（Mnova）在 NMR 数据处理和化学结构解析中的标准工作流、关键插件与操作要点，便于团队按"预测-比对-验证"三步法或 CASE 自动化路径落地。
>
> 信息来源：Mestrelab Research 官网功能说明、Reachsoft 中文代理页、相关学术综述与公开教程。最新功能请以官方手册为准。
>
> 文档整理日期：2026-07-22

## 一、软件概览

MestReNova（简称 Mnova）是 Mestrelab Research 推出的多功能实验室数据分析平台，覆盖 NMR、LC/GC/MS、UV/Vis、NIR/MIR、Raman、电子/振动光谱、qNMR 等场景。在 NMR 结构解析领域，Mnova 兼具三类核心能力：

1. **数据处理**：FID 导入、相位校正、基线校正、化学位移定标、积分、标峰、多重峰分析、谱图美化和报告生成。
2. **谱图预测**：由化合物结构预测 1H、13C、杂原子谱、HSQC 等理论谱图，并与实测谱比对归属（`NMRPredict` 插件）。
3. **结构解析（Computer-Assisted Structure Elucidation, CASE）**：从分子式 + 1D/2D NMR 数据从头产生候选结构，并对每个候选打分排序（`Structure Elucidation` 插件）。

辅助插件包括 `Verify`（结构验证打分）、`Database / DB`（结构/谱图检索）、`StereoFitter`（3D 构象/构型）、`Binding`（蛋白-配体 HSQC 滴定）、`SMA`（混合物定量）、`qNMR`（定量核磁）、`RM`（反应监测）、`Screen`（基于片段的筛选）和 `IUPAC Name`（结构命名）。

### 功能地图

| 功能域 | 主要能力 | 代表模块/插件 |
|---|---|---|
| 一维谱处理 | 相位、基线、定标、积分、标峰 | 内核 |
| 多重峰分析 | J 耦合常数、峰形归属 | 内核 + `Multiplet Analysis` |
| 二维谱 | COSY / HSQC / HMBC / NOESY | 内核 |
| 谱图预测 | 由结构预测 NMR | `NMRPredict` |
| 结构验证 | 候选结构与实测谱打分 | `Verify` |
| 库检索 | 谱图/结构数据库 | `Database (DB)` |
| 自动结构生成 | 分子式 + 2D NMR → 候选 | `Structure Elucidation (CASE)` |
| 3D 构象 | NOE 约束的构象/构型分析 | `StereoFitter` |
| 报告 | JSON / PDF / HTML 报告 | 内核 |

---

## 二、安装与许可

### 2.1 系统要求

- Windows 10/11、macOS、Linux 主流发行版均可运行。
- 内存建议 ≥ 8 GB，处理大型 2D 谱建议 ≥ 16 GB。
- 二维谱处理和 `Structure Elucidation` 模块对 CPU 有一定要求，建议多核。

### 2.2 安装与激活

Mnova 是商业软件，需要 license（单机 license 或浮动 license）。典型步骤：

1. 从 Mestrelab 官方下载安装包（当前主版本 14/15/16/17）。
2. 安装后启动 Mnova，会提示输入 Host ID（Windows 下可在 Mnova `Help → License` 中查看，也可用 `getmac` 或 `ipconfig /all` 获取网卡 MAC）。
3. 将 Host ID 发给 Mestrelab 或国内代理（如 Reachsoft）申请 `.lic` 文件。
4. 把 `.lic` 文件放到 Mnova 安装目录下的 `licenses` 子目录，或通过 `Help → Install License` 加载。
5. 重启 Mnova，`Help → About` 确认已激活的插件列表（`Structure Elucidation` / `Verify` / `NMRPredict` 等需单独购买）。

> 不同插件的 license 相互独立；本项目做结构推断时最关键的三个插件是 `NMRPredict`、`Verify` 和 `Structure Elucidation`。

### 2.3 版本历史要点

- Mnova 14：Qt 升级，UI 调整。
- Mnova 15.1：多模态分析平台的多方位升级。
- Mnova 16：强化 Python 集成，支持 ASM 格式。
- Mnova 17：迁移至 Qt 6，UI 焕新，内置 Python 包管理器，支持 JSON 数据格式。

Python 集成对自动化有意义：自 Mnova 16 起，可通过 `mnova` Python 包脚本化执行打开谱图、设置参数、批量导出 JSON 等任务。本项目暂未启用该集成，但可作为后续自动化方向。

---

## 三、数据导入与预处理

### 3.1 数据来源

Mnova 支持 Bruker、Varian/Agilent、JEOL、MestreNova 自有格式以及 JCAMP-DX、ASCII 等通用格式。对本项目对应的 Bruker 数据：

- 直接把含 `fid` 的文件夹拖入 Mnova 窗口，Mnova 自动完成 FID 读取、相位校正、基线校正、傅里叶变换全流程。
- 若文件夹中同时存在 `.mol` 文件，谱图旁会自动显示分子结构。
- 也可走菜单 `File → Open` 选择具体实验。

### 3.2 预处理关键步骤

**3.2.1 参考定标（Reference）**

菜单 `Analysis → Reference → Reference (R)`：

- 选 TMS 或已知溶剂峰；
- 单击谱图上的目标峰，输入其已知化学位移；
- 确认后整张谱按该偏移重新校准。

常用氘代溶剂残余峰参考值（务必显式记录原始值，不静默校正）：

| 溶剂 | ¹H (ppm) | ¹³C (ppm) |
|---|---|---|
| CDCl₃ | 7.26 | 77.16 |
| DMSO-d₆ | 2.50 | 39.52 |
| CD₃OD | 3.31 | 49.00 |
| D₂O | 4.79 | — |
| Acetone-d₆ | 2.05 | 29.84 |
| CD₂Cl₂ | 5.32 | 53.84 |

实际参考值应以本批次 NMR 实验的 `acqus`/`audita.txt` 中 `SOLVENT` 字段和文献标准值共同确认。

**3.2.2 相位校正**

Mnova 默认自动相位校正，必要时手动调整：

- `Processing → Phase Correction`，可零阶（p0）和一阶（p1）独立校正。
- 调整判据：基线平整、所有峰左右对称、TMS 峰两侧基线一致。

**3.2.3 基线校正**

`Processing → Baseline Correction`：

- 默认 Bernstein 多项式拟合；
- 可选其他算法（如 Whittaker、线性）；
- 对宽峰或基线严重倾斜的情况需要分段校正。

### 3.3 批量处理与模板

- 通过 `Processing → Template` 保存当前处理参数为模板，应用于其他谱图保持一致。
- 在谱图上右键 `Apply Processing Template` 可复用模板。
- 模板对长时间序列实验、多通道切换场景特别重要。

---

## 四、一维谱分析（¹H / ¹³C / DEPT）

### 4.1 ¹H 谱工作流

1. **积分**：菜单 `Analysis → Integration`：
   - `Auto`（应急用，但可能纳入杂质/边缘信号）；
   - 推荐手动：放大局部，逐峰拉基线、设积分区间，给出相对积分。
2. **多重峰分析**：菜单 `Analysis → Multiplet Analysis`：
   - `Automatic` 自动得到化学位移、积分、J 耦合常数、峰形归属（s/d/t/q/m 等）；
   - 也可手动指定多重峰中心和范围。
3. **标峰（Peak Picking）**：`Analysis → Peak Picking` 标注所有峰顶位置。
4. **谱图美化**：调整线宽、颜色、坐标轴字号、位移标注字体；将化学位移、积分值标在谱图上。

### 4.2 ¹³C 谱与 DEPT 推断碳级数

¹³C 谱本身不直接显示碳级数，需结合 DEPT：

| DEPT 类型 | CH | CH₂ | CH₃ | Cq |
|---|---|---|---|---|
| DEPT-45 | 正 | 正 | 正 | 消失 |
| DEPT-90 | 正 | 消失 | 消失 | 消失 |
| DEPT-135 | 正 | 负 | 正 | 消失 |

操作：

- 在 Mnova 中打开 ¹³C 和 DEPT 谱，`Stack → Stack Plot` 叠加；
- 逐一比对：正峰出现在 DEPT-135 → CH 或 CH₃；负峰 → CH₂；正峰不在 DEPT-135 → 进一步用 DEPT-90 区分 CH 与 CH₃；谱图中存在但 DEPT 全部消失 → Cq。

### 4.3 谱图归属与参数报告

- `View → Tables → Parameters` 插入实验参数（场强、溶剂、温度、采样次数等）；
- `View → Tables → Peaks` 插入峰位表；
- 一键 `Report` 按钮生成 PDF / HTML / JSON 报告。

---

## 五、二维谱关联（COSY / HSQC / HMBC / NOESY）

### 5.1 通用处理流程

1. `File → Open Multiple Files` 同时打开多个 2D 谱；
2. 对每个 2D 谱执行相位校正（2D 维度独立）、基线校正、参考定标；
3. 调整等高线（contour）水平、缩放、平滑；
4. 设置谱图显示窗口（芳香区、脂肪族区等）；
5. 在谱图上标注关键交叉峰。

### 5.2 各 2D 谱的核心用途

| 实验 | 主要提供信息 | 解谱作用 |
|---|---|---|
| COSY / TOCSY | ¹H–¹H 耦合相关 | 拼出邻碳 H 的耦合链 |
| HSQC / HSQC-DEPT | ¹H–¹³C 直接连接（¹J） | 给每个 H/C 配对 |
| HMBC | ¹H–¹³C 远程相关（²J, ³J） | 连接片段、确认取代位点 |
| NOESY / ROESY | 空间近邻（< 5 Å） | 立体构型、相对构型 |
| INADEQUATE | ¹³C–¹³C 直接连接 | 完整碳骨架 |

### 5.3 Mnova 二维谱增强功能

- 自动峰标注（`Analysis → Peak Picking → 2D`）；
- 多谱图同步缩放 / 联动平移；
- 2D 谱与一维谱联动（点击 2D 峰高亮对应 1D 位置）；
- 谱图数据库检索（把 2D 谱作为查询条件）。

---

## 六、结构解析核心：Structure Elucidation（CASE）

> 该模块需购买 `Structure Elucidation` 插件 license，是本项目目前缺失的能力。下文概述其工作流以作为功能对标参考。

### 6.1 输入要求

- 分子式（由 HRMS 或元素分析得到；Mnova 自动算 DoU）；
- 1H NMR；
- 13C NMR；
- HSQC（必须，用于关联 H/C）；
- HMBC（强烈推荐，用于连接片段）；
- COSY（可选，用于 H-H 邻接）；
- 任何可用的 DEPT-135/90。

### 6.2 工作流

1. `Mnova → Structure Elucidation → New Project`；
2. 输入分子式，选择元素范围（如允许 C/H/N/O/卤素）；
3. 导入各谱图（按格式要求填入通道、溶剂、参考）；
4. 设置 HMBC 距离窗口（默认 ²J + ³J，必要时排除 ¹J 干扰）；
5. 运行 CASE：Mnova 用专家系统 + 谱图一致性搜索生成所有满足 2D 相关性的候选结构；
6. 自动对每个候选进行谱图预测（`NMRPredict`）并与实测谱比对打分（类似 `Verify`）；
7. 候选列表按打分排序，导出归属表。

### 6.3 三种典型场景对照

| 场景 | 推荐工具组合 | 关键步骤 |
|---|---|---|
| 已知化合物，确认结构 | `Database` → `Verify` | 库检索 → 候选 → Verify 打分 |
| 几种可能结构待定 | `NMRPredict` + `Verify` | 对每个候选预测 1H/13C/HSQC → 比对 → 归属 |
| 完全未知结构 | `Structure Elucidation` | 分子式 + 1H/13C/HSQC/HMBC → 自动生成候选 → 排序 |

### 6.4 与传统人工解谱的关系

- CASE 生成的候选集合往往很大（>100 个），需要靠 `Verify` 排序；
- 高质量 HMBC 和 HSQC 是 CASE 收敛速度的关键；
- 互变异构、共振结构、立体异构仍需人工从打分前列候选中复核。

---

## 七、Verify 与候选打分

### 7.1 Verify 输入

- 候选结构（mol/SMILES）；
- 实测 1H 和/或 13C 和/或 HSQC；
- 可选 LC/MS 数据。

### 7.2 打分维度

- 化学位移 RMS（预测 vs 实测）；
- 峰形一致性（多重峰模式 vs 耦合常数预测）；
- HSQC 关联匹配率；
- HMBC 远程相关匹配率；
- LC/MS 分子量、碎片离子一致性。

### 7.3 输出

- 总分 + 分项得分；
- 归属表（每个 H/C 标在结构上）；
- 不一致项告警。

> 本项目 `evidence.score_candidates` 的打分思路与 Verify 类似，但维度更简化：coverage + RMS penalty + peak count diff + DEPT coverage bonus + H 总数一致性。

---

## 八、报告生成与可复现性

### 8.1 报告内容

- 实验参数（场强、溶剂、温度、采样参数）；
- 谱图（1D / 2D 可缩放窗口）；
- 峰位表、化学位移表、积分表、多重峰表；
- 归属表（每个 H/C 标在结构上）；
- 耦合常数 J；
- 候选结构及其打分（如使用 `Verify`）。

### 8.2 输出格式

- PDF（论文 SI 友好）；
- HTML（嵌入 base64 PNG，单文件可移植）；
- JSON（自动化下游使用）；
- CSV（数据表二次分析）。

### 8.3 可复现性建议

- 保存处理模板（`.mnp` 或 `.tpl`）；
- 记录参考定标、相位、基线参数；
- 保留原始 FID 与导出 JSON；
- 对 `Structure Elucidation` 项目，导出全部候选与打分快照。

---

## 九、与本项目（`NMR2Structure`）的能力对照

| 功能 | Mnova | 本项目增强后 | 差异 |
|---|---|---|---|
| Bruker peaklist 读取 | 内核 | `io_bruker.read_bruker_folder` | 本项目只读已处理文件，无相位/基线 |
| 1r 实谱读取 | 内核 | `io_processed.read_processed_folder` | 同上 |
| 多重峰 / J 常数分析 | `Multiplet Analysis` | `bridge_nmr2struct`（仅 DEPT 推断） | 本项目无 J/峰形 |
| 候选库检索 | `Database`（~26 万 13C / ~10 万 1H） | `retrieve.retrieve_candidates` + `retrieve_pubchem.retrieve_pubchem_candidates` | Mnova 闭源，本项目开源可定制；增加 PubChem 拓展 |
| 候选打分 | `Verify` | `evidence.score_candidates` + `determination.compute_unified_score` | 本项目多源融合打分 + 四级置信度 |
| 谱图预测 | `NMRPredict` | `nmr_predict.verify_candidate`（HOSE code 13C 预测） | 精度 MAE ~5 ppm（限 nmrshiftdb2 10K 建表） |
| **从头生成候选** | **`Structure Elucidation (CASE)`** | **`fragment_assembler.assemble_candidates`**（仅 1D，限 200 候选） | 本项目仅做 1D 片段组装，不需 2D |
| **官能团推断** | Mnova 内置 | `shift_ranges.infer_functional_groups` | 本项目覆盖 20+ 常见官能团 |
| **置信度输出** | Mnova 内置评分 | `determination.assess_determination`（determined/probable/ambiguous/no_match） | 本项目四级中文标签 + 证据链 |
| 3D 构象 | `StereoFitter` | 缺失 | 本项目只做 2D 推断 |
| 报告 | PDF/HTML/JSON | `report.write_report_{json,markdown,html}` | 格式对齐 + `meta.determination` 章节 |

> 如果要让本项目进一步增强结构推断能力，可考虑以下外部路径：
> - `external/chefnmr/`：基于 ML 的 1D/2D 谱图推断（需 GPU）；
> - `external/LSDos-3.5.3/`：基于 2D NMR 的传统 CASE 专家系统；
> - `external/NMR2Struct/`：基于深度学习的 HSQC → 结构预测（需 GPU）；
> - `external/lucy-ng/`：基于 HSQC + HMBC 的 CASE 工具，HOSE 表 928K 条可与本项目 `nmr_predict` 模块互换。
>
> 当前 CPU 路径增强后，本项目已覆盖：
> - 库检索 (nmrshiftdb2 + PubChem)
> - 谱图预测 (HOSE + nmrsim)
> - 候选打分 + 置信度
> - 有限从头生成 (片段组装)
> - 报告 + 证据链
> - 无 GPU / 无 ML 权重 即可运行完整闭环

---

## 十、解谱实操：以本项目 `clxh3` 为例

> 假设在装有 Mnova 14+ 的工作站上对 `E:\NMR2Structure\NMR\clxh3\` 中的谱图做手工复核。

### 10.1 数据导入

1. 打开 Mnova，把 `E:\NMR2Structure\NMR\clxh3\21`（1H）、`22`（C-13 或 DEPT-135）、`23`、`24` 四个实验文件夹依次拖入；
2. 窗口内并排显示 1H、13C、DEPT-135、DEPT-90（若有）。

### 10.2 一维谱分析

1. **定标**：检查 `audita.txt` 的 `SOLVENT` 字段，按上节溶剂表选参考峰；记录原始参考偏移量。
2. **积分**：手动对 1H 谱每个清晰可分的峰积分，记录相对积分比。
3. **多重峰分析**：自动 + 手动校正，记录每个峰的 δ、积分、J、峰形（s/d/t/q/m 等）。
4. **DEPT 推断碳级数**：与 13C 叠加，标注每个 13C 峰为 CH / CH₂ / CH₃ / Cq。
5. **导出**：保存为 `clxh3_1H.mnova`、`clxh3_13C.mnova`、`clxh3_dept135.mnova`，并导出 JSON。

### 10.3 二维谱关联（如有）

1. 打开 HSQC，确认每个 1H 峰都有对应 13C，反之亦然；
2. 打开 HMBC，做关键远程相关标注（2–3 键）；
3. 用 COSY 把 H-H 邻接链连起来；
4. 用 NOESY（若存在）确认立体。

### 10.4 候选打分

```powershell
# 用本项目自动检索 + 打分
cd E:\NMR2Structure\NMR2Structure
python -X utf8 -m nmr_elucidator `
    --folder "E:\NMR2Structure\NMR\clxh3" `
    --sdf "E:\NMR2Structure\data\nmrshiftdb2\nmrshiftdb2withsignals.sd" `
    --top-n 10 `
    --with-real-spectrum
```

将本项目打分前列候选与 Mnova `Verify` 结果交叉对比：

- 若 Top 1 重合 → 高置信；
- 若 Top 3 部分重合 → 取并集人工复核；
- 若完全不重合 → 检查打分维度、定标偏差、候选筛选阈值。

### 10.5 结果汇总

- 在 `E:\NMR2Structure\output\clxh3\` 中已有本项目生成的 `{json,md,html}` 报告与结构图；
- Mnova 端导出 `clxh3_verify.html` 与归属表 PNG；
- 把两份结果在 Markdown 中并列对照，标记共同命中与差异候选。

---

## 十一、常见问题与排错

| 现象 | 可能原因 | 处理建议 |
|---|---|---|
| 拖入文件夹无响应 | 文件夹不含 `fid` 或 `acqus` | 确认 Bruker 目录结构完整 |
| 参考定标后峰位偏移 | 选择了错误的溶剂 | 重新对照 `audita.txt` 的 SOLVENT |
| DEPT 谱无信号 | 参数不当或样品浓度低 | 检查 DEPT 参数、增益、增加 NS |
| CASE 候选数过多 | HMBC 窗口太宽或元素范围过大 | 收紧 HMBC 窗口（限定 ²J + ³J）、限定元素 |
| Verify 打分偏低 | 谱图参考定标不准 | 重定标后重跑 Verify |
| 二维谱交叉峰找不到 | 等高线水平太高 | 降低 contour level，逐级确认 |
| 积分偏小 | 自动基线校正把峰截了 | 关自动基线或改手动基线 |

---

## 十二、安全与合规

- Mnova license 与 Host ID 绑定，勿随意变更网卡或虚拟机 MAC；
- 输出报告若用于发表或商业交付，确认 Mnova 版本号、license 范围与第三方插件来源；
- 与外部数据库（如 PubChem、ChemSpider）联用时遵守各库的使用条款与速率限制。

---

## 十三、参考链接

- 官方功能总览：<https://www.reachsoft.com.cn/mnova>
- Mnova 中文代理页：<https://www.reachsoft.com.cn/productinfo/3528227.html>
- Mestrelab 官方功能说明：<https://www.go-soft.cn/news/cjwt/810.html>
- 解谱软件教程：<https://www.go-soft.cn/news/cjwt/1087.html>
- 核磁软件使用宝典：<https://www.molaid.com/baike/tool/569>
- 培训日程（1D/2D + 数据库 + 结构验证）：<https://m.sohu.com/a/867201037_610519>
- Mnova Tools for Structure Elucidation Using 1D and 2D NMR（彭琛讲座）：<https://ibook.antpedia.com/x/652627.html>
- 手把手教你用 Mnova + ChemDraw：<https://m.blog.csdn.net/weixin_29290947/article/details/159134830>
- 5 分钟搞懂核磁共振谱图解析（Chemdraw + Mnova）：<https://wenku.csdn.net/column/r95rv9h6oxp>
- DeepSAT（HSQC → 结构，文中比较 ACD/SE、Mnova SE、LSD）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC10406729/>
- HSQC Spectra Simulation and Matching：<https://pubs.acs.org/doi/10.1021/acs.jcim.3c01735>

---

## 结论

MestReNova 是 NMR 解谱和化学结构推断领域的工业标准工具，覆盖从谱图处理、谱图预测到 CASE 自动结构解析的完整链路。本项目（`NMR2Structure`）目前的 `io_bruker → retrieve → evidence → report` 链路在功能上接近 Mnova 的 `Database` + `Verify` 组合，但在 `Structure Elucidation`（从头生成候选）和 `NMRPredict`（预测谱图）两块能力上仍存在差距。建议在团队后续路线图中评估引入 `external/chefnmr/`、`external/LSDos-3.5.3/` 或 `external/NMR2Struct/` 等替代实现，以补齐"无库可查的新骨架"场景下的结构生成能力。