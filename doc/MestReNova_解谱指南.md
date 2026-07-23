# MestReNova 解谱指南

> 本文档面向本项目（`NMR2Structure`）解谱与结构推断工程背景，整理 MestReNova（Mnova）在 NMR 数据处理和化学结构解析中的标准工作流、关键插件与操作要点，便于团队按"预测-比对-验证"三步法或 CASE 自动化路径落地。
>
> 信息来源：Mestrelab Research 官网功能说明、Reachsoft 中文代理页、相关学术综述与公开教程。最新功能请以官方手册为准。
>
> 文档整理日期：2026-07-22（第五章 2D 解谱工作流补充完整; 新增零章 NMR 基础概念与术语）

## 零、NMR 基础概念与术语

> 本章是阅读后续章节的必备前提。化学家拿到 NMR 谱图后的所有判断——官能团识别、邻位连接推断、结构归属——都建立在**化学位移**这一基础概念之上。

### 0.1 化学位移 δ 与单位 ppm

**δ (delta)** = 化学位移（chemical shift）

NMR 谱图中原子核的共振位置是相对值,不是绝对频率。它表示原子核相对于参考物（标准品,通常为 **TMS,四甲基硅烷**,人为设定其位移为 0）的偏移量。

**ppm** = parts per million（百万分之一）,是化学位移的单位

之所以不用 Hz 而用 ppm,是因为不同场强的仪器（300 MHz / 600 MHz / 800 MHz）测出的同一原子核的共振频率（Hz）不同,但化学位移（ppm）相同,使不同仪器上测得的谱图可以直接比对。

**公式**

\[
\delta(\text{ppm}) = \frac{\nu_{\text{样}} - \nu_{\text{参}}}{\nu_{\text{谱仪}}} \times 10^6
\]

- ν_样 = 样品原子核的共振频率（Hz）
- ν_参 = 参考物 TMS 的共振频率（Hz,人为定为 0）
- ν_谱仪 = 仪器工作频率（Hz）

### 0.2 典型化学位移区段

化学家按 δ 区间"对号入座"地初步判断官能团。以下经验值涵盖 90% 以上的常见结构。

#### ¹H（δH）

| δH (ppm) | 典型归属 |
|---|---|
| 0 – 1 | 脂肪族饱和 C–H（CH₃, –CH₂–） |
| 1 – 2 | 烯丙位/苄位 H（与 C=C 或芳环间隔一个饱和碳的 H） |
| 2 – 3 | α-羰基/α-杂原子 C–H（–CO–CH₃, –N–CH₂– 等） |
| 3 – 4 | 与 O/N 直接相连的 C–H（–OCH₃, –OCH₂–, –NCH₂–） |
| 4 – 6 | 烯烃 sp² =C–H |
| 6.5 – 8.5 | 芳香 H |
| 9 – 10 | 醛基 H（–CHO） |
| 10 – 13 | 羧酸/酚 –OH（宽峰, 交换性强） |

#### ¹³C（δC）

| δC (ppm) | 典型归属 |
|---|---|
| 0 – 50 | 饱和 sp³ C（烷烃碳、季碳常出现在 30–50 区） |
| 50 – 90 | 与 O/N 直接相连的 sp³ C（–O–CH₃ ≈ 55–60, –N–CH₂– ≈ 40–60） |
| 100 – 150 | 烯烃/芳香 sp² C（烯烃 120–140, 芳环 115–135） |
| 160 – 180 | 羰基（酮、醛、酰胺、酯的 C=O 多在此区） |
| 180 – 210 | 羧酸/酯/酸酐 C=O（更高场） |

> **注意**：边界并非绝对,受取代基电负性、立体效应、氢键、温度、溶剂等因素影响可在 ±5 ppm 内漂移。化学家解谱时需结合 2D（HSQC/HMBC）确认归属,而不是孤立看 δ 数值。

### 0.3 化学家的"看图说话"

化学位移反映原子核周围的电子云密度（即屏蔽程度）：

- **低场（高 δ）** = 去屏蔽 = 电子云密度低 = 邻接吸电子基团
- **高场（低 δ）** = 屏蔽 = 电子云密度高 = 邻接给电子基团

化学家拿到一个未知样品时,基本顺序是：

1. 看 ¹H 谱的**积分总数** → 确认 H 总数与分子式匹配
2. 看每一组的 **δH 位置** → 判断官能团类型
3. 看**裂分模式（J 耦合）** → 判断邻位 H 个数
4. 看 **¹³C 谱峰数** → 确认分子对称性 / 不同 C 数
5. 交叉比对 2D（HSQC/COSY/HMBC）→ 完成归属

这就是后续 5.4 节解谱流程的化学直觉来源。**δ 是因, ppm 是量纲**;化学位移 δ 才是化学家解读结构的"密码本"。

### 0.4 常用参考峰速查

| 化合物 | δH (ppm) | δC (ppm) |
|---|---|---|
| TMS (内标) | 0.00 | 0.00 |
| CDCl₃ (氘代氯仿) | 7.26 | 77.00（t） |
| DMSO-d₆ (氘代 DMSO) | 2.50 | 39.50（7 重峰） |
| CD₃OD (氘代甲醇) | 3.31 | 49.00（7 重峰） |
| D₂O (重水) | 4.79 | — |
| 丙酮-d₆ | 2.05 | 29.80 / 206.00 |

> 化学家在归属时若发现"某个峰对不上任何官能团",先检查它是不是**溶剂残留峰**或**水峰**——这是新手最常踩的坑。

---

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

### 5.1 化学家的 2D NMR 解谱逻辑（“拼图”策略）

化学家使用 2D NMR 解谱的核心逻辑是“拼图”——先建立小块片段，再拼接成完整结构：

```
步骤 1: HSQC —— 建立 C-H “拼图块”
    每个交叉峰 = 一个直接 C-H 连接对
    → 得到所有 CH/CH₂/CH₃ 单元 (季碳无 HSQC 峰)
        ↓
步骤 2: COSY / TOCSY —— 拼出氢自旋系统
    从已知的 H 出发, 沿 COSY 交叉峰找到相邻 H
    → 得到 H-C-C-H 片段 (邻碳连接链)
        ↓
步骤 3: HMBC —— 连接片段、确认取代位点
    从已归属的 H 出发, 找 2-3 键外的 C
    → 跨越季碳/杂原子连接片段, 确定取代位置
        ↓
步骤 4: NOESY / ROESY —— 确定立体化学
    空间近邻 (< 5 Å) 的 H-H 相关
    → 相对构型、双键几何、构象
```

**关键原则**：
- HSQC 是起点：给每个 H 找到其直接连接的 C，建立“拼图块”
- COSY 是横向拼接：在同一自旋系统内连接相邻 H
- HMBC 是纵向连接：跨越季碳/杂原子把片段串联成完整骨架
- NOESY 是立体约束：确定平面结构后的 3D 排列

### 5.2 Mnova 中 2D 谱的处理流程

**5.2.1 导入与自动处理**

1. 把含 `ser`（或 `fid`）的 Bruker 2D 实验文件夹拖入 Mnova；
2. Mnova 自动执行：FT（二维傅里叶变换）、相位校正、基线校正；
3. 若 `File → Preferences → NMR → Import` 中勾选了“Auto attach traces”，已打开的 1D 谱会自动作为 2D 图的外围投影；
4. 也可手动拖放：将 1H 谱拖到 HSQC 上方→水平投影；将 13C 谱拖到左侧→垂直投影。

**5.2.2 参考定标（2D 谱）**

- 点击工具栏 `Reference` 按钮，先沿 f2（¹H 维）定标，再沿 f1（¹³C 维）定标；
- 若已加载 1H 和 13C 谱，可用 `Analysis → Reference → Absolute Reference` 一键同步定标所有 2D 谱；
- **重要**：确保 HSQC/HMBC/COSY 的 F2 (¹H) sr 值与 1H 谱一致，F1 (¹³C) sr 值与 13C 谱一致。

**5.2.3 等高线调整**

- 右键谱图 → `Contours` 调整等高线级别：
  - `Positive Base`：起始强度（越小显示越多峰，但噪声也越多）；
  - `Total Number`：等高线层数（典型 12–20）；
  - `Increment`：层间倍增因子（典型 1.2–1.5）；
- HMBC 远程相关峰弱，需降低 Positive Base 以显示弱峰；
- COSY 无需调相位（绝对值谱）；HSQC 相敏谱需调相（CH/CH₂/CH₃ 颜色区分）。

**5.2.4 峰检测与标注**

- `Analysis → Peak Picking → 2D`：自动检测交叉峰，标注 (f2, f1) 坐标；
- 手动标注：放大局部，双击交叉峰添加文本标注；
- `View → Tables → Peaks` 插入峰位表（包含所有交叉峰的 F2/F1 坐标和强度）。

### 5.3 各 2D 谱的核心用途与解读要点

| 实验 | 轴 | 提供信息 | 解谱作用 | 解读要点 |
|---|---|---|---|---|
| COSY / TOCSY | ¹H × ¹H | ³J(HH) 耦合相关 | 拼出邻碳 H 的耦合链 | 对角线两侧的对称交叉峰；取一个即可 |
| HSQC / HSQC-DEPT | ¹H × ¹³C | ¹J(CH) 直接连接 | 给每个 H/C 配对 | CH/CH₃ 正峰(红), CH₂ 负峰(蓝); 季碳无峰 |
| HMBC | ¹H × ¹³C | ²J, ³J 远程相关 | 连接片段、确认取代位点 | 峰较弱; ²J 通常比 ³J 强; 注意排除 ¹J 干扰 |
| NOESY / ROESY | ¹H × ¹H | 空间近邻 (< 5 Å) | 立体构型、相对构型 | 与分子内距离相关, 与键数无关 |
| INADEQUATE | ¹³C × ¹³C | ¹J(CC) 直接连接 | 完整碳骨架 | 灵敏度极低, 需大量样品 |

### 5.4 化学家的完整 2D 解谱实操步骤（以 Mnova 为例）

> 以下是化学家在 Mnova 中进行 2D NMR 全归属的标准操作流程。

**第一步：从 1D 谱提取已知信息**

1. 在 1H 谱上完成积分、多重峰分析、标峰，记录每个 H 的 δ、J、积分、峰形；
2. 在 13C 谱上标峰，结合 DEPT 确定每个 C 的级数 (CH/CH₂/CH₃/Cq)；
3. 根据化学位移范围初步判断官能团（芳香区 / 脂肪区 / 羰基区 / 活泼氢）。

**第二步：HSQC —— 建立 C-H 对应表**

1. 打开 HSQC 谱，将 1H 谱拖到上方作为水平投影，13C 谱拖到左侧作为垂直投影；
2. 从 1H 谱的每个已确认峰出发，向下（或向上）在 HSQC 中找对应交叉峰；
3. 记录每个交叉峰的 (δH, δC) 对，建立归属表：

| H 编号 | δH (ppm) | δC (ppm) | 碳级数 (DEPT) | 归属 |
|---|---|---|---|---|
| H-1 | 7.26 | 128.5 | CH | Ar-H |
| H-2 | 3.85 | 55.2 | CH₃ | OCH₃ |
| ... | ... | ... | ... | ... |

4. 注意：HSQC-DEPT (multiplicity-edited) 中 CH/CH₃ 为正峰（红色），CH₂ 为负峰（蓝色），可替代 DEPT 实验；
5. 季碳（Cq）在 HSQC 中无交叉峰，需从 13C 谱中确认哪些 C 没有 HSQC 对应。

**第三步：COSY —— 拼出氢自旋系统**

1. 打开 COSY 谱，将 1H 谱拖到上方和左侧作为双向投影；
2. 从已归属的 H 峰出发，在对角线两侧找交叉峰：
   - 交叉峰 (δHa, δHb) 表示 Ha 与 Hb 存在 ³J 耦合（相隔 3 个键，H-C-C-H）；
3. 沿耦合链逐步延伸：H-1 → H-2 → H-3 → ...，拼出完整的自旋系统；
4. 对于复杂分子，TOCSY 比 COSY 更有用：从任何一个分离良好的峰（如糖的端基氢）出发，可找到同一自旋系统内的所有 H；
5. 记录 COSY 相关对：

| 相关对 | δH (ppm) | 解读 |
|---|---|---|
| H-5 ↔ H-6 | 7.84 ↔ 7.52 | 邻位芳香 H, 确认取代模式 |
| H-2 ↔ H-3 | 3.85 ↔ 4.12 | 相邻 CH-CH |

**第四步：HMBC —— 连接片段、确认取代位点**

1. 打开 HMBC 谱，同样附加 1H 和 13C 投影；
2. 从已归属的 H 出发，找 2-3 键外的 C 相关：
   - 交叉峰 (δH, δC) 表示该 H 与该 C 相隔 2–3 个键；
3. **关键应用**：
   - 跨越季碳：如芳香 H 与季碳 C 的 HMBC 相关确认取代位点；
   - 连接片段：如 CH₃ 与相邻 C=O 的相关确认乙酰基连接位置；
   - 糖苷键：糖端基 H 与苷元 C 的 HMBC 相关确定糖苷键位置；
4. 注意：
   - HMBC 峰较弱，需降低 contour level 才能看到；
   - ²J 相关通常比 ³J 强；
   - 要排除 ¹J 干扰（某些 HMBC 序列未完全抑制 ¹J）；
   - 某些弱相关可能未检测到，不能因为“没看到”就认为“不存在”；
5. 记录关键 HMBC 相关：

| H (δH) | → C (δC) | 键数 | 解读 |
|---|---|---|---|
| H-2 (3.85) | C-1 (170.5) | ²J | OCH₃ 连在 C=O 上 (甲酯) |
| H-5 (7.84) | C-7 (135.2, Cq) | ³J | 确认取代位点 |
| CH₃ (2.10) | C=O (170.5) | ²J | 乙酰基 |

**第五步：NOESY / ROESY —— 确定立体化学**

1. 打开 NOESY/ROESY 谱，将 1H 谱作为投影；
2. 找空间近邻的 H-H 相关（与键数无关，仅与空间距离 < 5 Å 相关）；
3. 关键应用：
   - 双键几何：cis/trans 的 NOE 差异；
   - 相对构型：环己烷 axial/equatorial、糖环构象；
   - 构象偏好：优势构象的确认；
4. ROESY 适合中等分子量 (MW 500–1500)，NOESY 适合小分子 (< 500) 或大分子 (> 1500)。

**第六步：汇总归属表 + 结构确认**

1. 将以上所有 2D 相关信息汇总为完整归属表：

| 位置 | δC | 碳级数 | δH (J, 积分) | COSY | HSQC | HMBC | NOESY |
|---|---|---|---|---|---|---|---|
| C-1 | 170.5 | Cq | — | — | — | H-2, CH₃ | — |
| C-2 | 55.2 | CH₃ | 3.85 (s, 3H) | — | ✔ | C-1 | H-3 |
| ... | ... | ... | ... | ... | ... | ... | ... |

2. 核对：
   - 每个非季碳都应有 HSQC 对应；
   - 每个 H 的 COSY 相关应符合其峰形 (J 值)；
   - HMBC 相关应能串联所有片段成完整骨架；
   - 分子式、DoU 与结构一致。

### 5.5 Mnova 二维谱增强功能

- **自动峰标注**：`Analysis → Peak Picking → 2D`；
- **多谱图同步缩放 / 联动平移**：同时打开 HSQC + HMBC + COSY，缩放一个其他同步；
- **2D 与 1D 联动**：点击 2D 峰高亮对应 1D 位置；
- **谱图数据库检索**：把 2D 谱作为查询条件检索 Mnova DB（~26 万 ¹³C / ~10 万 ¹H）；
- **归属表自动编号**：在结构上标注原子编号，与谱峰表联动。

### 5.6 用 Mnova Verify 进行 2D 结构验证

当已有候选结构时，Mnova Verify 可自动验证：

1. 在 Mnova 中打开候选结构（mol/SMILES/绘制）；
2. 菜单 `Analysis → Verify` 或工具栏 Verify 按钮；
3. Verify 自动执行：
   - 用 NMRPredict 预测候选结构的 ¹H、¹³C、HSQC 谱图；
   - 将预测谱与实测 1D/2D 谱比对；
   - 检查 HSQC 关联匹配率、HMBC 远程相关匹配率；
   - 给出总分 (-1.0 ~ +1.0) + 分项得分；
4. 输出：
   - 归属表（每个 H/C 标在结构上）；
   - 不一致项告警（哪个峰预测偏差大）；
   - 多候选排序（分数从高到低）。

> **本项目对标**：`correlation_2d.score_candidate_2d()` 实现了类似的 HSQC/HMBC/COSY 匹配率计算，
> 但使用 nmrshiftdb2 的 per-atom 归属而非 NMRPredict 预测。

### 5.7 用 Mnova Structure Elucidation (CASE) 进行 2D 自动解谱

当结构完全未知时，CASE 模块从 2D NMR 数据自动生成候选：

1. 输入分子式（由 HRMS 得到）；
2. 导入 1H、13C、HSQC（必须）、HMBC（强烈推荐）、COSY（可选）；
3. CASE 引擎工作流：
   - 从 HSQC 提取所有 C-H 对（“拼图块”）；
   - 从 COSY 提取 H-H 邻接关系（片段内连接）；
   - 从 HMBC 提取 H-C 远程相关（片段间连接）；
   - 用专家系统枚举所有满足 2D 约束的候选结构；
   - 对每个候选用 NMRPredict 预测谱图并打分；
4. 输出排序候选列表 + 归属表。

> **本项目对标**：`fragment_assembler.assemble_candidates()` 实现了简化版片段组装（仅 1D），
> `correlation_2d.py` 实现了 2D 约束验证。完整 CASE 能力可对接 `external/lucy-ng/` 或 `external/LSDos-3.5.3/`。

### 5.8 2D 解谱中的常见陷阱

| 陷阱 | 说明 | 应对 |
|---|---|---|
| HMBC 中的 ¹J 残留 | 某些 HMBC 序列未完全抑制直接相关 | 与 HSQC 对照，排除已知的 ¹J 对 |
| 弱相关未检测到 | HMBC/NOESY 中某些弱峰被噪声淹没 | 不能因“没看到”就认为“不存在” |
| 假峰 (t1 noise) | 2D 谱中沿 f1 方向的条纹噪声 | 与真实交叉峰区分：调 contour level 后假峰消失 |
| 化学位移简并 | 两个 H 或 C 化学位移相同/极近 | HSQC 中表现为一个交叉峰对应两个 C；用 TOCSY 或更高场强区分 |
| 溶剂峰干扰 | 残余溶剂峰在 2D 谱中产生假相关 | 确认溶剂峰位置，排除其相关 |
| 旋转异构体 | 酰胺键等慢交换导致双倍峰 | 变温 NMR 或延长混合时间确认 |

### 5.9 重叠峰的详细分析流程（Mnova 实操）

> 重叠峰是 NMR 解谱的核心难点。化学家使用 2D NMR（尤其 HSQC）将 1D 中重叠的峰在第二维展开，再通过等高线图目视检查、积分分析和自旋系统追踪来分解重叠信号。

#### 5.9.1 重叠峰的识别

**1D 谱中的识别信号：**

| 信号 | 说明 |
|---|---|
| 积分异常 | 1H 积分 ≥ 4 通常意味着多个峰重叠（如 2×CH₂ 或 CH₃+CH） |
| 峰形异常 | 峰宽明显大于邻近峰，或峰形不对称 |
| 多重峰分析失败 | Mnova Multiplet Analysis 无法给出合理的 J 值 |
| 化学位移简并 | 多个 H 恰好落在相同/极近的 δH |

**2D HSQC 中的确认：**

1. 打开 HSQC 谱，将 1H 谱拖到上方作为水平投影；
2. 在 1H 谱上找到积分异常的区域（如 δH 3.5–4.0，积分 4H）；
3. 向下看 HSQC 中该 δH 范围内有几个交叉峰：
   - **若多个交叉峰在 ¹³C 维分开** → 确认是重叠峰，每个交叉峰对应一个 C-H 对；
   - **若只有一个交叉峰** → 可能是同碳等价 H（如 CH₂ 的两个 H）或真正的简并；
4. 调整等高线级别（降低 Positive Base）查看是否有弱峰被掩盖。

#### 5.9.2 Mnova 等高线图分析操作

**步骤 1：设置合适的等高线级别**

```
右键谱图 → 2D Spectrum → Contours:
  - Positive Base: 从最大峰强度的 5–10% 开始
  - Number of Levels: 12–20 层
  - Increment: 1.2–1.5（层间倍增因子）
  - Negative: 勾选（显示负峰，CH₂ 在相敏 HSQC 中为负）
```

**步骤 2：添加 1D 投影（关键！）**

```
1. 将 1H 谱拖到 HSQC 上方 → 水平投影
2. 将 13C 谱拖到 HSQC 左侧 → 垂直投影
3. 在投影图上标峰，与 2D 交叉峰对照
```

**步骤 3：局部放大检查重叠区**

```
1. 鼠标框选重叠区域 → 滚轮放大
2. 观察等高线形状：
   - 圆形/椭圆形 → 单个峰
   - 拉长/不规则 → 可能是多个峰叠加
   - 同心圆套叠 → 强峰掩盖了弱峰
3. 逐步降低 contour level，看是否有新峰出现
```

**步骤 4：沿 ¹³C 维切片**

```
1. 在重叠峰的 δH 位置画一条垂直线（F1 方向切片）
2. 观察切片中有几个 ¹³C 峰：
   - 每个 ¹³C 峰对应一个独立的 C-H 对
   - 峰面积比 ≈ H 个数比
```

#### 5.9.3 重叠峰的分解策略

**策略 A：利用 ¹³C 维分辨率（首选）**

- HSQC 的 ¹³C 维分散度是 ¹H 的 ~20 倍；
- 即使 δH 完全重叠，δC 通常能分开；
- 在 HSQC 中直接读取每个交叉峰的 (δH, δC)。

**策略 B：利用 TOCSY 追踪自旋系统**

- 从分离良好的峰（如糖的端基 H）出发；
- TOCSY 显示同一自旋系统内的所有 H；
- 通过自旋系统归属反推重叠峰的归属。

**策略 C：利用 HMBC 远程相关**

- 重叠峰的 H 可能与不同的季碳有 HMBC 相关；
- 通过 HMBC 相关模式区分不同片段。

**策略 D：变场强/变溶剂**

- 更高场强（600 → 800 MHz）增加化学位移分散度；
- 换氘代溶剂可能改变化学位移，解除简并。

#### 5.9.4 重叠峰的积分处理

| 情况 | 处理方法 |
|---|---|
| HSQC 中完全分开 | 分别积分每个交叉峰，H 个数 = 1 |
| HSQC 中部分重叠 | 沿 ¹³C 维切片，分峰拟合后积分 |
| 同碳等价 H（CH₂） | 积分 = 2H，HSQC 中为一个交叉峰 |
| 无法分解 | 记录为“重叠区”，H 个数 = 积分值，归属待定 |

#### 5.9.5 本项目对标实现

| Mnova 操作 | 本项目实现 |
|---|---|
| 1H 积分异常检测 | `interpret_2d._step1`: 积分 ≥ 4 标记为潜在重叠 |
| HSQC 等高线图 | `interpret_2d._draw_hsqc_contour_for_high_integral`: 为高积分 C-H 对绘制局部等高线图 |
| ¹³C 维分辨率检查 | 等高线图中可目视确认是否有多个交叉峰 |
| 重叠峰分解 | 待实现：沿 ¹³C 维切片 + 分峰拟合 |
| TOCSY 自旋系统追踪 | `interpret_2d._step2`: COSY 自旋系统建立（TOCSY 待扩展） |

> **注意**：当前实现仅检测并可视化重叠峰，尚未实现自动分峰拟合。化学家需查看 `step1_hsqc_ch_pairs/` 文件夹中的等高线图，人工确认重叠峰的分解结果。

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
- MestReNova 核磁软件使用介绍：<https://www.go-soft.cn/news/cjwt/2138.html>
- 解谱软件教程：<https://www.go-soft.cn/news/cjwt/1087.html>
- 核磁软件使用宝典：<https://www.molaid.com/baike/tool/569>
- 培训日程（1D/2D + 数据库 + 结构验证）：<https://m.sohu.com/a/867201037_610519>
- Mnova Tools for Structure Elucidation Using 1D and 2D NMR（彭琛讲座）：<https://ibook.antpedia.com/x/652627.html>
- 手把手教你用 Mnova + ChemDraw：<https://m.blog.csdn.net/weixin_29290947/article/details/159134830>
- 5 分钟搞懂核磁共振谱图解析（Chemdraw + Mnova）：<https://wenku.csdn.net/column/r95rv9h6oxp>
- DeepSAT（HSQC → 结构，文中比较 ACD/SE、Mnova SE、LSD）：<https://pmc.ncbi.nlm.nih.gov/articles/PMC10406729/>
- HSQC Spectra Simulation and Matching：<https://pubs.acs.org/doi/10.1021/acs.jcim.3c01735>
- NMR 综合解析分享（1H, 13C, DEPT, COSY, HSQC, HMBC）：<http://www.360doc.com/content/24/0425/19/68893998_1121448752.shtml>
- 有机化合物 NMR 谱图解析（2D 详解）：<https://zhidao.baidu.com/question/213028094745872765.html>
- Structure verification through computer-assisted spectral assignment (CCASA/LSD)：<https://m.zhangqiaokeyan.com/journal-foreign-detail/0704039582786.html>
- Mnova Verify 结构验证功能介绍：<https://m.instrument.com.cn/netshow/SH104663/C504841.html>

---

## 结论

MestReNova 是 NMR 解谱和化学结构推断领域的工业标准工具，覆盖从谱图处理、谱图预测到 CASE 自动结构解析的完整链路。本项目（`NMR2Structure`）目前的 `io_bruker → retrieve → evidence → report` 链路在功能上接近 Mnova 的 `Database` + `Verify` 组合，但在 `Structure Elucidation`（从头生成候选）和 `NMRPredict`（预测谱图）两块能力上仍存在差距。建议在团队后续路线图中评估引入 `external/chefnmr/`、`external/LSDos-3.5.3/` 或 `external/NMR2Struct/` 等替代实现，以补齐"无库可查的新骨架"场景下的结构生成能力。