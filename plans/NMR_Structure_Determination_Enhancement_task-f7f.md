# NMR 结构推断能力增强方案

## 背景与约束

- 目标：从 NMR 数据推断出真实结构（而非仅给出候选列表）
- 数据约束：仅有 1D NMR（1H / 13C / DEPT），无 2D 谱（HSQC/HMBC/COSY）
- 硬件约束：无 GPU，排除 ChefNMR / NMR2Struct 的 ML 推理
- 2D 谱是实验测量结果，不能从 1D 谱计算生成；但可以从已知结构预测 2D 谱（用于验证）

## 技术路线

采用"扩展候选生成 + 预测验证 + 置信度评估"三层架构：

```
Layer 1: 候选生成（扩大搜索空间）
  ├─ nmrshiftdb2 检索（已有）
  ├─ PubChem 分子式/子结构检索（新增）
  └─ 片段组装生成（新增，用于库外新骨架）

Layer 2: 预测验证（区分真实结构）
  ├─ HOSE code 13C 预测 vs 实测（新增，来自 lucy-ng）
  ├─ nmrsim 谱图模拟 vs 实测（已有，增强）
  └─ 1H 化学位移经验规则验证（新增）

Layer 3: 置信度评估（输出结论）
  ├─ 多源打分融合
  ├─ 置信度分级（determined / probable / ambiguous）
  └─ 证据链报告
```

---

## Phase 1: 集成 HOSE Code 13C 预测验证模块

**目标**：对每个候选结构预测 13C 化学位移，与实测对比，作为核心验证手段。

**来源**：`external/lucy-ng/src/lucy_ng/prediction/` 中的 `C13Predictor`

**实现**：
- 新建 `NMR2Structure/nmr_elucidator/nmr_predict.py`
- 从 lucy-ng 提取 HOSE code 生成器（`hose.py`）和预测器（`predictor.py`）
- 需要 HOSE 统计数据库（SQLite，~2.8 GB）或从 nmrshiftdb2 SDF 自建轻量版
- 接口：`predict_c13(smiles) -> list[(atom_idx, predicted_shift, confidence)]`
- 验证函数：`verify_candidate(smiles, observed_c13, c_tol) -> VerificationScore`

**关键文件**：
- 新建：`nmr_elucidator/nmr_predict.py`（HOSE 预测封装）
- 新建：`nmr_elucidator/hose_generator.py`（从 lucy-ng 提取 HOSE code 生成逻辑）
- 修改：`nmr_elucidator/evidence.py`（集成预测验证分数）

**备选方案**（若 HOSE 数据库太大）：
- 用 nmrshiftdb2 SDF 中已有的 13C 数据自建 HOSE 统计表（~50K 化合物，覆盖常见环境）
- 或用 RDKit 的经验化学位移规则（精度较低但零依赖）

---

## Phase 2: 扩展候选来源 — PubChem 检索

**目标**：当 nmrshiftdb2 无命中时，从 PubChem 按分子式/子结构检索候选。

**来源**：已安装的 `pubchempy` 库

**实现**：
- 新建 `NMR2Structure/nmr_elucidator/retrieve_pubchem.py`
- 按分子式检索（`get_compounds(formula, 'formula')`）
- 按 13C 特征位移推断的子结构做子结构检索
- 对 PubChem 候选用 HOSE 预测验证（Phase 1）
- 限制：单次最多返回 100 条，按分子量/LogP 预筛

**关键文件**：
- 新建：`nmr_elucidator/retrieve_pubchem.py`
- 修改：`nmr_elucidator/__main__.py`（添加 `--use-pubchem` 选项）

---

## Phase 3: 片段组装候选生成（库外新骨架）

**目标**：当数据库检索均无命中时，从 1D 约束出发组装候选结构。

**方法**：
1. 从 13C 化学位移范围推断官能团片段（羰基 160-220 ppm、芳香 110-160 ppm、脂肪 0-90 ppm 等）
2. 从 DEPT 推断碳级数约束（CH/CH2/CH3/Cq 数量）
3. 从 1H 积分推断氢总数
4. 从分子式（若已知）计算 DoU
5. 用 RDKit 枚举满足约束的结构异构体（限制重原子数 <= 15）
6. 对每个候选用 HOSE 预测验证

**实现**：
- 新建 `NMR2Structure/nmr_elucidator/fragment_assembler.py`
- 化学位移范围表（硬编码，来自文献）
- 片段库：常见官能团 SMARTS + 对应 13C 范围
- 组装逻辑：片段组合 → RDKit 连接 → 价态检查 → 分子式过滤
- 限制：组合爆炸控制（最多生成 500 个候选）

**关键文件**：
- 新建：`nmr_elucidator/fragment_assembler.py`
- 新建：`nmr_elucidator/shift_ranges.py`（化学位移范围表）

---

## Phase 4: 统一验证与置信度评估

**目标**：融合多源证据，输出置信度分级结论。

**验证维度**：
| 维度 | 来源 | 权重 |
|---|---|---|
| 13C 覆盖率 | 已有 evidence | 0.20 |
| 1H 覆盖率 | 已有 evidence | 0.15 |
| HOSE 13C 预测 RMS | Phase 1 | 0.30 |
| DEPT 碳级数一致性 | 已有 evidence | 0.15 |
| nmrsim cosine | 已有（增强） | 0.10 |
| 分子式/DoU 一致性 | 已有 formula_solver | 0.10 |

**置信度分级**：
- `determined`（确定）：Top-1 综合分 > 0.85 且与 Top-2 差距 > 0.15
- `probable`（很可能）：Top-1 综合分 > 0.70 且与 Top-2 差距 > 0.08
- `ambiguous`（不确定）：多个候选分数接近，需 2D NMR 进一步确认
- `no_match`（无匹配）：所有候选验证分数 < 0.50

**实现**：
- 新建 `NMR2Structure/nmr_elucidator/determination.py`
- 修改：`nmr_elucidator/evidence.py`（添加 HOSE 验证分数）
- 修改：`nmr_elucidator/report.py`（输出置信度结论）

**关键文件**：
- 新建：`nmr_elucidator/determination.py`（置信度评估引擎）
- 修改：`nmr_elucidator/evidence.py`
- 修改：`nmr_elucidator/report.py`

---

## Phase 5: CLI 与报告更新

**目标**：更新命令行接口和报告格式，体现"结构推断"而非"候选列表"。

**CLI 变更**：
```
python -X utf8 -m nmr_elucidator \
    --folder <...> --sdf <...> \
    --use-pubchem          # 新增：PubChem 扩展检索
    --use-hose             # 新增：HOSE 预测验证
    --use-fragment-assembly # 新增：片段组装（库外）
    --determine            # 新增：输出结构推断结论（含置信度）
```

**报告变更**：
- JSON 增加 `determination` 字段：`{conclusion, confidence, top_candidate, evidence_chain}`
- HTML 报告顶部显示结论（"推断结构为 XXX，置信度：确定/很可能/不确定"）
- Markdown 报告增加证据链摘要

**关键文件**：
- 修改：`nmr_elucidator/__main__.py`
- 修改：`nmr_elucidator/report.py`

---

## Phase 6: 测试与验证

- 用 `NMR/clxh3` 数据端到端测试
- 验证 HOSE 预测对已知结构的 RMS 误差
- 验证 PubChem 检索的命中率和延迟
- 验证置信度分级是否合理（Top-1 是否确实是正确结构）
- 更新 `tests/` 中相关测试

---

## 实施优先级

| 优先级 | Phase | 理由 |
|---|---|---|
| P0 | Phase 1 (HOSE 预测) | 核心验证能力，区分"候选"与"确定" |
| P1 | Phase 4 (置信度) | 输出"真实结构"结论的关键 |
| P2 | Phase 2 (PubChem) | 扩大候选来源，提高命中率 |
| P3 | Phase 5 (CLI/报告) | 用户体验 |
| P4 | Phase 3 (片段组装) | 最复杂，仅库外场景需要 |
| P5 | Phase 6 (测试) | 贯穿各阶段 |

---

## 依赖与风险

| 风险 | 影响 | 缓解 |
|---|---|---|
| HOSE 数据库 2.8 GB 下载慢 | Phase 1 阻塞 | 从 nmrshiftdb2 自建轻量 HOSE 表 |
| lucy-ng 依赖 Python 3.10+ 且用 Pydantic v2 | 集成复杂度 | 仅提取核心算法，不引入整个包 |
| PubChem API 速率限制 | Phase 2 延迟 | 缓存 + 限制并发 |
| 片段组装组合爆炸 | Phase 3 超时 | 严格限制重原子数和候选数上限 |
| 仅 1D 数据无法唯一确定结构 | 置信度上限 | 明确报告"需 2D NMR 确认"的情况 |

---

## 与 MestReNova 功能对照（增强后）

| Mnova 功能 | 本项目增强后 | 状态 |
|---|---|---|
| Database 检索 | nmrshiftdb2 + PubChem | 增强 |
| NMRPredict（结构→谱图） | HOSE code 13C 预测 + nmrsim | 新增 |
| Verify（候选验证打分） | 多源融合验证 + 置信度 | 新增 |
| Structure Elucidation (CASE) | 片段组装（有限，仅 1D） | 部分实现 |
| 置信度/结论输出 | determination 模块 | 新增 |
