# 归属表完全匹配(OK)迭代改进计划

## 当前状态

- 当前: **23/41 OK** (改进点 #22 已完成, #23 已实施但未改变 OK 数)
- 目标: **41/41 OK** (100%)
- 剩余: 16 个 DIFF 行 + 2 个未匹配行

## 已完成

- #22: 端基碳区(95-110ppm) CH2/CH3 强制 CH -> +1 OK (C05a)
- #23: CH2/CH3 仅 1 个 HSQC 峰降级 CH -> 已实施, C27b 在 3 峰组中未生效

## 剩余 DIFF 行分类

| 类别 | 程序类型 | 真实类型 | 行数 | 具体行 |
|------|---------|---------|------|--------|
| A | CH2 | CH | 5 | C44b<->25, C08<->5', C11c<->2', C51c<->8, C09a<->16 |
| B | CH | CH2 | 5 | C52<->24a, C40a<->15a, C41b<->15b, C38b<->1b, C39c<->1a |
| C | C | CH | 2 | Cq64<->6(121.5), Cq65<->3'(79.1) |
| D | CH | C | 2 | C37a<->10(37.2), C29a<->13(47.1) |
| E | CH2/CH3 | CH | 1 | C27b<->17(60.5) |
| F | C | CH3 | 1 | Cq66<->21(30.4) |
| 未匹配 | - | CH2 | 2 | 11(29.6), 24b(29.2) |

## 下一步改进 (从 #24 开始)

### #24: C->CH 修复 (HSQC 漏检)
- 文件: `assembly.py` 季碳行生成逻辑 (line 180-194)
- 分析: Cq64(121.0) 应为 CH(H=5.42, 烯碳), Cq65(78.6) 应为 CH(H=4.26, 糖碳)
- 方向: 将 DEPT 信息 (dept_summary) 传入 assemble_assignment_table, 在季碳行生成时检查 DEPT-90 是否有该 13C 峰; 若 DEPT-90 有峰则不是季碳, 创建 CH 行
- 预期: +2 OK
- 需先确认: DEPT-90 (exp 23) 峰表中 121.0/78.6 ppm 是否有峰

### #25: CH->C 修复 (假 HSQC 峰)
- 文件: `interpret_2d.py` 或 `assembly.py`
- 分析: C37a(37.7, H=2.04) 对应真实 10(37.2 C), C29a(46.7, H=2.60) 对应真实 13(47.1 C)
- 方向: 若 DEPT-90 和 DEPT-135 均无该 13C 峰 -> 强制为 Cq
- 预期: +2 OK

### #26: C->CH3 修复 (Cq66 30.0)
- 方向: 同 #24, HSQC 漏检; 用 DEPT-135 正峰 + 1H 峰表 (H~2.1) 推断 CH3
- 预期: +1 OK

### #27-#30: CH<->CH2 互换修复 (10 行)
- 类别 A (CH2->CH, 5行): 用 DEPT-90 确认 CH, 细化 13C 分组
- 类别 B (CH->CH2, 5行): 放宽 CH2 升级条件 (dH 阈值/CH 保护)

### #31+: 未匹配行 (11: 29.6 CH2, 24b: 29.2 CH2)
- 方向: 降低 HSQC 检峰阈值或从 1H/13C 峰表补充

## 执行流程 (每个改进点)

1. 修改代码 (单一改动)
2. 运行: `cd E:\NMR2Structure\NMR2Structure && python -m nmr_elucidator --folder "E:\NMR2Structure\NMR\clxh100"`
3. 读取 comparison_report.txt 统计 OK 行数
4. 记录到 `NMR2Structure/plans/改进点和效果.md`
5. 若 OK 行数提升则保留, 否则回滚并尝试下一方向

## 关键文件

- `nmr_elucidator/interpret_2d.py`: DEPT 解读 + CH3 检测 + HSQC 峰检测
- `nmr_elucidator/assembly.py`: CH2 升级逻辑 + 季碳行生成
- `nmr_elucidator/standard_assignment.py`: 匹配评估 (c_tol=2.5, h_tol=1.0)
- `nmr_elucidator/__main__.py`: CLI 参数 + 流程编排
- `NMR2Structure/plans/改进点和效果.md`: 改进记录
