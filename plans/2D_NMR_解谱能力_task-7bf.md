# 重叠峰程序化解谱实现

## Summary

按 5.9 节化学家工作流，将当前"仅检测+可视化"升级为**完整的程序化重叠峰分解**：

1. 识别重叠区（1H 积分 >= 4 的 C-H 对按 δH 聚类）
2. 对每个重叠区做 ¹³C 维切片（对标 Mnova 步骤 4）
3. 在切片中自动检峰，确定该区域实际包含几个独立 C-H 对
4. 分解重叠峰：为每个独立 ¹³C 峰分配正确的 H 个数
5. 输出分析报告 CSV + 带切片叠加的等高线图

---

## 修改文件: `nmr_elucidator/interpret_2d.py`

### 新增数据结构

```python
@dataclass
class OverlapRegion:
    """一个重叠区的分析结果."""
    h_ppm_center: float          # 重叠区中心 δH
    integral_total: float        # 该区域总积分
    n_peaks_in_13c_slice: int    # ¹³C 切片中检出的峰数
    c_peaks: List[float]         # 切片中各峰的 δC 位置
    c_peak_intensities: List[float]  # 各峰相对强度
    resolved_pairs: List[CHPair] # 分解后的 C-H 对
    status: str                  # 'resolved'/'partial'/'unresolvable'
    note: str                    # 分析说明
```

在 `Interpretation2D` 中新增字段：
```python
overlap_regions: List[OverlapRegion] = field(default_factory=list)
```

### 新增核心函数: `_step1b_resolve_overlapping`

```python
def _step1b_resolve_overlapping(
    hsqc_spectrum: Spectrum2D,
    ch_pairs: List[CHPair],
    integral_threshold: float = 4.0,
    h_slice_width: float = 0.15,   # ¹H 切片宽度 (ppm)
    c_min_distance: float = 2.0,   # ¹³C 峰最小间距 (ppm)
) -> List[OverlapRegion]:
    """对标 5.9 流程: 对重叠区做 ¹³C 切片 + 分峰检测 + 分解.

    算法:
    1. 筛选 h_integral >= threshold 的 C-H 对
    2. 按 δH 聚类 (间距 < 0.3 ppm 的归为同一重叠区)
    3. 对每个重叠区:
       a. 在 HSQC 数据中提取 δH ± h_slice_width 的窄带
       b. 沿 ¹H 维求和 → 得到 ¹³C 切片 (trace)
       c. 用 scipy.signal.find_peaks 检测切片中的峰
       d. 峰数 = 该区域实际 C-H 对数
       e. 每个 C-H 对的 H 个数 = total_integral / n_peaks
    """
```

**¹³C 切片算法细节**（对标 5.9.2 步骤 4）：

```python
# 提取窄带: δH_center ± h_slice_width
f2_mask = (f2_ppm >= h_center - h_slice_width) & (f2_ppm <= h_center + h_slice_width)
# 沿 ¹H 维 (axis=1) 求和 → ¹³C 投影
trace_13c = np.sum(data[:, f2_mask], axis=1)
# 在 trace_13c 中检峰
from scipy.signal import find_peaks
peaks_idx, props = find_peaks(trace_13c, distance=min_dist_pts, height=threshold)
```

### 修改 `_step1_hsqc_to_ch_pairs`

在现有逻辑后追加重叠峰分解：
- 对积分 >= 4 的 C-H 对，不再简单标记 `unknown`
- 调用 `_step1b_resolve_overlapping` 做 ¹³C 切片分析
- 用切片结果修正每个 C-H 对的 `h_integral` 和 `multiplicity`

### 修改 `write_step_files`

新增输出文件：
- `step1b_overlap_analysis.csv` — 每个重叠区的分析结果
- 等高线图中叠加 ¹³C 切片曲线（右侧子图）

### 修改等高线图 `_draw_hsqc_contour_for_high_integral`

在现有 Mnova 风格布局基础上，新增**右侧 ¹³C 切片子图**：
```
┌─────────────────────────┬──────┐
│     1H 投影 (上方)       │      │
├────┬────────────────────┤ ¹³C  │
│13C │                    │ 切片 │
│投影│   2D 等高线图       │ (带  │
│    │                    │ 检峰)│
└────┴────────────────────┴──────┘
```

切片图中标注检出的峰位置（红色三角）和峰数。

---

## 修改文件: `nmr_elucidator/__main__.py`

- `interpret_2d_sequential` 调用后，将 `overlap_regions` 写入 meta
- 控制台摘要新增重叠峰分解信息

---

## 修改文件: `nmr_elucidator/report.py`

Markdown/HTML 报告"2D 解谱过程"章节新增：
```
### 步骤 1b: 重叠峰分解 (N 个区域)
| 区域 | δH | 总积分 | ¹³C切片峰数 | 分解结果 | 状态 |
```

---

## 输出文件结构 (clxh100 示例)

```
output/clxh100/
├── step1_hsqc_ch_pairs.csv          (更新: 积分>=4的行标注分解结果)
├── step1b_overlap_analysis.csv      (新增: 重叠区分析详情)
├── step1_hsqc_ch_pairs/             (更新: 等高线图含切片)
│   ├── C-19_H-19_integral16.png     (含 ¹³C 切片 + 检峰标注)
│   └── ...
```

---

## 实施顺序

1. `interpret_2d.py`: 新增 `OverlapRegion` + `_step1b_resolve_overlapping`
2. `interpret_2d.py`: 修改 `_step1_hsqc_to_ch_pairs` 集成分解逻辑
3. `interpret_2d.py`: 修改等高线图添加 ¹³C 切片子图
4. `interpret_2d.py`: 修改 `write_step_files` 输出 `step1b_overlap_analysis.csv`
5. `__main__.py`: meta 新增 overlap_regions
6. `report.py`: 报告新增重叠峰分解章节
7. 验证: 运行 clxh100 确认分解结果正确

---

## 关键设计决策

| 决策 | 选择 | 理由 |
|---|---|---|
| 切片宽度 | ±0.15 ppm (¹H) | 太宽引入噪声，太窄丢信号；0.15 约为典型 ¹H 线宽的 3 倍 |
| 检峰最小间距 | 2.0 ppm (¹³C) | ¹³C 线宽通常 < 1 ppm，2 ppm 避免重复检峰 |
| 聚类容差 | 0.3 ppm (¹H) | 与现有 h_tol 一致 |
| 分解后 H 个数 | total_integral / n_peaks | 假设各峰贡献相等（化学家可修正） |
| 状态判定 | n_peaks >= round(integral) → resolved | 切片峰数能解释积分 → 完全分解 |
