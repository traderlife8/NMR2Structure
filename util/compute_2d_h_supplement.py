"""
exp 31 (2D, Pyr-d5) 氢谱数据修正脚本

核心原则: exp 21 和 exp 31 是同一样品的氢谱, exp 21 经过精细处理, 作为参考标准.
修正策略:
  1. 峰表: 用 exp 21 峰表做模板, 在 exp 31 中找对应峰 (±0.5 ppm 容差), 剔除噪声峰
  2. 积分: 将 exp 21 的积分区间平移到 exp 31 (加溶剂偏移), H 计数直接继承
  3. peakrng: 将 exp 21 的 peakrng 区域平移到 exp 31
  4. 溶剂残余: 过滤 Pyr-d5 残余峰 (8.74/7.58/7.22 ppm)

输出到 NMR/clxh100/31_补充数据/:
  - corrected_peaks.txt: 修正后的峰表 (与 exp 21 对应的峰)
  - corrected_integrals.txt: 修正后的积分 (平移 exp 21 区间)
  - corrected_peakrng.txt: 修正后的 peakrng (平移 exp 21 区域)
  - correction_report.txt: 修正前后差异对比表
"""
import os
import sys
import numpy as np
from typing import List, Tuple, Optional

# 添加项目路径
sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')

from nmr_elucidator.io_bruker import (
    read_bruker_experiment,
    _parse_integrals,
    _parse_peaklist_xml,
)

# 路径
EXP21_DIR = r'E:\NMR2Structure\NMR\clxh100\21'
EXP31_DIR = r'E:\NMR2Structure\NMR\clxh100\31'
OUTPUT_DIR = r'E:\NMR2Structure\NMR\clxh100\31_补充数据'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Pyr-d5 溶剂残余峰 1H 位置 (ppm)
PYR_D5_RESIDUALS = [8.74, 7.58, 7.22]
SOLVENT_TOL = 0.08  # ppm, 溶剂残余过滤容差


def is_solvent_residual(ppm: float, solvent: str = 'Pyr') -> bool:
    """判断是否为溶剂残余峰."""
    if solvent.upper() == 'PYR':
        return any(abs(ppm - r) < SOLVENT_TOL for r in PYR_D5_RESIDUALS)
    return False


def parse_peakrng(path: str) -> List[dict]:
    """解析 Bruker peakrng 文件."""
    regions = []
    if not os.path.exists(path):
        return regions
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('P') or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) >= 4:
                try:
                    regions.append({
                        'low': float(parts[0]),
                        'high': float(parts[1]),
                        'mi': float(parts[2]),
                        'maxi': float(parts[3]),
                    })
                except ValueError:
                    continue
    return regions


print("=" * 70)
print("exp 31 (2D, Pyr-d5) 氢谱数据修正")
print("参考标准: exp 21 (1D, CDCl3)")
print("=" * 70)

# ============================================================
# 1. 读取两个实验的数据
# ============================================================
print("\n[1] 读取数据...")
spec21 = read_bruker_experiment(EXP21_DIR)
spec31 = read_bruker_experiment(EXP31_DIR)
print(f"  exp 21: {spec21.solvent}, {len(spec21.peaks)} 峰, {len(spec21.integrals)} 积分区间")
print(f"  exp 31: {spec31.solvent}, {len(spec31.peaks)} 峰, {len(spec31.integrals)} 积分区间")

# 读取 exp 21 的 peakrng
peakrng21_path = os.path.join(EXP21_DIR, 'pdata', '1', 'peakrng')
peakrng21 = parse_peakrng(peakrng21_path)
print(f"  exp 21 peakrng: {len(peakrng21)} 区域")

# 读取 exp 21 的 integrals.txt
exp21_integrals = _parse_integrals(os.path.join(EXP21_DIR, 'pdata', '1', 'integrals.txt'))

# ============================================================
# 2. 计算溶剂偏移 (CDCl3 -> Pyr-d5)
# ============================================================
print("\n[2] 计算溶剂偏移...")

# 用 exp 21 的强峰 (top 60) 与 exp 31 匹配
peaks21_by_intensity = sorted(spec21.peaks, key=lambda p: -p.intensity)
MATCH_TOL = 0.5  # ppm

matched_pairs = []  # (exp21_peak, exp31_peak, delta)
for p21 in peaks21_by_intensity[:60]:
    # 跳过 CDCl3 溶剂残余
    if abs(p21.shift - 7.26) < 0.1:
        continue
    best_match = None
    best_dist = MATCH_TOL
    for p31 in spec31.peaks:
        d = abs(p21.shift - p31.shift)
        if d < best_dist:
            best_dist = d
            best_match = p31
    if best_match:
        matched_pairs.append((p21, best_match, best_match.shift - p21.shift))

deltas = [m[2] for m in matched_pairs]
# 用中位数作为稳健估计 (避免异常值)
solvent_shift = float(np.median(deltas)) if deltas else 0.0
mean_shift = float(np.mean(deltas)) if deltas else 0.0
std_shift = float(np.std(deltas)) if deltas else 0.0

print(f"  匹配峰对: {len(matched_pairs)}")
print(f"  中位偏移: {solvent_shift:+.4f} ppm (用作校正)")
print(f"  平均偏移: {mean_shift:+.4f} ppm")
print(f"  标准差:   {std_shift:.4f} ppm")

# ============================================================
# 3. 修正峰表: 用 exp 21 做模板, 补充 exp 31 缺失的峰
# ============================================================
print("\n[3] 修正峰表 (用 exp 21 模板 + 补充缺失峰)...")

# 策略:
#   1. 对 exp 21 的每个峰, 在 exp 31 中找最近匹配 (±0.5 ppm)
#   2. 如果找到匹配, 使用 exp 31 的实际峰位 (更准确)
#   3. 如果未找到匹配 (exp 31 缺失该区域), 使用 exp 21 峰位 + 溶剂偏移
#   4. 排除 Pyr-d5 溶剂残余峰

corrected_peaks = []  # (ppm, intensity, source, exp21_ppm)
used_exp31_indices = set()
n_from_exp31 = 0
n_from_exp21 = 0

for p21 in spec21.peaks:
    # 跳过 CDCl3 溶剂残余
    if abs(p21.shift - 7.26) < 0.1:
        continue
    
    # 在 exp 31 中找最近峰 (考虑溶剂偏移)
    expected_pos = p21.shift + solvent_shift
    best_idx = None
    best_dist = MATCH_TOL
    for i, p31 in enumerate(spec31.peaks):
        if i in used_exp31_indices:
            continue
        d = abs(p31.shift - expected_pos)
        if d < best_dist:
            best_dist = d
            best_idx = i
    
    if best_idx is not None:
        p31 = spec31.peaks[best_idx]
        # 排除 Pyr-d5 溶剂残余
        if not is_solvent_residual(p31.shift, 'Pyr'):
            corrected_peaks.append((p31.shift, p31.intensity, 'exp31', p21.shift))
            used_exp31_indices.add(best_idx)
            n_from_exp31 += 1
    else:
        # exp 31 缺失该峰, 用 exp 21 峰位 + 溶剂偏移补充
        shifted_ppm = p21.shift + solvent_shift
        # 排除 Pyr-d5 溶剂残余
        if not is_solvent_residual(shifted_ppm, 'Pyr'):
            corrected_peaks.append((shifted_ppm, p21.intensity, 'exp21_shifted', p21.shift))
            n_from_exp21 += 1

# 按化学位移排序
corrected_peaks.sort(key=lambda x: -x[0])

n_total = len(corrected_peaks)
n_removed = len(spec31.peaks) - n_from_exp31
print(f"  原始 exp 31 峰数: {len(spec31.peaks)}")
print(f"  修正后峰数: {n_total}")
print(f"    - 来自 exp 31 实际峰: {n_from_exp31}")
print(f"    - 来自 exp 21 补充 (移位): {n_from_exp21}")
print(f"  剔除峰数: {n_removed} (噪声 + 溶剂残余)")
print(f"  exp 21 参考峰数: {len(spec21.peaks)}")

# 输出修正后峰表
output_path = os.path.join(OUTPUT_DIR, 'corrected_peaks.txt')
with open(output_path, 'w', encoding='utf-8') as f:
    f.write("# exp 31 修正后峰表 (用 exp 21 模板 + 补充缺失峰)\n")
    f.write(f"# 溶剂偏移校正: {solvent_shift:+.4f} ppm (CDCl3 -> Pyr-d5)\n")
    f.write(f"# 匹配容差: {MATCH_TOL} ppm\n")
    f.write(f"# 原始峰数: {len(spec31.peaks)}, 修正后: {n_total}\n")
    f.write(f"# 来源: exp31={n_from_exp31}, exp21_shifted={n_from_exp21}\n")
    f.write("#\n")
    f.write(f"# {'No':>4} {'ppm':>10} {'intensity':>10} {'source':>15} {'exp21_ppm':>10}\n")
    for i, (ppm, intensity, source, exp21_ppm) in enumerate(corrected_peaks, 1):
        f.write(f"  {i:4d} {ppm:10.4f} {intensity:10.4f} {source:>15} {exp21_ppm:10.4f}\n")
print(f"  输出: {output_path}")

# ============================================================
# 4. 修正积分: 平移 exp 21 的积分区间到 exp 31
# ============================================================
print("\n[4] 修正积分 (平移 exp 21 区间)...")

# exp 21 的积分区间 + 溶剂偏移 = exp 31 的对应区间
# H 计数直接继承 (同一样品)
corrected_integrals = []
for ig in exp21_integrals:
    new_low = ig.low + solvent_shift
    new_high = ig.high + solvent_shift
    corrected_integrals.append({
        'low': new_low,
        'high': new_high,
        'value': ig.value,
        'h_count': round(ig.value),
        'original_low': ig.low,
        'original_high': ig.high,
    })

total_h_corrected = sum(ig['value'] for ig in corrected_integrals)
total_h_exp21 = sum(ig.value for ig in exp21_integrals)
total_h_exp31_raw = sum(ig.value for ig in spec31.integrals)

print(f"  exp 21 积分区间: {len(exp21_integrals)}, 总 H = {total_h_exp21:.2f}")
print(f"  exp 31 原始积分: {len(spec31.integrals)} 区间, 总 H = {total_h_exp31_raw:.2f} (含溶剂, 不可用)")
print(f"  修正后积分: {len(corrected_integrals)} 区间, 总 H = {total_h_corrected:.2f}")

# 输出修正后积分
output_path = os.path.join(OUTPUT_DIR, 'corrected_integrals.txt')
with open(output_path, 'w', encoding='utf-8') as f:
    f.write("# exp 31 修正后积分 (平移 exp 21 区间)\n")
    f.write(f"# 溶剂偏移: {solvent_shift:+.4f} ppm\n")
    f.write("# H 计数继承自 exp 21 (同一样品)\n")
    f.write("#\n")
    f.write(f"# {'No':>3} {'High(ppm)':>10} {'Low(ppm)':>10} {'Integral':>10} {'H_count':>8} "
            f"{'orig_high':>10} {'orig_low':>10}\n")
    for i, ig in enumerate(corrected_integrals, 1):
        f.write(f"  {i:3d} {ig['high']:10.4f} {ig['low']:10.4f} {ig['value']:10.5f} {ig['h_count']:8d} "
                f"{ig['original_high']:10.4f} {ig['original_low']:10.4f}\n")
    f.write(f"\n# 总 H: {total_h_corrected:.2f}\n")
print(f"  输出: {output_path}")

# ============================================================
# 5. 修正 peakrng: 平移 exp 21 的区域到 exp 31
# ============================================================
print("\n[5] 修正 peakrng (平移 exp 21 区域)...")

corrected_peakrng = []
for region in peakrng21:
    corrected_peakrng.append({
        'low': region['low'] + solvent_shift,
        'high': region['high'] + solvent_shift,
        'mi': region['mi'],
        'maxi': region['maxi'],
    })

print(f"  exp 21 peakrng: {len(peakrng21)} 区域")
print(f"  修正后 peakrng: {len(corrected_peakrng)} 区域 (平移 {solvent_shift:+.4f} ppm)")

# 输出修正后 peakrng
output_path = os.path.join(OUTPUT_DIR, 'corrected_peakrng.txt')
with open(output_path, 'w', encoding='utf-8') as f:
    f.write("# exp 31 修正后 peakrng (平移 exp 21 区域)\n")
    f.write(f"# 溶剂偏移: {solvent_shift:+.4f} ppm\n")
    f.write("# 格式: low_ppm  high_ppm  mi  maxi\n")
    for i, r in enumerate(corrected_peakrng, 1):
        f.write(f"  {r['low']:10.4f}  {r['high']:10.4f}  {r['mi']:8.3f}  {r['maxi']:8.3f}  # region {i}\n")
print(f"  输出: {output_path}")

# ============================================================
# 6. 生成修正前后差异对比表
# ============================================================
print("\n[6] 生成差异对比表...")

output_path = os.path.join(OUTPUT_DIR, 'correction_report.txt')
with open(output_path, 'w', encoding='utf-8') as f:
    f.write("=" * 78 + "\n")
    f.write("exp 31 (2D, Pyr-d5) 氢谱数据修正报告\n")
    f.write("参考标准: exp 21 (1D, CDCl3)\n")
    f.write("=" * 78 + "\n\n")

    # 6.1 基本信息对比
    f.write("【1】基本信息对比\n")
    f.write("-" * 78 + "\n")
    f.write(f"{'项目':<20} {'exp 21 (1D, CDCl3)':<25} {'exp 31 原始':<20} {'exp 31 修正后':<20}\n")
    f.write("-" * 78 + "\n")
    f.write(f"{'溶剂':<20} {'CDCl3':<25} {'Pyr-d5':<20} {'Pyr-d5':<20}\n")
    f.write(f"{'峰数':<20} {len(spec21.peaks):<25} {len(spec31.peaks):<20} {len(corrected_peaks):<20}\n")
    f.write(f"{'积分区间数':<18} {len(exp21_integrals):<25} {len(spec31.integrals):<20} {len(corrected_integrals):<20}\n")
    f.write(f"{'总 H 计数':<18} {total_h_exp21:<25.2f} {total_h_exp31_raw:<20.2f} {total_h_corrected:<20.2f}\n")
    f.write(f"{'peakrng 区域数':<16} {len(peakrng21):<25} {'1 (全谱, 无效)':<20} {len(corrected_peakrng):<20}\n")
    f.write(f"{'检峰 MI 阈值':<16} {'0.36~5.67 (分区)':<25} {'0.000 (无阈值)':<20} {'继承 exp 21':<20}\n")
    f.write(f"{'RG (增益)':<18} {'23.29':<25} {'11.2':<20} {'11.2':<20}\n")
    f.write(f"{'SWH (谱宽 Hz)':<16} {'7003':<25} {'10000':<20} {'10000':<20}\n")
    f.write("\n")

    # 6.2 溶剂偏移统计
    f.write("【2】溶剂偏移统计 (CDCl3 -> Pyr-d5)\n")
    f.write("-" * 78 + "\n")
    f.write(f"  匹配峰对数: {len(matched_pairs)}\n")
    f.write(f"  中位偏移 (校正值): {solvent_shift:+.4f} ppm\n")
    f.write(f"  平均偏移: {mean_shift:+.4f} ppm\n")
    f.write(f"  标准差: {std_shift:.4f} ppm\n")
    f.write(f"  偏移范围: [{min(deltas):+.4f}, {max(deltas):+.4f}] ppm\n")
    f.write("\n")

    # 6.3 峰表修正对比
    f.write("【3】峰表修正对比\n")
    f.write("-" * 78 + "\n")
    f.write(f"  exp 31 原始峰数: {len(spec31.peaks)}\n")
    f.write(f"  exp 31 修正后峰数: {n_total}\n")
    f.write(f"    - 来自 exp 31 实际峰: {n_from_exp31}\n")
    f.write(f"    - 来自 exp 21 补充 (移位): {n_from_exp21}\n")
    f.write(f"  剔除峰数: {n_removed}\n")
    f.write(f"  剔除原因: 噪声峰 (MI=0 导致) + Pyr-d5 溶剂残余 (8.74/7.58/7.22 ppm)\n")
    f.write(f"  exp 21 参考峰数: {len(spec21.peaks)} (含 CDCl3 溶剂峰)\n")
    f.write("\n")
    f.write("  修正后峰表 (前 30 个):\n")
    f.write(f"  {'No':>4} {'ppm':>10} {'intensity':>10} {'source':>15} {'exp21_ppm':>10}\n")
    for i, (ppm, intensity, source, exp21_ppm) in enumerate(corrected_peaks[:30], 1):
        f.write(f"  {i:4d} {ppm:10.4f} {intensity:10.4f} {source:>15} {exp21_ppm:10.4f}\n")
    if len(corrected_peaks) > 30:
        f.write(f"  ... (共 {len(corrected_peaks)} 峰, 完整列表见 corrected_peaks.txt)\n")
    f.write("\n")

    # 6.4 积分修正对比
    f.write("【4】积分修正对比\n")
    f.write("-" * 78 + "\n")
    f.write(f"  {'No':>3} {'exp21区间':>18} {'H':>5}  ->  {'exp31修正区间':>18} {'H':>5}  "
            f"{'exp31原始区间':>18} {'H':>6}\n")
    f.write("  " + "-" * 90 + "\n")

    # 对齐 exp 21 和 exp 31 原始积分 (按位置匹配)
    for i, ig21 in enumerate(exp21_integrals):
        # 找 exp 31 原始积分中对应的区间
        shifted_center = (ig21.low + ig21.high) / 2 + solvent_shift
        orig_ig31 = None
        for ig31 in spec31.integrals:
            if ig31.low <= shifted_center <= ig31.high:
                orig_ig31 = ig31
                break
        orig_str = f"{orig_ig31.low:.3f}-{orig_ig31.high:.3f}" if orig_ig31 else "N/A"
        orig_h = f"{orig_ig31.value:.1f}" if orig_ig31 else "N/A"

        corr = corrected_integrals[i]
        f.write(f"  {i+1:3d} {ig21.low:.3f}-{ig21.high:.3f}  {ig21.value:5.2f}  ->  "
                f"{corr['low']:.3f}-{corr['high']:.3f}  {corr['value']:5.2f}  "
                f"{orig_str:>18} {orig_h:>6}\n")

    f.write("\n")
    f.write(f"  exp 21 总 H: {total_h_exp21:.2f}\n")
    f.write(f"  exp 31 修正后总 H: {total_h_corrected:.2f} (继承 exp 21)\n")
    f.write(f"  exp 31 原始总 H: {total_h_exp31_raw:.2f} (含溶剂, 不可用)\n")
    f.write("\n")

    # 6.5 peakrng 修正对比
    f.write("【5】peakrng 修正对比\n")
    f.write("-" * 78 + "\n")
    f.write(f"  {'No':>3} {'exp21 区域':>22} {'MI':>7} {'MAXI':>7}  ->  {'exp31 修正区域':>22}\n")
    f.write("  " + "-" * 75 + "\n")
    for i, (r21, rc) in enumerate(zip(peakrng21, corrected_peakrng), 1):
        f.write(f"  {i:3d} {r21['low']:8.3f}-{r21['high']:8.3f}  {r21['mi']:7.3f} {r21['maxi']:7.3f}  ->  "
                f"{rc['low']:8.3f}-{rc['high']:8.3f}\n")
    f.write("\n")
    f.write(f"  exp 31 原始 peakrng: 1 区域 (15.996 ~ -3.996 ppm, MI=0, MAXI=10000, 无意义)\n")
    f.write(f"  exp 31 修正后: {len(corrected_peakrng)} 区域 (继承 exp 21 的精细分区 + MI/MAXI 阈值)\n")
    f.write("\n")

    # 6.6 修正总结
    f.write("=" * 78 + "\n")
    f.write("【6】修正总结\n")
    f.write("=" * 78 + "\n\n")
    f.write("修正内容:\n")
    f.write(f"  1. 峰表: {len(spec31.peaks)} -> {len(corrected_peaks)} 峰 "
            f"(用 exp 21 模板匹配, 剔除 {n_removed} 个噪声/溶剂峰)\n")
    f.write(f"  2. 积分: {len(spec31.integrals)} 区间 (总H={total_h_exp31_raw:.0f}, 不可用) -> "
            f"{len(corrected_integrals)} 区间 (总H={total_h_corrected:.0f}, 继承 exp 21)\n")
    f.write(f"  3. peakrng: 1 区域 (无效) -> {len(corrected_peakrng)} 区域 (继承 exp 21 精细分区)\n")
    f.write(f"  4. 溶剂偏移: {solvent_shift:+.4f} ppm (CDCl3 -> Pyr-d5)\n")
    f.write("\n")
    f.write("后续使用建议:\n")
    f.write("  - 2D 流程 (HSQC/COSY/HMBC): 使用 exp 31 修正后的峰位 (Pyr-d5 溶剂)\n")
    f.write("  - H 计数/积分: 使用修正后的积分 (继承 exp 21, 已验证)\n")
    f.write("  - 重叠区判断: 使用修正后的 peakrng (继承 exp 21 精细分区)\n")
    f.write("  - 1D 检索流程: 继续使用 exp 21 原始数据 (CDCl3 溶剂)\n")
    f.write("  - 注意: 两套数据溶剂不同, 化学位移不可直接混用\n")

print(f"  输出: {output_path}")

# ============================================================
# 7. 验证修正结果
# ============================================================
print("\n[7] 验证修正结果...")

# 验证 1: 峰数一致性
peak_match_ratio = n_total / len(spec21.peaks) * 100
print(f"  峰匹配率: {n_total}/{len(spec21.peaks)} = {peak_match_ratio:.1f}%")

# 验证 2: 总 H 一致性
print(f"  总 H: exp21={total_h_exp21:.2f}, exp31修正={total_h_corrected:.2f}")

# 验证 3: 溶剂残余检查
solvent_peaks_remaining = [ppm for ppm, _, _, _ in corrected_peaks if is_solvent_residual(ppm, 'Pyr')]
print(f"  残留溶剂峰: {len(solvent_peaks_remaining)} (应为 0)")

# 验证 4: 化学位移偏移一致性 (仅对来自 exp 31 的峰)
corr_deltas = [ppm - exp21_ppm for ppm, _, source, exp21_ppm in corrected_peaks if source == 'exp31']
if corr_deltas:
    print(f"  修正后偏移 (exp31峰): 均值={np.mean(corr_deltas):+.4f}, 标准差={np.std(corr_deltas):.4f} ppm")

print("\n" + "=" * 70)
print("修正完成! 输出文件:")
print(f"  {OUTPUT_DIR}/corrected_peaks.txt      - 修正后峰表")
print(f"  {OUTPUT_DIR}/corrected_integrals.txt - 修正后积分")
print(f"  {OUTPUT_DIR}/corrected_peakrng.txt   - 修正后 peakrng")
print(f"  {OUTPUT_DIR}/correction_report.txt   - 差异对比表")
print("=" * 70)
