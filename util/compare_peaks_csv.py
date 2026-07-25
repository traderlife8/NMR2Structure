"""比较 exp 21 和 31_修正的数据 的所有峰数值"""
import sys
import csv
import statistics
sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')
from nmr_elucidator.io_bruker import read_bruker_experiment

# 读取两个实验
exp21 = read_bruker_experiment(r'E:\NMR2Structure\NMR\clxh100\21')
exp31_corr = read_bruker_experiment(r'E:\NMR2Structure\NMR\clxh100\31_修正的数据')

# 溶剂偏移
SOLVENT_SHIFT = 0.0040
MATCH_TOL = 0.05  # 匹配容差 (ppm)

# 按 ppm 降序排列
peaks21 = sorted(exp21.peaks, key=lambda p: -p.shift)
peaks31 = sorted(exp31_corr.peaks, key=lambda p: -p.shift)

# 匹配峰
matched = []
used_31 = set()

for p21 in peaks21:
    expected_pos = p21.shift + SOLVENT_SHIFT
    best_match = None
    best_dist = MATCH_TOL
    for i, p31 in enumerate(peaks31):
        if i in used_31:
            continue
        dist = abs(p31.shift - expected_pos)
        if dist < best_dist:
            best_dist = dist
            best_match = (i, p31)
    
    if best_match:
        i, p31 = best_match
        used_31.add(i)
        delta = p31.shift - p21.shift
        matched.append({
            'exp21_ppm': p21.shift,
            'exp21_intensity': p21.intensity,
            'exp31_ppm': p31.shift,
            'exp31_intensity': p31.intensity,
            'delta_ppm': delta,
            'intensity_ratio': p31.intensity / p21.intensity if p21.intensity > 0.01 else 0,
            'status': 'matched'
        })
    else:
        matched.append({
            'exp21_ppm': p21.shift,
            'exp21_intensity': p21.intensity,
            'exp31_ppm': '',
            'exp31_intensity': '',
            'delta_ppm': '',
            'intensity_ratio': '',
            'status': 'exp21_only'
        })

# 添加 exp31 独有的峰
for i, p31 in enumerate(peaks31):
    if i not in used_31:
        matched.append({
            'exp21_ppm': '',
            'exp21_intensity': '',
            'exp31_ppm': p31.shift,
            'exp31_intensity': p31.intensity,
            'delta_ppm': '',
            'intensity_ratio': '',
            'status': 'exp31_only'
        })

# 按 exp21_ppm 降序排列 (空值排最后)
matched.sort(key=lambda x: -float(x['exp21_ppm']) if x['exp21_ppm'] != '' else -999)

# 输出 CSV
output_path = r'E:\NMR2Structure\NMR\clxh100\31_补充数据\峰值比较.csv'
with open(output_path, 'w', encoding='utf-8-sig', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['序号', 'exp21_ppm', 'exp21_强度', 'exp31_ppm', 'exp31_强度', 'delta_ppm', '强度比', '状态'])
    for i, row in enumerate(matched, 1):
        writer.writerow([
            i,
            f"{row['exp21_ppm']:.4f}" if row['exp21_ppm'] != '' else '',
            f"{row['exp21_intensity']:.2f}" if row['exp21_intensity'] != '' else '',
            f"{row['exp31_ppm']:.4f}" if row['exp31_ppm'] != '' else '',
            f"{row['exp31_intensity']:.2f}" if row['exp31_intensity'] != '' else '',
            f"{row['delta_ppm']:+.4f}" if row['delta_ppm'] != '' else '',
            f"{row['intensity_ratio']:.3f}" if row['intensity_ratio'] != '' else '',
            row['status']
        ])

# 统计
n_matched = sum(1 for r in matched if r['status'] == 'matched')
n_exp21_only = sum(1 for r in matched if r['status'] == 'exp21_only')
n_exp31_only = sum(1 for r in matched if r['status'] == 'exp31_only')

print(f'输出: {output_path}')
print(f'总行数: {len(matched)}')
print(f'  匹配峰: {n_matched}')
print(f'  exp21 独有: {n_exp21_only}')
print(f'  exp31 独有: {n_exp31_only}')

# 计算匹配峰的 delta 统计
deltas = [r['delta_ppm'] for r in matched if r['status'] == 'matched' and r['delta_ppm'] != '']
if deltas:
    print(f'\n匹配峰 delta 统计:')
    print(f'  均值: {statistics.mean(deltas):+.4f} ppm')
    print(f'  标准差: {statistics.stdev(deltas):.4f} ppm')
    print(f'  范围: {min(deltas):+.4f} ~ {max(deltas):+.4f} ppm')
