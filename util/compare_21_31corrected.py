"""比较 exp 21 和 31_修正的数据 的完整性和准确性"""
import sys
sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')
from nmr_elucidator.io_bruker import read_bruker_experiment

# 读取两个实验
exp21 = read_bruker_experiment(r'E:\NMR2Structure\NMR\clxh100\21')
exp31_corr = read_bruker_experiment(r'E:\NMR2Structure\NMR\clxh100\31_修正的数据')

print('=' * 70)
print('exp 21 vs 31_修正的数据 完整性对比')
print('=' * 70)

print(f'\n【1】基本信息')
print(f'  exp 21: 溶剂={exp21.solvent}, 峰数={len(exp21.peaks)}, 积分区间={len(exp21.integrals)}')
print(f'  31_修正: 溶剂={exp31_corr.solvent}, 峰数={len(exp31_corr.peaks)}, 积分区间={len(exp31_corr.integrals)}')

# 峰位范围
peaks21 = sorted([p.shift for p in exp21.peaks], reverse=True)
peaks31 = sorted([p.shift for p in exp31_corr.peaks], reverse=True)
print(f'\n【2】峰位范围')
print(f'  exp 21: {peaks21[0]:.4f} ~ {peaks21[-1]:.4f} ppm')
print(f'  31_修正: {peaks31[0]:.4f} ~ {peaks31[-1]:.4f} ppm')

# 检查缺失的峰区域
print(f'\n【3】缺失峰分析')
shift = 0.004  # 溶剂偏移
missing_peaks = []
for p21 in exp21.peaks:
    expected = p21.shift + shift
    found = any(abs(p31.shift - expected) < 0.5 for p31 in exp31_corr.peaks)
    if not found:
        missing_peaks.append(p21.shift)

print(f'  exp 21 总峰数: {len(exp21.peaks)}')
print(f'  31_修正 总峰数: {len(exp31_corr.peaks)}')
print(f'  缺失峰数: {len(missing_peaks)}')
if missing_peaks:
    missing_sorted = sorted(missing_peaks, reverse=True)
    print(f'  缺失峰位 (前20): {[round(p, 4) for p in missing_sorted[:20]]}')
    
    # 按区域统计
    aromatic = [p for p in missing_peaks if p > 6.0]
    sugar = [p for p in missing_peaks if 3.0 <= p <= 6.0]
    aliphatic = [p for p in missing_peaks if p < 3.0]
    print(f'  缺失分布: 芳香区(>6ppm)={len(aromatic)}, 糖区(3-6ppm)={len(sugar)}, 脂肪区(<3ppm)={len(aliphatic)}')

# 积分对比
print(f'\n【4】积分对比')
print(f'  exp 21 积分区间: {len(exp21.integrals)}')
for ig in exp21.integrals:
    print(f'    {ig.low:.3f}-{ig.high:.3f}: {ig.value:.3f} H')
print(f'  31_修正 积分区间: {len(exp31_corr.integrals)}')
for ig in exp31_corr.integrals:
    print(f'    {ig.low:.3f}-{ig.high:.3f}: {ig.value:.3f} H')

# 检查积分格式是否正确解析
print(f'\n【5】积分解析验证')
print(f'  exp 21 总 H: {sum(ig.value for ig in exp21.integrals):.2f}')
print(f'  31_修正 总 H: {sum(ig.value for ig in exp31_corr.integrals):.2f}')

# 检查 peaklist.xml 格式
print(f'\n【6】peaklist.xml 格式检查')
print(f'  exp 21 峰强度范围: {min(p.intensity for p in exp21.peaks):.2f} ~ {max(p.intensity for p in exp21.peaks):.2f}')
print(f'  31_修正 峰强度范围: {min(p.intensity for p in exp31_corr.peaks):.4f} ~ {max(p.intensity for p in exp31_corr.peaks):.4f}')
