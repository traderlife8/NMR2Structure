"""
创建 31_修正的数据 文件夹 (标准 Bruker 格式版)

原则: 不修改原始数据 (exp 21, exp 31 保持不动)
策略:
  1. 创建 31_修正的数据/pdata/1/ 结构
  2. 从 exp 31 复制原始谱图数据 (1r, procs 等) - 这些是真实数据不能改
  3. 用修正后的文件替换 peaklist.xml, integrals.txt, peakrng, intrng
  4. 所有输出文件格式与 exp 21 完全一致 (标准 Bruker 格式)

输出: NMR/clxh100/31_修正的数据/
"""
import os
import re
import shutil
import sys
from datetime import datetime
from typing import List, Tuple

sys.path.insert(0, r'E:\NMR2Structure\NMR2Structure')

# 路径
CLXH100_DIR = r'E:\NMR2Structure\NMR\clxh100'
EXP21_DIR = os.path.join(CLXH100_DIR, '21')
EXP31_DIR = os.path.join(CLXH100_DIR, '31')
SUPPLEMENT_DIR = os.path.join(CLXH100_DIR, '31_补充数据')
OUTPUT_DIR = os.path.join(CLXH100_DIR, '31_修正的数据')

# 溶剂偏移 (Pyr-d5 vs CDCl3)
SOLVENT_SHIFT = 0.0040

print("=" * 70)
print("创建 31_修正的数据 文件夹 (标准 Bruker 格式)")
print("=" * 70)

# ============================================================
# 1. 创建目录结构
# ============================================================
print("\n[1] 创建目录结构...")
pdata1_dir = os.path.join(OUTPUT_DIR, 'pdata', '1')
os.makedirs(pdata1_dir, exist_ok=True)
print(f"  创建: {pdata1_dir}")

# ============================================================
# 2. 从 exp 31 复制原始数据文件 (不修改)
# ============================================================
print("\n[2] 从 exp 31 复制原始数据文件...")

copy_files_31 = [
    # 根目录
    ('acqu', 'acqu'),
    ('acqus', 'acqus'),
    ('audita.txt', 'audita.txt'),
    ('fid', 'fid'),
    ('format.temp', 'format.temp'),
    ('pulseprogram', 'pulseprogram'),
    ('scon2', 'scon2'),
    ('shimvalues', 'shimvalues'),
    ('specpar', 'specpar'),
    ('uxnmr.info', 'uxnmr.info'),
    ('uxnmr.par', 'uxnmr.par'),
    # pdata/1 目录
    (os.path.join('pdata', '1', '1i'), os.path.join('pdata', '1', '1i')),
    (os.path.join('pdata', '1', '1r'), os.path.join('pdata', '1', '1r')),
    (os.path.join('pdata', '1', 'auditp.txt'), os.path.join('pdata', '1', 'auditp.txt')),
    (os.path.join('pdata', '1', 'outd'), os.path.join('pdata', '1', 'outd')),
    (os.path.join('pdata', '1', 'parm.txt'), os.path.join('pdata', '1', 'parm.txt')),
    (os.path.join('pdata', '1', 'peaks'), os.path.join('pdata', '1', 'peaks')),
    (os.path.join('pdata', '1', 'proc'), os.path.join('pdata', '1', 'proc')),
    (os.path.join('pdata', '1', 'procs'), os.path.join('pdata', '1', 'procs')),
    (os.path.join('pdata', '1', 'thumb.png'), os.path.join('pdata', '1', 'thumb.png')),
    (os.path.join('pdata', '1', 'title'), os.path.join('pdata', '1', 'title')),
]

for src_rel, dst_rel in copy_files_31:
    src = os.path.join(EXP31_DIR, src_rel)
    dst = os.path.join(OUTPUT_DIR, dst_rel)
    if os.path.exists(src):
        shutil.copy2(src, dst)
        print(f"  复制: {src_rel}")
    else:
        print(f"  跳过 (不存在): {src_rel}")

# ============================================================
# 3. 读取 exp 21 的 peakrng 获取区域信息
# ============================================================
print("\n[3] 读取 exp 21 peakrng 区域信息...")

def parse_exp21_peakrng(path: str) -> List[Tuple[float, float, float, float, int]]:
    """解析 exp 21 peakrng, 返回 [(low_ppm, high_ppm, mi, maxi, region_id), ...]"""
    regions = []
    if not os.path.isfile(path):
        return regions
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()
    # 匹配: low high mi maxi # for region N
    pattern = re.compile(
        r'^\s*([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s*#\s*for\s+region\s+(\d+)',
        re.MULTILINE,
    )
    for m in pattern.finditer(content):
        low_ppm = float(m.group(1))
        high_ppm = float(m.group(2))
        mi = float(m.group(3))
        maxi = float(m.group(4))
        region_id = int(m.group(5))
        regions.append((low_ppm, high_ppm, mi, maxi, region_id))
    regions.sort(key=lambda x: x[4])
    return regions

exp21_peakrng_path = os.path.join(EXP21_DIR, 'pdata', '1', 'peakrng')
exp21_regions = parse_exp21_peakrng(exp21_peakrng_path)
print(f"  exp 21 peakrng: {len(exp21_regions)} 区域")
for low, high, mi, maxi, rid in exp21_regions:
    print(f"    region {rid}: {low:.4f} - {high:.4f} ppm, MI={mi:.3f}, MAXI={maxi:.3f}")

# ============================================================
# 4. 读取修正后的峰表
# ============================================================
print("\n[4] 读取修正后的峰表...")

corrected_peaks_path = os.path.join(SUPPLEMENT_DIR, 'corrected_peaks.txt')
corrected_peaks = []  # [(ppm, intensity, source, exp21_ppm), ...]
if os.path.exists(corrected_peaks_path):
    with open(corrected_peaks_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    ppm = float(parts[1])
                    intensity = float(parts[2])
                    source = parts[3]
                    exp21_ppm = float(parts[4]) if len(parts) > 4 else ppm
                    corrected_peaks.append((ppm, intensity, source, exp21_ppm))
                except (ValueError, IndexError):
                    continue
print(f"  修正后峰数: {len(corrected_peaks)}")

# ============================================================
# 5. 生成标准 Bruker peaklist.xml (多区域格式)
# ============================================================
print("\n[5] 生成标准 Bruker peaklist.xml...")

# 按区域分组峰 (保留 source 信息)
peaks_by_region = {rid: [] for _, _, _, _, rid in exp21_regions}
for ppm, intensity, source, exp21_ppm in corrected_peaks:
    # 找到峰所属区域 (平移后的区域边界)
    for low, high, mi, maxi, rid in exp21_regions:
        # 平移区域边界
        shifted_low = low + SOLVENT_SHIFT
        shifted_high = high + SOLVENT_SHIFT
        region_min = min(shifted_low, shifted_high)
        region_max = max(shifted_low, shifted_high)
        if region_min <= ppm <= region_max:
            # source: 'exp31' -> 'measured', 'exp21_shifted' -> 'estimated'
            source_tag = 'measured' if source == 'exp31' else 'estimated'
            peaks_by_region[rid].append((ppm, intensity, source_tag))
            break

# 生成 XML
peaklist_path = os.path.join(pdata1_dir, 'peaklist.xml')
timestamp = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
with open(peaklist_path, 'w', encoding='utf-8') as f:
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write(f'<PeakList modified="{timestamp}">\n')
    
    for low, high, mi, maxi, rid in exp21_regions:
        # 平移区域边界
        shifted_low = low + SOLVENT_SHIFT
        shifted_high = high + SOLVENT_SHIFT
        
        f.write('  <PeakList1D>\n')
        f.write(f'    <PeakList1DHeader creator="correction" date="{timestamp}" '
                f'expNo="31" name="clxh100" owner="" procNo="1" '
                f'source="E:/NMR2Structure/NMR/clxh100/31_修正的数据">\n')
        f.write(f'      <PeakPickDetails>F1={shifted_low:.3f}ppm, F2={shifted_high:.3f}ppm, '
                f'MI={mi:.2f}cm, MAXI={maxi:.2f}cm, PC=1.000\n')
        f.write('      </PeakPickDetails>\n')
        f.write('    </PeakList1DHeader>\n')
        
        # 写入该区域的峰 (按 ppm 降序, 含 source 属性)
        region_peaks = sorted(peaks_by_region.get(rid, []), key=lambda x: -x[0])
        for ppm, intensity, source_tag in region_peaks:
            f.write(f'    <Peak1D F1="{ppm:.4f}" intensity="{intensity:.2f}" type="0" source="{source_tag}"/>\n')
        
        f.write('  </PeakList1D>\n')
    
    f.write('</PeakList>\n')

total_peaks_in_xml = sum(len(peaks) for peaks in peaks_by_region.values())
print(f"  生成: peaklist.xml ({total_peaks_in_xml} 峰, {len(exp21_regions)} 区域)")

# ============================================================
# 6. 读取 exp 21 积分并生成标准 Bruker integrals.txt
# ============================================================
print("\n[6] 生成标准 Bruker integrals.txt...")

def parse_exp21_integrals(path: str) -> List[Tuple[float, float, float]]:
    """解析 exp 21 integrals.txt, 返回 [(low_ppm, high_ppm, value), ...]"""
    regions = []
    if not os.path.isfile(path):
        return regions
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    # 格式: Number low high value
                    low = float(parts[1])
                    high = float(parts[2])
                    value = float(parts[3])
                    regions.append((low, high, value))
                except ValueError:
                    continue
    return regions

exp21_integrals_path = os.path.join(EXP21_DIR, 'pdata', '1', 'integrals.txt')
exp21_integrals = parse_exp21_integrals(exp21_integrals_path)
print(f"  exp 21 积分区间: {len(exp21_integrals)}")

# 生成标准 Bruker integrals.txt
integrals_path = os.path.join(pdata1_dir, 'integrals.txt')
with open(integrals_path, 'w', encoding='utf-8') as f:
    f.write('   Current data set:\n')
    f.write('   NAME =\tclxh100\tEXPNO =\t31\tPROCNO = 1\n')
    f.write('   DU   =\tE:/NMR2Structure/NMR/clxh100/31_修正的数据\n')
    f.write('\n')
    f.write('   Number   Integrated Region     Integral\n')
    for i, (low, high, value) in enumerate(exp21_integrals, 1):
        # 平移积分区间
        shifted_low = low + SOLVENT_SHIFT
        shifted_high = high + SOLVENT_SHIFT
        f.write(f'   {i:4d}      {shifted_low:.3f}      {shifted_high:.3f}         {value:.5f}\n')
    f.write('\n')

print(f"  生成: integrals.txt ({len(exp21_integrals)} 区间)")

# ============================================================
# 7. 生成标准 Bruker peakrng
# ============================================================
print("\n[7] 生成标准 Bruker peakrng...")

peakrng_path = os.path.join(pdata1_dir, 'peakrng')
with open(peakrng_path, 'w', encoding='utf-8') as f:
    f.write('P    #regions in PPM\n')
    f.write(f'# {timestamp} +0800  correction@NMR2Structure\n')
    f.write('# low field   high field   mi   maxi\n')
    for low, high, mi, maxi, rid in exp21_regions:
        # 平移区域边界
        shifted_low = low + SOLVENT_SHIFT
        shifted_high = high + SOLVENT_SHIFT
        f.write(f'  {shifted_low:.15f}  {shifted_high:.15f}  {mi:.3f}  {maxi:.3f}  # for region {rid}\n')

print(f"  生成: peakrng ({len(exp21_regions)} 区域)")

# ============================================================
# 8. 生成标准 Bruker intrng
# ============================================================
print("\n[8] 生成标准 Bruker intrng...")

intrng_path = os.path.join(pdata1_dir, 'intrng')
with open(intrng_path, 'w', encoding='utf-8') as f:
    f.write('A 1.0 #regions in PPM\n')
    f.write('# low field   high field  bias        slope\n')
    for i, (low, high, value) in enumerate(exp21_integrals, 1):
        # 平移积分区间
        shifted_low = low + SOLVENT_SHIFT
        shifted_high = high + SOLVENT_SHIFT
        f.write(f'  {shifted_low:.15f}  {shifted_high:.15f}  0.0  0.0  # for region {i}\n')

print(f"  生成: intrng ({len(exp21_integrals)} 区间)")

# ============================================================
# 9. 创建 README 说明文件
# ============================================================
print("\n[9] 创建 README 说明...")
readme_path = os.path.join(OUTPUT_DIR, 'README_修正说明.txt')
with open(readme_path, 'w', encoding='utf-8') as f:
    f.write("=" * 70 + "\n")
    f.write("31_修正的数据 说明\n")
    f.write("=" * 70 + "\n\n")
    f.write("本文件夹包含 exp 31 (2D, Pyr-d5) 氢谱的修正数据。\n\n")
    f.write("修正原则:\n")
    f.write("  - 原始数据 (exp 21, exp 31) 保持不动\n")
    f.write("  - 以 exp 21 (1D, CDCl3) 作为参考标准\n")
    f.write("  - 修正 exp 31 的峰表、积分、peakrng\n")
    f.write("  - 所有输出文件格式与 exp 21 完全一致 (标准 Bruker 格式)\n\n")
    f.write("修正内容:\n")
    f.write(f"  1. peaklist.xml: 用 exp 21 模板过滤, 329 -> {total_peaks_in_xml} 峰\n")
    f.write("  2. integrals.txt: 继承 exp 21 积分 (平移溶剂偏移 +0.004 ppm)\n")
    f.write("  3. peakrng: 继承 exp 21 的 13 区域精细分区\n")
    f.write("  4. intrng: 继承 exp 21 积分区间\n\n")
    f.write("数据来源:\n")
    f.write("  - 1r, fid, procs 等: 复制自 exp 31 (原始谱图数据)\n")
    f.write("  - peaklist.xml, integrals.txt, peakrng, intrng: 来自 exp 21 模板 + 修正\n\n")
    f.write("使用方法:\n")
    f.write("  - 2D 流程应读取本文件夹而非 exp 31\n")
    f.write("  - 1D 流程继续使用 exp 21\n")
    f.write("  - 注意: 溶剂为 Pyr-d5, 化学位移与 CDCl3 有 ~0.004 ppm 偏移\n")
print(f"  创建: README_修正说明.txt")

# ============================================================
# 10. 验证输出
# ============================================================
print("\n[10] 验证输出...")

# 验证 peaklist.xml 格式
with open(peaklist_path, 'r', encoding='utf-8') as f:
    content = f.read()
    n_peaklist1d = content.count('<PeakList1D>')
    n_peak1d = content.count('<Peak1D ')
print(f"  peaklist.xml: {n_peaklist1d} 区域, {n_peak1d} 峰")

# 验证 integrals.txt 格式
with open(integrals_path, 'r', encoding='utf-8') as f:
    lines = [l for l in f if l.strip() and not l.startswith('   Number')]
    n_integrals = len([l for l in lines if l.strip() and l.strip()[0].isdigit()])
print(f"  integrals.txt: {n_integrals} 区间")

# 验证 peakrng 格式
with open(peakrng_path, 'r', encoding='utf-8') as f:
    content = f.read()
    n_regions = content.count('# for region')
print(f"  peakrng: {n_regions} 区域")

# 验证 intrng 格式
with open(intrng_path, 'r', encoding='utf-8') as f:
    content = f.read()
    n_regions = content.count('# for region')
print(f"  intrng: {n_regions} 区间")

print("\n" + "=" * 70)
print("完成! 修正数据文件夹已创建 (标准 Bruker 格式):")
print(f"  {OUTPUT_DIR}")
print("=" * 70)
