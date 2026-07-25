from nmr_elucidator.io_bruker_2d import read_bruker_2d_experiment, read_hmbc_f1_projection
import numpy as np

hm = read_bruker_2d_experiment(r'E:\NMR2Structure\NMR\clxh100\27')
hmbc_proj = read_hmbc_f1_projection(r'E:\NMR2Structure\NMR\clxh100\27')

# Distribution of projection peak SNR
print(f'Total projection peaks: {len(hmbc_proj)}')
print('SNR distribution:')
import collections
bins = collections.Counter()
for p in hmbc_proj:
    snr = p.snr or 0
    if snr < 5: b = '<5'
    elif snr < 20: b = '5-20'
    elif snr < 100: b = '20-100'
    elif snr < 1000: b = '100-1000'
    else: b = '>=1000'
    bins[b] += 1
for b in ['<5', '5-20', '20-100', '100-1000', '>=1000']:
    print(f'  {b}: {bins[b]} peaks')

# Show top 20 projection peaks by intensity
print('\nTop 20 projection peaks by intensity (with SNR):')
for p in hmbc_proj[:20]:
    print(f'  ppm={p.shift:7.2f} snr={p.snr:.1f} intensity={p.intensity:.1f}')

# 25-150 ppm range filtering
print('\n13C range (10-200 ppm) only:')
in_range = [p for p in hmbc_proj if 10 < p.shift < 200]
print(f'  {len(in_range)} peaks')

# Strict: snr >= 100 AND in 13C range
strict = [p for p in hmbc_proj if 10 < p.shift < 200 and (p.snr or 0) >= 100]
print(f'  SNR >= 100: {len(strict)} peaks')

stricter = [p for p in hmbc_proj if 10 < p.shift < 200 and (p.snr or 0) >= 1000]
print(f'  SNR >= 1000: {len(stricter)} peaks')
