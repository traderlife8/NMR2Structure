import sys, importlib

# Force fresh import
for m in list(sys.modules.keys()):
    if 'nmr_elucidator' in m:
        del sys.modules[m]

from nmr_elucidator.io_bruker import read_bruker_folder
from nmr_elucidator.io_bruker_2d import read_bruker_2d_folder
from nmr_elucidator.interpret_2d import interpret_2d_sequential

folder = r'E:\NMR2Structure\NMR\clxh100'
specs = read_bruker_folder(folder)
specs_2d = read_bruker_2d_folder(folder)
h_spec_pyr = next((s for s in specs if s.solvent == 'Pyr'), None)
pyr_ints = []
if h_spec_pyr and h_spec_pyr.integrals:
    for r in h_spec_pyr.integrals:
        center = (r.low + r.high) / 2.0
        pyr_ints.append((center, r.value))
print(f'Pyr integrals count: {len(pyr_ints)}')
print(f'Sample: {pyr_ints[:3]}')

# Also try CDCl3
h_spec_cd = next((s for s in specs if s.solvent == 'CDCl3'), None)
cd_ints = []
if h_spec_cd and h_spec_cd.integrals:
    for r in h_spec_cd.integrals:
        center = (r.low + r.high) / 2.0
        cd_ints.append((center, r.value))
print(f'CDCl3 integrals count: {len(cd_ints)}')

# Try with PYR integrals
result = interpret_2d_sequential(
    spectra_2d=specs_2d,
    carbon_shifts_1d=[c for s in specs if s.nucleus == '13C' for c in [p.shift for p in s.peaks]],
    proton_shifts_1d=[p.shift for s in specs if s.nucleus == '1H' for p in s.peaks],
    dept_summary=[],
    proton_integrals=pyr_ints,
    multiplets=[],
    solvent='Pyr',
)
n_ch3_pyr = sum(1 for cp in result.ch_pairs if cp.multiplicity == 'CH3')
print(f'PYR integrals -> Total: {len(result.ch_pairs)}, CH3: {n_ch3_pyr}')
for cp in result.ch_pairs:
    if cp.multiplicity == 'CH3' or (0.5 <= cp.h_ppm <= 1.4):
        print(f'  h_ppm={cp.h_ppm:.3f} c_ppm={cp.c_ppm:.2f} mult={cp.multiplicity}')

# Try with CD integrals
result2 = interpret_2d_sequential(
    spectra_2d=specs_2d,
    carbon_shifts_1d=[c for s in specs if s.nucleus == '13C' for c in [p.shift for p in s.peaks]],
    proton_shifts_1d=[p.shift for s in specs if s.nucleus == '1H' for p in s.peaks],
    dept_summary=[],
    proton_integrals=cd_ints,
    multiplets=[],
    solvent='CDCl3',
)
n_ch3_cd = sum(1 for cp in result2.ch_pairs if cp.multiplicity == 'CH3')
print(f'\nCDCl3 integrals -> Total: {len(result2.ch_pairs)}, CH3: {n_ch3_cd}')
for cp in result2.ch_pairs:
    if cp.multiplicity == 'CH3' or (0.5 <= cp.h_ppm <= 1.4):
        print(f'  h_ppm={cp.h_ppm:.3f} c_ppm={cp.c_ppm:.2f} mult={cp.multiplicity}')
