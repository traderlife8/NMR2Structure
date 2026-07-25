import sys
for m in list(sys.modules.keys()):
    if 'nmr_elucidator' in m: del sys.modules[m]
from nmr_elucidator.io_bruker import read_bruker_folder
from nmr_elucidator.io_bruker_2d import read_bruker_2d_folder
from nmr_elucidator.interpret_2d import interpret_2d_sequential
folder = r'E:\NMR2Structure\NMR\clxh100'
specs = read_bruker_folder(folder)
spectra_2d = read_bruker_2d_folder(folder)
h_spec_pyr = next((s for s in specs if s.solvent == 'Pyr'), None)
pyr_ints = []
if h_spec_pyr and h_spec_pyr.integrals:
    for r in h_spec_pyr.integrals:
        center = (r.low + r.high) / 2.0
        pyr_ints.append((center, r.value))

result = interpret_2d_sequential(
    spectra_2d=spectra_2d,
    carbon_shifts_1d=[c for s in specs if s.nucleus == '13C' for c in [p.shift for p in s.peaks]],
    proton_shifts_1d=[p.shift for s in specs if s.nucleus == '1H' for p in s.peaks],
    dept_summary=[],
    proton_integrals=pyr_ints,
    multiplets=[],
    solvent='CDCl3',
    min_snr=5.0,
)
print(f'Total CHPairs: {len(result.ch_pairs)}')
for cp in sorted(result.ch_pairs, key=lambda x: x.c_ppm, reverse=True):
    h = cp.h_ppm
    c = cp.c_ppm
    print(f'  c={c:6.2f} h={h:6.3f} mult={cp.multiplicity}')
