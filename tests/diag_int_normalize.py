from nmr_elucidator.io_bruker import read_bruker_experiment
spec = read_bruker_experiment(r'E:\NMR2Structure\NMR\clxh100\31')
ints = sorted(spec.integrals, key=lambda x: x.value)
print('Pyr integrals normalized to smallest:')
min_v = min(i.value for i in ints if i.value > 100)
for ig in ints:
    if ig.value > 100:
        ratio = ig.value / min_v
        n_h = round(ratio)
        if n_h == 1:
            mult = 'CH'
        elif n_h == 2:
            mult = 'CH2'
        elif n_h == 3:
            mult = 'CH3'
        else:
            mult = str(n_h)
        print(f'  {ig.low:.2f}-{ig.high:.2f}: value={ig.value:.0f}, ratio={ratio:.2f}, mult={mult}')
