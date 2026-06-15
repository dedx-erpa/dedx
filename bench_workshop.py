"""
bench_workshop.py -- run the model on the charged-particle stopping-power
workshop benchmark cases (Grabowski et al. 2020; Stanek et al. 2024) and
tabulate / plot the alpha-particle stopping vs velocity.

The workshops use an alpha (Z=2) projectile and compare ELECTRONIC stopping vs
alpha velocity vp; average-atom models are reported to agree beyond the peak.

Model dirs (alpha = zp 2; Stanek 2024 Table I single-species cases):
  H1  : python dedx.py --zt=1  --zp=2 --d=1.0  --t=2   ... --od=/tmp/bench_H1
  C1  : python dedx.py --zt=6  --zp=2 --d=10.0 --t=2   ... --od=/tmp/bench_C1
  Be1 : python dedx.py --zt=4  --zp=2 --d=1.84 --t=4.4 ... --od=/tmp/bench_Be1
  Al1 : python dedx.py --zt=13 --zp=2 --d=2.7  --t=1   ... --od=/tmp/bench_Al1
  Cu1 : python dedx.py --zt=29 --zp=2 --d=8.96 --t=1   ... --od=/tmp/bench_Cu1
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HARTREE_EV = 27.211386
AMU_ME = 1822.888486
M_ALPHA = 4.0026

CASES = [  # label, dir, color
    ('H1  (H, 1 g/cc, 2 eV)',    '/tmp/bench_H1',  'C0'),
    ('C1  (C, 10 g/cc, 2 eV)',   '/tmp/bench_C1',  'C1'),
    ('Be1 (Be, 1.84 g/cc, 4.4 eV)', '/tmp/bench_Be1', 'C2'),
    ('Al1 (Al, 2.7 g/cc, 1 eV)', '/tmp/bench_Al1', 'C3'),
    ('Cu1 (Cu, 8.96 g/cc, 1 eV)', '/tmp/bench_Cu1', 'C4'),
]


def vp_au(E_amu):
    """alpha velocity [a.u.] from energy per amu [MeV]."""
    E_alpha_ha = E_amu * M_ALPHA * 1e6 / HARTREE_EV
    return np.sqrt(2.0 * E_alpha_ha / (M_ALPHA * AMU_ME))


fig, ax = plt.subplots(figsize=(7.5, 5.5))
rows = []
for lab, od, c in CASES:
    try:
        d = np.loadtxt('%s/dedx_nuc.dat' % od, comments='#')
    except Exception:
        continue
    E, de, dn, dt = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    vp = vp_au(E)
    ax.plot(vp, de, c + '-', lw=2, label=lab)
    # peak electronic stopping and its velocity; value at vp=2 a.u. reference
    ip = np.argmax(de)
    de_ref = np.interp(2.0, vp, de)
    rows.append((lab, de[ip], vp[ip], de_ref))

ax.set_xlabel('alpha velocity vp (a.u.)')
ax.set_ylabel(r'electronic dE/dx ($10^{-15}$ eV cm$^2$/atom)')
ax.set_title('Workshop benchmark: alpha electronic stopping (Stanek 2024 Table I cases)')
ax.set_xlim(0, 12)
ax.legend(fontsize=9); ax.grid(alpha=.3)
fig.tight_layout()
fig.savefig('exp_bench.pdf'); fig.savefig('exp_bench.png', dpi=130)
print('wrote exp_bench.pdf / .png\n')

print('Table BENCH -- model alpha electronic stopping (10^-15 eV cm^2/atom):')
print('  case                          peak Se   @vp(a.u.)   Se(vp=2 a.u.)')
for lab, pk, vpk, ref in rows:
    print('  %-28s  %6.2f    %5.2f       %6.2f' % (lab, pk, vpk, ref))
