"""
exp_stanek.py -- overlay our eRPA-AA alpha electronic stopping on the TD-DFT-MD
reference points of Stanek 2024 (2nd charged-particle transport workshop), Fig 12,
for the single-species cases H1, C1, Al1.

Stanek Fig 12 plots dEe/dx (eV/cm) vs alpha velocity vp (cm/s).  Our model gives
the stopping cross section in 1e-15 eV cm^2/atom; convert with the case number
density n: dEe/dx[eV/cm] = eps * 1e-15 * n.

Note (Mehlhorn): TD-DFT-MD treats valence electrons explicitly and ADDS an
estimated deep-core contribution; the eRPA-AA includes all bound+free electrons
self-consistently (with DFT-like local-field corrections).  See the radial /
core-depth analysis (stopping_depth.py) for where the AA contributions arise.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NA = 6.022e23
VP = lambda E_amu: 1.389e9 * np.sqrt(E_amu)          # cm/s from MeV/amu

CASES = [  # label, bench dir, A, rho(g/cc), TD-DFT-MD csv
    ('(a) H, 1 g/cc, 2 eV',  '/tmp/bench_H1',  1.008, 1.0,  'data/refs/Stanek_fig12_H.csv'),
    ('(b) C, 10 g/cc, 2 eV', '/tmp/bench_C1',  12.011, 10.0, 'data/refs/Stanek_fig12_C.csv'),
    ('(c) Al, 2.7 g/cc, 1 eV', '/tmp/bench_Al1', 26.98, 2.7, 'data/refs/Stanek_fig12_Al.csv'),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
for ax, (lab, od, A, rho, csv) in zip(axes, CASES):
    d = np.loadtxt('%s/dedx_nuc.dat' % od, comments='#')
    E, de = d[:, 0], d[:, 1]                          # E/amu, electronic
    n = rho / A * NA                                  # atoms/cm^3
    vp = VP(E)
    dedx_cm = de * 1e-15 * n                          # eV/cm
    md = np.loadtxt(csv, delimiter=',')
    ax.loglog(vp, dedx_cm, 'b-', lw=2, label='eRPA-AA (this work)')
    ax.loglog(md[:, 0], md[:, 1], 'ks', ms=5, mfc='0.3', label='TD-DFT-MD (Stanek)')
    ax.set_xlim(1e7, 1e10); ax.set_ylim(3e8, 2e10)
    ax.set_xlabel(r'$v_p$ (cm/s)'); ax.set_title(lab)
    ax.grid(alpha=.3, which='both')
axes[0].set_ylabel(r'$dE_e/dx$ (eV/cm)')
axes[0].legend(fontsize=9, loc='lower center')
fig.suptitle('Stanek 2024 Fig 12: eRPA-AA vs TD-DFT-MD alpha electronic stopping')
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('exp_stanek.pdf'); fig.savefig('exp_stanek.png', dpi=130)
print('wrote exp_stanek.pdf / .png\n')

print(' case   peak(model)  peak(TD-DFT-MD)  model/data@peak')
for lab, od, A, rho, csv in CASES:
    d = np.loadtxt('%s/dedx_nuc.dat' % od, comments='#')
    n = rho / A * NA
    pk_m = (d[:, 1] * 1e-15 * n).max()
    pk_d = np.loadtxt(csv, delimiter=',')[:, 1].max()
    print('  %-22s %.2e   %.2e     %.2f' % (lab, pk_m, pk_d, pk_m / pk_d))
