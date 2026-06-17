"""
exp_stanek_C_AA.py -- carbon (10 g/cc, 2 eV) detail panel from Stanek 2024 Fig 12:
overlay our eRPA-AA on BOTH the TD-DFT-MD reference points and Stanek's own two
average-atom (AA) curves.

Units note: the TD-DFT-MD points are digitized vs projectile velocity vp [cm/s];
the two AA curves are digitized vs energy per nucleon [eV/u].  Converting the AA
x-axis with vp = 1.389e9*sqrt(E[MeV/u]) aligns all three peaks (our model 6.6e8,
TD-DFT-MD 6.2e8, Stanek AA 5.7e8 cm/s), confirming the unit assignment.

Takeaway: our independent eRPA-AA lands between Stanek's two AA curves and, like
both of them, sits ~20-30% above the TD-DFT-MD data near the Bragg peak.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NA = 6.022e23
VP = lambda E_amu: 1.389e9 * np.sqrt(E_amu)           # cm/s from MeV/amu
VP_eVu = lambda x: 1.389e9 * np.sqrt(x / 1e6)         # cm/s from eV/u

# our eRPA-AA model
d = np.loadtxt('/tmp/bench_C1/dedx_nuc.dat', comments='#')
n = 10.0 / 12.011 * NA
vp_m = VP(d[:, 0]); dedx_m = d[:, 1] * 1e-15 * n

td = np.loadtxt('data/refs/Stanek_fig12_C.csv', delimiter=',')
aa_s = np.loadtxt('data/refs/Stanek_fig12_C_AA_solid.csv', delimiter=',')
aa_d = np.loadtxt('data/refs/Stanek_fig12_C_AA_dotdash.csv', delimiter=',')

fig, ax = plt.subplots(figsize=(7.5, 6))
ax.loglog(vp_m, dedx_m, 'b-', lw=2.5, label='eRPA-AA (this work)')
ax.loglog(VP_eVu(aa_s[:, 0]), aa_s[:, 1], 'g-', lw=1.8, label='Stanek AA (solid)')
ax.loglog(VP_eVu(aa_d[:, 0]), aa_d[:, 1], 'g-.', lw=1.8, label='Stanek AA (dot-dash)')
ax.loglog(td[:, 0], td[:, 1], 'ks', ms=6, mfc='0.3', label='TD-DFT-MD (Stanek)')
ax.set_xlim(8e7, 2.2e9); ax.set_ylim(3e8, 2e10)
ax.set_xlabel(r'$v_p$ (cm/s)'); ax.set_ylabel(r'$dE_e/dx$ (eV/cm)')
ax.set_title('Stanek Fig 12 (b): C, 10 g/cc, 2 eV')
ax.legend(fontsize=10, loc='lower center'); ax.grid(alpha=.3, which='both')
fig.tight_layout()
fig.savefig('exp_stanek_C_AA.pdf'); fig.savefig('exp_stanek_C_AA.png', dpi=130)
print('wrote exp_stanek_C_AA.pdf / .png\n')

# peak comparison + ratio to TD-DFT-MD
pk = {'eRPA-AA (ours)': dedx_m.max(),
      'Stanek AA solid': aa_s[:, 1].max(),
      'Stanek AA dotdash': aa_d[:, 1].max(),
      'TD-DFT-MD': td[:, 1].max()}
print(' %-20s  peak dEdx (eV/cm)   /TD-DFT-MD' % 'model')
for k, v in pk.items():
    print('  %-20s  %.3e        %.2f' % (k, v, v / pk['TD-DFT-MD']))
