"""
exp_stanek_C_AA.py -- carbon (10 g/cc, 2 eV) detail panel from Stanek 2024 Fig 12:
overlay our eRPA-AA on BOTH the TD-DFT-MD reference points and Stanek's own two
average-atom (AA) curves.

All three datasets are on the figure's vp [cm/s] axis (the AA curves were
re-digitized directly against the vp axis; an earlier digitization had picked up
the secondary energy/u axis).  No transform is applied.

Takeaway: our independent eRPA-AA peak (1.59e10) lands between Stanek's two AA
curves (solid 1.53e10, dot-dash 1.69e10), all ~19-32% above the TD-DFT-MD data at
the Bragg peak; the three average-atom high-velocity tails converge on the shared
Bethe limit.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NA = 6.022e23
VP = lambda E_amu: 1.389e9 * np.sqrt(E_amu)           # cm/s from MeV/amu

# our eRPA-AA model
d = np.loadtxt('/tmp/bench_C1/dedx_nuc.dat', comments='#')
n = 10.0 / 12.011 * NA
vp_m = VP(d[:, 0]); dedx_m = d[:, 1] * 1e-15 * n

td = np.loadtxt('data/refs/Stanek_fig12_C.csv', delimiter=',')
aa_s = np.loadtxt('data/refs/Stanek_fig12_C_AA_solid.csv', delimiter=',')
aa_d = np.loadtxt('data/refs/Stanek_fig12_C_AA_dotdash.csv', delimiter=',')

fig, ax = plt.subplots(figsize=(7.5, 6))
ax.loglog(vp_m, dedx_m, 'b-', lw=2.5, label='eRPA-AA (this work)')
ax.loglog(aa_s[:, 0], aa_s[:, 1], 'g-', lw=1.8, label='Stanek AA (solid)')
ax.loglog(aa_d[:, 0], aa_d[:, 1], 'g-.', lw=1.8, label='Stanek AA (dot-dash)')
ax.loglog(td[:, 0], td[:, 1], 'ks', ms=6, mfc='0.3', label='TD-DFT-MD (Stanek)')
ax.set_xlim(1e7, 1e10); ax.set_ylim(3e8, 3e10)
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
