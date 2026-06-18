"""
range_proj.py -- total-stopping (electronic + Gordon-Kim nuclear) proton ranges
in C, Al, Ag, Au: CSDA pathlength and projected (practical) range, vs PSTAR.

Both ranges that the standard tables (Ziegler/SRIM, Janni, PSTAR, IAEA) quote are
now produced by the package (dedx_nuc.dat proj_range/detour columns).  The CSDA
range includes the nuclear energy loss; the projected range adds the nuclear
multiple-angle scattering (transport cross section sigma_tr = 4 pi S).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ELEMS = [('C', 6), ('Al', 13), ('Ag', 47), ('Au', 79)]
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, (sym, Z) in zip(axes.ravel(), ELEMS):
    d = np.loadtxt('/tmp/rng_%s/dedx_nuc.dat' % sym, comments='#')   # E,dEe,dEn,dEtot,range,proj,detour
    E, csda, proj = d[:, 0], d[:, 4], d[:, 5]                        # mg/cm^2
    ax.loglog(E, csda, 'b-', lw=2, label='CSDA pathlength (this work)')
    ax.loglog(E, proj, 'r-', lw=2, label='projected range (this work)')
    p = np.loadtxt('data/%s/pstar.txt' % sym, skiprows=8)
    m = (p[:, 0] > E.min()) & (p[:, 0] < E.max())
    ax.loglog(p[m, 0], p[m, 4] * 1e3, 'bo', mfc='none', ms=3.5, markeredgewidth=0.7,
              linestyle='none', label='PSTAR CSDA')
    ax.loglog(p[m, 0], p[m, 5] * 1e3, 'rs', mfc='none', ms=3.5, markeredgewidth=0.7,
              linestyle='none', label='PSTAR projected')
    ax.set_title(sym); ax.grid(alpha=.3, which='both')
    ax.set_xlabel('proton energy (MeV)'); ax.set_ylabel('range (mg/cm$^2$)')
axes[0, 0].legend(fontsize=8, loc='upper left')
fig.suptitle('Total (electronic + nuclear) proton range: CSDA pathlength and projected range vs PSTAR',
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig('range_proj.pdf'); fig.savefig('range_proj.png', dpi=130)
print('wrote range_proj.pdf / .png')
