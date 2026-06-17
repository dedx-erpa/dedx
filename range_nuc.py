"""
range_nuc.py -- proton CSDA range in C, Al, Ag, Au: electronic-only vs
electronic+nuclear (Gordon-Kim), compared with PSTAR.

This is the nuclear-section companion to Fig. 7 of the main dedx paper (which
showed the eRPA electronic-only range vs PSTAR).  PSTAR's CSDA range includes
nuclear stopping, so at low proton energy the electronic-only eRPA range runs
slightly long; adding the GK nuclear term shortens it toward PSTAR.

Run after generating /tmp/rng_<sym> with:
  python dedx.py --zt=Z --zp=1 --d=rho --t=0.025 --aa=2 --nuc=1 --npot=gk \
                 --emin=0.001 --emax=100 --mep=100 --od=/tmp/rng_<sym>
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ELEMS = [('C', 6), ('Al', 13), ('Ag', 47), ('Au', 79)]

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
for ax, (sym, Z) in zip(axes.ravel(), ELEMS):
    de = np.loadtxt('/tmp/rng_%s/dedx.dat' % sym, comments='#')      # E, dEdx_e, range_e
    dn = np.loadtxt('/tmp/rng_%s/dedx_nuc.dat' % sym, comments='#')  # E, dEe, dEn, dEtot, range_tot
    E = de[:, 0]
    rng_e, rng_t = de[:, 2], dn[:, 4]                                # mg/cm^2 (precomputed)
    p = np.loadtxt('data/%s/pstar.txt' % sym, skiprows=8)            # col4 = CSDA range g/cm^2
    # interpolate PSTAR CSDA range onto the eRPA energy grid (log-log)
    ps = 10 ** np.interp(np.log10(E), np.log10(p[:, 0]), np.log10(p[:, 4] * 1e3))
    m = (E >= p[:, 0].min()) & (E <= p[:, 0].max())
    ax.semilogx(E[m], (rng_e / ps)[m], 'b--', lw=1.8, label='electronic only')
    ax.semilogx(E[m], (rng_t / ps)[m], 'r-', lw=2.0, label='electronic + nuclear (GK)')
    ax.axhline(1.0, color='k', lw=0.8, ls=':')
    ax.set_title('%s' % sym); ax.grid(alpha=.3, which='both')
    ax.set_ylim(0.85, 1.8)
    ax.set_xlabel('Energy (MeV)'); ax.set_ylabel('range / PSTAR CSDA')
axes[0, 0].legend(fontsize=9, loc='upper right')
fig.suptitle('Proton range relative to PSTAR: the GK nuclear term removes the\n'
             'low-energy electronic-only overshoot (C, Al, Ag, Au)', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('range_nuc.pdf'); fig.savefig('range_nuc.png', dpi=130)
print('wrote range_nuc.pdf / .png\n')

print(' elem  E(MeV)   range_e   range_tot  PSTAR   e-err%%  tot-err%%')
for sym, Z in ELEMS:
    de = np.loadtxt('/tmp/rng_%s/dedx.dat' % sym, comments='#')
    dn = np.loadtxt('/tmp/rng_%s/dedx_nuc.dat' % sym, comments='#')
    E = de[:, 0]
    p = np.loadtxt('data/%s/pstar.txt' % sym, skiprows=8)
    for Eq in [1e-3, 1e-2]:
        ie = np.argmin(np.abs(E - Eq))
        re, rt = de[ie, 2], dn[ie, 4]
        ps = np.interp(Eq, p[:, 0], p[:, 4] * 1e3)
        print('  %-3s  %.0e   %7.4f   %7.4f   %7.4f  %+5.0f   %+5.0f' %
              (sym, Eq, re, rt, ps, 100 * (re / ps - 1), 100 * (rt / ps - 1)))
