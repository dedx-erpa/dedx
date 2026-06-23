# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
exp_compare.py -- compare the total eRPA stopping model against recent
warm-dense-matter / plasma stopping-power measurements.

Datasets (see docs/recent_stopping_experiments.md):
  - Malko 2022 (Nat. Commun. 13, 2893): ~0.5 MeV protons in warm dense carbon,
    Te ~ 7.5 eV, Z* <= 4.  Data digitized in malko_fig.txt.
  - Cayzac 2017 (Nat. Commun. 8, 15693): 0.586 MeV/u N ions in fully-ionized
    carbon plasma, ne ~ 5e20, Te ~ 150 eV.  Headline numbers only.

Generate the model dirs first (carbon target; nuclear/GK included):
  python dedx.py --zt=6 --zp=1 --d=2.0    --t=0.025 --aa=2 --nuc=1 --npot=gk --od=/tmp/malko_solid   --emin=5e-2 --emax=2 --mep=40
  python dedx.py --zt=6 --zp=1 --d=2.0    --t=7.5   --aa=2 --nuc=1 --npot=gk --od=/tmp/malko_wdm     --emin=5e-2 --emax=2 --mep=40
  python dedx.py --zt=6 --zp=7 --d=1.7e-3 --t=150   --aa=2 --nuc=1 --npot=gk --od=/tmp/cayzac_plasma --emin=0.2  --emax=2 --mep=40
  python dedx.py --zt=6 --zp=7 --d=2.0    --t=0.025 --aa=2 --nuc=1 --npot=gk --od=/tmp/cayzac_solid  --emin=0.2  --emax=2 --mep=40
Then: python exp_compare.py  -> exp_compare.pdf/png and a printed comparison table.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

# 10^-15 eV cm^2/atom  ->  keV/(ug/cm^2) for a carbon target (A=12)
A_C = 12.0
CONV = 1e-18 / (A_C * 1.67e-24 * 1e6)


def tot(od):
    """E/AMU [MeV], total stopping [10^-15 eV cm^2/atom] from dedx_nuc.dat."""
    d = np.loadtxt('%s/dedx_nuc.dat' % od, comments='#')
    return d[:, 0], d[:, 3]          # E, dEdx_tot (electronic + GK nuclear)


def at(E0, E, y):
    return np.interp(E0, E, y)


fig, (axM, axC) = plt.subplots(1, 2, figsize=(13, 5))

# ---------- Malko 2022: proton in warm dense carbon ----------
Es, ys = tot('/tmp/malko_solid'); ys = ys * CONV
Ew, yw = tot('/tmp/malko_wdm');   yw = yw * CONV
axM.plot(Es, ys, 'b--', lw=2, label='model: solid C (Te=0.025 eV)')
axM.plot(Ew, yw, 'r-', lw=2, label='model: warm dense C (Te=7.5 eV)')
# Malko data points (digitized): d[0]=keV, paired value/value+err
d = np.loadtxt('malko_fig.txt', unpack=1)
mx = d[0][::2] / 1e3
mv = d[1][::2]
merr = d[1][1::2] - d[1][::2]
axM.errorbar(mx, mv, yerr=merr, fmt='ko', capsize=3, ms=5, label='Malko 2022 (expt/DFT)')
axM.set_xlim(0, 1.0); axM.set_ylim(0, 1.2)
axM.set_xlabel('proton energy (MeV)')
axM.set_ylabel(r'dE/dx (keV/($\mu$g/cm$^2$))')
axM.set_title('Malko 2022: proton in warm dense carbon')
axM.legend(fontsize=8); axM.grid(alpha=.3)

# ---------- Cayzac 2017: N ions in carbon plasma ----------
Ep, yp = tot('/tmp/cayzac_plasma')          # N in C plasma, Te=150 eV
Es2, ys2 = tot('/tmp/cayzac_solid')         # N in solid C
axC.loglog(Ep, yp, 'r-', lw=2, label='model: C plasma (Te=150 eV)')
axC.loglog(Es2, ys2, 'b--', lw=2, label='model: solid C')
E_N = 0.586                                  # MeV/u (probe energy)
rp = at(E_N, Ep, yp); rs = at(E_N, Es2, ys2)
axC.axvline(E_N, color='0.5', ls=':', label='0.586 MeV/u (N probe)')
axC.set_xlabel('N energy (MeV/u)')
axC.set_ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
axC.set_title('Cayzac 2017: N ions in carbon plasma')
axC.legend(fontsize=8); axC.grid(alpha=.3, which='both')

fig.suptitle('Total eRPA stopping (electronic + GK nuclear) vs recent WDM/plasma measurements')
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('exp_compare.pdf'); fig.savefig('exp_compare.png', dpi=130)
print('wrote exp_compare.pdf / .png')

# ---------- printed comparison ----------
print('\n=== Malko 2022 (proton in warm dense C, Te~7.5 eV) ===')
print(' proton E      model solid   model WDM    WDM/solid    Malko pt')
for E0 in [0.226, 0.401, 0.506, 0.625]:
    ms = at(E0, Es, ys); mw = at(E0, Ew, yw)
    md = at(E0, mx, mv) if E0 >= mx.min() and E0 <= mx.max() else float('nan')
    print('  %.3f MeV   %.3f        %.3f       %.2f         %.3f' % (E0, ms, mw, mw / ms, md))
print(' [units keV/(ug/cm2)]  measured ratio WDM/solid (Malko) ~ 0.74-0.87')

print('\n=== Cayzac 2017 (N 0.586 MeV/u, C plasma Te=150 eV vs solid) ===')
print('  model dE/dx plasma = %.3e, solid = %.3e  ->  plasma/solid = %.2f' %
      (rp, rs, rp / rs))
print('  measured: plasma energy loss enhanced up to ~1.5x vs solid')
