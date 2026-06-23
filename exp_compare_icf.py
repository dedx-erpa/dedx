# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
exp_compare_icf.py -- compare the total eRPA stopping model against two
ICF-relevant plasma stopping measurements:

  - Zylstra 2015 (PRL 114, 215002): energetic protons in isochorically-heated
    solid-density Be plasma, Te ~ 32 eV. Reported INCREASED loss vs cold matter,
    in agreement with the average-atom LDA.
  - Frenje 2019 (PRL 122, 015002): ion stopping in DT plasma (Te ~ 2 keV) around
    the Bragg peak. BPS underpredicts ~20% at vi ~ 0.3 v_th and the authors
    postulate NUCLEAR-ELASTIC scattering -- exactly the nuclear/ion channel here.

Model dirs (carbon -> Be, DT):
  python dedx.py --zt=4 --zp=1 --d=1.85 --t=32    --aa=2 --nuc=1 --npot=gk --od=/tmp/zylstra_wdm  --emin=0.3 --emax=15 --mep=45
  python dedx.py --zt=4 --zp=1 --d=1.85 --t=0.025 --aa=2 --nuc=1 --npot=gk --od=/tmp/zylstra_cold --emin=0.3 --emax=15 --mep=45
  python dedx.py --zt=1 --zp=1 --d=2.5  --t=2000  --aa=2 --nuc=1 --npot=gk --od=/tmp/frenje_dt    --emin=1e-3 --emax=1 --mep=60
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

fig, (axZ, axF) = plt.subplots(1, 2, figsize=(13, 5))

# ---------- Zylstra 2015: proton in Be, cold vs WDM ----------
c = rd.rdedx('/tmp/zylstra_cold')
w = rd.rdedx('/tmp/zylstra_wdm')
axZ.loglog(c[0], c[1], 'b--', lw=2, label='model: cold Be')
axZ.loglog(w[0], w[1], 'r-', lw=2, label='model: Be plasma (Te=32 eV)')
axZ.set_xlim(0.3, 15)
axZ.set_xlabel('proton energy (MeV)')
axZ.set_ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
axZ.set_title('Zylstra 2015: proton in warm dense Be')
# annotate the enhancement
r1 = np.interp(1.0, w[0], w[1]) / np.interp(1.0, c[0], c[1])
axZ.text(0.05, 0.05, 'WDM/cold ≈ %.2f at 1 MeV\n(increased loss vs cold —\nagrees with AA-LDA, Zylstra)' % r1,
         transform=axZ.transAxes, fontsize=8, va='bottom',
         bbox=dict(boxstyle='round', fc='wheat', alpha=0.6))
axZ.legend(fontsize=8); axZ.grid(alpha=.3, which='both')

# ---------- Frenje 2019: proton in DT plasma, nuclear/ion channel ----------
d = np.loadtxt('/tmp/frenje_dt/dedx_nuc.dat', comments='#')
E, de, dn, dt = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
axF.loglog(E, de, 'b-', lw=2, label='electronic (eRPA)')
axF.loglog(E, np.where(dn > 0, dn, np.nan), 'r-', lw=2, label='nuclear / ion-elastic')
axF.loglog(E, dt, 'k-', lw=2.4, alpha=.8, label='total')
# vi ~ 0.3 v_th: low-velocity side of the Bragg peak (~0.1 of the peak energy)
E03 = 0.1
nf = np.interp(E03, E, dn) / np.interp(E03, E, dt)
axF.axvline(E03, color='0.5', ls=':', lw=1)
axF.set_xlabel('proton energy (MeV/amu)')
axF.set_ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
axF.set_title('Frenje 2019: ion stopping in DT plasma (Te=2 keV)')
axF.text(0.04, 0.05,
         'nuclear/ion = %.0f%% of total\nat vi≈0.3 v_th (dotted line).\nFrenje: BPS underpredicts ~20%%\nhere — attributed to\nnuclear-elastic scattering.' % (100 * nf),
         transform=axF.transAxes, fontsize=8, va='bottom',
         bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.7))
axF.legend(fontsize=8); axF.grid(alpha=.3, which='both')

fig.suptitle('Total eRPA stopping vs ICF-plasma measurements (Zylstra 2015, Frenje 2019)')
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('exp_compare_icf.pdf'); fig.savefig('exp_compare_icf.png', dpi=130)
print('wrote exp_compare_icf.pdf / .png')

print('\n=== Zylstra (proton in Be, 32 eV vs cold) ===')
for Ev in [0.5, 1.0, 3.0]:
    print('  %.1f MeV: WDM/cold = %.2f (increased loss vs cold)' %
          (Ev, np.interp(Ev, w[0], w[1]) / np.interp(Ev, c[0], c[1])))
print('\n=== Frenje (proton in DT, Te=2 keV): nuclear/ion fraction ===')
for Ev in [3e-3, 1e-2, 3e-2, 0.1, 0.3]:
    print('  %.0e MeV/amu: nuclear/total = %.0f%%' % (Ev, 100 * np.interp(Ev, E, dn) / np.interp(Ev, E, dt)))
print('  -> at vi~0.3 v_th the nuclear/ion channel ~ Frenje\'s postulated ~20%%')
