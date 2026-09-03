# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""degeneracy_panels.py -- alpha stopping in equimolar DT at T=1 keV vs density
(electron degeneracy).  Three panels: rho = 10 (Theta~15), 1000 (Theta~0.7),
5000 g/cc (Theta~0.24).  Compares e-RPA/LDA (this work, Bloch-fixed) against the
classical Brown-Preston-Singleton (BPS) and Li-Petrasso-Zylstra (LPZ) models of
Borscz et al.  e-RPA carries Fermi degeneracy via the Maynard-Deutsch RPA; the
classical BPS does not, so e-RPA falls below it as the fuel becomes degenerate
(Pauli blocking).

Regenerate the e-RPA inputs first (H surrogate at the DT electron density):
  for pair in 10:4.01 1000:403.2 5000:2016; do DT=${pair%:*}; H=${pair#*:};
    python dedx.py --zp=2 --zt=1 --d=$H --t=1000 --aa=2 --mloss=24 \
       --emin=0.05 --emax=20 --od=Deg_$DT; done
Then:  python degeneracy_panels.py
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/Users/mehlhorn/Codes/260624 Borscz dEdx/extracted/dEdx')
from dEdx import plasma_properties, dEdx

QM = 4.0026                      # He mass number
AMU_E = 5.030 / 2.0             # DT amu per electron
CONV = 1e-21 * (6.022e23 / AMU_E) / 1e3   # dedx col1 -> MeV cm^2/mg (DT mass)
PANELS = [(10.0, r'$\Theta\approx15$'), (1000.0, r'$\Theta\approx0.7$'),
          (5000.0, r'$\Theta\approx0.24$')]

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), sharey=True)
for ax, (rho, thlab) in zip(axes, PANELS):
    # e-RPA electronic stopping (Bloch-fixed)
    d = np.loadtxt('Deg_%d/dedx.dat' % int(rho), comments='#')
    E_e = QM * d[:, 0]                       # total alpha energy [MeV]
    S_e = d[:, 1] * CONV                     # MeV cm^2/mg
    m = E_e <= 20.0
    ax.loglog(E_e[m], S_e[m], 'b-', lw=2.4, label='e-RPA/LDA (this work)')
    # BPS + LPZ (Borscz)
    pl = plasma_properties(m=[2, 3], Z=[1, 1], d=rho, Ti=1.0, Te=1.0, w=[1, 1])
    Eg = np.geomspace(0.2, 20.0, 22)
    Sb = -dEdx(pl, mp=4, Zp=2, Ep=Eg, model='BrownPrestonSingleton') / (rho * 1e3)
    El = np.geomspace(0.2, 20.0, 60)
    Sl = -dEdx(pl, mp=4, Zp=2, Ep=El, model='LiPetrassoZylstra') / (rho * 1e3)
    ax.loglog(Eg, Sb, 'r--', lw=1.8, label='Brown--Preston--Singleton')
    ax.loglog(El, Sl, 'g-.', lw=1.8, label='Li--Petrasso--Zylstra')
    ax.axvline(3.5, color='0.6', lw=0.7, ls=':')
    ax.set_title(r'DT %g g/cc  (%s)' % (rho, thlab))
    ax.set_xlabel('alpha energy (MeV)')
    ax.grid(alpha=0.3, which='both')
axes[0].set_ylabel(r'stopping power (MeV cm$^2$/mg)')
axes[0].legend(fontsize=8.5, loc='lower left')
fig.tight_layout()
fig.savefig('degeneracy_panels.png', dpi=140)
fig.savefig('degeneracy_panels.pdf')
print('wrote degeneracy_panels.png / .pdf')
# report eRPA/BPS ratio at 3.5 MeV per panel
for rho, _ in PANELS:
    d = np.loadtxt('Deg_%d/dedx.dat' % int(rho), comments='#')
    se = np.interp(3.5, QM * d[:, 0], d[:, 1] * CONV)
    pl = plasma_properties(m=[2, 3], Z=[1, 1], d=rho, Ti=1.0, Te=1.0, w=[1, 1])
    sb = float(-dEdx(pl, mp=4, Zp=2, Ep=np.array([3.5]), model='BrownPrestonSingleton')[0] / (rho * 1e3))
    print('  rho=%5g: eRPA/BPS at 3.5 MeV = %.2f' % (rho, se / sb))
