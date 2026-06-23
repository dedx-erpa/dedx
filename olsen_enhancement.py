# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
olsen_enhancement.py -- connect the e-RPA/LDA model to the Olsen 1985 data.

Olsen et al. (J. Appl. Phys. 58, 2958, 1985) measured the stopping-power
enhancement of 1.6 MeV protons in beam-heated Al and Ni foils: ~100% (Al) and
~50% (Ni) at peak intensity, with inferred average ionizations Zbar=8-10 (Al)
and 9-12 (Ni) at peak electron temperatures of 48 eV (Al) and 42 eV (Ni).

The enhancement is strongly density dependent, so we do NOT guess the expanded
foil density.  Instead we fix the temperature at the measured peak value and
scan density; the average-atom ionization balance then gives Zbar(rho), and the
density consistent with the measured Zbar is the physically correct one (Saha
constraint).  We plot the enhancement -- the reduction in the 1.6 MeV proton
range relative to cold matter -- along each measured isotherm vs Zbar, with the
measured enhancements overlaid.

Run (fac env):  python olsen_enhancement.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from multiprocessing import Pool

EPROBE = 1.6                                   # MeV proton
# (name, Z, peak T eV, density scan, colour, marker, measured: Zlo,Zhi,Elo,Ehi)
CASES = [('Al', 13, 48, [0.002, 0.005, 0.01, 0.02, 0.04, 0.08],
          'tab:blue', 'o', (8, 10, 90, 110)),
         ('Ni', 28, 42, [0.005, 0.01, 0.02, 0.04, 0.08, 0.15],
          'tab:green', 's', (9, 12, 40, 60))]
WORK = '/tmp/olsenh'


def run(job):
    name, zt, T, rho = job
    od = '%s/%s_%g_%g' % (WORK, name, T, rho)
    os.system('python dedx.py --zt=%d --zp=1 --d=%g --t=%g --aa=2 --nuc=0 '
              '--emin=0.02 --emax=3 --mep=40 --od=%s >/dev/null 2>&1'
              % (zt, rho, T, od))


def read(od):
    zb = np.nan
    for ln in open('%s/dedx.dat' % od):
        if 'zbar' in ln:
            zb = float(ln.split('=')[1]); break
    d = np.loadtxt('%s/dedx.dat' % od)
    R = np.exp(np.interp(np.log(EPROBE), np.log(d[:, 0]), np.log(d[:, 2])))
    return zb, R


if __name__ == '__main__':
    os.makedirs(WORK, exist_ok=True)
    jobs = []
    for name, zt, T, rhos, _, _, _ in CASES:
        for rho in rhos:
            jobs.append((name, zt, T, rho))
            jobs.append((name, zt, 0.025, rho))    # cold reference, same density
    with Pool(14) as p:
        p.map(run, jobs)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for name, zt, T, rhos, col, mk, meas in CASES:
        Z, E = [], []
        for rho in rhos:
            zb, Rh = read('%s/%s_%g_%g' % (WORK, name, T, rho))
            _, Rc = read('%s/%s_%g_%g' % (WORK, name, 0.025, rho))
            Z.append(zb); E.append(100.0 * (Rc / Rh - 1.0))
        o = np.argsort(Z); Z = np.array(Z)[o]; E = np.array(E)[o]
        ax.plot(Z, E, mk + '-', color=col, label='e-RPA/LDA %s, T = %g eV' % (name, T))
        zlo, zhi, elo, ehi = meas
        ax.add_patch(Rectangle((zlo, elo), zhi - zlo, ehi - elo,
                               color=col, alpha=0.18, ec=col, lw=1.2))
        ax.text((zlo + zhi) / 2, ehi + 4, 'Olsen %s (meas.)' % name,
                color=col, fontsize=8.5, ha='center', fontweight='bold')
        print('%s (T=%g eV): Zbar %.1f-%.1f, enhancement %.0f-%.0f%%'
              % (name, T, Z.min(), Z.max(), E.min(), E.max()))

    ax.set_xlabel(r'average ionization $\bar{Z}$')
    ax.set_ylabel('stopping-power enhancement (%)')
    ax.set_title('Olsen 1985: proton stopping enhancement at the measured temperature')
    ax.set_xlim(5, 13); ax.set_ylim(0, 130)
    ax.legend(loc='upper left', fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('olsen_enhancement.png', dpi=150)
    fig.savefig('olsen_enhancement.pdf')
    print('wrote olsen_enhancement.png/.pdf')
