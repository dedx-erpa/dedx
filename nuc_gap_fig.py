"""
nuc_gap_fig.py -- the 4-panel Levine-Louie band-gap comparison (LiF, SiO2,
Al2O3, H2O) with the nuclear (ion-ion) stopping power added in.

For each insulator it plots:
  - gapless (kappa=0) electronic + GK nuclear          (blue)
  - kappa=0.9  electronic + GK nuclear                 (red, solid)
  - kappa=0.9  electronic only                         (red, dotted; shows the
                                                         nuclear increment)
  - PSTAR                                              (black dashed)
  - experiment (IAEA)                                  (open circles)

The average-atom density is gap-independent, so the nuclear stopping (computed
once per material with the Gordon-Kim potential, number-weighted over the
constituent atoms) is identical for the gapless and kappa=0.9 curves.

Run in the `fac` conda env:  python nuc_gap_fig.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pfac import fac, aa

import rd
import nuclear as N

MATERIALS = [('LiF',   2.64, 14.0),
             ('SiO2',  2.65,  9.0),
             ('Al2O3', 3.97,  8.8),
             ('H2O',   1.00,  7.0)]   # (formula, density g/cc, E_g eV for title)
KAPPA = 0.7                           # best electronic calibration (eRPA band gap)
TE = 0.025                            # cold


def kappa_electronic(m, d, t=TE):
    """Run (or reuse) the gap-corrected electronic stopping at the best kappa."""
    od = '/tmp/gapk%02d_%s' % (round(KAPPA * 10), m)
    if not os.path.exists('%s/dedx.dat' % od):
        cmd = ('python dedx.py --fc=%s --d=%g --t=%g --aa=2 '
               '--mloss=324 --kgap=%g --od=%s' % (m, d, t, KAPPA, od))
        print('running:', cmd)
        os.system(cmd)
    y = rd.rdedx(od)
    return y[0], y[1]                  # E/AMU [MeV], dedx_e [1e-15 eV cm^2/atom]


def nuclear_compound(m, E, te=TE):
    """Number-weighted Gordon-Kim nuclear stopping for compound `m` on grid E,
    using the per-component average-atom densities under data/<m>/<sym>."""
    zc, wc = aa.zw4c(m)
    eps = np.zeros_like(E)
    wtot = 0.0
    for z, w in zip(zc, wc):
        sod = 'data/%s/%s' % (m, fac.ATOMICSYMBOL[z])
        h = rd.rdedx(sod, header='')
        e = np.atleast_1d(N.nuclear_stopping(E, 1, z, fac.ATOMICMASS[1],
                                             fac.ATOMICMASS[z], float(h['rs']),
                                             te, 0.0, float(h['zbar']),
                                             potential='gk', od=sod))
        eps += e * w
        wtot += w
    return eps / wtot


def read_exp(m):
    """IAEA experimental points (replicating rd_mod.cmp_dedx_mod's reader)."""
    for fd in [m + '.txt', m.lower() + '.txt', m.upper() + '.txt']:
        p = 'data/iaea/' + fd
        if not os.path.exists(p):
            continue
        try:
            r = np.loadtxt(p, unpack=1, usecols=(0, 1), dtype=str)
            x = np.array([rd.tofloat(v) for v in r[0]])
            y = np.array([rd.tofloat(v) for v in r[1]])
            if (x[0] == 0 and y[0] == 0 and x[1] == 0) or max(y) < 0.1:
                y = y * fac.ATOMICMASS[8] * 1.67   # rough; rarely hit for these
            else:
                x = x / 1e3                         # keV -> MeV
            return x, y
        except Exception:
            pass
    return None, None


def main():
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(r'Total proton stopping: eRPA electronic ($\kappa=%.1f$ band gap) '
                 r'+ Gordon-Kim nuclear vs PSTAR' % KAPPA, fontsize=14)

    for ax, (m, d, eg) in zip(axes.flat, MATERIALS):
        # best-kappa electronic
        Ek, e_k = kappa_electronic(m, d)
        nuc_k = nuclear_compound(m, Ek)

        ax.semilogx(Ek, e_k + nuc_k, '-', color='tab:red', lw=2.4,
                    label=r'$\kappa$=%.1f electronic + nuclear (total)' % KAPPA)
        ax.semilogx(Ek, e_k, ':', color='tab:red', lw=1.5,
                    label=r'$\kappa$=%.1f electronic only' % KAPPA)
        ax.semilogx(Ek, nuc_k, '-', color='tab:green', lw=1.3, alpha=0.8,
                    label='nuclear (Gordon-Kim)')

        if os.path.exists('data/%s/pstar.txt' % m):
            r = rd.rpstar('data/%s' % m)
            ax.semilogx(r[0], r[1], 'k--', lw=2, label='PSTAR')

        xe, ye = read_exp(m)
        if xe is not None:
            ax.semilogx(xe, ye, 'o', mfc='none', mec='0.4', ms=4,
                        mew=0.6, ls='none', label='Exp.')

        ax.set_xlim(1e-3, 2.0)
        ax.set_xlabel('Proton E (MeV)')
        ax.set_ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
        ax.set_title('%s (E_g~%g eV)' % (m, eg))
        ax.legend(fontsize=9)
        ax.grid(True, which='both', alpha=0.25)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('nuc_gap_4panel.pdf')
    fig.savefig('nuc_gap_4panel.png', dpi=130)
    print('wrote nuc_gap_4panel.pdf / .png')


if __name__ == '__main__':
    main()
