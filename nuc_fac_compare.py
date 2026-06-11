"""
nuc_fac_compare.py -- compare the FAC-based implementation of the Faussurier
nuclear-stopping model against the original paper (which used the SCAALP
average-atom model), for the paper's canonical case: proton in aluminum.

Left panel  (cf. Faussurier Fig. 1): cold solid-density Al, the three pair
potentials vs NIST/ZBL.
Right panel (cf. Faussurier Fig. 2): Gordon-Kim nuclear stopping vs proton
energy at Te = cold, 10, 100, 1000 eV, with the FAC mean ionization Zbar(Te).

Reuses the average-atom runs in data/Al (cold) and /tmp/nucsweep_Al_* (hot,
created by `python nuc_test.py sweep`).  Run: python nuc_fac_compare.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import rd
import nuclear as N
from nuclear_gk import GordonKim

ZP, ZT, M0, MT = 1, 13, 1.008, 26.98


def stop_with(U, E, m0=M0, mt=MT):
    """Ti=0 nuclear stopping [1e-15 eV cm^2/atom] for an explicit potential U."""
    m0a, mta = m0 * N.AMU_ME, mt * N.AMU_ME
    Ela = np.atleast_1d(E) * m0 * 1e6 / N.HARTREE_EV
    ec = (mta / (m0a + mta)) * Ela
    Sf = N.make_S_interp(U, ec.min(), ec.max())
    return np.array([N.eps_n_Ti0(Ela[i], m0a, mta, Sf) * N.EPS_UNIT
                     for i in range(len(Ela))])


def zbl(E, Z1=1, Z2=13, M1=1.008, M2=26.98):
    Ek = np.atleast_1d(E) * 1e3
    a = Z1**0.23 + Z2**0.23
    eps = 32.53 * M2 * Ek / (Z1 * Z2 * (M1 + M2) * a)
    sn = np.where(eps <= 30, np.log1p(1.1383 * eps) /
                  (2 * (eps + 0.01321 * eps**0.21226 + 0.19593 * np.sqrt(eps))),
                  np.log(eps) / (2 * eps))
    return 8.462 * Z1 * Z2 * M1 * sn / ((M1 + M2) * a)


def main():
    E = np.geomspace(1e-3, 2.0, 60)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

    # ---- Left: cold Al, three potentials vs NIST/ZBL (Fig 1) ----
    h = rd.rdedx('data/Al', header='')
    rs, zbar = float(h['rs']), float(h['zbar'])
    axL.loglog(E, zbl(E), 'k-', lw=2.4, label='NIST / ZBL')
    ne = zbar * N.Nion(rs)
    for p, s, lab in [('ionsphere', 'b--', 'ion sphere'),
                      ('yukawa', 'g-.', 'Yukawa')]:
        y = np.atleast_1d(N.nuclear_stopping(E, ZP, ZT, M0, MT, rs, 0.025, 0.0,
                                             zbar, potential=p, od='data/Al'))
        axL.loglog(E, y, s, lw=1.8, label=lab)
    # Gordon-Kim: electrostatic-only vs full (electrostatic + kinetic/exch/corr)
    gk_es = GordonKim(ZP, zbar, ne, od='data/Al', rs=rs, full=False)
    gk_full = GordonKim(ZP, zbar, ne, od='data/Al', rs=rs, full=True, scale=True)
    axL.loglog(E, stop_with(gk_es, E), 'r:', lw=1.6, label='GK electrostatic only')
    axL.loglog(E, stop_with(gk_full, E), 'r-', lw=2.2,
               label='GK full (+ kinetic/exch/corr)')
    axL.set_xlabel('proton energy (MeV)')
    axL.set_ylabel(r'nuclear stopping ($10^{-15}$ eV cm$^2$/atom)')
    axL.set_title('Cold Al: FAC implementation vs NIST (cf. Fig. 1)')
    axL.set_ylim(1e-3, 1)
    axL.legend(); axL.grid(True, which='both', alpha=0.3)

    # ---- Right: GK nuclear vs E at several Te (Fig 2) ----
    cases = [('cold (0.025 eV)', 0.025, 'data/Al')]
    for t in (10, 100, 1000):
        od = '/tmp/nucsweep_Al_%g' % t
        if os.path.exists(od + '/dedx.dat'):
            cases.append(('%g eV' % t, float(t), od))
    for name, te, od in cases:
        hh = rd.rdedx(od, header='')
        rs2, zb2 = float(hh['rs']), float(hh['zbar'])
        y = np.atleast_1d(N.nuclear_stopping(E, ZP, ZT, M0, MT, rs2, te, 0.0,
                                             zb2, potential='gk', od=od))
        axR.loglog(E, y, lw=2, label=r'%s, $\bar Z$=%.1f' % (name, zb2))
    axR.set_xlabel('proton energy (MeV)')
    axR.set_ylabel(r'GK nuclear stopping ($10^{-15}$ eV cm$^2$/atom)')
    axR.set_title('Al: nuclear growth with Te (FAC Zbar; cf. Fig. 2)')
    axR.legend(); axR.grid(True, which='both', alpha=0.3)

    fig.tight_layout()
    fig.savefig('nuc_fac_compare.pdf')
    fig.savefig('nuc_fac_compare.png', dpi=130)
    print('wrote nuc_fac_compare.pdf / .png')


if __name__ == '__main__':
    main()
