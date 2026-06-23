# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
nuc_ion_verify.py -- verify the finite-Ti ion-ion stopping of our model against
the Fokker-Planck slowing-down (the projectile-ion / target-ion drag), i.e. the
"ion stopping" that the low-T "nuclear" stopping becomes at finite temperature
(cf. Mehlhorn, J. Appl. Phys. 52, 6522, 1981; Faussurier 2013 Eq. 6).

Two checks:
  1. THRESHOLD (energy gain -> loss).  The analytic Fokker-Planck single-particle
     threshold solves   erf(x)/(x erf'(x)) = 1 + m2/m1 ,  x=v/v_th2,
     v_th2=sqrt(2 kT2/m2).  For He on a fully-ionized H plasma at Te=Ti=100 eV
     this gives E* = 131 eV, matching Faussurier's MD-validated 132.7 eV.
  2. FULL CURVE.  The analytic FP ion-ion stopping
        -dE/dx = (1/2)(4 pi Z1^2 Z2^2 e^4 n2 lnL / (m2 v^2)) [erf(x) - (1+m2/m1) x erf'(x)]
     vs our Maxwellian model (nuclear.eps_n_Maxwell) using (a) the same Coulomb-log
     cross section -- which isolates and validates our slowing-down machinery --
     and (b) the physical screened (Yukawa) interaction.

Run (fac env):  python nuc_ion_verify.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import brentq

import nuclear as N

ERFP = lambda x: (2.0 / np.sqrt(np.pi)) * np.exp(-x * x)     # erf'(x)


def fp_threshold(m1, m2, T2_eV):
    """Analytic FP single-particle energy-gain threshold E* [eV]."""
    rhs = 1.0 + m2 / m1
    x = brentq(lambda x: erf(x) / (x * ERFP(x)) - rhs, 1e-4, 8.0)
    return (m1 / m2) * x**2 * T2_eV, x


def coulomb_lnL(zp, zt, m1, m2, v, ne, T_eV):
    """Velocity-dependent ION-ION Coulomb logarithm ln(lambda_D / b_min).
    b_90 uses the ion-ion relative velocity (projectile + target-ion thermal),
    NOT the electron thermal velocity."""
    m1a, m2a = m1 * N.AMU_ME, m2 * N.AMU_ME
    mu = m1a * m2a / (m1a + m2a)
    vth2 = 2.0 * (T_eV / N.HARTREE_EV) / m2a          # field-ion thermal speed^2
    lamD = np.sqrt((T_eV / N.HARTREE_EV) / (4.0 * np.pi * ne))   # Debye length
    vrel2 = v**2 + vth2                                # ion-ion relative velocity^2
    b90 = zp * zt / (mu * vrel2)                       # classical 90-deg impact param
    return np.maximum(1.0, np.log(lamD / b90))


def fp_dedx(E_eV, zp, zt, m1, m2, T2_eV, ne):
    """Analytic FP ion-ion stopping cross section [1e-15 eV cm^2/atom], with a
    velocity-dependent Coulomb log."""
    m1a, m2a = m1 * N.AMU_ME, m2 * N.AMU_ME
    E = np.atleast_1d(E_eV) / N.HARTREE_EV
    v = np.sqrt(2.0 * E / m1a)
    vth2 = np.sqrt(2.0 * (T2_eV / N.HARTREE_EV) / m2a)
    x = v / vth2
    lnL = coulomb_lnL(zp, zt, m1, m2, v, ne, T2_eV)
    eps_ha = 0.5 * 4.0 * np.pi * zp**2 * zt**2 * lnL / (m2a * v**2) \
        * (erf(x) - (1.0 + m2a / m1a) * x * ERFP(x))
    return eps_ha * N.EPS_UNIT


def model_dedx(E_eV, zp, zt, m1, m2, T2_eV, U, ti_eV=None, Scoul=None):
    """Our Maxwellian model ion stopping [1e-15 eV cm^2/atom] over E_eV.
    If Scoul is given, use that Coulomb-log S(Ec) instead of the potential U."""
    if ti_eV is None:
        ti_eV = T2_eV
    m1a, m2a = m1 * N.AMU_ME, m2 * N.AMU_ME
    ti = ti_eV / N.HARTREE_EV
    mu = m1a * m2a / (m1a + m2a)
    out = np.empty(len(E_eV))
    for i, Eev in enumerate(E_eV):
        El = Eev / N.HARTREE_EV
        V0 = np.sqrt(2.0 * El / m1a)
        wth = np.sqrt(ti / m2a)
        g0 = max(1e-6, abs(V0 - 5 * wth))
        g1 = np.sqrt((5 * wth)**2 + (V0 + 5 * wth)**2)
        Sf = Scoul if Scoul is not None else \
            N.make_S_interp(U, 0.5 * mu * g0**2, 0.5 * mu * g1**2)
        out[i] = N.eps_n_Maxwell(El, m1a, m2a, Sf, ti) * N.EPS_UNIT
    return out


def main():
    import rd
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # ---- He on fully-ionized H plasma, Te=Ti=100 eV (paper E*=132.7 eV) ----
    ni = 1e24 * (0.529177e-8)**3
    zbar = 1.0
    ne = zbar * ni
    E = np.geomspace(20, 3000, 50)
    U = N.build_potential('yukawa', 2, zbar, ne, 100.0 / N.HARTREE_EV)
    ax = axes[0]
    ax.semilogx(E, model_dedx(E, 2, 1, 4.0, 1.008, 100.0, U),
                'g-', lw=2, label='our model (screened)')
    ax.semilogx(E, fp_dedx(E, 2, 1, 4.0, 1.008, 100.0, ne), 'k:', lw=2,
                label='analytic Fokker-Planck')
    ax.axhline(0, color='0.6', lw=.7)
    ax.axvline(132.7, color='r', ls='--', lw=1, label='paper E*=132.7 eV')
    ax.set_title('He on H plasma, Te=Ti=100 eV')
    ax.set_xlabel('He energy (eV)'); ax.set_ylabel(r'ion stopping ($10^{-15}$ eV cm$^2$/atom)')
    ax.set_ylim(-6, 6)
    ax.legend(fontsize=8); ax.grid(alpha=.3, which='both')

    # ---- p on ionized Al at Te=Ti = 100 and 1000 eV ----
    for ax, (te, od) in zip(axes[1:], [(100, '/tmp/nucsweep_Al_100'),
                                       (1000, '/tmp/nucsweep_Al_1000')]):
        h = rd.rdedx(od, header=''); rsA = float(h['rs']); zb = float(h['zbar'])
        neA = zb * N.Nion(rsA)
        Uy = N.build_potential('yukawa', 1, zb, neA, te / N.HARTREE_EV)
        Ep = np.geomspace(20, 5000, 44)
        ym = model_dedx(Ep, 1, zb, 1.008, 26.98, te, Uy)
        ya = fp_dedx(Ep, 1, zb, 1.008, 26.98, te, neA)
        ax.semilogx(Ep, ym, 'g-', lw=2, label='our model (screened)')
        ax.semilogx(Ep, ya, 'k:', lw=2, label='analytic Fokker-Planck')
        ax.axhline(0, color='0.6', lw=.7)
        ax.set_title(r'p on Al, Te=Ti=%g eV ($\bar Z$=%.1f)' % (te, zb))
        ax.set_xlabel('proton energy (eV)')
        top = max(ym.max(), ya.max()) * 1.3
        ax.set_ylim(-top, top)
        ax.legend(fontsize=8); ax.grid(alpha=.3, which='both')

    fig.suptitle('Ion stopping: our finite-Ti model vs Fokker-Planck projectile-ion/target-ion drag')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('nuc_ion_verify.pdf'); fig.savefig('nuc_ion_verify.png', dpi=130)
    print('wrote nuc_ion_verify.pdf / .png')

    # printed threshold summary
    print('\nthreshold (energy gain -> loss):')
    for lbl, m1, m2, T in [('He-on-H', 4.0, 1.008, 100.0),
                           ('p-on-Al', 1.008, 26.98, 100.0),
                           ('p-on-Al', 1.008, 26.98, 1000.0)]:
        e, x = fp_threshold(m1, m2, T)
        print('  %-8s T=%5g eV:  analytic FP E*=%.0f eV (x*=%.2f)' % (lbl, T, e, x))


if __name__ == '__main__':
    main()
