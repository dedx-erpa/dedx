# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
alpha_dt_verify.py -- verify alpha-particle (3.5 MeV He-4) stopping and range in
a DT plasma, and the electron/ion deposition split, vs temperature.

This is the canonical ICF self-heating calculation (Fraley et al., Phys. Fluids
17, 474, 1974; Mehlhorn, J. Appl. Phys. 52, 6522, 1981).  Electronic and ion
stopping are the same Fokker-Planck/Chandrasekhar drag with the field being the
electrons or the (D,T) ions respectively -- the framework verified in
nuc_ion_verify.py (and reproduced by our screened nuclear model).

Outputs:
  - alpha range rho_R (g/cm^2) vs Te for several DT densities  (cf. left panel)
  - ion/electron deposition fraction vs Te, crossover near 30 keV (right panel)

Run (fac env):  python alpha_dt_verify.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erf

HARTREE_EV = 27.211386
BOHR_CM = 0.529177210903e-8
AMU_ME = 1822.888486
NA = 6.022140760e23
ERFP = lambda x: (2.0 / np.sqrt(np.pi)) * np.exp(-x * x)

# alpha projectile
E0_eV = 3.52e6
ZP, MP = 2.0, 4.0026
# DT field: 50/50 D and T (both Z=1), fully ionized
DT = [(1.0, 2.0141), (1.0, 3.0160)]      # (Z, mass amu) for D, T
MDT = 0.5 * (2.0141 + 3.0160)            # mean ion mass amu


def coulomb_lnL(zp, zf, mp, mf, v, ne_au, Te_eV, Tf_eV):
    """Velocity-dependent Coulomb log ln(lambda_D/b_min) for projectile-field."""
    mpa, mfa = mp * AMU_ME, mf * AMU_ME
    mu = mpa * mfa / (mpa + mfa)
    vthf2 = 2.0 * (Tf_eV / HARTREE_EV) / mfa            # field thermal speed^2
    lamD = np.sqrt((Te_eV / HARTREE_EV) / (4.0 * np.pi * ne_au))  # electron Debye
    vrel2 = v**2 + vthf2
    bcl = abs(zp * zf) / (mu * vrel2)                    # classical 90-deg
    bqm = 1.0 / (2.0 * mu * np.sqrt(vrel2))              # de Broglie /2
    bmin = np.maximum(bcl, bqm)
    return np.maximum(1.0, np.log(lamD / bmin))


def S_field(E_eV, zp, mp, zf, mf, nf_au, Tf_eV, ne_au, Te_eV):
    """Fokker-Planck stopping cross section on one field species
    [Hartree/Bohr], i.e. -dE/dx contribution."""
    mpa, mfa = mp * AMU_ME, mf * AMU_ME
    E = np.atleast_1d(E_eV) / HARTREE_EV
    v = np.sqrt(2.0 * E / mpa)
    vthf = np.sqrt(2.0 * (Tf_eV / HARTREE_EV) / mfa)
    x = v / vthf
    lnL = coulomb_lnL(zp, zf, mp, mf, v, ne_au, Te_eV, Tf_eV)
    # Chandrasekhar friction: 4 pi zp^2 zf^2 nf lnL / (mf v^2) * [erf(x) - (1+mf/mp) x erf'(x)].
    # (A spurious 1/2 here previously halved the stopping and doubled the range --
    #  e.g. the 3.5 MeV alpha in DT came out ~0.89 instead of ~0.44 g/cm^2, matching BPS.)
    return 4.0 * np.pi * zp**2 * zf**2 * nf_au * lnL / (mfa * v**2) \
        * (erf(x) - (1.0 + mfa / mpa) * x * ERFP(x))


def stopping(E_eV, rho_gcc, Te_eV):
    """Total, electronic, ion stopping [Hartree/Bohr] for the alpha at energy E."""
    ni_cc = rho_gcc * NA / MDT                  # total ion density 1/cc
    ne_cc = ni_cc                               # Z=1 -> ne = ni
    cc2bohr3 = BOHR_CM**3
    ne = ne_cc * cc2bohr3
    Se = S_field(E_eV, ZP, MP, 1.0, 1.0 / AMU_ME, ne, Te_eV, ne, Te_eV)  # electrons
    Si = np.zeros_like(np.atleast_1d(E_eV), dtype=float)
    for zf, mf in DT:
        nf = 0.5 * ni_cc * cc2bohr3             # half D, half T
        Si = Si + S_field(E_eV, ZP, MP, zf, mf, nf, Te_eV, ne, Te_eV)
    return np.atleast_1d(Se), np.atleast_1d(Si)    # signed (ion term -> gain near thermal)


def integrate_alpha(rho_gcc, Te_eV, nE=600):
    """Slow the alpha from E0 down to thermalization (where the net stopping
    first vanishes, ~ (3/2)kTe).  Return rho_R [g/cm^2] and electron fraction."""
    E = np.linspace(E0_eV, 0.3 * Te_eV, nE)      # down toward thermal
    Se, Si = stopping(E, rho_gcc, Te_eV)         # signed
    Stot = Se + Si
    # slow only while net stopping is positive; stop at thermalization crossing
    bad = np.where(Stot <= 0.0)[0]
    end = bad[0] if len(bad) else nE
    E, Se, Si, Stot = E[:end], Se[:end], Si[:end], Stot[:end]
    if len(E) < 3:
        return np.nan, np.nan
    dE_ha = -np.gradient(E) / HARTREE_EV
    rng_bohr = np.sum(dE_ha / Stot)
    rho_R = rho_gcc * rng_bohr * BOHR_CM         # g/cm^2
    Edep = np.sum(dE_ha)                          # energy actually deposited
    efrac = np.sum((Se / Stot) * dE_ha) / Edep
    return rho_R, efrac


def main():
    Te = np.geomspace(1.0, 300.0, 30) * 1e3     # eV  (1 - 300 keV)
    dens = [0.213, 10.0, 50.0, 100.0]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
    for rho in dens:
        rr = np.array([integrate_alpha(rho, te)[0] for te in Te])
        axL.loglog(Te / 1e3, rr, label='DT %g g/cc' % rho)
    axL.axhline(0.5, color='r', ls=':', lw=1, label=r'~0.5 g/cm$^2$')
    axL.set_xlabel('Te (keV)'); axL.set_ylabel(r'$\rho_R$ (g/cm$^2$)')
    axL.set_title('Alpha range in DT'); axL.legend(fontsize=8)
    axL.grid(True, which='both', alpha=.3)

    # deposition fractions vs Te (use a representative density)
    for rho, c in [(0.213, 'C0'), (100.0, 'C3')]:
        ef = np.array([integrate_alpha(rho, te)[1] for te in Te])
        axR.semilogx(Te / 1e3, ef, c + '-', label='efract %g g/cc' % rho)
        axR.semilogx(Te / 1e3, 1 - ef, c + '--', label='ifract %g g/cc' % rho)
    axR.axvline(30, color='k', ls=':', lw=1, label='30 keV')
    axR.axhline(0.5, color='0.6', lw=.6)
    axR.set_xlabel('Te (keV)'); axR.set_ylabel('Fraction')
    axR.set_title('Alpha deposition: electron vs ion'); axR.legend(fontsize=8)
    axR.grid(True, which='both', alpha=.3); axR.set_ylim(0, 1)

    fig.suptitle('Alpha (3.5 MeV He-4) stopping in DT plasma -- our model')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('alpha_dt_verify.pdf'); fig.savefig('alpha_dt_verify.png', dpi=130)
    print('wrote alpha_dt_verify.pdf / .png')

    # crossover temperature (efract = 0.5)
    print('\n electron/ion deposition crossover (efract=ifract=0.5):')
    for rho in dens:
        ef = np.array([integrate_alpha(rho, te)[1] for te in Te])
        i = np.where(np.diff(np.sign(ef - 0.5)))[0]
        if len(i):
            j = i[0]
            tc = np.interp(0.5, [ef[j + 1], ef[j]], [Te[j + 1] / 1e3, Te[j] / 1e3])
            rr = integrate_alpha(rho, 30e3)[0]
            print('  DT %6.3f g/cc:  crossover Te = %5.1f keV, rho_R(30keV)=%.3f g/cm2'
                  % (rho, tc, rr))


if __name__ == '__main__':
    main()
