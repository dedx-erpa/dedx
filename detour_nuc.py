# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
detour_nuc.py -- the nuclear scattering model also gives the pathlength-vs-
projected-range detour (and the end-of-range angular straggling).

The deflection integral the model already evaluates for nuclear stopping,
S(Ec) = int sin^2(theta/2) p dp, is exactly the transport (momentum-transfer)
cross section up to a constant:  sigma_tr = int (1-cos theta) dsigma = 4 pi S.
sigma_tr drives the projectile's multiple-angle scattering, hence the difference
between the CSDA pathlength and the projected (practical) range.

Multiple-scattering / diffusion estimate (proton on heavy Al: theta_lab ~ theta_cm):
  d<theta^2>/ds = 2 N sigma_tr(E)
  R - R_proj    = (1/2) int_0^R <theta^2>(s) ds
              = (N_A/A) int_0^E0 (1/S(E)) [ int_E^E0 sigma_tr/S dE'' ] dE   [areal units]
Detour factor = R_proj / R_CSDA, compared here with PSTAR for protons in Al.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, interpolate
import nuclear as nuc

NA, A, Z2, ZBAR, RS = 6.022e23, 26.98, 13.0, 2.9748752, 2.99005     # cold Al
BOHR_CM = 0.529177210903e-8
HARTREE_EV = 27.211386245988
AMU_ME = 1822.888486209
CONV = NA / A * 1e-21                                                # 1e-15 eVcm2/atom -> MeV/(g/cm2)

# total stopping S(E) [MeV/(g/cm2)] and CSDA range (proton in Al, e+nuclear)
d = np.loadtxt('/tmp/rng_Al/dedx_nuc.dat', comments='#')
E = d[:, 0]
S = d[:, 3] * CONV
S_of_E = interpolate.interp1d(E, S, bounds_error=False, fill_value=(S[0], S[-1]))
R = integrate.cumulative_trapezoid(1.0 / S, E, initial=0.0)          # g/cm2

# nuclear transport cross section sigma_tr(E) = 4 pi S(Ec)
m0_au = 1.007 * AMU_ME
mt_au = A * AMU_ME
ne = ZBAR * nuc.Nion(RS)
U = nuc.build_potential('gk', 1.0, ZBAR, ne, 0.025 / HARTREE_EV, od='/tmp/rng_Al', rs=RS)
E_lab = E * 1.007 * 1e6 / HARTREE_EV                                  # Hartree
Ec = (mt_au / (m0_au + mt_au)) * E_lab                               # CM energy, Hartree
Sfunc = nuc.make_S_interp(U, Ec.min(), Ec.max())
S_imp = np.array([Sfunc(ec) for ec in Ec])                          # Bohr^2 (int sin^2 p dp)
sigma_tr_cm2 = 4.0 * np.pi * S_imp * BOHR_CM**2                      # cm^2
sigtr_of_E = interpolate.interp1d(E, sigma_tr_cm2, bounds_error=False,
                                  fill_value=(sigma_tr_cm2[0], sigma_tr_cm2[-1]))

# detour: R_proj/R for a grid of incident energies
E0_grid = np.array([1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0])
detour = []
for E0 in E0_grid:
    Eint = np.linspace(1e-4, E0, 3000)
    Si = S_of_E(Eint); stri = sigtr_of_E(Eint)
    # inner(E) = int_E^E0 sigma_tr/S dE''   (cumulative from top)
    integrand = stri / Si
    inner = np.concatenate([[0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(Eint))])
    inner = inner[-1] - inner                                        # int from E to E0
    dR = (NA / A) * np.trapezoid(inner / Si, Eint)                  # g/cm2 deficit
    R0 = np.interp(E0, E, R)
    detour.append(max(0.0, 1.0 - dR / R0))
detour = np.array(detour)

# PSTAR detour
p = np.loadtxt('data/Al/pstar.txt', skiprows=8)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
axL.loglog(E, sigma_tr_cm2, 'b-', lw=2)
axL.set_xlabel('proton energy (MeV)'); axL.set_ylabel(r'nuclear $\sigma_{tr}=4\pi S$ (cm$^2$)')
axL.set_title('Nuclear transport cross section (same integral as nuclear stopping)')
axL.grid(alpha=.3, which='both')

axR.axvspan(8e-4, 0.1, color='0.9', zorder=0)
axR.semilogx(E0_grid, detour, 'bo-', lw=2, label='from nuclear model ($\\sigma_{tr}=4\\pi S$)')
axR.semilogx(p[:, 0], p[:, 6], 'k--', lw=1.5, label='PSTAR detour factor')
axR.set_xlabel('proton energy (MeV)'); axR.set_ylabel('projected / CSDA range')
axR.set_title('Pathlength-to-projected-range detour: model vs PSTAR')
axR.set_ylim(0, 1.05); axR.set_xlim(8e-4, 20)
axR.text(0.0011, 0.08, 'small-angle diffusion\napprox. breaks down\n(end-of-range single\nlarge-angle scattering)',
         fontsize=7.5, color='0.4', va='bottom')
axR.legend(fontsize=9, loc='lower right'); axR.grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig('detour_nuc.pdf'); fig.savefig('detour_nuc.png', dpi=130)
print('wrote detour_nuc.pdf / .png\n')
print(' E0(MeV)   detour(model)   detour(PSTAR)')
for E0, dt in zip(E0_grid, detour):
    ps = np.interp(E0, p[:, 0], p[:, 6])
    print('  %.0e     %.3f          %.3f' % (E0, dt, ps))
