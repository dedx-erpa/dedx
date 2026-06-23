# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
radial_contrib.py -- radial origin of the eRPA-LDA proton stopping in cold Al.

For a proton at three energies, decompose the LDA stopping integral
    dE/dx = integral_0^rs  [4 pi r^2 rho(r)] * L(rho(r), v)  dr
by radius and show, in one figure: the radial electron density 4 pi r^2 rho(r),
the loss function L(r), the integrand 4 pi r^2 rho L, and the cumulative
fraction int_0^r 4 pi r^2 rho L dr.  A red dot marks the median (50%) stopping
depth.  This replaces the earlier single-energy figures and matches the APS-2024
talk version.

Run after:
  python dedx.py --zt=13 --zp=1 --d=2.7 --t=0.025 --aa=2 --nuc=0 --mout=1 \
                 --emin=0.105 --emax=113 --mep=49 --od=/tmp/alrad2
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

d, r = rd.rcdedx('/tmp/alrad2')
rad = np.asarray(d[0])            # radius (atomic units)
dens = np.asarray(d[1])          # 4 pi r^2 rho(r)  (total, bound+free)
slices = [0, 4, 8]               # 0.105, 1.076, 11.03 MeV
norm = lambda y: y / np.max(y[np.isfinite(y)])

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
for ax, i in zip(axes, slices):
    E = r[0, i, 0] * 1e3          # keV
    integ = np.asarray(r[1, i, :])
    cum = np.asarray(r[2, i, :])  # already normalized 0..1
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.where(dens > 1e-12 * dens.max(), integ / dens, 0.0)
    ax.plot(rad, norm(dens),  color='tab:blue',   lw=1.8, label=r'$4\pi r^2\rho(r)$')
    ax.plot(rad, norm(L),     color='tab:orange', lw=1.8, label=r'$L(r)$')
    ax.plot(rad, norm(integ), color='tab:green',  lw=1.8, label=r'$4\pi r^2\rho(r)L(r)$')
    ax.plot(rad, cum,         color='tab:red',    lw=1.8, label=r'$\int_0^r 4\pi r^2\rho L\,dr$')
    r50 = np.interp(0.5, cum, rad)
    ax.plot([r50], [0.5], 'o', color='tab:red', ms=9, zorder=5)
    ax.axvline(r50, color='tab:red', ls='--', lw=0.8, alpha=0.6)
    ax.set_xlim(0, 3.0); ax.set_ylim(0, 1.05)
    ax.set_xlabel('radius (atomic unit)')
    ax.set_title('E = %.2f keV' % E)
    ax.grid(alpha=0.3)
axes[0].set_ylabel('arb. unit')
axes[0].legend(fontsize=8, loc='upper center')
fig.suptitle('Radial origin of proton stopping in cold Al (red dot = median stopping depth)')
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('radial_contrib.png', dpi=140); fig.savefig('radial_contrib.pdf')
print('wrote radial_contrib.png/.pdf')
for i in slices:
    print('  E=%.3f MeV  median depth r50=%.3f a.u.'
          % (r[0, i, 0], np.interp(0.5, r[2, i, :], rad)))
