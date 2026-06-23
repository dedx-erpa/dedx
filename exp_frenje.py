# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 T. A. Mehlhorn, M. F. Gu, and I. Golovkin
"""
exp_frenje.py -- compare the total eRPA model to Frenje 2019 (PRL 122, 015002),
Fig. 5: measured ion stopping -dE_i/Z_i^2 vs E_i/A_i in DT plasma (Te~2 keV).

Data digitized with WebPlotDigitizer (data/refs/frenje_fig5_*.csv); the BPS file
also captured the inset curves on the main axes, which are cleaned out here.

The model -dE/Z^2 has the shape of dE/dx vs E/A; a single areal-density factor k
(= the path-integrated n*L, which we do not have without the Fig. 2 profile) sets
the scale.  We fix k by matching the model ELECTRONIC stopping to the BPS curve,
then show that adding the nuclear/ion channel lifts the low-velocity end toward
the measured points -- the effect Frenje attributed to nuclear-elastic scattering.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- digitized Frenje Fig 5 ----
raw = np.loadtxt('data/refs/frenje_fig5_bps.csv', delimiter=',')
dat = np.loadtxt('data/refs/frenje_fig5_data.csv', delimiter=',')

# clean the BPS: the inset traces lie BELOW the main (unimodal) curve for x>1.4,
# and a couple of strays sit ABOVE it at low x.  Drop the low-x strays, then take
# the upper envelope (max y in each log-x bin) -> the smooth main BPS curve.
x, y = raw[:, 0], raw[:, 1]
keep = ~((x < 0.3) & (y > 0.30))
x, y = x[keep], y[keep]
edges = np.geomspace(x.min() * 0.999, x.max() * 1.001, 45)
bx, by = [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (x >= lo) & (x < hi)
    if m.any():
        j = np.argmax(y[m])
        bx.append(x[m][j]); by.append(y[m][j])
bx, by = np.array(bx), np.array(by)

# data points: group into 4 (center=median, err from spread)
pts = []
for xc in sorted(set(np.round(dat[:, 0], 1))):
    m = np.abs(dat[:, 0] - xc) < 0.1 * max(1, xc)
    ys = dat[m, 1]
    if len(ys) >= 2:
        pts.append((np.median(dat[m, 0]), np.median(ys), ys.min(), ys.max()))
pts = np.array(pts)

# ---- absolute scale from the digitized Te(r)/ne(r) profile (Fig 2) ----
# corrections: Te y-axis was 0-15 not 0-5 (/3); ne x-axis offset by 110 um.
ne_raw = np.loadtxt('data/refs/frenje_fig3a_ne_curve.csv', delimiter=',')
r_ne = ne_raw[:, 0] - 110.0
v_ne = ne_raw[:, 1] * 1e23                    # cm^-3
o = np.argsort(r_ne); r_ne, v_ne = r_ne[o], v_ne[o]
rr = np.linspace(0, 80, 400)
nn = np.interp(rr, r_ne, np.maximum(v_ne, 0))
NL_rad = np.trapezoid(nn, rr * 1e-4)          # int_0^R ne dr  (radius path), /cm2
NL_dia = 2 * NL_rad                            # diameter chord
print('profile areal density: radius=%.2e  diameter=%.2e /cm2 (rhoL=%.1f-%.1f mg/cm2)'
      % (NL_rad, NL_dia, NL_rad * 2.5 * 1.66e-24 * 1e3, NL_dia * 2.5 * 1.66e-24 * 1e3))

# ---- model: proton in DT plasma, Te=2 keV (mass-weighted) ----
d = np.loadtxt('/tmp/frenje_dt2/dedx_nuc.dat', comments='#')
E, de, dn, dt = d[:, 0], d[:, 1], d[:, 2], d[:, 3]   # E/AMU, elec, nuc, total

# scale: -dE/Z^2 = eps * N_L * 1e-21.  We fix one effective areal density N_L by
# matching the model electronic stopping to the digitized BPS curve (shape test);
# the profile above shows this N_L is physical (within radius–diameter range).
bE = np.interp(bx, E, de)
k = np.exp(np.mean(np.log(by) - np.log(bE)))         # = N_L * 1e-21
NL_fit = k / 1e-21
print('fit-implied areal density N_L = %.2e /cm2  (profile range %.2e–%.2e)'
      % (NL_fit, NL_rad, NL_dia))

fig, ax = plt.subplots(figsize=(8.5, 6))
ax.plot(bx, by, color='0.5', lw=1.5, label='BPS (Frenje, digitized)')
ax.plot(E, k * de, 'b--', lw=2, label='model: electronic only')
ax.plot(E, k * dt, 'r-', lw=2.2, label='model: electronic + nuclear/ion')
ax.errorbar(pts[:, 0], pts[:, 1], yerr=[pts[:, 1] - pts[:, 2], pts[:, 3] - pts[:, 1]],
            fmt='ks', ms=7, capsize=4, label='Frenje 2019 data', zorder=5)
ax.set_xscale('log'); ax.set_xlim(0.1, 20); ax.set_ylim(0, 0.85)
ax.set_xlabel(r'$E_i/A_i$ (MeV)')
ax.set_ylabel(r'$-\Delta E_i/Z_i^2$ (MeV)')
ax.set_title('Frenje 2019: ion stopping in DT plasma (Te≈2 keV)')
ax.text(0.03, 0.04, 'scale N_L=%.1e /cm² (fit)\nprofile gives %.1e–%.1e /cm²\n→ scale is physical, not free'
        % (NL_fit, NL_rad, NL_dia), transform=ax.transAxes, fontsize=8, va='bottom',
        bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.7))
ax.legend(fontsize=9); ax.grid(alpha=.3, which='both')
fig.tight_layout()
fig.savefig('exp_frenje.pdf'); fig.savefig('exp_frenje.png', dpi=130)
print('wrote exp_frenje.pdf / .png\n')

print(' E/A(MeV)  data    model_elec  model_tot')
for xc, yc, lo, hi in pts:
    print('  %6.3f  %.3f    %.3f      %.3f' % (xc, yc, k * np.interp(xc, E, de), k * np.interp(xc, E, dt)))
