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

# ---- model: proton in DT plasma, Te=2 keV ----
d = np.loadtxt('/tmp/frenje_dt2/dedx_nuc.dat', comments='#')
E, de, dn, dt = d[:, 0], d[:, 1], d[:, 2], d[:, 3]   # E/AMU, elec, nuc, total

# calibrate k so k*electronic best matches the BPS curve (log least squares)
bE = np.interp(bx, E, de)                    # model electronic at BPS x-values
k = np.exp(np.mean(np.log(by) - np.log(bE)))
print('areal-density calibration k = %.4g (set by model-electronic = BPS)' % k)

fig, ax = plt.subplots(figsize=(8.5, 6))
ax.errorbar(pts[:, 0], pts[:, 1], yerr=[pts[:, 1] - pts[:, 2], pts[:, 3] - pts[:, 1]],
            fmt='ks', ms=7, capsize=4, label='Frenje 2019 data', zorder=5)
ax.plot(bx, by, color='0.5', lw=1.5, label='BPS (Frenje, digitized)')
ax.plot(E, k * de, 'b--', lw=2, label='model: electronic only')
ax.plot(E, k * dt, 'r-', lw=2.2, label='model: electronic + nuclear/ion')
ax.set_xscale('log'); ax.set_xlim(0.1, 20); ax.set_ylim(0, 0.85)
ax.set_xlabel(r'$E_i/A_i$ (MeV)')
ax.set_ylabel(r'$-\Delta E_i/Z_i^2$ (MeV)')
ax.set_title('Frenje 2019: ion stopping in DT plasma (Te≈2 keV)')
ax.legend(fontsize=9); ax.grid(alpha=.3, which='both')
fig.tight_layout()
fig.savefig('exp_frenje.pdf'); fig.savefig('exp_frenje.png', dpi=130)
print('wrote exp_frenje.pdf / .png\n')

print(' E/A(MeV)  data   model_elec  model_tot   (k-scaled)')
for xc, yc, lo, hi in pts:
    me = k * np.interp(xc, E, de); mt = k * np.interp(xc, E, dt)
    print('  %6.3f  %.3f   %.3f      %.3f' % (xc, yc, me, mt))
