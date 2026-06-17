"""
stopping_depth.py -- which electron shells supply the eRPA-AA electronic stopping,
as a function of projectile velocity?

The average-atom model integrates the stopping over the FULL bound+free radial
electron density, so dE/dx can be decomposed by radius:

    dE/dx = integral_0^rs  [ 4 pi r^2 rho(r) ] * L(r,v)  dr

This is something a TD-DFT-MD calculation cannot do directly: it propagates the
valence electrons explicitly and ADDS an estimated deep-core contribution.

PHYSICAL RESULT (Al, 2.7 g/cc, 1 eV):  the stopping is dominated by the loosely
bound valence/conduction electrons at the Wigner-Seitz edge at low velocity
(median depth r50 ~ 2.4 a.u.); as the projectile speeds up the median stopping
depth moves steadily INWARD -- the core re-engages -- reaching ~1.0 a.u. by
1 MeV/u.  In the Stanek alpha velocity range (Bragg peak and above, vp ~ 3-10e8
cm/s) the median depth is the valence/free + outer-L region, NOT the deep K-shell;
so the AA-vs-TD-DFT-MD gap near the peak is a valence/free-electron treatment
difference, not deep-core over-counting.

KNOWN DIAGNOSTIC CAVEAT (investigated June 2026): the radial-distribution output
(dedx --mout=1, arrays rdedx/cdedx) carries a spurious inner contribution at the
SINGLE LOWEST energy grid point (ie=1, e0=emin).  It is an edge artifact of the
mout diagnostic only:
  * it tracks emin exactly (same physical velocity is clean as the 2nd grid point
    but artifacted as the 1st), so it is not a physical velocity threshold;
  * it is NOT the binding correction (icb=0 leaves it unchanged) and NOT the
    stopping-table resolution (the table grid is fixed, 100..1e-4 MeV/u);
  * the INTEGRATED stopping in dedx.dat is unaffected -- matched-velocity totals
    are grid-independent -- so all stopping curves / experimental comparisons stand.
We therefore set emin well below the analysis range and DROP the first slice.

Run after:  python dedx.py --zt=13 --zp=2 --d=2.7 --t=1 --aa=2 --nuc=0 --mout=1 \
                           --od=/tmp/aldepth --emin=0.004 --emax=5 --mep=145
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

VP = lambda E_amu: 1.389e9 * np.sqrt(E_amu)          # cm/s from MeV/amu (alpha)
d, r = rd.rcdedx('/tmp/aldepth')
rad = d[0]
rho4pir2 = d[1]
M = np.concatenate([[0], np.cumsum(0.5 * (rho4pir2[1:] + rho4pir2[:-1]) * np.diff(rad))])
Z = M[-1]
sel = list(range(1, r.shape[1]))                     # DROP slice 0 (ie=1 edge artifact)
Evals = np.array([r[0, i, 0] for i in sel])
vp = VP(Evals)


def quantile_depth(i, q):
    cum = r[2, i, :] / r[2, i, -1]
    return np.interp(q, cum, rad)


fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# left: median stopping depth (with 25-75% band) vs velocity, shell bands
r50 = np.array([quantile_depth(i, 0.5) for i in sel])
r25 = np.array([quantile_depth(i, 0.25) for i in sel])
r75 = np.array([quantile_depth(i, 0.75) for i in sel])
for y0, y1, lab in [(0.03, 0.15, 'K'), (0.3, 0.8, 'L'), (1.0, 2.99, 'valence/free')]:
    axL.axhspan(y0, y1, color='0.93', zorder=0)
    axL.text(1.1e8, np.sqrt(y0 * y1), lab, fontsize=8, color='0.4', va='center')
axL.fill_between(vp, r25, r75, color='tab:blue', alpha=0.2, label='25-75% range')
axL.plot(vp, r50, 'o-', color='tab:blue', lw=2, label='median depth')
axL.axvspan(3e8, 1.1e9, color='tab:orange', alpha=0.12)
axL.text(5.7e8, 0.06, 'Stanek\nvp range', fontsize=8, color='tab:orange', ha='center')
axL.set_xscale('log'); axL.set_xlim(1e8, 2.5e9); axL.set_ylim(0, 3.0)
axL.set_xlabel('alpha velocity vp (cm/s)')
axL.set_ylabel('stopping depth (atomic units)')
axL.set_title('Al: median stopping depth moves inward with velocity')
axL.legend(fontsize=8, loc='upper right'); axL.grid(alpha=.3, which='both')

# right: cumulative fraction vs radius for a few velocities
pick = [sel[0], sel[2], sel[4], sel[6], sel[8]]
cmap = plt.cm.plasma(np.linspace(0.05, 0.85, len(pick)))
for c, i in zip(cmap, pick):
    cum = r[2, i, :] / r[2, i, -1]
    axR.plot(rad, cum, color=c, lw=2, label='%.0f keV/u' % (r[0, i, 0] * 1e3))
axR.set_xscale('log'); axR.set_xlim(0.05, 3.0); axR.set_ylim(0, 1.02)
axR.set_xlabel('radius (atomic units)')
axR.set_ylabel('cumulative fraction of $dE_e/dx$')
axR.set_title('Cumulative stopping vs depth')
axR.legend(fontsize=8, loc='upper left'); axR.grid(alpha=.3, which='both')

fig.suptitle('eRPA-AA stopping depth in Al (valence-dominated at low v; core re-engages when fast)')
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('stopping_depth.pdf'); fig.savefig('stopping_depth.png', dpi=130)
print('wrote stopping_depth.pdf / .png  (first slice ie=1 dropped)\n')
print(' vp(cm/s)   E(keV/u)   r(50%% of dEdx)  shell')
for i in sel:
    cum = r[2, i, :] / r[2, i, -1]
    r50i = np.interp(0.5, cum, rad)
    sh = 'L/inner' if r50i < 1.0 else 'valence/free'
    print('  %.2e   %6.1f      %.3f a.u.    %s' % (VP(r[0, i, 0]), r[0, i, 0] * 1e3, r50i, sh))
