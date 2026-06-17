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

(The radial-output corruption at coarse energy grids -- a mre=mep/nre truncation
bug that let high-energy slices overwrite low-energy radial columns -- was fixed
in dedx.f; the decomposition is now clean and grid-independent for all slices.)

Run after:  python dedx.py --zt=13 --zp=2 --d=2.7 --t=1 --aa=2 --nuc=0 --mout=1 \
                           --od=/tmp/aldepth --emin=0.008 --emax=5 --mep=49
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

VP = lambda E_amu: 1.389e9 * np.sqrt(E_amu)          # cm/s from MeV/amu (alpha)
d, r = rd.rcdedx('/tmp/aldepth')
rad = d[0]
nE = r.shape[1]
Evals = np.array([r[0, i, 0] for i in range(nE)])
vp = VP(Evals)


def quantile_depth(i, q):
    cum = r[2, i, :] / r[2, i, -1]
    return np.interp(q, cum, rad)


fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# left: median stopping depth (with 25-75% band) vs velocity, shell bands
r50 = np.array([quantile_depth(i, 0.5) for i in range(nE)])
r25 = np.array([quantile_depth(i, 0.25) for i in range(nE)])
r75 = np.array([quantile_depth(i, 0.75) for i in range(nE)])
for y0, y1, lab in [(0.03, 0.15, 'K'), (0.3, 0.8, 'L'), (1.0, 2.99, 'valence/free')]:
    axL.axhspan(y0, y1, color='0.93', zorder=0)
    axL.text(1.05e8, np.sqrt(y0 * y1), lab, fontsize=8, color='0.4', va='center')
axL.fill_between(vp, r25, r75, color='tab:blue', alpha=0.2, label='25-75% range')
axL.plot(vp, r50, 'o-', color='tab:blue', lw=2, label='median depth')
axL.axvspan(3e8, 1.1e9, color='tab:orange', alpha=0.12)
axL.text(5.7e8, 0.06, 'Stanek\nvp range', fontsize=8, color='tab:orange', ha='center')
axL.set_xscale('log'); axL.set_xlim(1e8, 1.6e9); axL.set_ylim(0, 3.0)
axL.set_xlabel('alpha velocity vp (cm/s)')
axL.set_ylabel('stopping depth (atomic units)')
axL.set_title('Al: median stopping depth moves inward with velocity')
axL.legend(fontsize=8, loc='upper right'); axL.grid(alpha=.3, which='both')

# right: cumulative fraction vs radius for a few velocities
pick = [0, 2, 4, 6, 9]
cmap = plt.cm.plasma(np.linspace(0.05, 0.85, len(pick)))
for c, i in zip(cmap, pick):
    cum = r[2, i, :] / r[2, i, -1]
    axR.plot(rad, cum, color=c, lw=2, label='%.0f keV/u' % (Evals[i] * 1e3))
axR.set_xscale('log'); axR.set_xlim(0.05, 3.0); axR.set_ylim(0, 1.02)
axR.set_xlabel('radius (atomic units)')
axR.set_ylabel('cumulative fraction of $dE_e/dx$')
axR.set_title('Cumulative stopping vs depth')
axR.legend(fontsize=8, loc='upper left'); axR.grid(alpha=.3, which='both')

fig.suptitle('eRPA-AA stopping depth in Al (valence-dominated at low v; core re-engages when fast)')
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('stopping_depth.pdf'); fig.savefig('stopping_depth.png', dpi=130)
print('wrote stopping_depth.pdf / .png\n')
print(' vp(cm/s)   E(keV/u)   r(50%% of dEdx)  shell')
for i in range(nE):
    sh = 'L/inner' if r50[i] < 1.0 else 'valence/free'
    print('  %.2e   %6.1f      %.3f a.u.    %s' % (vp[i], Evals[i] * 1e3, r50[i], sh))
