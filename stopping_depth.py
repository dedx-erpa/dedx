"""
stopping_depth.py -- where (which electron shells) do the eRPA-AA electronic
stopping contributions come from, as a function of projectile velocity?

The average-atom model integrates the stopping over the FULL bound+free radial
electron density, so dE/dx can be decomposed by radius:

    dE/dx = integral_0^rs  [ 4 pi r^2 rho(r) ] * L(r,v)  dr

where the bracket is the electron count per shell (d[1], integrates to Z) and
L is the local stopping number.  This is something a TD-DFT-MD calculation cannot
do directly: it propagates the valence electrons explicitly and ADDS an estimated
deep-core contribution.  Here we show the shell-resolved origin of the AA stopping.

Finding (Al, 2.7 g/cc, 1 eV): the dominant shell SHIFTS with velocity --
  - slow alpha   (<~0.03 MeV/u): bound L-shell (r ~ 0.5 a.u.)
  - near Bragg peak (~0.05-0.2): free/conduction electrons at the WS edge (r ~ rs)
  - fast alpha   (>~0.4):        core re-engages (r ~ 0.6 a.u.)
so the contribution is genuinely shell-resolved and non-monotonic in velocity.

Run after:  python dedx.py --zt=13 --zp=2 --d=2.7 --t=1 --aa=2 --nuc=0 --mout=1 \
                           --od=/tmp/aldepth --emin=0.01 --emax=5 --mep=49
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rd

VP = lambda E_amu: 1.389e9 * np.sqrt(E_amu)          # cm/s from MeV/amu (alpha)
d, r = rd.rcdedx('/tmp/aldepth')                     # d=[radius,4pi r^2 rho,...]; r=[E,diff,cumul]
rad = d[0]
rho4pir2 = d[1]                                       # electrons per unit radius (integrates to Z)
M = np.concatenate([[0], np.cumsum(0.5 * (rho4pir2[1:] + rho4pir2[:-1]) * np.diff(rad))])
Z = M[-1]
nE = r.shape[1]
Evals = np.array([r[0, i, 0] for i in range(nE)])     # MeV/amu
vp = VP(Evals)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# Al shell radii (approx, a.u.) for orientation
for x0, x1, lab, col in [(0.03, 0.15, 'K', '0.85'),
                         (0.3, 0.8, 'L', '0.92'),
                         (1.0, 2.99, 'valence/free', '0.97')]:
    axL.axvspan(x0, x1, color=col, zorder=0)
    axL.text(np.sqrt(x0 * x1), 0.93, lab, ha='center', fontsize=8, color='0.4')

# normalized electron density (per unit radius)
axL.plot(rad, rho4pir2 / rho4pir2.max(), 'k--', lw=1.2, label=r'$4\pi r^2\rho(r)$ (density)')

# differential stopping contribution per unit radius, several velocities
sel = [0, 2, 4, 6, 9]
cmap = plt.cm.plasma(np.linspace(0.05, 0.85, len(sel)))
for c, i in zip(cmap, sel):
    integ = r[1, i, :] / np.maximum(rad, 1e-12)       # 4 pi r^2 rho L
    axL.plot(rad, integ / integ.max(), color=c, lw=2,
             label='%.0f keV/u (vp=%.1e)' % (Evals[i] * 1e3, vp[i]))
axL.set_xscale('log'); axL.set_xlim(0.02, 3.0); axL.set_ylim(0, 1.02)
axL.set_xlabel('radius (atomic units)')
axL.set_ylabel('normalized contribution / density')
axL.set_title('Al: shell-resolved origin of $dE_e/dx$')
axL.legend(fontsize=7.5, loc='upper right'); axL.grid(alpha=.3, which='both')

# right: cumulative fraction vs radius (monotonic, honest) for the same velocities
for c, i in zip(cmap, sel):
    cum = r[2, i, :] / r[2, i, -1]
    axR.plot(rad, cum, color=c, lw=2, label='%.0f keV/u' % (Evals[i] * 1e3))
axR.set_xscale('log'); axR.set_xlim(0.02, 3.0); axR.set_ylim(0, 1.02)
axR.set_xlabel('radius (atomic units)')
axR.set_ylabel('cumulative fraction of $dE_e/dx$')
axR.set_title('Cumulative stopping vs depth')
axR.legend(fontsize=8, loc='lower right'); axR.grid(alpha=.3, which='both')

fig.suptitle('eRPA-AA stopping depth in Al: the contributing shell shifts with velocity')
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('stopping_depth.pdf'); fig.savefig('stopping_depth.png', dpi=130)
print('wrote stopping_depth.pdf / .png\n')

# half-contribution radius (median depth) per velocity
print(' vp(cm/s)   E(keV/u)   r(50%% of dEdx)  shell')
for i in range(nE):
    cum = r[2, i, :] / r[2, i, -1]
    r50 = np.interp(0.5, cum, rad)
    sh = 'L-bound' if r50 < 0.9 else ('valence/free' if r50 > 1.0 else 'mixed')
    print('  %.2e   %6.1f      %.3f a.u.    %s' % (vp[i], Evals[i] * 1e3, r50, sh))
