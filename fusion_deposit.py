"""
fusion_deposit.py -- draft "fusion isotropic source" panel.

In a fusion hot spot the alphas are not a collimated monoenergetic beam: they are
born isotropically (4 pi), uniformly through the burning volume, with a thermal
energy spread.  A 3-D Monte-Carlo transport of 3.5 MeV alphas (range R_alpha)
through a uniform DT sphere of radius rho R shows that:

  (left)  the deposited energy density is volumetric -- flat in the interior,
          rounded over one alpha range at the edge -- with NO Bragg peak; and
  (right) the self-heated (retained) energy fraction depends only on rho R / R_alpha,
          i.e. on the INTEGRATED range, which the eRPA model gives accurately.

So in the fusion context the differential Bragg shape -- and with it the fine
details of the stopping model near the peak -- is irrelevant; only the range matters.

Stopping: 3.5 MeV alpha in a DT plasma (~50 g/cc-equiv, 5 keV), eRPA total
(electronic + Gordon-Kim nuclear).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, interpolate

CONV = 6.022e23 / 2.5 * 1e-21          # DT mean A~2.5 -> MeV/(g/cm2) (illustrative)
d = np.loadtxt('/tmp/alpha_dt/dedx_nuc.dat', comments='#')
Eu, S = d[:, 0], d[:, 3] * CONV
E0 = 0.875                              # 3.5 MeV / 4 amu
R = integrate.cumulative_trapezoid(1.0 / S, Eu, initial=0.0)
R_of_E = interpolate.interp1d(Eu, R, fill_value=(0, R[-1]), bounds_error=False)
E_of_R = interpolate.interp1d(R, Eu, fill_value=(0, Eu[-1]), bounds_error=False)
Ra = float(R_of_E(E0))                  # alpha range (g/cm2)

rng = np.random.default_rng(0)
NS = 200                                # path steps
sgrid = np.linspace(0, Ra, NS)
# energy deposited in each path step (Bragg), per alpha, in units of E0
Ed = E_of_R(np.clip(Ra - sgrid, 0, None))            # residual energy along path
dep_step = -np.diff(Ed, append=0.0)                  # energy lost per step (MeV/amu)


def mc(a, N=120000):
    """3-D MC: alphas born uniform-in-volume, isotropic, range Ra, in a sphere
    of radius a.  Returns (radial deposition density, retained fraction)."""
    rb = a * rng.random(N) ** (1.0 / 3.0)                 # birth radius (uniform in vol)
    cb = 2 * rng.random(N) - 1                            # cos of birth polar (axis = ray)
    # birth.direction dot product: place ray along z, birth at (rho, theta) -> b.u = rb*cb
    bu = rb * cb
    rbin = np.linspace(0, a, 60)
    rcen = 0.5 * (rbin[1:] + rbin[:-1])
    vol = (4.0 / 3.0) * np.pi * (rbin[1:]**3 - rbin[:-1]**3)
    dens = np.zeros(len(rcen)); ret = 0.0
    for k in range(NS - 1):
        s = sgrid[k]
        pos = np.sqrt(rb**2 + 2 * bu * s + s * s)        # |b + s u|
        inside = pos <= a
        ret += dep_step[k] * inside.sum()
        idx = np.clip(np.searchsorted(rbin, pos[inside]) - 1, 0, len(rcen) - 1)
        np.add.at(dens, idx, dep_step[k])
    ret /= (E0 * N)
    return rcen, dens / vol / N, ret


fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

for fac, lab, c in [(0.5, r'$\rho R = R_\alpha$', 'tab:green'),
                    (1.5, r'$\rho R = 3\,R_\alpha$', 'tab:orange'),
                    (5.0, r'$\rho R = 10\,R_\alpha$', 'tab:red')]:
    rc, dens, _ = mc(fac * Ra)
    axL.plot(rc / Ra, dens / dens.max(), color=c, lw=2, label=lab)
axL.set_xlabel(r'radius / alpha range'); axL.set_ylabel('deposited energy density (norm.)')
axL.set_title('Hot-spot alpha deposition is volumetric — no Bragg peak')
axL.legend(fontsize=9); axL.grid(alpha=.3)

taus = np.linspace(0.1, 12, 40)
fr = np.array([mc(t * Ra, N=60000)[2] for t in taus])
axR.plot(taus, fr, 'b-', lw=2.5)
axR.axhline(0.5, color='0.6', ls=':'); axR.axvline(1, color='0.6', ls=':')
axR.set_xlabel(r'hot-spot size  $\rho R / R_\alpha$'); axR.set_ylabel('retained (self-heated) alpha fraction')
axR.set_title('Self-heating is set by the integrated range, not the Bragg shape')
axR.set_ylim(0, 1.02); axR.grid(alpha=.3)

fig.tight_layout(); fig.savefig('fusion_deposit.pdf'); fig.savefig('fusion_deposit.png', dpi=130)
print('wrote fusion_deposit.pdf / .png')
print('alpha range R_alpha = %.4f g/cm2 (3.5 MeV alpha in DT, ~50 g/cc-equiv, 5 keV)' % Ra)
print(' rhoR/Ra   retained fraction (3-D MC)')
for t in [0.5, 1, 2, 5, 10]:
    print('  %5.1f      %.3f' % (t, mc(t * Ra, N=80000)[2]))
