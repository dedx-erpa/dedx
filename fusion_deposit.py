"""
fusion_deposit.py -- "fusion isotropic source" panel.

In a fusion hot spot the alphas are not a collimated monoenergetic beam: they are
born isotropically (4 pi), uniformly through the burning volume, with a thermal
energy spread.  A 3-D Monte-Carlo transport of 3.5 MeV alphas (CSDA range
R_alpha) through a uniform DT sphere of areal radius rho R, using the eRPA total
(electronic + Gordon-Kim nuclear) Bragg deposition profile, shows that:

  (left)  the deposited energy density is volumetric -- flat in the interior,
          rounded over one alpha range at the edge -- with NO Bragg peak; and
  (right) the self-heated (retained) energy fraction depends only on rho R / R_alpha,
          i.e. on the INTEGRATED range, which the eRPA model gives accurately.

So in the fusion context the differential Bragg shape -- and with it the fine
details of the stopping model near the peak -- is irrelevant; only the range matters.
The eRPA self-heating curve reaches ~0.68 at rho R = R_alpha, consistent with the
classic ICF result (Krokhin & Rozanov 1973; Atzeni & Meyer-ter-Vehn 2004).  A
constant-stopping (uniform-deposition) reference agrees to within a few percent,
confirming that the self-heating fraction is set by the integrated range and is
insensitive to the detailed shape of the deposition profile.

Stopping: 3.5 MeV alpha in DT, rho = 100 g/cm^3, T = 10 keV (R_alpha = 0.087 g/cm^2),
from the eRPA total model (electronic + Gordon-Kim nuclear).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, interpolate

# eRPA alpha-in-DT stopping (H at the DT electron density, T=10 keV); areal range
# in true DT mass units via mass/electron ratio 2.5/1.008.
CONV = 6.022e23 / 1.008 * 1e-21
d = np.loadtxt('/tmp/alpha_dt10/dedx_nuc.dat', comments='#')
Eu, S = d[:, 0], d[:, 3] * CONV
E0 = 0.875                              # 3.5 MeV / 4 amu
R = integrate.cumulative_trapezoid(1.0 / S, Eu, initial=0.0)
R_of_E = interpolate.interp1d(Eu, R, fill_value=(0, R[-1]), bounds_error=False)
E_of_R = interpolate.interp1d(R, Eu, fill_value=(0, Eu[-1]), bounds_error=False)
RaDT = float(R_of_E(E0)) * 2.5 / 1.008    # DT areal range (g/cm2), for the label

rng = np.random.default_rng(0)
NS = 240
sgrid = np.linspace(0, float(R_of_E(E0)), NS)             # work in (H-proxy) range units
Ra = sgrid[-1]
Ed = E_of_R(np.clip(Ra - sgrid, 0, None))                # residual energy along path
dep_bragg = -np.diff(Ed, append=0.0)                     # eRPA Bragg deposition per step
dep_unif = np.full(NS, E0 / NS)                          # constant-stopping reference


def mc(a, dep, N=200000, prof=False):
    """3-D MC: alphas born uniform-in-volume, isotropic, range Ra, sphere radius a.
    Returns retained fraction (and, if prof, the radial deposited-energy density)."""
    rb = a * rng.random(N) ** (1.0 / 3.0)
    cb = 2 * rng.random(N) - 1
    bu = rb * cb
    ret = 0.0
    if prof:
        rbin = np.linspace(0.02 * a, a, 50)
        rcen = 0.5 * (rbin[1:] + rbin[:-1])
        vol = (4.0 / 3.0) * np.pi * (rbin[1:]**3 - rbin[:-1]**3)
        dens = np.zeros(len(rcen))
    for k in range(NS - 1):
        pos = np.sqrt(rb**2 + 2 * bu * sgrid[k] + sgrid[k]**2)
        inside = pos <= a
        ret += dep[k] * inside.sum()
        if prof:
            idx = np.clip(np.searchsorted(rbin, pos[inside]) - 1, 0, len(rcen) - 1)
            np.add.at(dens, idx, dep[k])
    ret /= (E0 * N)
    return (ret, rcen, dens / vol / N) if prof else ret


fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

for fac, lab, c in [(0.5, r'$\rho R = R_\alpha$', 'tab:green'),
                    (1.5, r'$\rho R = 3\,R_\alpha$', 'tab:orange'),
                    (5.0, r'$\rho R = 10\,R_\alpha$', 'tab:red')]:
    _, rc, dens = mc(fac * Ra, dep_bragg, N=400000, prof=True)
    axL.plot(rc / Ra, dens / dens.max(), color=c, lw=2, label=lab)
axL.set_xlabel(r'radius / alpha range'); axL.set_ylabel('deposited energy density (norm.)')
axL.set_title('Hot-spot alpha deposition is volumetric — no Bragg peak')
axL.legend(fontsize=9); axL.grid(alpha=.3); axL.set_ylim(0, 1.05)

taus = np.linspace(0.1, 12, 44)
fb = np.array([mc(t * Ra, dep_bragg, N=120000) for t in taus])
fu = np.array([mc(t * Ra, dep_unif, N=120000) for t in taus])
axR.plot(taus, fb, 'b-', lw=2.5, label='eRPA Bragg deposition (this work)')
axR.plot(taus, fu, color='0.5', lw=1.8, ls='--', label='constant-stopping reference')
axR.axhline(0.5, color='0.7', ls=':'); axR.axvline(1, color='0.7', ls=':')
axR.set_xlabel(r'hot-spot size  $\rho R / R_\alpha$'); axR.set_ylabel('retained (self-heated) alpha fraction')
axR.set_title('Self-heating is set by the integrated range, not the Bragg shape')
axR.set_ylim(0, 1.02); axR.legend(fontsize=9, loc='lower right'); axR.grid(alpha=.3)

fig.tight_layout(); fig.savefig('fusion_deposit.pdf'); fig.savefig('fusion_deposit.png', dpi=130)
print('wrote fusion_deposit.pdf / .png')
print('DT alpha range R_alpha = %.4f g/cm2 (3.5 MeV alpha, DT 100 g/cc, 10 keV)' % RaDT)
print(' rhoR/Ra   retained (eRPA Bragg)   retained (const-stopping)')
for t in [0.5, 1, 2, 5, 10]:
    print('  %5.1f       %.3f                  %.3f' %
          (t, mc(t * Ra, dep_bragg, N=160000), mc(t * Ra, dep_unif, N=160000)))
