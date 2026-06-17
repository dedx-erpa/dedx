"""
divergence_spread.py -- how finite beam divergence (and the combination of
energy spread + divergence) smears the energy-deposition profile.

A ray launched at polar angle theta from the beam axis reaches axial depth
z = s*cos(theta) after pathlength s, and deposits dE/dz = S(E(s))/cos(theta).
Its Bragg peak therefore lands at a shallower axial depth R*cos(theta) and the
ensemble over a divergence distribution smears the peak -- the geometric analog
of energy spread.  We model the divergence as a 2-D Gaussian angular
distribution (Rayleigh polar angle, rms half-angle theta_rms).

Regimes:
  linac / RF accel.  : <~0.1 deg, sub-% energy spread       -> essentially sharp
  LWFA / fast-ignition proton beams : ~5-20 deg, 5-30% energy spread -> washed out
  fusion-born ions   : isotropic (4 pi), thermal energy spread, distributed
                       birth radii -> deposition is fully volumetric (text).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, interpolate

NA, A, RHO = 6.022e23, 26.98, 2.70
CONV = NA / A * 1e-21
E0 = 10.0

d = np.loadtxt('/tmp/rng_Al/dedx_nuc.dat', comments='#')
E = d[:, 0]; S = d[:, 3] * CONV
R = integrate.cumulative_trapezoid(1.0 / S, E, initial=0.0)
R_of_E = interpolate.interp1d(E, R, bounds_error=False, fill_value=(0.0, R[-1]))
E_of_R = interpolate.interp1d(R, E, bounds_error=False, fill_value=(0.0, E[-1]))
S_of_E = interpolate.interp1d(E, S, bounds_error=False, fill_value=0.0)
R0 = float(R_of_E(E0))


def dep_axial(E0v, z, cth):
    """deposition per axial depth at axial z for a ray of energy E0v, cos(theta)=cth."""
    s = z / cth
    Rres = R_of_E(E0v) - s
    return np.where(Rres > 0, S_of_E(E_of_R(np.clip(Rres, 0, None))) / cth, 0.0)


def smear(E0v, sE, thrms_deg, z, nth=120, nE=81):
    """deposition vs axial depth z, averaged over Gaussian energy spread sE and
    Rayleigh divergence with rms half-angle thrms_deg."""
    # divergence quadrature (Rayleigh: P(th) ~ th exp(-th^2/thrms^2))
    thr = np.radians(max(thrms_deg, 1e-4))
    th = np.linspace(0, 3.2 * thr, nth)
    wth = th * np.exp(-(th / thr)**2) if thrms_deg > 1e-3 else None
    cth = np.cos(th)
    # energy quadrature (Gaussian)
    if sE > 0:
        g = np.linspace(-4, 4, nE); wE = np.exp(-0.5 * g * g); Es = E0v * (1 + sE * g)
    else:
        g = [0.0]; wE = np.array([1.0]); Es = np.array([E0v])
    out = np.zeros_like(z); wsum = 0.0
    for wEi, Ei in zip(wE, Es):
        Ei = max(Ei, 1e-3)
        if thrms_deg <= 1e-3:
            out += wEi * dep_axial(Ei, z, 1.0); wsum += wEi
        else:
            for wt, ct in zip(wth, cth):
                out += wEi * wt * dep_axial(Ei, z, ct); wsum += wEi * wt
    return out / wsum


z = np.linspace(0, R0 * 1.12, 1200)
zmm = z / RHO * 10.0

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# (a) pure divergence, monoenergetic
divs = [0.0, 5.0, 10.0, 20.0]
labs = ['0° (linac / collimated)', '5° divergence', '10° (LWFA)', '20° (fast-ignition / TNSA)']
cmap = plt.cm.plasma(np.linspace(0.05, 0.8, len(divs)))
for c, dv, lab in zip(cmap, divs, labs):
    axL.plot(zmm, smear(E0, 0.0, dv, z), color=c, lw=2, label=lab)
axL.set_xlabel('axial depth in Al (mm)'); axL.set_ylabel('deposition  dE/dz  (MeV/(g/cm$^2$))')
axL.set_title('10 MeV protons: beam divergence smears the Bragg peak')
axL.legend(fontsize=8.5); axL.grid(alpha=.3)

# (b) realistic source combinations: energy spread + divergence
axR.plot(zmm, smear(E0, 0.0, 0.0, z), 'k-', lw=2, label='ideal: monoenergetic, 0°')
axR.plot(zmm, smear(E0, 0.005, 0.1, z), color='tab:green', lw=2, label='linac: 0.5% , 0.1°')
axR.plot(zmm, smear(E0, 0.10, 10.0, z), color='tab:red', lw=2.5, label='LWFA / fast-ignition: 10% , 10°')
axR.plot(zmm, smear(E0, 0.30, 20.0, z), color='tab:purple', lw=2, ls='--', label='broad TNSA: 30% , 20°')
axR.set_xlabel('axial depth in Al (mm)'); axR.set_ylabel('deposition  dE/dz  (MeV/(g/cm$^2$))')
axR.set_title('Energy spread + divergence compound the smearing')
axR.legend(fontsize=8.5); axR.grid(alpha=.3)

fig.tight_layout(); fig.savefig('divergence_spread.pdf'); fig.savefig('divergence_spread.png', dpi=130)
print('wrote divergence_spread.pdf / .png\n')


def stats(prof):
    pk = prof.max(); ip = prof.argmax()
    half = pk / 2; ab = np.where(prof > half)[0]
    return pk, zmm[ip], (zmm[ab[-1]] - zmm[ab[0]] if len(ab) > 1 else 0)


print(' divergence (mono)   peak   peak_depth(mm)  FWHM(mm)')
for dv, lab in zip(divs, labs):
    pk, pd, fw = stats(smear(E0, 0.0, dv, z))
    print('  %-22s %6.1f   %.4f       %.4f' % (lab, pk, pd, fw))
print('\n combined source       peak   peak_depth(mm)  FWHM(mm)')
for sE, dv, lab in [(0, 0, 'ideal'), (0.005, 0.1, 'linac'), (0.10, 10, 'LWFA'), (0.30, 20, 'TNSA')]:
    pk, pd, fw = stats(smear(E0, sE, dv, z))
    print('  %-20s %6.1f   %.4f       %.4f' % (lab, pk, pd, fw))
