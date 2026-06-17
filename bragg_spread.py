"""
bragg_spread.py -- how finite ion-source energy spread washes out the Bragg
peak and erases the differences between stopping-power models.

Thesis (closing discussion of the dedx paper): electronic stopping is a
non-resonant, Fokker-Planck-type sum of many small-angle collisions, so an
accurate average-atom model (eRPA) already gives the right energy loss and
range.  The fine differences between stopping models (eRPA vs RPA vs TD-DFT,
~10-20% near the Bragg peak) are smaller than the smearing of the deposition
profile produced by realistic source energy spread -- a few % for thermal
fusion products or linac beams, 10-30% for laser-wakefield/TNSA beams.

We build the depth-deposition (Bragg) curve from the total (electronic +
Gordon-Kim nuclear) eRPA stopping for 10 MeV protons in cold Al, convolve it
over a Gaussian energy distribution of several relative widths, and overlay a
representative +-15% model band (the documented eRPA-RPA / AA-TD-DFT spread).
We also compute the intrinsic Bohr range straggling, which sets the Bragg
width floor even for a perfectly monoenergetic beam.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, interpolate

NA, A, Z2, RHO = 6.022e23, 26.98, 13.0, 2.70           # Al
CONV = NA / A * 1e-21                                    # 1e-15 eVcm2/atom -> MeV/(g/cm2)
E0 = 10.0                                                # MeV proton

d = np.loadtxt('/tmp/rng_Al/dedx_nuc.dat', comments='#')
E = d[:, 0]
S = d[:, 3] * CONV                                       # total stopping, MeV/(g/cm2)
# range R(E) = int_0^E dE'/S  (g/cm2)
R = integrate.cumulative_trapezoid(1.0 / S, E, initial=0.0)
R_of_E = interpolate.interp1d(E, R, bounds_error=False, fill_value=(0.0, R[-1]))
E_of_R = interpolate.interp1d(R, E, bounds_error=False, fill_value=(0.0, E[-1]))
S_of_E = interpolate.interp1d(E, S, bounds_error=False, fill_value=0.0)


def bragg(E0v, x):
    """deposition S(E_res) [MeV/(g/cm2)] at depth x [g/cm2] for mono energy E0v."""
    Rres = R_of_E(E0v) - x
    return np.where(Rres > 0, S_of_E(E_of_R(np.clip(Rres, 0, None))), 0.0)


def ensemble(mu, sfrac, x, ng=241):
    """Bragg curve averaged over a Gaussian energy distribution, width sfrac*mu."""
    if sfrac <= 0:
        return bragg(mu, x)
    g = np.linspace(-4, 4, ng)
    w = np.exp(-0.5 * g * g)
    Es = mu * (1 + sfrac * g)
    out = np.zeros_like(x)
    for wi, Ei in zip(w, Es):
        out += wi * bragg(max(Ei, 1e-3), x)
    return out / w.sum()


# depth grid (convert to mm via Al density)
R0 = float(R_of_E(E0))
x = np.linspace(0, R0 * 1.12, 1000)                     # g/cm2
xmm = x / RHO * 10.0                                    # mm

# intrinsic Bohr range straggling: dOmega2/dxi = 0.1569 Z1^2 (Z2/A) [MeV^2/(g/cm2)]
dOm2 = 0.1569 * 1.0 * (Z2 / A)
EE = np.linspace(1e-3, E0, 4000)
SS = S_of_E(EE)
sigR = np.sqrt(np.trapezoid(dOm2 / SS**3, EE))          # g/cm2
sigR_mm = sigR / RHO * 10.0

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# left: deposition vs depth at several energy spreads
spreads = [0.0, 0.01, 0.05, 0.10, 0.20]
labels = ['monoenergetic', '1% (thermal / linac)', '5%', '10% (LWFA)', '20% (LWFA/TNSA)']
cmap = plt.cm.viridis(np.linspace(0.1, 0.85, len(spreads)))
for c, s, lab in zip(cmap, spreads, labels):
    axL.plot(xmm, ensemble(E0, s, x), color=c, lw=2, label=lab)
axL.axvspan((R0 - sigR) / RHO * 10, (R0 + sigR) / RHO * 10, color='0.85', zorder=0)
axL.set_xlabel('depth in Al (mm)')
axL.set_ylabel('deposition  dE/dx  (MeV/(g/cm$^2$))')
axL.set_title('10 MeV protons: source energy spread washes out the Bragg peak')
axL.legend(fontsize=8.5); axL.grid(alpha=.3)
axL.text(0.02, 0.02, 'grey band = intrinsic Bohr range straggling (%.0f µm)' % (sigR_mm * 1e3),
         transform=axL.transAxes, fontsize=8, color='0.4')

# right: model-difference washout. eRPA vs a +15% peak-region model.
Emask = E < 0.5                                          # Bragg-peak region
S2 = S.copy(); S2[Emask] *= 1.15
S2_of_E = interpolate.interp1d(E, S2, bounds_error=False, fill_value=0.0)
R2 = integrate.cumulative_trapezoid(1.0 / S2, E, initial=0.0)
R2_of_E = interpolate.interp1d(E, R2, bounds_error=False, fill_value=(0, R2[-1]))
E2_of_R = interpolate.interp1d(R2, E, bounds_error=False, fill_value=(0, E[-1]))


def bragg2(E0v, x):
    Rres = R2_of_E(E0v) - x
    return np.where(Rres > 0, S2_of_E(E2_of_R(np.clip(Rres, 0, None))), 0.0)


def ens2(mu, sfrac, x, ng=241):
    g = np.linspace(-4, 4, ng); w = np.exp(-0.5 * g * g); Es = mu * (1 + sfrac * g)
    out = np.zeros_like(x)
    for wi, Ei in zip(w, Es):
        out += wi * bragg2(max(Ei, 1e-3), x)
    return out / w.sum()


axR.plot(xmm, ensemble(E0, 0, x), 'b-', lw=2, label='eRPA, monoenergetic')
axR.plot(xmm, bragg2(E0, x), 'r--', lw=2, label='+15% peak model, monoenergetic')
axR.plot(xmm, ensemble(E0, 0.10, x), 'b-', lw=3.5, alpha=0.45, label='eRPA, 10% spread')
axR.plot(xmm, ens2(E0, 0.10, x), 'r--', lw=2, label='+15% peak model, 10% spread')
axR.set_xlabel('depth in Al (mm)')
axR.set_ylabel('deposition  dE/dx  (MeV/(g/cm$^2$))')
axR.set_title('A 15% stopping-model difference vanishes at 10% energy spread')
axR.legend(fontsize=8.5); axR.grid(alpha=.3)
axR.set_xlim(0.4, 0.72)

fig.tight_layout()
fig.savefig('bragg_spread.pdf'); fig.savefig('bragg_spread.png', dpi=130)
print('wrote bragg_spread.pdf / .png\n')

# quantify peak height and FWHM vs spread; model peak difference vs spread
def peak_stats(prof):
    pk = prof.max(); ipk = prof.argmax()
    half = pk / 2
    above = np.where(prof > half)[0]
    fwhm = (xmm[above[-1]] - xmm[above[0]]) if len(above) > 1 else 0
    return pk, xmm[ipk], fwhm

print(' spread   peak(MeV/g/cm2)  peak_depth(mm)  FWHM(mm)')
for s, lab in zip(spreads, labels):
    pk, pd, fw = peak_stats(ensemble(E0, s, x))
    print('  %-20s %7.1f       %.4f      %.4f' % (lab, pk, pd, fw))
print('\nintrinsic Bohr range straggling sigma_R = %.1f um (%.2f%% of range)' %
      (sigR_mm * 1e3, 100 * sigR / R0))
print('CSDA range R0 = %.4f g/cm2 = %.4f mm' % (R0, R0 / RHO * 10))
# model peak-deposition difference at 0% and 10% spread
for s in [0.0, 0.10]:
    p1 = ensemble(E0, s, x).max(); p2 = ens2(E0, s, x).max()
    print('  peak deposition: eRPA=%.1f  +15%%model=%.1f  -> differ %.1f%% at %.0f%% spread' %
          (p1, p2, 100 * (p2 / p1 - 1), s * 100))
