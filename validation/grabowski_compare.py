# SPDX-License-Identifier: GPL-3.0-or-later
"""
Grabowski et al. (PRL 111, 215002, 2013) strong-coupling stopping comparison.

STATUS: framework + our-side (e-RPA) curves. HELD FOR FUTURE -- not part of the
PoP paper. Purpose: exercise the strong-collision (Zwicknagel) correction, which
is inactive in every weakly-coupled benchmark in the paper (Cayzac Gamma~0.01).

What this does:
  * runs the e-RPA pipeline (proton in H) at several electron-coupling values
    Gamma_e, with the strong-collision correction OFF (mloss=1, bare RPA) and
    ON (mloss=24, +LFC+strong-collision+Bloch);
  * plots dE/dx vs v_p/v_th,e for each Gamma, so the strong-collision shift is
    visible, with an empty overlay slot for Grabowski's digitized classical-MD
    points (data/refs/grabowski_md_*.csv, WebPlotDigitizer, as for Frenje).

IMPORTANT physical caveat (write this into any future figure caption):
  Grabowski's MD is a CLASSICAL plasma. For a real electron plasma, strong
  coupling (Gamma_e > 1) and non-degeneracy (Theta = T/E_F > 1) barely coexist --
  electrons degenerate before they strongly couple -- so the comparison is
  apples-to-apples only at Gamma_e < 1 (Theta > 1). At Gamma_e >~ 1 the e-RPA
  point is degenerate (Theta < 1) while Grabowski's is classical; it is then an
  indicative test of the correction's direction/magnitude, not a strict match.
  Theta is computed and printed/annotated for every point so this is explicit.
"""
import os, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY_ENV = sys.executable
HARTREE_EV = 27.211386
A0_CM = 0.529177e-8

# (density g/cc, Te eV) chosen to span Gamma_e ~ 0.3, 1, 3 for hydrogen (zbar~1)
CONDITIONS = [
    dict(d=1.0,  t=66.0, tag='g03'),
    dict(d=10.0, t=43.0, tag='g10'),
    dict(d=50.0, t=20.0, tag='g30'),
]

def nion_per_bohr3(d_gcc, A=1.008):
    return d_gcc / A * 6.02214076e23 * A0_CM**3

def gamma_theta(d, t):
    n = nion_per_bohr3(d)                 # zbar~1 for H plasma
    rs = (3.0/(4*np.pi*n))**(1/3)         # a_e in Bohr
    Te = t/HARTREE_EV
    Gamma = 1.0/(rs*Te)
    Ef = 0.5*(3*np.pi**2*n)**(2/3)        # Hartree
    Theta = Te/Ef
    return Gamma, Theta, rs

def run(cond, mloss):
    od = os.path.join(REPO, 'validation', 'gr_%s_m%d' % (cond['tag'], mloss))
    cmd = [PY_ENV, 'dedx.py', '--zp=1', '--zt=1',
           '--d=%g' % cond['d'], '--t=%g' % cond['t'],
           '--aa=2', '--mloss=%d' % mloss,
           '--emin=1e-3', '--emax=10', '--mep=60', '--od=' + od]
    subprocess.run(cmd, cwd=REPO, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return np.loadtxt(os.path.join(od, 'dedx.dat'), comments='#')

def main():
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
    for ax, cond in zip(axes, CONDITIONS):
        G, Th, rs = gamma_theta(cond['d'], cond['t'])
        m1 = run(cond, 1)     # bare RPA
        m24 = run(cond, 24)   # + LFC + strong-collision + Bloch
        vthe = np.sqrt(2*cond['t']/511000.0)          # /c
        vp = np.sqrt(2*m1[:, 0]/931.0)                # /c, E/A in MeV
        x = vp/vthe
        ax.plot(x, m1[:, 1],  '--', lw=1.6, label='e-RPA (bare RPA)')
        ax.plot(x, m24[:, 1], '-',  lw=1.8, label='e-RPA (+strong-collision)')
        # --- overlay slot for Grabowski classical MD (digitize into this CSV) ---
        md = os.path.join(REPO, 'data', 'refs', 'grabowski_md_%s.csv' % cond['tag'])
        if os.path.exists(md):
            g = np.loadtxt(md, delimiter=',')
            ax.plot(g[:, 0], g[:, 1], 'ko', ms=5, label='Grabowski MD (2013)')
        cls = 'classical' if Th > 1 else 'DEGENERATE'
        ax.set_title(r'$\Gamma_e\approx%.1f$, $\Theta\approx%.2f$ (%s)' % (G, Th, cls),
                     fontsize=10)
        ax.set_xlabel(r'$v_p/v_{th,e}$'); ax.set_xlim(0, x.max())
        ax.grid(alpha=0.3)
        if ax is axes[0]:
            ax.set_ylabel(r'$dE/dx$  [$10^{-15}$ eV cm$^2$/atom]')
        ax.legend(fontsize=8)
        print('%s: Gamma_e=%.2f  Theta=%.2f  rs=%.3f  %s'
              % (cond['tag'], G, Th, rs, cls))
    fig.suptitle('Strong-collision (Zwicknagel) correction vs coupling '
                 '-- e-RPA curves; drop in Grabowski MD to complete', fontsize=11)
    fig.tight_layout()
    out = os.path.join(REPO, 'validation', 'grabowski_compare.png')
    fig.savefig(out, dpi=140); print('wrote', out)

if __name__ == '__main__':
    main()
