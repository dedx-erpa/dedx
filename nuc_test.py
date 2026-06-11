"""
nuc_test.py -- validation / example driver for the nuclear stopping model
(nuclear.py), reproducing the key results of Faussurier, Blancard & Gauthier,
Phys. Plasmas 20, 012705 (2013).

Run (in the `fac` conda env):
    python nuc_test.py fig1      # cold Al: 3 potentials vs NIST/ZBL nuclear (Fig 1)
    python nuc_test.py sweep     # Te = 1/10/100/1000 eV electronic+nuclear (Figs 2/3)
    python nuc_test.py all

fig1 reuses the existing data/Al average-atom density (no AA run needed).
sweep runs dedx.py at several temperatures (each runs the average-atom model).
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import nuclear as N


# ---------------------------------------------------------------------------
# ZBL universal nuclear stopping -- the basis of the NIST/ICRU-49 nuclear
# stopping column; returned directly in 10^-15 eV cm^2/atom.
# ---------------------------------------------------------------------------
def zbl_nuclear(E_MeV, Z1, Z2, M1, M2):
    Ek = np.asarray(E_MeV, dtype=float) * 1e3                 # keV
    a = Z1**0.23 + Z2**0.23
    eps = 32.53 * M2 * Ek / (Z1 * Z2 * (M1 + M2) * a)
    sn = np.where(eps <= 30.0,
                  np.log1p(1.1383 * eps) /
                  (2.0 * (eps + 0.01321 * eps**0.21226 + 0.19593 * np.sqrt(eps))),
                  np.log(eps) / (2.0 * eps))
    return 8.462 * Z1 * Z2 * M1 * sn / ((M1 + M2) * a)        # 1e-15 eV cm^2/atom


# ---------------------------------------------------------------------------
# Fig 1: cold solid-density Al, proton -- three potentials vs NIST/ZBL nuclear
# ---------------------------------------------------------------------------
def fig1(od='data/Al', zp=1, zt=13, m0=1.008, mt=26.98):
    import rd
    h = rd.rdedx(od, header='')
    rs, te, zbar = float(h['rs']), float(h['Te']), float(h['zbar'])
    E = np.geomspace(3e-4, 1.0, 60)                          # MeV/amu

    pots = {'ionsphere': 'ion sphere', 'yukawa': 'Yukawa', 'gk': 'Gordon-Kim'}
    eps = {p: np.atleast_1d(N.nuclear_stopping(E, zp, zt, m0, mt, rs, te, 0.0,
                                               zbar, potential=p, od=od))
           for p in pots}
    ref = zbl_nuclear(E, zp, zt, m0, mt)

    # save a reference table
    os.makedirs('data/refs', exist_ok=True)
    np.savetxt('data/refs/Al_nuclear_zbl.dat',
               np.column_stack([E, ref]),
               header='proton-in-Al nuclear stopping (ZBL/NIST)\n'
                      'E/amu[MeV]   Sn[1e-15 eV cm^2/atom]')

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(E, ref, 'k-', lw=2, label='NIST / ZBL')
    sty = {'ionsphere': 'b--', 'yukawa': 'g-.', 'gk': 'r:'}
    for p in pots:
        ax.loglog(E, eps[p], sty[p], lw=1.8, label=pots[p])
    ax.set_xlabel('proton energy (MeV)')
    ax.set_ylabel(r'nuclear stopping ($10^{-15}$ eV cm$^2$/atom)')
    ax.set_title('Proton in cold solid-density Al (cf. Faussurier Fig. 1)')
    ax.legend(); ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout(); fig.savefig('nuc_fig1_coldAl.pdf')
    print('wrote nuc_fig1_coldAl.pdf')

    # quantitative summary
    print('\n E(MeV)     ZBL      GK     GK/ZBL   ionsph   yukawa')
    for i in (0, 15, 30, 45, 59):
        print('%.2e  %.4f  %.4f   %.2f   %.4f  %.4f' %
              (E[i], ref[i], eps['gk'][i], eps['gk'][i] / ref[i],
               eps['ionsphere'][i], eps['yukawa'][i]))


# ---------------------------------------------------------------------------
# Figs 2/3: Te = 1/10/100/1000 eV, electronic + nuclear, and the range trend
# ---------------------------------------------------------------------------
def sweep(zt=13, d=2.7, tes=(1.0, 10.0, 100.0, 1000.0)):
    import rd
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ranges_e, ranges_t = [], []
    for te in tes:
        od = '/tmp/nucsweep_Al_%g' % te
        if not os.path.exists('%s/dedx_nuc.dat' % od):
            cmd = ('python dedx.py --zt=%d --d=%g --t=%g --aa=2 --nuc=1 '
                   '--npot=gk --od=%s --emin=1e-4 --emax=10 --mep=80' %
                   (zt, d, te, od))
            print('running:', cmd)
            os.system(cmd)
        dn = np.loadtxt('%s/dedx_nuc.dat' % od, comments='#')
        E, dedx_e, dedx_n, dtot, rng = dn[:, 0], dn[:, 1], dn[:, 2], dn[:, 3], dn[:, 4]
        ax1.loglog(E, dedx_e, '-', label='e, %g eV' % te)
        ax1.loglog(E, dedx_n, '--', label='n, %g eV' % te)
        # range at a few energies: with vs without nuclear
        re = rd.int_range(np.array((E, dedx_e)), m=26.98)
        ranges_e.append((te, np.interp(np.log10(1e-3), np.log10(E), re)))
        ranges_t.append((te, np.interp(np.log10(1e-3), np.log10(E), rng)))
    ax1.set_xlabel('E (MeV)'); ax1.set_ylabel(r'stopping ($10^{-15}$ eV cm$^2$/atom)')
    ax1.set_title('Proton in Al: electronic (solid) + nuclear/GK (dashed)')
    ax1.legend(fontsize=7); ax1.grid(True, which='both', alpha=0.3)

    re = np.array([r[1] for r in ranges_e]); rt = np.array([r[1] for r in ranges_t])
    ax2.semilogx(tes, re, 'ko-', label='electronic only')
    ax2.semilogx(tes, rt, 'rs--', label='with nuclear (GK)')
    ax2.set_xlabel('Te (eV)'); ax2.set_ylabel('proton range at 1 keV (mg/cm$^2$)')
    ax2.set_title('Range reduction from nuclear stopping (cf. Fig. 3)')
    ax2.legend(); ax2.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig('nuc_sweep_Al.pdf')
    print('wrote nuc_sweep_Al.pdf')
    print('\n Te(eV)   range_e   range_e+n   reduction')
    for i, te in enumerate(tes):
        print('%7g  %.4e  %.4e   %5.1f%%' %
              (te, re[i], rt[i], 100 * (1 - rt[i] / re[i])))


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'fig1'
    if what in ('fig1', 'all'):
        fig1()
    if what in ('sweep', 'all'):
        sweep()
