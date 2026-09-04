# SPDX-License-Identifier: GPL-3.0-or-later
"""
V&V diagnostic: does the ionic-channel stopping recover the Rutherford 1/v^2
high-velocity limit?

For a fast projectile the elastic ion-ion stopping must fall as 1/v^2, i.e.
dEdx_n * Etot should be ~constant with energy (up to the slowly-growing Coulomb
log).  This is the single fastest check that you have chosen the right ionic pair
potential for your target:

  * ion-sphere (correct for a fully-ionized plasma): dEdx_n * Etot ~ flat.
  * Gordon-Kim on a fully-ionized plasma (WRONG regime): dEdx_n * Etot climbs
    ~30x over a decade of energy, because the neutral-atom exchange-correlation
    overlap adds a fixed large-impact-parameter feature that does not shrink with
    velocity.  That is the bug the runtime guard warns about.

This script uses the analytic ion-sphere potential directly (no FAC pipeline
needed), so it runs in seconds.  Example (alpha in a DT plasma at the DT electron
density, Te=Ti=10 keV):

    python examples/check_rutherford_limit.py --zp 2 --zt 1 --d 40.3 --t 10000

Expect the last column (dEdx_n * Etot) to be roughly flat; a factor >~3 climb
over the tabulated span is the neutral-atom artifact -- switch to --npot=ionsphere
(and, for DT, add --imass=2.014,3.016).
"""
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nuclear as N


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--zp', type=float, default=2.0, help='projectile charge (2=alpha)')
    ap.add_argument('--zt', type=int, default=1, help='target Z (1=H surrogate for DT)')
    ap.add_argument('--mt', type=float, default=None,
                    help='target ion mass amu (default: atomic mass of zt; use 2.5 for DT)')
    ap.add_argument('--mp', type=float, default=4.0026, help='projectile mass amu')
    ap.add_argument('--d', type=float, default=40.3, help='mass density g/cc (H-surrogate)')
    ap.add_argument('--t', type=float, default=10000.0, help='temperature eV (Te=Ti)')
    a = ap.parse_args()

    n_ion = a.d / 1.008 * 6.02214076e23 * (0.529177e-8) ** 3   # 1/Bohr^3 (H surrogate)
    rs = (3.0 / (4 * np.pi * n_ion)) ** (1.0 / 3.0)
    mt = a.mt if a.mt is not None else 1.008 * a.zt
    E = np.geomspace(0.05, 100.0, 12)                          # MeV/u
    eps = N.nuclear_stopping(E, a.zp, a.zt, a.mp, mt, rs, a.t, a.t, float(a.zt),
                             potential='ionsphere', od=None)
    Etot = E * a.mp
    prod = np.atleast_1d(eps) * Etot
    print('ion-sphere ionic stopping, zp=%g zt=%d d=%g g/cc T=%g eV  (rs=%.3f Bohr, mt=%.3f amu)'
          % (a.zp, a.zt, a.d, a.t, rs, mt))
    print(' E/A(MeV/u)   dEdx_n[1e-15 eVcm2/atom]   dEdx_n*Etot  (flat => Rutherford 1/v^2)')
    for i in range(len(E)):
        print('  %8.3f      %.4e                 %.4e' % (E[i], eps[i], prod[i]))
    climb = prod[-1] / prod[np.argmin(np.abs(E - 1.0))]
    verdict = 'OK (Rutherford recovered)' if climb < 3.0 else 'ARTIFACT -- wrong potential for this regime'
    print('\n dEdx_n*Etot climb from 1 -> %.0f MeV/u: %.1fx   => %s' % (E[-1], climb, verdict))


if __name__ == '__main__':
    main()
