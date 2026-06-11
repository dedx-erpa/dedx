"""
nuclear_gk.py -- finite-temperature Gordon-Kim effective potential for the
nuclear stopping model (Faussurier, Blancard & Gauthier 2013).

This is the most sophisticated of the three pair potentials and the default for
the total stopping power.  It builds the projectile-target interaction from the
real average-atom electron density (FAC), so -- unlike the analytic ion-sphere /
Yukawa potentials which screen only with the free electrons -- it includes the
bound electrons and reproduces the cold-matter NIST nuclear stopping (Fig. 1 of
the paper).

Full Gordon-Kim electron-gas model (Gordon & Kim, J. Chem. Phys. 56, 3122,
1972).  The interaction of two atoms/ions with FROZEN spherical electron
densities rho_A (target) and rho_B (projectile) is the energy difference between
the combined and separated systems:

    U_GK(R) = U_es(R) + dT(R) + dEx(R) + dEc(R)

with, for each density functional f in {kinetic, exchange, correlation},

    dF(R) = \int [ f(rho_A+rho_B) - f(rho_A) - f(rho_B) ] d^3r          (overlap)

evaluated over the projectile electron cloud (where the integrand is nonzero).
The functionals are the standard local electron-gas forms (Hartree a.u.):

    kinetic   (Thomas-Fermi):  t(rho)  = C_F rho^(5/3),  C_F = (3/10)(3 pi^2)^(2/3)
    exchange  (Dirac):         ex(rho) = -C_x rho^(4/3),  C_x = (3/4)(3/pi)^(1/3)
    correlation (PW92 LDA):    ec(rho) = rho * eps_c^PW92(rs)

Including the kinetic/exchange/correlation overlap (the Pauli repulsion and
exchange-correlation) is the same step, in spirit, as going from RPA to eRPA in
the electronic stopping: keep the full electron-gas physics, not just the
electrostatics.

Electrostatics (with the projectile electron cloud, so the bare nuclei emerge at
small R):

    U_es(R) = Z0'*phi_t(R) - rho_B0 * \int_{cloud} phi_t(r_A) d^3s

phi_t is the screened potential of the NEUTRAL target atom built from rho_A
(integrates to Zt over the Wigner-Seitz cell):
    phi_t(s) = (Zt - M(s))/s - \int_s^{R_WS} rho_A_4pir2/r dr,   M(s)=enclosed e-.

The projectile (charge state Z0') is a uniform electron cloud of Z0' electrons
with radius D0 set by the free density at the boundary, (4/3) pi D0^3 qe(R_WS)
= Z0' (Faussurier).  rho_B0 = 3 Z0'/(4 pi D0^3).

All quantities in Hartree atomic units (Bohr, Hartree).
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
import rd


def _lens_volume(p, q, R):
    """Overlap (lens) volume of two spheres of radii p, q with centers R apart."""
    if R >= p + q:
        return 0.0
    if R <= abs(p - q):
        return (4.0 / 3.0) * np.pi * min(p, q)**3
    return (np.pi * (p + q - R)**2 *
            (R**2 + 2.0 * R * (p + q) - 3.0 * (p - q)**2) / (12.0 * R))

# electron-gas functional constants (Hartree a.u.)
C_F = 0.3 * (3.0 * np.pi**2)**(2.0 / 3.0)      # Thomas-Fermi kinetic, ~2.8712
C_X = 0.75 * (3.0 / np.pi)**(1.0 / 3.0)        # Dirac exchange,        ~0.7386


def _ec_pw92(rho):
    """PW92 LDA correlation energy density [Hartree/Bohr^3], spin-unpolarized."""
    rho = np.maximum(rho, 1e-30)
    rs = (3.0 / (4.0 * np.pi * rho))**(1.0 / 3.0)
    A, a1 = 0.031091, 0.21370
    b1, b2, b3, b4 = 7.5957, 3.5876, 1.6382, 0.49294
    sq = np.sqrt(rs)
    Q0 = -2.0 * A * (1.0 + a1 * rs)
    Q1 = 2.0 * A * (b1 * sq + b2 * rs + b3 * rs * sq + b4 * rs * rs)
    eps_c = Q0 * np.log1p(1.0 / Q1)
    return rho * eps_c


def _t_tf(rho):
    return C_F * np.maximum(rho, 0.0)**(5.0 / 3.0)


def _ex_dirac(rho):
    return -C_X * np.maximum(rho, 0.0)**(4.0 / 3.0)


def _egas(rho):
    """Total local electron-gas energy density (kinetic + exchange + correlation)."""
    return _t_tf(rho) + _ex_dirac(rho) + _ec_pw92(rho)


class GordonKim:
    def __init__(self, z0p, zbar, ne, od=None, rs=None, full=True, scale=True,
                 muffin=True, nR=80, ns=40, na=32):
        """full=True   -> full Gordon-Kim (electrostatic + kinetic/exchange/corr)
           full=False  -> electrostatic term only (for comparison)
           scale=True  -> peanut-molecule volume-conservation scaling (Faussurier/
                          Stein): expand both densities by a(R) so the overlapping
                          pair conserves total volume, removing the frozen-density
                          over-counting of the kinetic (Pauli) repulsion.
           muffin=True -> treat the target electron density beyond the cell as the
                          constant interstitial free-electron sea n0=qe(R_WS)
                          (consistent with the eRPA muffin-tin), instead of an
                          isolated neutral atom (rho_A->0). Uses the background-
                          subtracted (embedded) overlap so U still ->0 at large R;
                          implies scale=False (background and scaling don't mix)."""
        if od is None:
            raise ValueError('GordonKim needs od= (dir with rho.functions)')
        d = rd.rden(od)                       # [r, rho0, rhof, rhoe, vzt]
        r = d[0]
        rho_tot = d[1].copy()                 # 4 pi r^2 rho_total (integrates to Zt)
        rho_free = d[2]                       # 4 pi r^2 rho_free  (integrates to Zbar)
        if rs is None:
            rs = r[-1]
        self.rs = rs
        self.z0p = z0p
        self.full = full
        self.muffin = muffin
        self.scale = scale and not muffin     # background and scaling don't mix

        m = r <= rs * (1.0 + 1e-9)
        r, rho_tot, rho_free = r[m], rho_tot[m], rho_free[m]

        # enclosed electrons and total Zt
        M = cumulative_trapezoid(rho_tot, r, initial=0.0)
        Zt = M[-1]
        self.Zt = Zt
        self.qchar = z0p * Zt

        # neutral-atom electrostatic potential phi_t(s)
        g = rho_tot / r
        Gcum = cumulative_trapezoid(g, r, initial=0.0)
        T = Gcum[-1] - Gcum                    # \int_s^{rs} rho_tot/r dr
        phi = (Zt - M) / r - T
        phi = phi - phi[-1] * (r / rs)         # enforce phi(rs)=0
        self._phi_interp = interp1d(r, phi, kind='cubic',
                                    bounds_error=False, fill_value=(phi[0], 0.0))

        # target electron number density rho_A(r) = rho_tot/(4 pi r^2)
        rhoA = rho_tot / (4.0 * np.pi * r**2)
        self._rhoA_interp = interp1d(r, rhoA, kind='linear',
                                     bounds_error=False,
                                     fill_value=(rhoA[0], rhoA[-1]))
        self._rmin = r[0]
        self._rmaxA = rs
        self._rhoA_ws = rhoA[-1]               # interstitial free density
        self._n0 = self._rhoA_ws if muffin else 0.0   # density beyond the cell

        # projectile uniform electron cloud
        qe_ws = rho_free[-1] / (4.0 * np.pi * rs**2)        # qe(R_WS)
        self.D0 = (3.0 * z0p / (4.0 * np.pi * qe_ws))**(1.0 / 3.0) if qe_ws > 0 else 0.0
        self.rho_B0 = 3.0 * z0p / (4.0 * np.pi * self.D0**3) if self.D0 > 0 else 0.0
        self.rcut = rs + self.D0

        # precompute U_GK(R) on a grid and build an interpolator
        self._build_UR(nR, ns, na)

    # -- target electron density; ->0 (isolated atom) or n0 (muffin-tin) beyond cell --
    def _rhoA(self, rA):
        return np.where(rA > self._rmaxA, self._n0, self._rhoA_interp(rA))

    def _avol(self, R):
        """Volume-conservation scale factor a(R) >= 1: expand both spheres so the
        overlapping pair conserves total volume V_WS + V_D0 (Faussurier/Stein)."""
        r1, r2 = self.rs, self.D0
        v0 = r1**3 + r2**3
        def f(a):
            return (a**3 - 1.0) * v0 - (3.0 / (4.0 * np.pi)) * _lens_volume(a * r1, a * r2, R)
        if f(1.0) >= -1e-12:                 # negligible overlap
            return 1.0
        ahi = 2.0
        while f(ahi) < 0.0 and ahi < 20.0:
            ahi *= 1.5
        return brentq(f, 1.0, ahi, rtol=1e-8)

    def _build_UR(self, nR, ns, na):
        """Precompute ONLY the smooth kinetic/exchange/correlation overlap
        U_xc(R); the (singular) electrostatic Z0'*phi_t(R) is added analytically
        at call time from the fine log-grid phi_t interpolator."""
        D0, rhoB0 = self.D0, self.rho_B0
        Rgrid = np.linspace(0.01, self.rcut, nR)
        Uxc = np.zeros(nR)

        if self.full and D0 > 0:
            gx, gw = np.polynomial.legendre.leggauss(ns)
            sx = 0.5 * (gx + 1.0)                              # nodes on (0,1)
            mu, wmu = np.polynomial.legendre.leggauss(na)      # cos(alpha)
            SX, MU = np.meshgrid(sx, mu, indexing='ij')
            WSX, WMU = np.meshgrid(0.5 * gw, wmu, indexing='ij')
            n0 = self._n0
            egas_n0 = float(_egas(np.array(n0))) if n0 > 0 else 0.0
            for i, R in enumerate(Rgrid):
                a = self._avol(R) if self.scale else 1.0
                D0e = a * D0                                    # scaled cloud radius
                rhoBe = rhoB0 / a**3                            # scaled cloud density
                # projectile self-energy reference: in vacuum f(rho_B); in a
                # muffin-tin medium the background-subtracted f(n0+rho_B)-f(n0)
                egas_B = float(_egas(np.array(n0 + rhoBe))) - egas_n0
                S = SX * D0e
                jac = 2.0 * np.pi * S**2 * (WSX * D0e) * WMU    # d^3s = 2pi s^2 ds dmu
                rA = np.sqrt(R**2 + S**2 - 2.0 * R * S * MU)    # |R zhat - s|
                rhoA = self._rhoA(rA / a) / a**3                # scaled target density
                dF = (_egas(rhoA + rhoBe) - _egas(rhoA)) - egas_B
                Uxc[i] = np.sum(jac * dF)

        self._Rgrid = Rgrid
        self._Uxc_vals = Uxc
        self._Uxc_interp = interp1d(Rgrid, Uxc, kind='cubic',
                                    bounds_error=False, fill_value=(Uxc[0], 0.0))

    def _Uxc(self, R):
        if not self.full or self.D0 <= 0 or R >= self.rcut:
            return 0.0
        return float(self._Uxc_interp(R))

    def __call__(self, R):
        Rr = np.atleast_1d(np.asarray(R, dtype=float))
        out = np.empty_like(Rr)
        for i, x in enumerate(Rr):
            if x >= self.rcut:
                out[i] = 0.0
            else:
                out[i] = self.z0p * float(self._phi_interp(x)) + self._Uxc(x)
        return out if out.size > 1 else out[0]
