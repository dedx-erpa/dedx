"""
nuclear.py -- nuclear (elastic ion-ion) stopping power for cold matter -> WDM.

Implements the classical two-body model of Faussurier, Blancard & Gauthier,
Phys. Plasmas 20, 012705 (2013), as a companion to the eRPA *electronic*
stopping in this package (the electronic dE/dx|e here is the Wang-Mehlhorn-
MacFarlane 1998 model the paper adopts as its Ref. 6).

The nuclear stopping cross section per target atom is returned in the SAME
units as the electronic column written by dedx.f / dedx.py:
    10^-15 eV cm^2 / atom
so the electronic and nuclear contributions add directly.

Three effective screened pair potentials are provided (Eqs. 4, 5 and the
finite-T Gordon-Kim construction), and both the Ti=0 form (Eq. 10) and the
finite ion-temperature Maxwellian average (Eq. 6).

All internal physics is done in Hartree atomic units (length: Bohr, energy:
Hartree, mass: electron mass m_e), then converted on output.

References to "Eq. (n)" below are to Faussurier et al. (2013).
"""

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# physical constants / unit conversions
# ---------------------------------------------------------------------------
HARTREE_EV = 27.211386245988      # 1 Hartree in eV
BOHR_CM    = 0.529177210903e-8    # Bohr radius in cm
BOHR2_CM2  = BOHR_CM**2           # Bohr^2 in cm^2  (~2.8003e-17)
AMU_ME     = 1822.888486209       # 1 atomic mass unit in electron masses

# Hartree*Bohr^2  ->  (10^-15 eV cm^2)
#   eps[10^-15 eV cm^2] = eps[Ha Bohr^2] * HARTREE_EV * BOHR2_CM2 / 1e-15
EPS_UNIT = HARTREE_EV * BOHR2_CM2 / 1e-15     # ~0.76208

# fixed Gauss-Legendre nodes on (0,1) for the (singularity-removed) orbit
# integral; nodes are interior so the integrand is evaluated away from the
# t=0 / t=1 endpoints.
_GL_N = 64
_gx, _gw = np.polynomial.legendre.leggauss(_GL_N)
_GL_T = 0.5 * (_gx + 1.0)        # nodes on (0,1)
_GL_W = 0.5 * _gw               # weights on (0,1)


# ---------------------------------------------------------------------------
# plasma helpers
# ---------------------------------------------------------------------------
def Nion(rs):
    """Ion (target-atom) number density [1/Bohr^3] from Wigner-Seitz radius rs [Bohr]."""
    return 3.0 / (4.0 * np.pi * rs**3)


def screening_k(ne, te_ha):
    """Thomas-Fermi screening wavenumber k [1/Bohr] for a free-electron gas.

    Finite-temperature interpolation that reduces to the Debye wavenumber at
    high T and the degenerate Thomas-Fermi wavenumber at T->0:
        k^2 = 4 pi ne / sqrt( (kB Te)^2 + (2 EF/3)^2 ),   EF = (1/2)(3 pi^2 ne)^(2/3)
    ne in 1/Bohr^3, te_ha = kB Te in Hartree.
    """
    ef = 0.5 * (3.0 * np.pi**2 * ne)**(2.0 / 3.0)
    denom = np.sqrt(te_ha**2 + (2.0 * ef / 3.0)**2)
    return np.sqrt(4.0 * np.pi * ne / denom)


# ---------------------------------------------------------------------------
# effective pair potentials  U(R) [Hartree], R [Bohr]
# Each is exposed as a small class with __call__(R) and a `rcut` attribute
# (R beyond which U==0, used to bound the impact-parameter integral).
# ---------------------------------------------------------------------------
class IonSphere:
    """Eq. (4): ion-sphere screened potential.

    U(R) = (Z0p*Zbar / R) [1 - 3R/(2 Rn) + R^3/(2 Rn^3)],  R <= Rn ; else 0
    Rn = [3 (Z0p + Zbar) / (4 pi Ne)]^(1/3),  Ne = Zbar * Ni.
    Charges are the *ion* charges: Z0p (projectile charge state) and Zbar
    (target mean ionization).
    """
    def __init__(self, z0p, zbar, ne):
        self.q = z0p * zbar
        self.qchar = z0p * zbar          # small-R Coulomb scale (sets close-collision p)
        self.Rn = (3.0 * (z0p + zbar) / (4.0 * np.pi * ne))**(1.0 / 3.0)
        self.rcut = self.Rn

    def __call__(self, R):
        x = R / self.Rn
        u = (self.q / R) * (1.0 - 1.5 * x + 0.5 * x**3)
        return np.where(R <= self.Rn, u, 0.0)


class Yukawa:
    """Eq. (5): Yukawa / Debye-Hueckel screened potential.

    U(R) = (Z0p*Zbar / R) (1 - kR/2) exp(-kR)
    k = Thomas-Fermi screening wavenumber.
    """
    def __init__(self, z0p, zbar, k):
        self.q = z0p * zbar
        self.qchar = z0p * zbar
        self.k = k
        # cut where exp(-kR) has decayed well past any contribution
        self.rcut = 40.0 / k

    def __call__(self, R):
        kR = self.k * R
        return (self.q / R) * (1.0 - 0.5 * kR) * np.exp(-kR)


class Coulomb:
    """Bare (unscreened) Coulomb -- for Rutherford verification only."""
    def __init__(self, q):
        self.q = q
        self.qchar = q
        self.rcut = np.inf

    def __call__(self, R):
        return self.q / R


# ---------------------------------------------------------------------------
# classical deflection angle  (Eqs. 7-9)
# ---------------------------------------------------------------------------
def _Rmin(p, Ec, U):
    """Largest root of f(R) = 1 - p^2/R^2 - U(R)/Ec = 0  (distance of closest approach)."""
    def f(R):
        return 1.0 - (p / R)**2 - U(R) / Ec
    # For p beyond the potential range there is no scattering: turning point = p.
    if p >= U.rcut:
        return p
    # f(p+) = -U(p)/Ec < 0 for a repulsive core; f(R->large) -> 1 > 0.
    lo = p * (1.0 + 1e-12)
    # start from a finite upper bracket (U.rcut may be inf for bare Coulomb)
    hi = p * 1.0 if not np.isfinite(U.rcut) else max(U.rcut, p)
    # ensure a sign change; expand hi if needed
    flo = f(lo)
    if flo >= 0.0:
        # potential already negligible at lo: essentially no deflection
        return p
    fhi = f(hi)
    ntry = 0
    while fhi < 0.0 and ntry < 60:
        hi *= 1.5
        fhi = f(hi)
        ntry += 1
    return brentq(f, lo, hi, rtol=1e-10, maxiter=200)


def deflection(p, Ec, U):
    """Center-of-mass deflection angle theta(p) [rad] for potential U at relative
    energy Ec (Eq. 7), via the s=Rmin/R, s=1-t^2 substitution that removes the
    turning-point inverse-sqrt singularity."""
    if p >= U.rcut:
        return 0.0
    Rm = _Rmin(p, Ec, U)
    if Rm <= 0.0:
        return 0.0
    pr = p / Rm
    t = _GL_T
    s = 1.0 - t * t
    R = Rm / s
    g = 1.0 - (pr * s)**2 - U(R) / Ec
    g = np.maximum(g, 1e-300)
    integ = np.sum(_GL_W * (2.0 * t) / np.sqrt(g))
    theta = np.pi - 2.0 * pr * integ
    return theta


def _sin2_halftheta(p, Ec, U):
    th = deflection(p, Ec, U)
    s = np.sin(0.5 * th)
    return s * s


# ---------------------------------------------------------------------------
# impact-parameter integral  S(Ec) = \int_0^pmax sin^2(theta/2) p dp   [Bohr^2]
# (the velocity/energy-independent geometric factor; depends on Ec via U/Ec)
# ---------------------------------------------------------------------------
def S_of_Ec(Ec, U, npp=300):
    """Transport-like impact-parameter integral  S = \int_0^pmax sin^2(theta/2) p dp
    [Bohr^2].

    The integrand peaks near the characteristic impact parameter p_c = qchar/(2 Ec)
    (where theta ~ 90 deg).  We use a log grid spanning well below p_c up to the
    screening cutoff so both the close-collision core and the screened tail are
    resolved at any energy; the p<p_lo remainder is added analytically assuming
    sin^2(theta/2)->1 there (head-on limit)."""
    pmax = U.rcut
    if not np.isfinite(pmax):
        pmax = 50.0
    pc = U.qchar / (2.0 * Ec) if Ec > 0 else pmax
    p_lo = min(1e-4 * pmax, 1e-3 * pc)
    p_lo = max(p_lo, 1e-12)
    if p_lo >= pmax:
        return 0.5 * pmax * pmax
    p = np.geomspace(p_lo, pmax, npp)
    integrand = np.array([_sin2_halftheta(pp, Ec, U) * pp for pp in p])
    S = np.trapezoid(integrand, p)
    # p < p_lo: head-on collisions, sin^2(theta/2) ~ 1  ->  contributes p_lo^2/2
    S += 0.5 * p_lo * p_lo
    return S


def make_S_interp(U, ec_min, ec_max, n=80, npp=300):
    """Precompute S(Ec) on a log-Ec grid once and return a fast interpolator
    S(Ec) [Bohr^2].  S is computed via S_of_Ec; outside [ec_min, ec_max] it is
    held constant (S saturates at low Ec and is tiny at high Ec)."""
    ec_min = max(ec_min, 1e-12)
    if ec_max <= ec_min:
        ec_max = ec_min * 10.0
    ecs = np.geomspace(ec_min, ec_max, n)
    svals = np.array([S_of_Ec(ec, U, npp) for ec in ecs])
    f = interp1d(np.log(ecs), svals, kind='cubic', bounds_error=False,
                 fill_value=(svals[0], svals[-1]))
    return lambda ec: float(f(np.log(max(ec, 1e-300))))


# ---------------------------------------------------------------------------
# Ti = 0 nuclear stopping cross section  (Eq. 10)
# ---------------------------------------------------------------------------
def eps_n_Ti0(E_lab_ha, m0_au, mt_au, Sfunc):
    """Per-atom nuclear stopping cross section at zero ion temperature.

    Eq. (10): eps_n = 2 pi * gamma * E * \int sin^2(theta/2) p dp
    with gamma = 4 m0 mt /(m0+mt)^2 and E the projectile lab kinetic energy.
    `Sfunc(Ec)` returns the precomputed impact-parameter integral [Bohr^2].
    Returns Hartree*Bohr^2 (convert with EPS_UNIT for output units).
    """
    gamma = 4.0 * m0_au * mt_au / (m0_au + mt_au)**2
    Ec = (mt_au / (m0_au + mt_au)) * E_lab_ha     # Eq. (8) with w=0
    return 2.0 * np.pi * gamma * E_lab_ha * Sfunc(Ec)


# ---------------------------------------------------------------------------
# finite-Ti nuclear stopping cross section  (Eq. 6)
# ---------------------------------------------------------------------------
def eps_n_Maxwell(E_lab_ha, m0_au, mt_au, Sfunc, ti_ha, nwpar=41, nwperp=21):
    """Per-atom nuclear stopping cross section, Maxwellian-averaged over the
    target-ion velocity distribution (Eq. 6).

    Reduces the 3-D velocity integral to (w_par, w_perp) by azimuthal symmetry
    about the projectile direction.  `Sfunc(Ec)` is the precomputed impact-
    parameter integral.  Returns Hartree*Bohr^2.  Can be negative below the
    energy-gain threshold E* (the projectile gains energy from a hot ion bath),
    as discussed around Eqs. (11)-(12).
    """
    V0 = np.sqrt(2.0 * E_lab_ha / m0_au)
    if ti_ha <= 0.0:
        return eps_n_Ti0(E_lab_ha, m0_au, mt_au, Sfunc)
    mu = m0_au * mt_au / (m0_au + mt_au)
    wth = np.sqrt(ti_ha / mt_au)                  # thermal speed scale

    # velocity grids: target ions are Maxwellian about REST (w=0), spread wth
    wpar = np.linspace(-5.0 * wth, 5.0 * wth, nwpar)
    wperp = np.linspace(0.0, 5.0 * wth, nwperp)

    norm = (mt_au / (2.0 * np.pi * ti_ha))**1.5
    pref = 4.0 * np.pi**2          # (2pi scattering azimuth) * (2pi velocity azimuth)
    cfac = 2.0 * m0_au * mt_au / (m0_au + mt_au)**2

    # 2-D trapezoid over (wpar, wperp) with the w_perp Jacobian
    WP, WPP = np.meshgrid(wpar, wperp, indexing='ij')
    gx = V0 - WP
    gmag = np.sqrt(gx**2 + WPP**2)
    f3 = norm * np.exp(-mt_au * (WP**2 + WPP**2) / (2.0 * ti_ha))
    # g . (m0 V0 + mt w) = (V0-wpar)(m0 V0 + mt wpar) - mt wperp^2
    gdotP = gx * (m0_au * V0 + mt_au * WP) - mt_au * WPP**2
    W = cfac * (gmag / V0) * gdotP
    Smat = np.vectorize(lambda g: Sfunc(0.5 * mu * g * g))(gmag)
    integ = f3 * Smat * W * WPP        # w_perp Jacobian
    inner = np.trapezoid(integ, wperp, axis=1)
    return pref * np.trapezoid(inner, wpar)


# ---------------------------------------------------------------------------
# top-level driver
# ---------------------------------------------------------------------------
def build_potential(name, z0p, zbar, ne, te_ha, od=None, rs=None, gk_muffin=True):
    """Construct one of the three effective potentials by name."""
    name = name.lower()
    if name in ('ionsphere', 'is', 'ion-sphere'):
        return IonSphere(z0p, zbar, ne)
    if name in ('yukawa', 'yk', 'debye'):
        return Yukawa(z0p, zbar, screening_k(ne, te_ha))
    if name in ('gk', 'gordonkim', 'gordon-kim'):
        from nuclear_gk import GordonKim
        return GordonKim(z0p, zbar, ne, od=od, rs=rs, muffin=gk_muffin)
    raise ValueError('unknown potential: %s' % name)


def nuclear_stopping(E_grid_MeVamu, zp, zt, m0_amu, mt_amu, rs, te_eV, ti_eV,
                     zbar, potential='ionsphere', od=None, gk_muffin=True):
    """Nuclear stopping cross section on an energy grid.

    Parameters
    ----------
    E_grid_MeVamu : array, projectile energy per AMU [MeV] (the dedx.dat x-axis)
    zp, zt        : projectile / target nuclear charge
    m0_amu,mt_amu : projectile / target mass [amu]
    rs            : Wigner-Seitz radius [Bohr]
    te_eV, ti_eV  : electron / ion temperature [eV]  (ti_eV<=0 -> Ti=0 form)
    zbar          : target mean ionization
    potential     : 'ionsphere' | 'yukawa' | 'gk'
    od            : output dir (needed for 'gk' to read rho.functions)

    Returns eps_n [10^-15 eV cm^2 / atom], same array shape as E_grid_MeVamu.
    """
    E = np.atleast_1d(np.asarray(E_grid_MeVamu, dtype=float))
    m0_au = m0_amu * AMU_ME
    mt_au = mt_amu * AMU_ME
    ne = zbar * Nion(rs)
    te_ha = te_eV / HARTREE_EV
    ti_ha = ti_eV / HARTREE_EV

    z0p = float(zp)                    # projectile charge state (fully stripped)
    U = build_potential(potential, z0p, zbar, ne, te_ha, od=od, rs=rs,
                        gk_muffin=gk_muffin)

    # lab kinetic energies [Hartree] on the grid
    E_lab = E * m0_amu * 1e6 / HARTREE_EV
    mu = m0_au * mt_au / (m0_au + mt_au)
    V0 = np.sqrt(2.0 * E_lab / m0_au)

    # determine the full center-of-mass-energy span, then precompute S(Ec) once
    if ti_eV > 0.0:
        wth = np.sqrt(ti_ha / mt_au)
        gmin = np.maximum(1e-6, np.abs(V0 - 5.0 * wth)).min()
        gmax = np.sqrt((5.0 * wth)**2 + (V0 + 5.0 * wth)**2).max()
        ec_min = 0.5 * mu * gmin**2
        ec_max = 0.5 * mu * gmax**2
    else:
        ecs = (mt_au / (m0_au + mt_au)) * E_lab
        ec_min, ec_max = ecs.min(), ecs.max()
    Sfunc = make_S_interp(U, ec_min, ec_max)

    out = np.zeros_like(E)
    for i in range(len(E)):
        if ti_eV > 0.0:
            eps = eps_n_Maxwell(E_lab[i], m0_au, mt_au, Sfunc, ti_ha)
        else:
            eps = eps_n_Ti0(E_lab[i], m0_au, mt_au, Sfunc)
        out[i] = eps * EPS_UNIT
    return out if out.size > 1 else out[0]


# ---------------------------------------------------------------------------
# output: augment a finished dedx.dat with nuclear columns -> dedx_nuc.dat
# ---------------------------------------------------------------------------
_POT_ORDER = ['gk', 'yukawa', 'ionsphere']     # preference for the "total" column


def write_nuclear(od, zp, te_eV, ti_eV, potlist, species, fout='dedx_nuc.dat',
                  gk_muffin=True):
    """Compute nuclear stopping on the grid of od/dedx.dat and write dedx_nuc.dat.

    Parameters
    ----------
    od       : output dir holding the (electronic) dedx.dat to augment
    zp       : projectile nuclear charge
    te_eV    : electron temperature [eV]
    ti_eV    : ion temperature [eV] (<=0 -> Ti=0 form)
    potlist  : list of potentials, subset of ['ionsphere','yukawa','gk']
    species  : list of dicts {'zt':int, 'w':float, 'sod':dir} -- one per target
               atomic species; 'sod' is the dir with that species' dedx.dat and
               rho.functions (== od for a single element).  Nuclear contributions
               are number-weighted across species, mirroring the electronic part.
    """
    from pfac import fac
    import rd
    data = rd.rdedx(od)                      # [E/AMU, dedx_e, range, (barkas)]
    E = data[0]
    dedx_e = data[1]
    m0_amu = fac.ATOMICMASS[int(zp)]

    cols = {p: np.zeros_like(E) for p in potlist}
    wtot = 0.0
    for sp in species:
        h = rd.rdedx(sp['sod'], header='')
        zt = int(sp['zt'])
        mt_amu = fac.ATOMICMASS[zt]
        for p in potlist:
            eps = nuclear_stopping(E, zp, zt, m0_amu, mt_amu, float(h['rs']),
                                   te_eV, ti_eV, float(h['zbar']),
                                   potential=p, od=sp['sod'], gk_muffin=gk_muffin)
            cols[p] += np.atleast_1d(eps) * sp['w']
        wtot += sp['w']
    for p in potlist:
        cols[p] /= wtot
        # A stopping power is an energy-LOSS rate, so floor at 0.  The finite-Ti
        # Maxwellian average (Eq. 6) goes negative below the energy-gain threshold
        # E* (the projectile is slower than the thermal ions and is heated by the
        # bath -- it has thermalized, so there is no net stopping there).  That
        # signed energy-exchange is real physics but a different quantity; it is
        # kept in eps_n_Maxwell (used by alpha_dt_verify for the thermalization
        # zero-crossing) and only floored here for the stopping-power table.
        cols[p] = np.maximum(cols[p], 0.0)

    # choose the potential used for the total/range column
    chosen = next((p for p in _POT_ORDER if p in potlist), potlist[0])
    total = dedx_e + cols[chosen]
    # range with the same target-mass convention used by dedx.py
    am = sum(fac.ATOMICMASS[int(sp['zt'])] * sp['w'] for sp in species) / wtot
    rng = rd.int_range(np.array((E, total)), m=am)

    path = '%s/%s' % (od, fout)
    with open(path, 'w') as f:
        f.write('#    zp = %10.4E\n' % zp)
        f.write('#    Te = %15.8E\n' % te_eV)
        f.write('#    Ti = %15.8E\n' % ti_eV)
        f.write('# npot = %s\n' % ','.join(potlist))
        f.write('# total = electronic + %s\n' % chosen)
        hdr = '#   E/AMU(MeV)        dEdx_e'
        for p in potlist:
            hdr += '   dEdx_n[%s]' % p
        hdr += '     dEdx_tot         range\n'
        f.write(hdr)
        for i in range(len(E)):
            f.write('%15.8E %13.6E' % (E[i], dedx_e[i]))
            for p in potlist:
                f.write(' %13.6E' % cols[p][i])
            f.write(' %13.6E %13.6E\n' % (total[i], rng[i]))
    return path


# ---------------------------------------------------------------------------
# plotting: electronic / nuclear(ionic) / total from a dedx_nuc.dat
# ---------------------------------------------------------------------------
def plot_nuclear(od, fin='dedx_nuc.dat', fout='dedx_nuc.pdf', show_range=True):
    """Plot electronic, nuclear/ionic, and total stopping (and the total range)
    from od/dedx_nuc.dat.  Returns the saved figure path."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    path = '%s/%s' % (od, fin)
    meta, potlist = {}, []
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if 'npot' in line:
                potlist = [s.strip() for s in line.split('=')[1].split(',') if s.strip()]
            elif '=' in line:
                k = line[1:].split('=')[0].strip()
                try:
                    meta[k] = float(line.split('=')[1])
                except ValueError:
                    meta[k] = line.split('=')[1].strip()

    d = np.loadtxt(path, comments='#')
    E = d[:, 0]
    dedx_e = d[:, 1]
    nuc = {p: d[:, 2 + i] for i, p in enumerate(potlist)}
    dedx_tot = d[:, 2 + len(potlist)]
    rng = d[:, 3 + len(potlist)]

    # material label from the electronic dedx.dat header, if present
    mat = ''
    try:
        import rd
        h = rd.rdedx(od, header='')
        zt = ','.join('%d' % int(round(z)) for z in np.atleast_1d(h['zt']))
        mat = 'Zt=%s, ' % zt
    except Exception:
        pass
    cond = '%srho-dep, Te=%g eV, Ti=%g eV, Zp=%g' % (
        mat, meta.get('Te', 0), meta.get('Ti', 0), meta.get('zp', 1))

    ncol = 2 if show_range else 1
    fig, axes = plt.subplots(1, ncol, figsize=(6.2 * ncol, 4.6))
    ax = axes[0] if show_range else axes

    # nuclear columns are floored at 0 (stopping = energy loss); mask the zeros
    # so the curves drop out cleanly on the log axis rather than implying a value.
    pos = lambda y: np.where(y > 0, y, np.nan)
    ax.loglog(E, dedx_e, 'C0-', lw=2, label='electronic (eRPA)')
    styles = ['C3-', 'C2--', 'C1-.', 'C4:']
    for i, p in enumerate(potlist):
        ax.loglog(E, pos(nuc[p]), styles[i % len(styles)], lw=1.7,
                  label='nuclear/ionic [%s]' % p)
    ax.loglog(E, pos(dedx_tot), 'k-', lw=2.4, alpha=.8,
              label='total (e + %s)' % meta.get('total', 'nuc').split('+')[-1].strip())
    ax.set_xlabel('E / AMU (MeV)')
    ax.set_ylabel(r'dE/dx ($10^{-15}$ eV cm$^2$/atom)')
    ax.set_title('Stopping power')
    ax.legend(fontsize=8); ax.grid(True, which='both', alpha=.3)

    if show_range:
        axr = axes[1]
        axr.loglog(E, rng, 'k-', lw=2)
        axr.set_xlabel('E / AMU (MeV)')
        axr.set_ylabel('range (mg/cm$^2$)')
        axr.set_title('Total range')
        axr.grid(True, which='both', alpha=.3)

    fig.suptitle(cond, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = '%s/%s' % (od, fout)
    fig.savefig(out)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# self-test: Rutherford analytic check of the scattering integrator
# ---------------------------------------------------------------------------
def _rutherford_check():
    """Verify deflection() against the analytic Rutherford relation
    sin^2(theta/2) = 1 / (1 + (2 Ec p / q)^2) for bare Coulomb U=q/R."""
    q = 1.0
    Ec = 0.5
    U = Coulomb(q)
    print('Rutherford check (q=%.1f, Ec=%.2f):' % (q, Ec))
    print('   p       sin2_num     sin2_exact     rel.err')
    ok = True
    for p in (0.5, 1.0, 2.0, 5.0, 10.0):
        num = _sin2_halftheta(p, Ec, U)
        exact = 1.0 / (1.0 + (2.0 * Ec * p / q)**2)
        rel = abs(num - exact) / exact
        ok = ok and (rel < 1e-4)
        print('  %5.1f   %.6e   %.6e   %.2e' % (p, num, exact, rel))
    print('  PASS' if ok else '  FAIL')
    return ok


if __name__ == '__main__':
    _rutherford_check()
