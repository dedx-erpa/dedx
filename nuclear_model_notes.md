# Nuclear stopping power model (cold matter → WDM)

Implementation notes for `nuclear.py` / `nuclear_gk.py`, which add the nuclear
(elastic ion–ion) stopping power `dE/dx|n` to the eRPA *electronic* stopping
`dE/dx|e` already produced by `dedx.f` / `dedx.py`.

Model: Faussurier, Blancard & Gauthier, *Phys. Plasmas* **20**, 012705 (2013).
The paper's electronic component is the Wang–Mehlhorn–MacFarlane (1998) eRPA
model — i.e. exactly what this package computes — so adding the nuclear
component completes a self-consistent **total** stopping power from cold solids
through warm/hot dense matter.

All internal physics is in Hartree atomic units (Bohr, Hartree, electron mass);
outputs are converted to the package convention `10⁻¹⁵ eV cm²/atom` so that
electronic and nuclear contributions add directly.

## Equations

Per-target-atom nuclear stopping cross section:

- **Ti = 0** (Eq. 10):   `ε_n(E) = 2π γ E ∫₀^{p_max} sin²(θ/2) p dp`,
  `γ = 4 M₀Mₜ/(M₀+Mₜ)²`, `E` = projectile lab kinetic energy.
- **Finite Ti** (Eq. 6): average the same integrand over a 3-D Maxwellian
  target-ion velocity `w`, reduced by azimuthal symmetry to a 2-D quadrature
  over `(w∥, w⊥)`. Below a threshold `E*` in hot plasma the net ion-ion energy
  exchange `ε_n` becomes a **gain** (the projectile is slower than the thermal
  ions and is heated toward equilibrium — it has thermalized) — Eqs. (11)–(12),
  Faussurier Fig. 4. This is the net energy-*exchange* rate, not a stopping power:
  a stopping power is an energy-*loss* rate (≥ 0), so `write_nuclear` floors the
  table at 0 below `E*` (total = electronic there). The signed exchange is kept in
  `eps_n_Maxwell` and used by `alpha_dt_verify` for the thermalization zero-
  crossing.
- **Deflection** (Eq. 7): `θ(p) = π − 2(p/R_min)∫₀¹ 2t dt/√g(1−t²)`,
  `g(s)=1 − p²s²/R_min² − U(R_min/s)/E_c`, `s=R_min/R`. The `s=1−t²` substitution
  removes the turning-point √ singularity. `R_min` (Eq. 9) is the largest root of
  `1 − p²/R² − U(R)/E_c = 0` (`brentq`); the orbit integral uses 64-pt
  Gauss–Legendre. `E_c = ½μ|V₀−w|²` (Eq. 8).

The impact-parameter integral `S(E_c) = ∫ sin²(θ/2) p dp` is the only expensive
piece; it is precomputed once on a log-`E_c` grid per potential and splined, so
both the Ti=0 and Maxwell paths are fast.

## Three effective pair potentials `U(R)`

1. **Ion-sphere** (Eq. 4):
   `U = (Z₀'Z̄/R)[1 − 3R/2Rₙ + R³/2Rₙ³]`, `R≤Rₙ`; `Rₙ=[3(Z₀'+Z̄)/4πNₑ]^{1/3}`,
   `Nₑ=Z̄Nᵢ`, `Nᵢ=3/4πrs³`. Free-electron screening only.
2. **Yukawa** (Eq. 5): `U = (Z₀'Z̄/R)(1−kR/2)e^{−kR}`, `k` = finite-T Thomas–Fermi
   screening wavenumber `k²=4πNₑ/√((kTe)²+(2E_F/3)²)` (Debye at high T, degenerate
   TF at T→0). Free-electron screening only.
3. **Gordon–Kim** (`nuclear_gk.py`, the recommended/default model): the full
   electron-gas pair potential `U = U_es + ΔT + ΔE_x + ΔE_c`.
   - `U_es = Z₀'·φ_t(R)`, `φ_t` = screened potential of the **neutral** target
     atom from the real average-atom total electron density `ρ_t(r)` (the `rho0`
     column of `rho.functions`, integrating to `Zₜ` over the cell). Recovers
     `R·U → Z₀'Zₜ` (bare nuclei) at small R, `→0` outside the cell.
   - `ΔT, ΔE_x, ΔE_c` = the kinetic (Thomas–Fermi `C_F ρ^{5/3}`), exchange
     (Dirac `−C_x ρ^{4/3}`) and correlation (PW92 LDA) **overlap** energies of
     the frozen target density and the projectile electron cloud (a uniform
     sphere of `Z₀'` electrons, radius `D₀` from `(4π/3)D₀³ qe(R_WS)=Z₀'`),
     `ΔF(R)=∫[f(ρ_A+ρ_B)−f(ρ_A)−f(ρ_B)]d³r`. Including these is, in spirit, the
     same step as RPA→eRPA for the electronic part: keep the full electron-gas
     physics, not just the electrostatics.
   - **Volume-conservation scaling** (Faussurier/Stein peanut-molecule): both
     densities are expanded by `a(R)≥1` so the overlapping pair conserves total
     volume, removing the frozen-density over-counting of the Pauli (kinetic)
     repulsion. (`scale=True` default; `a(R)∈[1,1.12]` for cold Al.)
   Result for cold proton-Al: GK lands within ~7% of NIST/ZBL — the kinetic
   overlap adds ~5–8% repulsion over the electrostatic-only term, the scaling
   trims ~1–2%, and the residual is within ZBL's own ~10% fit uncertainty.
   `full=False` (electrostatic only) and `scale=False` are available for
   comparison.

### Target boundary: muffin-tin (default) vs isolated atom

The GK overlap needs the target electron density beyond the Wigner-Seitz cell.
Two choices:
- **isolated neutral atom**: `ρ_A → 0` for `r > R_WS` (the atom in vacuum).
- **muffin-tin** (`muffin=True`, default): `ρ_A → n0 = qe(R_WS)`, the constant
  interstitial free-electron sea, consistent with how the eRPA *electronic*
  stopping treats the medium. To avoid the projectile cloud spuriously
  interacting with the infinite background (which would leave `U` non-zero at
  large R), the overlap is **background-subtracted** (embedded):
  `ΔF = ∫[f(ρ_A+ρ_B) − f(ρ_A) − (f(n0+ρ_B) − f(n0))]d³r`, which → 0 at large R.
  (Background and volume-scaling don't compose, so `muffin` implies `scale=False`.)

Effect on the nuclear stopping (proton-Al):
- **Cold** (`qe(R_WS)=0.022/Bohr³`): muffin vs isolated differ by **≤0.5%** — the
  interstitial sea is too dilute to matter; the NIST match is unchanged.
- **Hot WDM** (`Te=1000 eV`, `qe(R_WS)=0.112/Bohr³`): muffin gives **~5% lower**
  nuclear stopping at 1 keV — the denser sea pre-fills the overlap region, so the
  marginal Pauli repulsion when the projectile cloud meets the target is reduced.

The muffin-tin is the default for self-consistency with the eRPA electronic
treatment; it matters only in hot dense matter and there gives the (slightly
softer, physically appropriate) embedded-in-medium result.

Charge convention: `Z₀'` = projectile charge state (taken fully stripped = `zp`;
velocity-dependent effective charge is a future refinement), `Z̄` = target mean
ionization (Eqs. 4/5), `Zₜ` = target nuclear charge (the GK small-R limit). This
convention was **validated**, not assumed (below).

## Validation

- **Scattering core**: `python nuclear.py` checks `deflection()` against the
  analytic Rutherford relation `sin²(θ/2)=1/(1+(2E_c p/q)²)` for bare Coulomb —
  agreement < 1e-6.
- **Cold NIST (Fig. 1)**: `python nuc_test.py fig1` — Gordon–Kim reproduces the
  NIST/ZBL proton-in-Al nuclear stopping to within ~5% over keV→100 keV, while
  ion-sphere and Yukawa fall below it, exactly as in the paper. This confirms the
  charge/screening convention. Plot: `nuc_fig1_coldAl.pdf`; reference table:
  `data/refs/Al_nuclear_zbl.dat`.
- **Finite-Ti threshold**: `ε_n` changes sign (energy gain) below `E*` at high
  Ti, consistent with the Butler–Buckingham behavior of the paper.
- **WDM trend (Figs. 2/3)**: `python nuc_test.py sweep` — at fixed density the
  nuclear stopping grows with Te (reduced screening, higher ionization) and cuts
  the low-energy proton range. Plot: `nuc_sweep_Al.pdf`.

## FAC implementation vs the original paper (SCAALP)

The paper builds the electron density and screening from the SCAALP average-atom
model; this package uses FAC (`pfac`). Comparison for the paper's canonical case,
proton in Al (`python nuc_fac_compare.py`):

**Agrees with the paper / NIST:**
- **Gordon–Kim** reproduces the NIST/ZBL cold nuclear stopping (FAC-GK ~5% high;
  paper-GK "good agreement"). This is robust because GK uses the *total* electron
  density, which integrates to Z (the neutral-atom constraint) independent of the
  AA model — so SCAALP-GK and FAC-GK both land on NIST.
- Gross ionization: FAC `Z̄(Al)` = 2.97 (cold) → 3.2 (10 eV) → 7.9 (100 eV) →
  12.8 (1000 eV, "nearly stripped"), matching the paper's SCAALP statements.
- ion-sphere and Yukawa fall below NIST at cold; nuclear grows with Te and only
  matters at low energy (the paper's Figs 1–2 conclusions).

**Differs from the paper:**
- **ion-sphere vs Yukawa ordering at room temperature is reversed**: FAC gives
  Yukawa > ion-sphere; the paper has Yukawa < ion-sphere. This traces to the
  screening wavenumber — our finite-T Thomas–Fermi interpolation built from the
  FAC free-electron density differs from the paper's specific TF model — and to
  `Z̄`, which both analytic potentials use as the interacting charge.
- **FAC `Z̄(Te)` is shell-resolved** (sharp M-/L-shell ionization thresholds)
  vs SCAALP's smoother TF-like ionization, so the `Z̄`-based potentials
  (ion-sphere, Yukawa) differ most at intermediate T. At 10 keV, cold→1000 eV the
  nuclear stopping grows ×1.7 (GK), ×13 (ion-sphere), ×23 (Yukawa) — the analytic
  potentials are far more sensitive to the AA model than GK.
- Our GK is the **full** electron-gas potential (electrostatic + TF kinetic +
  Dirac exchange + PW92 correlation, with volume-conservation scaling); it lands
  ~7% above NIST/ZBL, consistent with the known frozen-density over-repulsion of
  Gordon–Kim and within ZBL's own ~10% fit uncertainty.

**Bottom line:** for the recommended **GK** potential the FAC implementation
agrees with the paper and NIST; the FAC-vs-SCAALP differences show up in the
`Z̄`-based analytic potentials (ion-sphere/Yukawa), which are exactly the
model-sensitive ones — another reason to use GK for the total stopping power.

## Verification: finite-T ion stopping vs the Fokker-Planck ion-ion drag

At low T the ion term is the classical "nuclear" stopping; at finite T it is the
projectile-ion / target-ion drag of the Fokker-Planck slowing-down (cf. Mehlhorn,
J. Appl. Phys. 52, 6522, 1981). `nuc_ion_verify.py` checks our finite-Ti model
(`eps_n_Maxwell`, Faussurier Eq. 6) against the analytic FP drag for protons on
ionized aluminum and He on a hydrogen plasma.

- **Slowing-down machinery is exact.** Fed the pure Coulomb-log cross section,
  `eps_n_Maxwell` reproduces the analytic FP single-particle energy-gain
  threshold `erf(x)/(x erf'(x)) = 1+m2/m1` (x=v/v_th,ion) for both mass-ratio
  regimes: He-on-H 133.8 vs analytic 131.0 (vs Faussurier's MD-validated
  **132.7 eV**); p-on-Al 164 vs analytic 155. This validates the Eq.-6
  Maxwellian average as the projectile-ion/target-ion FP term.
- **Fast/loss regime matches in magnitude.** The high-velocity limit reproduces
  the textbook `−dE/dx = 4π Z1²Z2²e⁴ n2 lnΛ/(m2 v²)` with an effective
  lnΛ ≈ 4–5; the model and the velocity-dependent-lnΛ FP overlay in the loss
  regime for He-on-H and p-on-Al.
- **Threshold sensitivity (heavy field).** For p-on-Al the screened-potential
  threshold (~0.5×(3/2)kT_i) sits above the fixed-lnΛ analytic (~0.1×(3/2)kT_i):
  the energy-gain regime is governed by slow collisions, where the screened
  cross section (finite, no lnΛ divergence) differs from the fixed Coulomb-log.
  The screened treatment is the more rigorous; the analytic is recovered when the
  field is light/fast (He-on-H). GK, ion-sphere, and Yukawa all give the same
  screened threshold (~740 eV at Te=Ti=1000 eV), confirming this is a cross-
  section-treatment effect, not a potential-specific artifact.

Conclusion: the model reproduces the Fokker-Planck ion stopping where it is
benchmarked (He-on-H, MD-validated; and the universal fast limit), and provides a
self-consistent screened-cross-section ion drag through the cold→WDM transition.

## Application: alpha self-heating in DT (3.5 MeV He-4)

`alpha_dt_verify.py` applies the model to the canonical ICF case -- a 3.5 MeV
fusion alpha slowing in a DT plasma -- where electronic and ion stopping are the
same Fokker-Planck drag with the field being the electrons or the (D,T) ions.
The alpha slowing-down is integrated from 3.5 MeV to thermalization (net stopping
-> 0 near (3/2)kTe) for several densities and Te = 1-300 keV.

Reproduces the known result (cf. Fraley 1974; Mehlhorn 1981):
- **Electron/ion deposition crossover at ~30 keV.** efract falls 0.97 -> 0.05 and
  ifract rises 0.03 -> 0.95 across 1-300 keV, crossing at **28.9-32.6 keV** for
  rho = 0.213-100 g/cc (alpha energy goes mostly to ions above ~30 keV, because
  the electronic stopping drops ~ Te^-3/2 as the alpha becomes slow vs the hot
  electrons while the ion drag stays ~ flat).
- **Range rho_R ~ 0.5 g/cm^2** scale, rising from ~0.03 g/cm^2 (1 keV) through
  ~0.5 g/cm^2 (~10 keV) with the correct density ordering -- within ~2x of the
  published curves at the extremes (Coulomb-log-level treatment).

The ion term here is the same projectile-ion/target-ion drag our screened
nuclear model reproduces (verified in nuc_ion_verify.py); the electronic term can
be refined with the eRPA dielectric (dief.tabdedx at keV) for a fully
package-native calculation, but the FP-level result already reproduces the
30 keV crossover and the rho_R ~ 0.5 g/cm^2 range.

## Known limitation

The analytic **Yukawa** potential overestimates the nuclear stopping at high
projectile energy because of its long screened-attractive tail (`U<0` for
`R>2/k`). Nuclear stopping is negligible compared to electronic there, so this
does not affect totals; **Gordon–Kim is the recommended potential** and the one
used for the `dEdx_tot` column.
