# Draft sections — nuclear/ion stopping, total stopping power, and verification

Prose to slot into the manuscript, written to match the *implemented* code and the
verification results from this work. `[NOTE]` markers flag editorial decisions for Tom
(e.g. whether this lives in the gap-correction paper as new sections, or a companion paper).
Companion to `gap_correction_paper_draft.md` (electronic band-gap correction).

Validation numbers below are reproduced by the drivers `nuc_test.py`, `nuc_fac_compare.py`,
`nuc_ion_verify.py`, and `alpha_dt_verify.py`; equations and implementation are documented
in `nuclear_model_notes.md`.

---

## A. Motivation — a total stopping power from cold solids to fusion plasmas

The eRPA average-atom model gives the *electronic* stopping power. The energy a projectile
loses to the target *ions* is a separate channel that the electronic model does not capture:
at low temperature it is the classical elastic ion–ion ("nuclear") stopping familiar from
cold-matter range tables; at finite temperature it becomes the projectile-ion/target-ion
drag of the Fokker–Planck slowing-down (the ion term of Mehlhorn 1981). We add this channel
so that the package returns a *total* stopping power — electronic plus nuclear/ionic —
that is consistent from cold solids, through warm and hot dense matter, and on to fully
ionized fusion plasmas. The two contributions are returned in the same units
(10⁻¹⁵ eV cm²/atom) and add directly.

`[NOTE]` Framing choice for Tom: present this as "completing the total stopping power" of
the existing eRPA model, or as a standalone nuclear/ion-stopping capability. The text below
is written for the former.

## B. Method — nuclear/ion stopping

We follow the classical two-body model of Faussurier, Blancard & Gauthier (Phys. Plasmas
**20**, 012705, 2013), whose electronic component is the same eRPA model used here (their
Ref. 6 is Wang–Mehlhorn–MacFarlane 1998), so the two channels are mutually consistent.

**Scattering integral.** The per-atom nuclear stopping cross section is obtained from the
classical deflection of the projectile in an effective screened ion–ion pair potential
U(R), integrated over impact parameter. At zero ion temperature,
ε_n(E) = 2π γ E ∫ sin²(θ/2) p dp with γ = 4M₀Mₜ/(M₀+Mₜ)². The orbit/turning-point integral
is evaluated with a substitution that removes the inverse-sqrt singularity; the
implementation reproduces the analytic Rutherford cross section to < 10⁻⁶.

**Effective pair potentials.** Three potentials are provided. Two are analytic and screen
with the *free* electrons only: an ion-sphere form and a Yukawa/Debye–Hückel form with a
finite-temperature Thomas–Fermi screening length. The recommended, default potential is the
finite-temperature **Gordon–Kim** electron-gas potential, built from the *real* average-atom
electron density (bound + free):

U_GK(R) = U_es(R) + ΔT(R) + ΔE_x(R) + ΔE_c(R),

the electrostatic interaction of the neutral target atom with the projectile plus the
kinetic (Thomas–Fermi), exchange (Dirac), and correlation (PW92 LDA) *overlap* energies of
the frozen densities, with the Faussurier/Stein volume-conservation ("peanut-molecule")
scaling. Because U_es is built from the full average-atom density, R·U → Z₀'Zₜ (the bare
nuclear product) at small R, so the Gordon–Kim potential recovers the unscreened
nucleus–nucleus interaction at close approach — which is why it reproduces the cold-matter
reference data while the free-electron-screened analytic potentials fall below it.

`[NOTE]` Including the kinetic/exchange/correlation overlap (not just the electrostatic
term) is the same modelling choice, in spirit, as going from RPA to eRPA in the electronic
stopping: retain the full electron-gas physics. Worth one sentence to that effect.

**Average-atom density boundary (eRPA-consistent).** The Gordon–Kim overlap requires the
target density beyond the Wigner–Seitz cell. We use the muffin-tin interstitial
free-electron sea n₀ = qe(R_WS), consistent with the eRPA electronic treatment, with a
background-subtracted (embedded) overlap so that U → 0 at large R. This differs from an
isolated-neutral-atom truncation by < 0.5 % in cold matter and by up to ~5 % in hot dense
matter (where the interstitial density is appreciable).

**Finite ion temperature.** The target ions carry a Maxwellian velocity distribution, over
which the energy exchange is averaged. Below a threshold E* the projectile is slower than
the thermal ions and the *net* ion–ion energy exchange becomes a gain — the projectile is
heated toward equilibrium (it thermalizes). This signed energy *exchange* is physical, but
it is distinct from a stopping *power*, which is an energy-*loss* rate (≥ 0); in the
stopping-power output we therefore floor the nuclear term at zero below E* (so the total
returns to the electronic value where the projectile has thermalized). The signed exchange
is retained internally and is what bounds the alpha thermalization range in §E.

## C. Correction — low-velocity electronic stopping (LSS/Chandrasekhar limit)

`[NOTE]` This is a correction to the electronic model; place it wherever the electronic
stopping is described, not necessarily in the nuclear section.

Below the lowest tabulated projectile energy (≈ 1 keV/amu) the table-based electronic
modes previously clamped the RPA stopping number to the grid-edge value, producing a
spurious dE/dx ∝ 1/E rise in warm/hot matter (the curve turned *up* at low energy). The
correct low-velocity limit is the Lindhard/LSS friction L ∝ v³, i.e. dE/dx ∝ v ∝ √E. We
impose this scaling below the table edge; the change is confined to E < 1 keV/amu (every
validated energy is unchanged, e.g. cold-Al protons vs PSTAR are byte-identical), and the
warm-matter low-energy slope d(ln dE/dx)/d(ln E) goes from −1.0 to the correct +0.50,
matching the analytic stopping-number path.

## D. Verification — nuclear stopping and the Fokker–Planck ion drag

**Cold-matter reference (Fig. nuc_fig1_coldAl).** For protons in cold solid-density
aluminum the Gordon–Kim nuclear stopping reproduces the NIST/ICRU-49 (ZBL) nuclear stopping
to within ~7 % across keV–MeV (the residual is consistent with the known frozen-density
over-repulsion of Gordon–Kim and within ZBL's own ~10 % fit uncertainty). The analytic
ion-sphere and Yukawa potentials fall below the reference, as in the source paper, because
they screen with the free electrons only.

**Proton range vs PSTAR (Fig. range_nuc).** Because the PSTAR CSDA range includes nuclear
stopping, it is a parameter-free cold-target test of the *total* (electronic + nuclear)
model. For protons in C, Al, Ag, and Au the electronic-only range overshoots the PSTAR CSDA
range at the lowest energies — by +66 % (C), +30 % (Al), +23 % (Ag), and +6 % (Au) at 1 keV
— because the nuclear channel is absent. Adding the Gordon–Kim nuclear term brings C, Al,
and Ag into agreement with PSTAR to within a few percent, and slightly overcorrects gold
(≈ −10 %), the same frozen-density over-repulsion seen in the cold-Al nuclear comparison.
This is an independent confirmation that the nuclear contribution carries the correct
magnitude and energy dependence at the level that controls ranges.

**FAC vs. SCAALP (Fig. nuc_fac_compare).** The original model used the SCAALP average-atom
model; here the densities come from FAC. For the recommended Gordon–Kim potential the two
agree (and both land on NIST), because the potential is built from the *total* electron
density, which integrates to Z independent of the average-atom model. The free-electron
mean ionization Z̄(Al) from FAC (2.97 cold → 12.8 at 1 keV, "nearly stripped") matches the
SCAALP behaviour. The model-dependence shows up only in the Z̄-based analytic potentials.

**Fokker–Planck ion drag (Fig. nuc_ion_verify).** The finite-T ion stopping is the
projectile-ion/target-ion term of the Fokker–Planck slowing-down. Driven by the same
Maxwellian average and a Coulomb-log cross section, the model reproduces the analytic
single-particle energy-gain threshold erf(x)/(x erf'(x)) = 1 + mₜ/m₀ in both mass-ratio
regimes: for He in a fully-ionized H plasma at Tₑ = Tᵢ = 100 eV it gives E* = 134 eV,
matching the MD-validated value of 132.7 eV in the source paper; for protons in ionized
aluminum at 1 keV it gives E* = 164 eV (analytic 155 eV). The high-velocity magnitude
reproduces the textbook ion–ion stopping 4πZ₀²Zₜ²e⁴nᵢ lnΛ/(mₜv²) with an effective Coulomb
logarithm lnΛ ≈ 4–5.

## E. Application — alpha self-heating in DT (Fig. alpha_dt_verify)

As an end-to-end test of the total stopping power we treat a 3.5 MeV fusion alpha slowing
in a DT plasma, integrating the slowing-down from 3.5 MeV to thermalization (the net
stopping crosses zero near (3/2)kTₑ). The calculation reproduces the established ICF result
(Fraley 1974; Mehlhorn 1981): the alpha deposits energy predominantly to electrons at low
temperature and predominantly to ions above ~30 keV, with the electron/ion deposition
fractions crossing at **28.9–32.6 keV** for densities 0.213–100 g/cc, and an alpha range
ρ_R of order **0.5 g/cm²**. The crossover arises because the electronic stopping falls as
~Tₑ^(−3/2) when the alpha becomes slow relative to the hot electrons, while the ion drag
stays roughly flat (the alpha remains fast relative to the DT ions).

## F. Comparison with recent plasma / WDM stopping measurements

Beyond the cold (NIST/PSTAR/IAEA) and analytic (Fokker–Planck) references, the total
model is compared to four recent, well-characterized WDM/plasma measurements (conditions
and citations in `docs/recent_stopping_experiments.md`; drivers `exp_compare.py` and
`exp_compare_icf.py`).

**Frenje 2019 — ion stopping in DT plasma (Tₑ≈2 keV): the nuclear/ion channel.** Frenje et
al. found that the electronic (Brown–Preston–Singleton) stopping underpredicts the measured
ion stopping near the Bragg peak and attributed the gap to nuclear-elastic scattering and
possibly coupled ion modes — the projectile-ion/target-ion channel added here. We digitized
their Fig. 5 (measured −ΔEᵢ/Zᵢ² vs Eᵢ/Aᵢ and the BPS curve) and overlaid the model for a
proton in DT plasma at Tₑ=2 keV, with a single areal-density scale fixed by matching the
model electronic stopping to BPS (the path-integrated profile, their Fig. 2, was not
digitized). The model + nuclear reproduces the measured stopping **within error at three of
the four probe energies** (0.61, 2.4, 15.7 MeV/u). At the lowest velocity (0.19 MeV/u) the
nuclear/ion channel lifts the stopping ~12% above the electronic value — in the direction of
the measured low-velocity excess Frenje attributed to nuclear-elastic scattering — and the
nuclear fraction grows toward lower velocity (≈23% of the total at 0.1 MeV/u, below the
measured range). A residual gap remains at the lowest point, consistent with Frenje's
additional suggestion of coupled ion modes (not modeled here). (Fig. exp_frenje.)

The single areal-density scale used in the overlay is **not a free parameter**: digitizing
their Fig. 2 Tₑ(r)/nₑ(r) profiles gives a path-integrated electron areal density of
3.6–7.1×10²¹ cm⁻² (radius to diameter chord; ρL ≈ 15–30 mg/cm²), and the scale that aligns
the model electronic stopping with BPS corresponds to 5.5×10²¹ cm⁻² — squarely within that
range. The model magnitude is therefore consistent with the measured plasma profile.

`[NOTE]` This is the strongest experimental motivation for the nuclear/ion addition. A fully
first-principles −ΔE per probe would integrate the slowing-down over the (shell-peaked)
profile with Frenje's emission/transport geometry; the shape comparison + the profile-
confirmed scale already make the case without that machinery.

**Cayzac 2017 — N ions in carbon plasma (Tₑ≈150 eV).** At full ionization the model
reproduces the measured plasma stopping enhancement: dE/dx(plasma)/dE/dx(solid) = 1.62 at
0.586 MeV/u vs the measured enhancement of up to ~1.5. (Fig. exp_compare, right.)

**Malko 2022 — proton in warm dense carbon (Tₑ≈7.5 eV).** The model reproduces the
measured/DFT stopping magnitude across 0.2–0.65 MeV (e.g. 0.40 MeV: 0.44 vs 0.44; 0.51 MeV:
0.38 vs 0.39 keV/(µg/cm²)). `[NOTE]` The measured ~13–26% reduction vs solid is small in the
model at this temperature (carbon barely ionizes further, Z̄ 2.7→3.0) — as in classical and
TD-DFT treatments; we note it rather than claim it. (Fig. exp_compare, left.)

**Zylstra 2015 — proton in warm dense Be (Tₑ≈32 eV): a parameter-free absolute test.** For
14.7 MeV D³He protons through the 94 mg/cm² Be foil (areal density taken from the paper, no
fitting), the model reproduces the measured energy downshift: ΔE = 2.87 (cold) and 2.97 MeV
(warm) vs measured 2.72 and 2.85–2.98 MeV (≈5% high on cold), and it matches Zylstra's own
AA-LDA values (2.80, 2.87) — the framework they found consistent with the data — to within
~3%. The model's warm/cold enhancement (+3.3%) agrees with Zylstra's AA-LDA (+2.5%), both at
the low end of the measured spread. Because the foil is flat with a known areal density,
this is a genuine parameter-free absolute comparison. (Fig. exp_zylstra.)

**Community benchmark — charged-particle stopping workshops.** The model is run on the
standardized cases of the two code-comparison workshops (Grabowski et al. 2020; Stanek et
al. 2024): alpha-particle stopping vs velocity for H, C, Be, Al, Cu at the defined ρ, T.
The model predictions (alpha electronic stopping, 10⁻¹⁵ eV cm²/atom) are:

| Case | target, ρ, T | Z̄ | peak Se | vp at peak (a.u.) | Se at vp=2 a.u. |
|---|---|---|---|---|---|
| H1  | H, 1 g/cc, 2 eV     | 0.79 | 14.9 | 2.1 | 14.7 |
| C1  | C, 10 g/cc, 2 eV    | 4.00 | 31.7 | 3.0 | 25.0 |
| Be1 | Be, 1.84 g/cc, 4.4 eV | 1.97 | 45.3 | 1.7 | 44.1 |
| Al1 | Al, 2.7 g/cc, 1 eV  | 2.97 | 87.1 | 1.5 | 79.9 |
| Cu1 | Cu, 8.96 g/cc, 1 eV | 4.46 | 101.2 | 2.1 | 100.7 |

Overlaying the eRPA-AA electronic stopping directly on the Stanek Fig. 12 TD-DFT-MD
reference points (Fig. exp_stanek; H, C, Al) reproduces the behavior Stanek et al. report
for average-atom models: the curves converge with the data beyond the stopping peak and sit
above the data near and below the peak, by 24–55 % at the peak (H 1.30, C 1.24, Al 1.55;
ratio model/data). The deviation grows with the size of the bound-electron complement,
largest for Al.

A sharper test is available for carbon, where Stanek Fig. 12 also reports two of Stanek's
*own* average-atom curves. Plotted on the same velocity axis (Fig. exp_stanek_C_AA), the
eRPA-AA peak (1.59×10¹⁰ eV/cm) lands *between* Stanek's two AA curves (solid 1.53×10¹⁰,
dot-dash 1.69×10¹⁰), and all three average-atom models sit 19–32 % above the TD-DFT-MD data
at the peak. The three average-atom high-velocity tails converge on the shared Bethe limit.
Our independent eRPA-AA therefore reproduces the workshop average-atom family across the full
velocity range, confirming that the ~25 % excess over TD-DFT-MD near the peak is a property
of the average-atom treatment as a class, not of our particular implementation.

**Where the average-atom stopping comes from — shell-resolved depth.** Because the eRPA-AA
integrates dE/dx over the full bound+free radial electron density, the contribution can be
decomposed by radius — equivalently by enclosed electron number, core to valence (Fig.
exp_stanek_depth). For Al this directly addresses the average-atom-vs-TD-DFT-MD comparison:
TD-DFT-MD propagates the valence electrons explicitly and *adds an estimated deep-core
contribution*, whereas the eRPA-AA includes all electrons self-consistently with DFT-like
local-field corrections. The decomposition shows the median stopping depth is the loosely
bound valence/conduction electrons at the Wigner–Seitz edge at low velocity (r₅₀ ≈ 2.4 a.u.)
and moves steadily *inward* as the projectile speeds up — the core re-engaging — reaching
≈1.0 a.u. only by 1 MeV/u. Crucially, throughout the Stanek velocity range (the Bragg peak
and above, where the average-atom excess appears) the median depth lies in the valence/free
+ outer-L region, **not** the deep K-shell. The average-atom-vs-TD-DFT-MD difference near the
peak is therefore attributable to the treatment of the valence/free electrons (where the
methods' local-field and exchange–correlation choices differ), not to an over-counting of
deep-core electrons — refining the common assumption that average-atom models run high near
the peak because they add inner-shell contributions a valence-only treatment omits.

## G. Figures

- **Fig. nuc_fig1_coldAl** (`nuc_fig1_coldAl.pdf`): proton in cold Al — ion-sphere, Yukawa,
  Gordon–Kim vs NIST/ZBL nuclear stopping.
- **Fig. nuc_fac_compare** (`nuc_fac_compare.pdf`): FAC implementation vs the paper — cold-Al
  three potentials vs NIST, and the Te-dependence with FAC Z̄.
- **Fig. nuc_ion_verify** (`nuc_ion_verify.pdf`): finite-T ion stopping vs the analytic
  Fokker–Planck drag for He–H and p–Al, with the energy-gain thresholds.
- **Fig. alpha_dt_verify** (`alpha_dt_verify.pdf`): alpha range ρ_R(Te) and the electron/ion
  deposition crossover at ~30 keV.
- **Fig. nuc_gap_4panel** (`nuc_gap_4panel.pdf`, optional): total = best-κ electronic +
  Gordon–Kim nuclear vs PSTAR/experiment for LiF, SiO₂, Al₂O₃, H₂O — shows the nuclear
  contribution filling part of the low-energy deficit.
- **Fig. exp_stanek** (`exp_stanek.pdf`): eRPA-AA alpha electronic stopping overlaid on the
  Stanek 2024 Fig. 12 TD-DFT-MD reference points for H, C, Al.
- **Fig. exp_stanek_C_AA** (`exp_stanek_C_AA.pdf`): carbon detail — eRPA-AA vs Stanek's two
  own average-atom curves and the TD-DFT-MD data on the labeled vp[cm/s] axis.
- **Fig. exp_stanek_depth** (`stopping_depth.pdf`): shell-resolved origin of the Al eRPA-AA
  stopping — median stopping depth vs velocity (valence-dominated at low v, moving inward as
  the projectile speeds up) and cumulative stopping fraction vs radius at several velocities.
- **Fig. range_nuc** (`range_nuc.pdf`): proton CSDA range in C, Al, Ag, Au relative to PSTAR
  — electronic-only vs electronic + Gordon–Kim nuclear; the nuclear term removes the
  low-energy overshoot of the electronic-only range.

## H. Scope and limitations (recommend including)

- The Gordon–Kim potential is a frozen-density electron-gas model; its ~7 % offset from the
  ZBL cold-nuclear reference is the expected sign and magnitude of frozen-density
  over-repulsion. The analytic ion-sphere/Yukawa potentials are provided for comparison and
  for fully-ionized conditions; the Yukawa form overestimates at high projectile energy
  (long screened-attractive tail), where the nuclear stopping is in any case negligible.
- The finite-T ion drag is the screened-potential generalization of the Fokker–Planck ion
  term; in the slow-projectile/heavy-field regime it differs from a fixed-lnΛ treatment
  because the screened cross section regularizes the Coulomb-log cutoff self-consistently.
- The projectile charge state is taken as fully stripped (exact for the proton/alpha cases
  here); a velocity-dependent effective charge is a natural extension.
