# Proposed paper wording — band-gap correction & §II.4 reconciliation

Draft prose to slot into the manuscript, written to match the *implemented* code and the
validation results from this work. `[NOTE]` markers flag editorial decisions for Tom.
All deviations are vs the experimental (IAEA) / PSTAR references through the Bragg-peak
region (~30–300 keV protons).

---

## A. §II.4 reconciliation (plasmon shift / local-field term)

**Issue.** The current §II.4 writes the shifted plasma frequency as
ω̃ₚ² = ωₚ² + γ²ωₚ² with a fixed γ² = 0.75. The code no longer uses a single γ²; the
bound/valence "binding" correction is applied through the energy-dependent parameters
`xepa`, `xepb` (defaults 1.0, 1.1) in `vlhfit` (and `xepc/xepd/xepe` in the strong-collision
and Barkas cutoffs).

**Proposed text (sketch — confirm physics intent before finalizing):**
> The bound-electron contribution to the dielectric response is incorporated through a
> local binding-energy shift rather than a single global plasmon-broadening constant. The
> shift is applied to the average-atom binding-energy density [c10 column] with an
> energy-dependent weight that turns the correction off above a cutoff set by the local
> binding scale (parameters `xepa`, `xepb`); this replaces the earlier fixed γ² = 0.75
> form, which is recovered approximately in the low-energy limit.

`[NOTE]` I'm not fully certain of the exact analytic relationship between the old γ² form
and the `xepa/xepb` weighting — this paragraph needs your confirmation. If γ²=0.75 is
still what you want to publish as the canonical model, then the code defaults should be
documented as the operative values instead.

---

## B. New subsection — Band-gap correction for cold insulators (Levine–Louie)

`[NOTE]` This is scoped deliberately as a **cold-insulator** correction, consistent with
what the validation supports. Dial the emphasis up/down to taste.

### B.1 Motivation
For cold ionic insulators the average-atom eRPA-LDA model overshoots the measured
electronic stopping through the Bragg peak (Fig. 17). The cause is that the model's
valence electrons are given a gapless free-electron-gas response: low-frequency
excitations that are forbidden by the material's band gap are spuriously allowed,
inflating the stopping at low projectile velocity.

### B.2 Method
We impose a Levine–Louie band gap on the dielectric response: the energy-loss function
Im[-1/ε] is zeroed below a threshold E_g and the oscillator strength is shifted upward,
preserving the f-sum rule by construction. Rather than a fixed material gap, the
threshold is tied to a *local* average-atom binding frequency,

    E_g(r) = κ · ω_b(r),

where ω_b(r) is obtained from the average-atom geometric-mean binding energy
[Eq. (7); in the code, ω_b(r) = c10(r)/r² from the .den output], and κ is a single
dimensionless calibration constant. Because ω_b derives from the AA occupations, it
evolves continuously with temperature and density. The correction is implemented by
tabulating the gapped loss function on a grid of gap frequencies and interpolating at the
local E_g(r) during the radial stopping integral.

### B.3 Results
A single κ ≈ 0.9 (κ < 1, as expected since the orbital binding ω_b exceeds the optical
gap) removes the cold overshoot for the materials where the LDA overshoots:

| target | gapless dev. | κ=0.9 dev. | optical gap |
|--------|-------------:|-----------:|------------:|
| LiF    | +16 %        | −3 %       | ~14 eV |
| SiO₂   | +12 %        | +2 %       | ~9 eV  |

The gap effect scales with the material gap (largest for LiF), emerging automatically from
ω_b with no per-material tuning (Fig. [calibration_trio]).

`[NOTE — honesty point, recommend keeping]` For several other wide-gap insulators
(Al₂O₃, CaF₂, KCl) the **gapless** model already agrees with the data to within a few
percent; there the correction is small or mildly over-corrects. This is the expected and
correct behavior — the gap only matters where sub-threshold excitation is being spuriously
allowed — but it means a single κ is an approximation, and the experimental spread
(±10–20 % between datasets) precludes a tighter calibration. We therefore present κ ≈ 0.9
as a representative value, not a precisely determined constant.

### B.4 Scope and limitations (recommend including)
- **Cold condensed regime only.** The correction targets the cold-insulator overshoot and
  is not active for the plasma/WDM targets that are the model's primary application; for
  those the band-gap concept dissolves as the material ionizes.
- **Liquid water excluded.** Water's stopping is anomalous (breakdown of Bragg additivity
  by 10–20 % near the peak, hydrogen-bond-modified dielectric response, phase effects) and
  lies outside an average-atom, Bragg-additive framework; we exclude it from the gap
  calibration and note it as a separate, known limitation.
- **Self-extinction in WDM is incomplete in the present form.** Tying E_g to the
  geometric-mean binding ω_b = c10/r² does not fully vanish at *partial* ionization
  (the surviving deep-bound electrons keep ω_b large), so the correction does not yet
  smoothly disappear across the warm-dense regime. A bound-fraction-weighted threshold
  (e.g. ω_b ∝ c10/(c06+c07)) is the natural extension and is left to future work.
  `[NOTE]` Include this only if you discuss the WDM limit; otherwise B.4 bullet 1 suffices.

---

## C. Figures
- **Fig. [calibration_trio]:** `data/LiF/dedx_calibration_trio_k09.png` — LiF/SiO₂/Al₂O₃,
  gapless vs κ=0.9 vs reference.
- **(optional) counterexample:** `data/CaF2/dedx_caf2_gap.png` (gapless already matches),
  `data/KCl/dedx_kcl_gap.png`.
- **(only if WDM discussed):** `data/LiF/wdm_fadeout_mechanism.png`,
  `data/LiF/wdm_stopping_panels.png`.
