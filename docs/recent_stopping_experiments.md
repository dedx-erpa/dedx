# Recent plasma / warm-dense-matter ion-stopping experiments — comparison resource

Candidate modern measurements for benchmarking the total eRPA stopping model
(electronic + nuclear/ion), beyond the cold IAEA/PSTAR references. Conditions for
Cayzac 2017 and Malko 2022 are taken from the papers (open-access); other entries
should have their citations/conditions verified before publication.

`[NOTE]` Both headline WDM experiments report **energy loss through a foil** at a
narrow projectile energy (plasma vs solid), not a wide dE/dx(E) curve — so the
natural comparison metric is the energy loss (or the plasma/solid ratio) at the
stated conditions, which is what the comparisons in this repo compute.

## Summary table

| Experiment | Projectile (energy) | Target / conditions | Regime | Measured result | Models supported / ruled out |
|---|---|---|---|---|---|
| **Cayzac 2017** (Nat. Commun. 8, 15693; GSI PHELIX+nhelix / UNILAC) | N ions, 0.586 MeV/u (degraded from 3.6) | C plasma, 100 µg/cm²; nₑ≈5×10²⁰ cm⁻³, Tₑ≈150 eV (fully ionized) | near Bragg peak, hot ionized | plasma energy loss enhanced up to ~50% vs solid (e.g. 1.14 vs 0.83 MeV) | disproved several standard models; supported strong ion–electron collision treatments |
| **Malko 2022** (Nat. Commun. 13, 2893; CLPU VEGA-II) | proton, 498±4 keV | warm dense C, 130 µg/cm²; Tₑ≈7.5 eV (mass-wt), Z*≤4, ρ≳0.1 g/cc, Γ≈0.1–2 | low velocity, below Bragg peak | ΔE_WDM = 39.4±5 keV vs ΔE_solid = 49±5 keV (13–26% **lower**) | matches TD-KS-DFT (51 keV); classical models too high (54–55 keV) |
| **Zylstra 2015** (PRL 114, 215002; OMEGA) | energetic protons | isochorically-heated **solid-density Be** plasma, Tₑ≈32 eV; Γ≈0.3, kTₑ/E_F≈2 (WDM) | around Bragg peak | **increased** energy loss vs cold matter (reduced mean ionization potential) | **agrees with average-atom LDA** and ad-hoc free+bound models — first WDM test |
| **Frenje 2019** (PRL 122, 015002; OMEGA) | DT-α / fusion ions, scan through vᵢ/v_th | **DT plasma**, Tₑ=1.4–2.8 keV, nₑ=4–8×10²³ cm⁻³ | low-v → Bragg peak → high-v | BPS best near peak, but **underpredicts ~20% at vᵢ≈0.3v_th** | BPS/Li–Petrasso; **postulates nuclear-elastic scattering** (= our nuclear/ion term) to close the gap |
| **GSI partially-ionized Ar / H** `[verify cite]` | protons | Ar (partially ionized) and H (fully ionized); Tₑ~few×100 eV, nₑ 10²⁰–10²¹ | medium energy | bound-electron contribution in partially-ionized plasma | free- vs bound-electron treatments |

## Community benchmark cases (code-comparison workshops)

The charged-particle transport workshops (Grabowski et al., HEDP 2020, "first"; Stanek
et al. 2024, "second") define standardized plasma conditions where many codes computed
stopping power — a community benchmark the model can be run against directly.

- **First workshop (Grabowski 2020):** pure **H**, pure **C**, and equimolar **CH**, on the
  grid ρ = 0.1, 1, 10, 100 g/cm³ × Tₑ = 0.2, 2, 20, 200, 2000 eV (60 cases; Z̄ from More TF).
  Proton/light-ion stopping vs velocity. The low-density high-T H case ≈ the Omega
  stopping experiments; the low-density moderate-T C case ≈ heated-Be Omega shots.
- **Second workshop (Stanek 2024), Table I (Priority 1 + subset):**

  | Case | species | n (cm⁻³) | ρ (g/cc) | T (eV) |
  |---|---|---|---|---|
  | H1 | H | 5.98×10²³ | 1 | 2 |
  | C1 | C | 5.01×10²³ | 10 | 2 |
  | CH1 | CH | 4.63×10²² (each) | 1 | 2 |
  | Al1 | Al | 6.03×10²² | 2.7 | 1 |
  | Cu1 | Cu | 8.49×10²² | 8.96 | 1 |
  | HCu1 | HCu | 1.68×10²² | 1.8 | 1 |
  | Be1 | Be | 1.23×10²³ | 1.84 | 4.4 |
  | CH2 | CH | 4.16×10²² | 0.9 | 7.8 |

  Finding: order-of-magnitude spread between approaches at low T, ~factor-5 among
  first-principles models — i.e. large model uncertainty exactly where the AA-eRPA +
  nuclear model is positioned. `[NOTE]` Consensus values live in the workshop figures;
  the model can be run on every case, with digitized consensus bands overlaid.

## Model-comparison status (this repo)

| Experiment | Data in repo? | Model run | Metric |
|---|---|---|---|
| Malko 2022 | yes (`malko_fig.txt`) | proton in C, solid vs Tₑ≈7.5 eV | model matches data magnitude (`exp_compare.py`) |
| Cayzac 2017 | headline numbers only | N (Z=7) in C plasma, Tₑ=150 eV vs solid | model plasma/solid = 1.62 vs measured ~1.5 (`exp_compare.py`) |
| Zylstra 2015 | headline numbers | proton in Be, Tₑ=32 eV vs cold | model gives ~5% increased loss; agrees with AA-LDA finding (`exp_compare_icf.py`) |
| Frenje 2019 | headline numbers | proton in DT, Tₑ=2 keV | **nuclear/ion = 23% of total at vᵢ≈0.3vₜₕ — matches the ~20% Frenje attributes to nuclear-elastic scattering** (`exp_compare_icf.py`) |

`[NOTE]` For Cayzac and the ICF-plasma sets we have the conditions but not the
digitized data points; the model predictions are produced here, and the
experimental points can be overlaid once digitized (or pulled from the papers'
supplementary data).

## Why these matter for the model

- **Malko (cold-to-warm, low velocity)** tests the electronic stopping below the
  Bragg peak in the WDM regime — exactly where the band-gap correction and the
  finite-T electronic response operate. The measured *reduction* vs solid is the
  key discriminator (classical models go the wrong way).
- **Cayzac (hot, near Bragg peak)** tests the model at full ionization where the
  ion-stopping (nuclear/ionic) channel and the strong ion–electron collision
  terms matter most.
- Together they bracket the WDM regime the model targets, complementing the cold
  IAEA/PSTAR/NIST validation already in the package.
