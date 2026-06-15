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
| **Zylstra 2015** (PRL 114, 215002; OMEGA) `[verify cite]` | protons / D³He α | ICF-relevant HED plasma, Tₑ~keV | around/above Bragg peak | charged-particle stopping in dense plasma | tested classical/QM plasma stopping |
| **Frenje 2015** (PRL 115, 205001; OMEGA) `[verify cite]` | protons | HED plasma, Tₑ~0.5–4 keV | around Bragg peak | ion stopping near the Bragg peak | Li–Petrasso / BPS-type models |
| **"Enhanced/reduced stopping regimes" 2018** (Sci. Rep. 8, 14586) `[verify cite]` | protons | hot plasma | velocity-dependent | enhanced and reduced stopping regimes vs velocity | finite-T stopping models |
| **GSI partially-ionized Ar / H** `[verify cite]` | protons | Ar (partially ionized) and H (fully ionized); Tₑ~few×100 eV, nₑ 10²⁰–10²¹ | medium energy | bound-electron contribution in partially-ionized plasma | free- vs bound-electron treatments |

## Model-comparison status (this repo)

| Experiment | Data in repo? | Model run | Metric |
|---|---|---|---|
| Malko 2022 | yes (`malko_fig.txt`) | proton in C, solid vs Tₑ≈7.5 eV | dE/dx [keV/(µg/cm²)] and ΔE reduction |
| Cayzac 2017 | headline numbers only | N (Z=7) in C plasma, Tₑ=150 eV, ρ≈1.7e-3 vs solid | ΔE plasma/solid enhancement |
| Zylstra/Frenje | not yet | proton/α in HED DT or D³He plasma | dE/dx near Bragg peak |

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
