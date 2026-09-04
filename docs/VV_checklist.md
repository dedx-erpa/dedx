# V&V checklist and gotchas

A short list of the checks and pitfalls that matter when you apply dedx-erpa to a
new target. Runnable anchors live in [`../examples/`](../examples/).

## 1. Pick the ionic pair potential by regime (the #1 pitfall)

| target | `--npot` | why |
|---|---|---|
| cold / warm dense matter (bound electrons) | `gk` (Gordon-Kim) | validated vs NIST/ZBL; neutral-atom model is correct |
| fully-ionized plasma (bare ions) | `ionsphere` | recovers the Rutherford 1/v² high-velocity limit |
| — (avoid) | `yukawa` | Debye length ≫ ion-sphere radius at solid density → over-screens |

Using `gk` on a fully-ionized plasma over-predicts the ionic stopping ~30–75× at
high projectile velocity (a fixed attractive well at the cell radius that does not
shrink with velocity). The code emits a **RuntimeWarning** if you do this, and the
`total`/range column auto-selects `ionsphere` for a (nearly) fully-ionized target.

## 2. The Rutherford-limit check

For a fast projectile the ionic term must fall as 1/v², i.e. **`dEdx_n × Etot`
should be ~flat** with energy (up to the slowly-growing Coulomb log). A factor
≳ 3 climb over a decade means the wrong potential for the regime. Run:

```bash
python examples/check_rutherford_limit.py --zp 2 --zt 1 --d <g/cc> --t <eV>
```

## 3. Ion masses in a plasma — use the real isotopes

The elastic ion-ion term goes as **1/m_target**. A Z=1 electron-density-matched
hydrogen surrogate (proton mass 1.008) over-weights a DT ion channel by ~2.2×.
Pass the true masses with `--imass=2.014,3.016 --iwt=1,1`; this overrides the
ionic-channel mass only — the electronic AA, ne, and range normalization stay on
the surrogate.

## 4. Range units and the H-surrogate conversion

- The range columns are the **physical (per-ion) range in mg/cm²**, not per-nucleon.
- With a Z=1 hydrogen surrogate at a matched electron density, the range is on the
  **surrogate** mass basis; multiply by `(A_target/A_H per electron)` for the real
  areal density — e.g. **× 2.515/1.008** for DT (`--d=40.3` H = DT 100 g/cc,
  `--d=171.6` H = DT 429 g/cc).
- Physical-unit sanity: `Etot[MeV]` = `E/A × mass number`; the file gives R(E) with
  both in real units.

## 5. Sanity numbers to expect

- A fast light ion's ionic channel is a **small fraction** of the electronic
  stopping (few % at birth), rising toward the ion thermal velocity.
- The **detour** (projected/CSDA) is ~0.8–1.0 for a fast light ion; a value near
  0.5 for a 3.5 MeV alpha signals the ionic-potential artifact.
- 3.5 MeV alpha in DT 100 g/cc, 10 keV (ionsphere + real D,T masses):
  ρR_α ≈ **0.51 g/cm²** (BPS 0.45, Zylstra-MD 0.46, classic Fraley/Atzeni ~0.3).

## 6. Corrections and their regimes

- **Bloch (Z⁴)** is a fast-projectile (Bethe-regime) correction; it is switched off
  for sub-thermal projectiles (v < v_th) so it does not over-correct a slow alpha.
- **Strong-collision (Zwicknagel) + local-field** activate only at strong coupling
  (Γ ≳ 1); they are inactive in weakly-coupled plasmas (e.g. Cayzac, Γ ≈ 0.01),
  where the model reduces to bare RPA.
- **Barkas** is suppressed at low projectile energy and negligible in hot plasma.

## 7. What NOT to trust quantitatively

- **`alpha_dt_verify.py`** is a rough analytic (Chandrasekhar) cross-check only. A
  spurious ×½ in it was corrected (ranges were ~2× long); still, prefer
  `dedx.py --nuc` for production numbers.
