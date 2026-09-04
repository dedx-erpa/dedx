# Worked examples & verification

Three runnable cases that anchor the model in the two regimes it spans and show
the one choice new users most often get wrong (the ionic pair potential).  Run
them with the FAC-env python that has `pfac` (see `../STUDIO_SETUP.md`); the
diagnostic in #3 needs no FAC pipeline and runs in seconds.

## 1. Cold-matter anchor — proton in solid Al vs NIST PSTAR

```bash
bash examples/cold_proton_Al.sh          # electronic + Gordon-Kim nuclear
```

Gordon-Kim (`--npot=gk`) is the *correct* potential here: cold matter has bound
electrons.  **Expected:** the range column reproduces the PSTAR CSDA range to a
few percent (the paper validates C/Al/Ag/Au in Fig. 18) — e.g. a 1 MeV proton in
Al gives **~3.9 mg/cm²**, matching PSTAR.

## 2. Plasma anchor — 3.5 MeV alpha in a DT hot spot (done right)

```bash
bash examples/alpha_DT_plasma.sh
```

A fully-ionized plasma: use **`--npot=ionsphere`** (bare ions, correct 1/v²
limit) and the real **`--imass=2.014,3.016 --iwt=1,1`** D,T ion masses.
**Expected:** ion share ~6% at the 3.5 MeV birth energy; multiply the range
column by **2.515/1.008** for the real DT areal density → **ρR_α ≈ 0.51 g/cm²**,
consistent with BPS (0.45) and the Zylstra MD fit (0.46).  The script also shows
the *wrong* way (`--npot=gk` on the plasma), which trips the runtime guard.

## 3. Potential-choice diagnostic — the Rutherford-limit check

```bash
python examples/check_rutherford_limit.py --zp 2 --zt 1 --d 40.3 --t 10000
```

The fastest sanity check on your ionic potential: for a fast projectile the ionic
term must fall as 1/v², so `dEdx_n × Etot` should be ~flat (up to the slow Coulomb
log).  **Expected:** ion-sphere gives a ~1.7× climb over 0.05→100 MeV/u (just the
log) → "Rutherford recovered."  Gordon-Kim on the same fully-ionized plasma would
climb ~30× — the neutral-atom artifact the guard warns about.

## Quick V&V checklist

See [`../docs/VV_checklist.md`](../docs/VV_checklist.md) for the full list of
diagnostics and gotchas (potential-by-regime, the H-surrogate ×2.515 conversion,
`--imass`, the `dEdx_n×Etot` flatness test, and what *not* to trust).

## Strong-coupling (future)

`../validation/grabowski_compare.py` is a framework to compare the strong-collision
(Zwicknagel) correction against Grabowski (2013) classical MD; it needs digitized
MD points to complete (see its docstring).  Held for future validation.
