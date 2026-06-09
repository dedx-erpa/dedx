# Band-gap (Levine–Louie) correction for insulators — implementation plan

## Context / why
The eRPA-LDA model overshoots PSTAR/experiment through the stopping peak for ionic
insulators (LiF, SiO2, H2O, Al2O3 — see Sec. IV / Fig. 17 of the draft). Cause: the
shipped `t##.dat` loss tables are generated with **no gap** (`wi = 0`), so the
average-atom model's valence/quasi-free electrons get free-electron-gas response with
no excitation threshold. A band gap should zero the loss below `E_g` and shift the
spectrum up (Levine–Louie). Verified in isolation: imposing a 14 eV gap on a uniform
gas (rs≈2) zeroes Im[-1/ε] below 14 eV, rigidly shifts the plasmon peak up, and
suppresses the stopping number at low projectile velocity while converging to the
no-gap result by v≈2 a.u. (~95 keV protons) — i.e. exactly the region of the overshoot.

The gap machinery already exists in `dief.py` (`gapcorr`, the `wi`/`iu` threading,
`u = sqrt(u²−iu²)` shift) but has been dormant since the first commit (Oct 2022),
never wired into production. This plan finishes it.

It is **orthogonal** to the `vzt`/`lindt` local-effective-temperature path (`icb=2`,
`mloss ≥ 200`), which is a bound-electron *characteristic-velocity* correction, not a
gap. Do not modify that path.

## Key design insight
The gap lives in the **loss table**, not in `dedx.f`. `rloss` just reads tabulated loss
numbers and has no notion of a gap. So generating tables with `wi = E_g` bakes the gap
in and `dedx.f` reads them transparently. This enables a near-zero-Fortran Stage 1.

## Units
- `wi` in `dief.py` is a gap **frequency in atomic units (Hartree)**, since ω[a.u.]=E[Ha].
- `wi[Ha] = E_g[eV] / 27.2114`. LiF optical gap ≈ 14 eV → `wi ≈ 0.514`.

---

## Environment / building FAC (pfac) — DO THIS FIRST
The whole pipeline (`dedx.py` → `from pfac import aa`) needs FAC's `pfac` Python package
installed. Everything downstream — table generation, the `dedx` runs, and any Claude Code
session — must run inside the env where `pfac` is installed, or imports fail.

### The working build recipe (macOS, Apple Silicon, miniforge)
Modern Python (3.12+) removed stdlib `distutils`, and modern setuptools refactored the
vendored compiler classes, which breaks FAC's `setup.py` (`MyCompiler(...)` called with an
arg the new base class rejects). The fix is a dedicated Python 3.11 env using the genuine
stdlib `distutils`:

```
conda create -n fac python=3.11 numpy
conda activate fac
conda install -c conda-forge compilers          # clang/clang++/gfortran for this arch
export SETUPTOOLS_USE_DISTUTILS=stdlib           # use real stdlib distutils (3.11 only)
```

Apple Silicon also needs up-to-date autotools helper scripts — the bundled `config.sub`
doesn't recognize `arm64-apple-darwin`:

```
conda install -c conda-forge gnuconfig
cp "$CONDA_PREFIX/share/gnuconfig/config.sub"   .
cp "$CONDA_PREFIX/share/gnuconfig/config.guess" .
```

Then build IN ORDER (the order matters — `make pfac` generates headers like `init.h` that
`python/fac.c` needs; running `make install-pfac` alone fails with "init.h not found"):

```
cd ~/Codes/fac
./configure
make clean        # if reconfigured / switching arch
make              # core C library (libfac.a)
make pfac         # build the python extension (generates headers)
make install-pfac # install pfac into the env's site-packages
```

### Verify
```
python -c "from pfac import aa, consts, fac; print('pfac ok')"
```

### Notes
- `conda activate fac` is required before every pipeline / Claude Code session.
- Stage 2 edits `wden` in `python/aa.py` — pure Python, so no recompile; just
  `make install-pfac` again (or re-copy `aa.py` into site-packages).
- The `distutils is deprecated` warnings on 3.11 are harmless.
- Staying on Python 3.13 instead: pin `pip install "setuptools<75"` (no stdlib distutils
  fallback there, so it's more fragile — the 3.11 env is the reliable path).
- Upstream: FAC's `setup.py` `MyCompiler` call is a genuine modern-Python incompatibility
  worth reporting to Ming Feng Gu.

---

## Stage 0 — recon (Claude Code: read before editing)
Confirm these in the current repo; the plan assumes them but verify:
1. **Table-generation entry point.** Find how `t##.dat` / `tct.dat` are produced.
   Check `dief.py` `__main__` / `tabdedx` / `savdedx`, and the batch scripts
   (`b.py`, `ba.sh`, `b.sh`, `bd.py`). Identify exactly where `wi` is (or can be) passed
   into `dedxmp(...)` and how the temperature grid `ts = 10**linspace(-1,4,21)` maps to
   the `t00..t20`/`tct` filenames. (In the version reviewed, `dedxmp(es, ds, ts, wi=...)`
   threads `wi → dedx → dedx1k → rpa`, and `dief.py __main__` took `wi` as a CLI arg.)
2. **Loss-table format.** `rloss` reads `ttti, rhoji, ek, cl1..cl6` (3 index + 6 loss
   columns). Confirm `savdedx` writes 9 columns: `log10(T), log10(N), E, r0, r1, ra0,
   ra1, rb0, rxi`. Gapped tables must use the identical format so `rloss` is unchanged.
3. **`floss` selection.** In `dedx.py`, `floss`/`floss0` are chosen from `opts.t`
   against the `ts` grid and written into `dedx.inp` (`odir.inp`). Stage 1/2 will point
   these at gapped table files.
4. **PSTAR overlay.** Find the script that plots eRPA vs PSTAR vs experiment for a
   target (used for Fig. 17). Likely uses `rd.py` (`rd.rdedx`) + the figure scripts.
   We reuse it to compare gapless vs gapped vs PSTAR.

---

## Stage 1 — minimal validation (gapped tables, binding correction OFF)
Goal: see whether the gap pulls the LiF curve onto PSTAR, with the least code change.

1. **Generate gapped cold tables.** Re-run the table generator with `wi = 14/27.2114`
   for the cold temperature(s) needed by LiF at T=0.025 eV (the `tct.dat` / `t00.dat`
   range). Name them distinctly, e.g. `tct_g14.dat`, `t00_g14.dat`. Everything else
   (density grid, energy grid, columns) identical to the standard tables.
2. **Run LiF with the gapped tables and binding correction off.** LiF AA densities are
   per-element (`data/LiF/Li`, `data/LiF/F`). For each element's `dedx` run, set
   `floss`/`floss0` to the gapped tables and use `mloss` with third digit = 1 so
   `icb = 0` (no `rhoe` bump — avoids double-counting gap + soft correction).
   - Practically: either pass the gapped filenames via the `--floss/--floss0` options,
     or temporarily swap the table the driver selects.
3. **Compare.** Overlay three curves vs energy for LiF:
   - current eRPA (gapless table + `icb=1` `rhoe` bump),
   - new gapped (gapped table, `icb=0`),
   - PSTAR / experiment.
   Success = the gapped curve drops from the overshoot toward PSTAR through the peak
   (~3e-3 to ~1e-1 MeV). Expect possible **over-suppression at the very lowest energies**
   (the demo gave ratio ~0.015 at 2 keV for a uniform full gap) — that's expected and
   motivates Stage 2; nuclear stopping (not yet in the model) also fills in there.

Note: applying the gap to the whole local density is OK for LiF (≈no free electrons) and
self-limiting for cores (local ω_p ≫ 14 eV; cores matter only in the Bethe regime where
the gap washes out). For F the gap is physical (2p valence below the gap); for Li⁺ it is
nearly inert (1s² core only) — the uniform treatment barely touches Li, fine for Stage 1.

---

## Stage 2 — gap tied to the AA binding frequency ω_b(r) (ionization-aware)
Goal: replace the soft `rhoe` density-bump with a Levine–Louie gap whose threshold is a
smooth *local* property derived from the average-atom binding frequency, so it (a) cures
the cold-insulator overshoot, (b) tracks ionization automatically, and (c) vanishes in
the fully-ionized/plasma limit — all within the single unified total-density convolution
(no hard bound/free split, no discontinuity).

**Why ω_b(r), not a fixed E_g.** A constant cold band gap is wrong once the material
ionizes: the surviving bound electrons belong to higher charge states with *larger*
excitation thresholds, and the band-gap concept itself dissolves into a plasma. The AA
binding frequency ω_b(r) (draft Eq. 7, ln ω_b = Σ f_i ln I_i, with I_i, f_i from FAC)
already evolves correctly and continuously with T and ρ, so the gap inherits the right
behavior for free and disappears smoothly as the bound population ionizes.

### Sourcing ω_b(r) — RESOLVED: it's already in the AA output
The `.den` column legend settles this. `dedx.f`'s `rhoe` = `.den` c10 ("avg binding
energy density"); the AA output also has c06 ("bound density"). The local average
binding energy is their ratio:

    ω_b(r) = c10(r) / c06(r)      (energy density ÷ density = energy)

No new FAC physics — FAC already computes the binding energy (also the scalar header
`ub`). Cleanest plumbing: add ~2 lines to `wden` (python/aa.py) to emit `d[10]/d[6]`
(guard c06→0) as a **7th array**, and read it in `dedx.f` as `wb(k)`. (`wden` already
loads all `.den` columns into `d`, so this is local and tiny.)

Robustness: κ is a free scale constant, so only the radial *shape* of ω_b(r) must be
right — large where bound, → 0 where free, tracking ionization via the AA occupations.
c10/c06 has that shape regardless of the exact internal weighting or a constant
normalization (which κ absorbs in the cold-PSTAR fit). Do NOT use the earlier
`sqrt(4π·rhoe/fpr)` form — c10 is an energy density, so divide, don't sqrt.

Units check (do at runtime, not on paper): `efm`≈11.45 reads as eV (Al EF), but `dedx.f`
applies `*27.21` to c15, so columns may mix eV/Hartree. Verify the recovered ω_b(r) comes
out as sensible binding energies (tens of eV) before trusting the absolute scale; κ and
the unit factor are degenerate, so fix units first, then calibrate κ.

Note: `vzt` = c15 "effective temperature" — this is the `lindt`/icb=2 local-temperature
mechanism, separate from the gap. Leave it untouched.

### Implementation
1. `dedx.py`: add `--kgap` (the ω_b→gap scale κ; default 0 = off) and, for Stage-1-style
   runs, keep a fixed-gap override `--egap` (eV). Write both into `dedxinp`.
2. `dedx.f`:
   - Add `kgap`/`egap` to the namelist and `xkgap`/`xegap` commons (defaults 0).
   - New mode parallel to the existing `icb` structure: `mloss ≥ 300 → icb = 3`
     (test ≥300 *before* ≥200; then `mloss = mod(mloss,100)`). In `icb=3`, **disable the
     `rhoe` density-bump** (no double-counting) and instead set a per-radius gap
     `E_g(r) = xkgap · ω_b(r)` (or `xegap` if the fixed override is used).
   - Evaluate the local loss from a **gapped** table at `E_g(r)`. Because `E_g(r)` varies
     continuously, this replaces the single `vlhtab` lookup with a 3-way interpolation in
     `(E, ρ, E_g)`.
3. `dief.py`: generate a gapped-table **family** — the standard tables plus copies at a
   grid of `wi` (e.g. {0, 0.1, 0.2, 0.4, 0.8} Ha) for the needed temperatures, identical
   9-column format so `rloss` only needs a gap index added.

### Calibration
κ is the single knob, analogous to γ². ω_b (orbital binding ≈ photoionization threshold)
exceeds the optical band gap, so expect κ < 1. Fix κ by matching cold-insulator PSTAR;
the target is that one κ lands the whole cold set. Ionization/WDM tracking then follows
automatically from ω_b(T,ρ) with no extra parameter — which is the unified-model payoff.

---

## Stage 3 — calibrate and roll out
1. Sweep `kgap` (κ) for LiF against PSTAR (CLI-tunable, no recompile).
2. Test whether the **same** κ lands the other insulators. Reference optical gaps for
   sanity: LiF (~14 eV), SiO2 (~9 eV), Al2O3 (~8.8 eV), H2O (~7 eV), MgF2 (~12 eV),
   CaF2 (~12 eV), NaCl (~8.5 eV). Since the gap is `κ·ω_b(r)`, check that the implied
   effective gaps order sensibly across materials.
3. If one κ works across the cold set, that is the paper result: "sub-gap electronic loss
   is suppressed via a Levine–Louie gapped dielectric with the threshold set by the
   average-atom binding frequency" — f-sum-rule preserving, single calibrated constant.
4. WDM bonus (no extra work): rerun a heated case (e.g. LiF or carbon at finite T) and
   confirm the gap effect fades smoothly as ω_b(T,ρ) evolves and the bound population
   ionizes — the unified cold→WDM→plasma behavior the model is supposed to have.

---

## Risks / open questions
- ω_b(r) source is RESOLVED: c10/c06 from the `.den` output, emitted via a ~2-line
  `wden` edit. Remaining runtime check is units (eV vs Hartree), not availability.
- Low-energy behavior: tying the gap to ω_b(r) and the bound population should self-limit
  the over-suppression a *uniform* gap showed at the lowest energies; any residual
  shortfall at very low E is expected to be filled by nuclear stopping (future module).
- Confirm the gapped Lindhard still satisfies the f-sum rule numerically (Levine–Louie
  preserves it by construction; check the gapless→gapped oscillator-strength balance).
- Make sure gap mode and the `icb=1` `rhoe` correction are mutually exclusive (never both).
- Leave the `vzt`/`lindt` (`icb=2`, c15 "effective temperature") path untouched.

## Paper reconciliation (do in Chat, not code)
- Sec. II.4 currently describes ω̃ₚ² = ωₚ² + γ²ωₚ² with γ²=0.75; the code now uses the
  energy-dependent `xepa`/`xepb` form (defaults 1.0 / 1.1). Rewrite to match the code, or
  decide which is canonical.
- Add a subsection for the band-gap correction (Levine–Louie), with the LiF and
  multi-insulator validation figures.
