# Running the corrected dedx-erpa (with FAC/pfac) — setup for Claude on the Studio

Three things the Studio Claude needs: (1) the corrected code, (2) a working FAC/pfac
Python it can actually call, (3) the run command with the nuclear/ionic term.

## 1. Get the corrected code onto the Studio

The Bloch fix, the physical-range/unit-label output, `nuclear.py`, and the new
`degeneracy_panels.py` live in this repo. From the machine that has the fixes:

```bash
cd ~/Codes/dedx-erpa
git add -A
git commit -m "Bloch sub-thermal switch; physical range + unit labels; degeneracy generator; package scaffold"
git push origin master        # remote 'origin' = github.com/dedx-erpa/dedx  (or: git push fork master)
```

On the **Studio**:

```bash
cd ~/Codes/dedx-erpa           # or wherever the Studio's checkout is; else: git clone <the repo>
git pull
make                           # RECOMPILE dedx.f -- the Fortran binary is machine-specific;
                               # the one committed for another Mac will not run here
```

## 2. Give Claude access to pfac

Claude's Bash tool runs in a **non-interactive shell where `conda activate` does not
work** (conda isn't initialized). The fix that worked all session: **call the fac
env's Python by its full path.** Find it:

```bash
ls ~/miniforge3/envs/fac/bin/python 2>/dev/null \
  || ls ~/miniconda3/envs/fac/bin/python 2>/dev/null \
  || ls ~/anaconda3/envs/fac/bin/python 2>/dev/null \
  || conda info --envs
```

Verify pfac imports with THAT interpreter (not the default `python3`):

```bash
~/miniforge3/envs/fac/bin/python -c "import pfac.fac; print('pfac OK')"
```

Then **tell Claude explicitly**: *"Use `~/miniforge3/envs/fac/bin/python` (the fac
conda env) to run dedx.py — the default python does not have pfac."*

If `import pfac.fac` fails, FAC/pfac isn't installed in that env. Install FAC from the
**master** branch of https://github.com/flexible-atomic-code/fac and its `pfac`
Python interface into the env (earlier FAC releases are not compatible).

## 3. Run with the nuclear/ionic term

```bash
~/miniforge3/envs/fac/bin/python dedx.py \
    --zp=2 --zt=1 --d=171.6 --t=10000 --ti=10000 \
    --aa=2 --nuc=1 --npot=gk --od=AlphaDT_hotspot
```

| flag | meaning |
|---|---|
| `--zp=2` | projectile nuclear charge (2 = alpha) |
| `--zt=1` | target Z (1 = hydrogen surrogate for DT at the DT electron density) |
| `--d`, `--t` | target mass density (g/cc) and electron temperature (eV) |
| `--ti` | ion temperature (eV); defaults to `--t` if omitted |
| `--aa=2` | run the average-atom model to (re)generate the density; 1 reuse, 0 read |
| `--nuc=1` | compute nuclear/ionic stopping -> writes `dedx_nuc.dat` |
| `--npot=gk` | Gordon-Kim potential (the validated default) |
| `--od` | output directory |

Outputs in `--od`:
- `dedx.dat` -- electronic only: `E/A[MeV/u]  dEdx[1e-15 eVcm2/atom]  Rphys[mg/cm2]`
- `dedx_nuc.dat` -- with the ion channel: columns
  `E/A(MeV)  dEdx_e  dEdx_n[gk]  dEdx_tot  range  proj_range  detour`

## Important notes (post-correction)

- **Range is now PHYSICAL (per-ion), in mg/cm2.** The old per-nucleon convention is
  gone -- do NOT multiply the alpha range by 4 anymore. (A 3.5 MeV alpha in cold Al now
  reads ~3.4 mg/cm2 = 12.7 um directly.)
- **H-surrogate mass conversion.** The range column is in the *target* (H-surrogate)
  mass. For a real DT areal density multiply by (DT amu/electron)/(H amu/electron) =
  2.515/1.008 ~= 2.5. (`--d=171.6` H matches DT 429 g/cc; `--d=40.3` matches DT 100 g/cc.)
- **Use the total, not electronic-only, for the alpha range in a hot plasma.** At
  ~10 keV the alpha is sub-thermal to the electrons but fast relative to the ions, so
  `dEdx_tot` (electronic + ionic) is what you integrate. Result at DT 100 g/cc, 10 keV:
  rho_R ~ 0.29 g/cm2 -- the classic Fraley/Atzeni/Lindl value.
- **Do not rely on `alpha_dt_verify.py` for quantitative numbers.** It carries a
  spurious factor of 0.5 in `S_field` (halves the stopping -> ranges ~2x too long). It
  is a rough analytic cross-check only; use `dedx.py --nuc` for the real answer.
- **Faussurier ionic caveat.** The finite-Ti ionic term over-predicts the ionic
  stopping at high projectile velocity (converges to BPS/LPZ near the ion thermal
  velocity); it does not bias the low-velocity-weighted integrated range.
