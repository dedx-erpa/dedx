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
    --aa=2 --nuc=1 --npot=ionsphere --od=AlphaDT_hotspot
```

| flag | meaning |
|---|---|
| `--zp=2` | projectile nuclear charge (2 = alpha) |
| `--zt=1` | target Z (1 = hydrogen surrogate for DT at the DT electron density) |
| `--d`, `--t` | target mass density (g/cc) and electron temperature (eV) |
| `--ti` | ion temperature (eV); defaults to `--t` if omitted |
| `--aa=2` | run the average-atom model to (re)generate the density; 1 reuse, 0 read |
| `--nuc=1` | compute nuclear/ionic stopping -> writes `dedx_nuc.dat` |
| `--npot=ionsphere` | **plasma** ionic channel (screened Coulomb at the ion-sphere scale) |
| `--od` | output directory |

> **IMPORTANT — potential choice by regime.** `--npot=gk` (Gordon-Kim) is a
> *cold-matter / warm-dense-matter* potential: it models scattering off a
> **neutral atom** (bound electrons + exchange-correlation overlap). For a
> **fully-ionized hot plasma** (the alpha-hotspot case) it is the wrong physics
> and grossly over-predicts the ionic stopping at high projectile velocity — its
> forced-zero + xc attractive well near the cell radius `rs` produces a fixed
> large-impact-parameter scattering feature that does not shrink with velocity,
> so the ionic term fails to fall as the Rutherford 1/v^2 and comes out
> ~30-75x too large at the 3.5 MeV alpha birth energy (77% of the total stopping
> instead of ~15%). Use **`--npot=ionsphere`** for the ionized plasma; reserve
> `--npot=gk` for cold/WDM targets with bound electrons (where it is validated
> against NIST/ZBL). `--npot=yukawa` (Debye) over-screens at solid density
> because the 10 keV Debye length is tens of Bohr, far larger than `rs`.

Outputs in `--od`:
- `dedx.dat` -- electronic only:
  `E/A[MeV/u]  dEdx[1e-15eVcm2/atom]  Rphys[mg/cm2]  Etot[MeV]`
- `dedx_nuc.dat` -- with the ion channel:
  `E/A(MeV)  dEdx_e  dEdx_n[gk]  dEdx_tot  range  proj_range  detour  Etot[MeV]`

The last column, `Etot[MeV]`, is the total projectile kinetic energy (= E/A x the
projectile mass number), so the file gives R(E) directly in physical units (plot
the physical range vs `Etot`, both real). `E/A` is retained for the velocity-based
comparisons.

## Important notes (post-correction)

- **Range is now PHYSICAL (per-ion), in mg/cm2.** The old per-nucleon convention is
  gone -- do NOT multiply the alpha range by 4 anymore. (A 3.5 MeV alpha in cold Al now
  reads ~3.4 mg/cm2 = 12.7 um directly.)
- **H-surrogate mass conversion.** The range column is in the *target* (H-surrogate)
  mass. For a real DT areal density multiply by (DT amu/electron)/(H amu/electron) =
  2.515/1.008 ~= 2.5. (`--d=171.6` H matches DT 429 g/cc; `--d=40.3` matches DT 100 g/cc.)
- **Use the total, not electronic-only, for the alpha range in a hot plasma** --
  but with the **ion-sphere** ionic channel (see the potential-choice box above).
  At ~10 keV the alpha is sub-thermal to the electrons but fast relative to the ions,
  so `dEdx_tot` (electronic + ionic) is what you integrate. With `--npot=ionsphere`
  the 3.5 MeV alpha in DT 429 g/cc, 10 keV gives **rho_R ~ 0.5 g/cm2**, consistent
  with BPS (0.45) and the Zylstra MD fit (0.46) and somewhat above the classic
  Fraley/Atzeni value (~0.3). The earlier "rho_R ~ 0.29 = classic" number came from
  `--npot=gk` and is a **cold-matter-potential artifact** -- do not use it.
- **Do not rely on `alpha_dt_verify.py` for quantitative numbers.** It carries a
  spurious factor of 0.5 in `S_field` (halves the stopping -> ranges ~2x too long). It
  is a rough analytic cross-check only; use `dedx.py --nuc` for the real answer.
- **Gordon-Kim (`--npot=gk`) high-velocity artifact.** In a hot plasma the GK ionic
  stopping does not recover the 1/v^2 Rutherford limit (its S(Ec) falls as ~1/Ec
  instead of 1/Ec^2, because ~99% of the impact-parameter integral comes from the
  fixed attractive-well feature at `p ~ rs`, not from the velocity-shrinking
  close-collision core). Diagnostic: `dEdx_n[gk] * Etot` should be ~flat at high E
  (Rutherford) but instead climbs ~30x. Correct with `--npot=ionsphere`, whose
  `S(Ec)*Ec^2` is flat up to the physical Coulomb log.
