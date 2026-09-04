# dedx-erpa — project guide for Claude

Wide-range ion **stopping power**: electronic (RPA + LDA via a FAC average-atom
density) plus a nuclear/ionic (elastic ion–ion) term, continuous from cold solids
through warm dense matter to fusion plasmas. Companion to the Physics of Plasmas
paper (POP26-AR-01190).

## Running it
- **Use the conda env that has `pfac`**, calling its Python by full path — `conda
  activate` does not work in non-interactive shells (e.g.
  `<conda-base>/envs/<env>/bin/python dedx.py ...`; `STUDIO_SETUP.md` shows how to
  locate it).
- **Rebuild the Fortran binary** after editing `dedx.f`: `make` (the committed binary
  is machine-specific and will not run elsewhere).
- Full setup notes: `STUDIO_SETUP.md`. Runnable V&V: `examples/`. Gotchas:
  `docs/VV_checklist.md`.

## Key physics rule — ionic pair potential by regime (don't get this wrong)
- `--npot=ionsphere` for a **fully-ionized plasma** (bare ions; correct 1/v² limit).
- `--npot=gk` (Gordon–Kim) for **cold / warm dense matter** (bound electrons).
- Using `gk` on an ionized plasma over-predicts the ionic term ~30–75× at high v;
  a runtime guard warns, and the total column auto-selects `ionsphere` for ionized
  targets. For a DT H-surrogate, multiply the range column by **2.515/1.008** for the
  real DT areal density, and pass real isotope masses with `--imass=2.014,3.016
  --iwt=1,1`.

## Current focus
**Z_eff — effective projectile charge overhaul.** See `docs/zeff_design.md`. Aim: one
continuous, *pluggable* Z_eff from cold/WDM to plasma, consistent across channels.
Model choice is **open** — candidates to survey before hardwiring: Brandt–Kitagawa
(cold), Peter & Meyer-ter-Vehn and Barriga-Carrasco (plasma), Ziegler (legacy).
Today: electronic uses a Ziegler fit (Z>2); ionic uses full Z.

## Before the next paper
**Complete the strong-collision (Zwicknagel) V&V.** `validation/grabowski_compare.py`
is a framework to compare the correction against Grabowski (2013) classical MD but
still needs digitized MD points (see its docstring). The current PoP paper only
*scopes* the correction (inactive at Γ≪1); any new paper that leans on the
strong-coupling regime should show it validated against MD first.

## Repo
CI (`.github/workflows/ci.yml`) runs pytest + `ruff check dedx_erpa tests` — keep
both green. Maintainer remotes and push workflow: see `STUDIO_SETUP.md`.
