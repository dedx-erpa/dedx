# Effective projectile charge (Z_eff) — design and build plan

**Status:** design / not started (Sep 2026). Startup description for a dedicated
work thread. Goal decided with T. Mehlhorn: build a *unified* Z_eff spanning
cold/WDM to hot plasma, consistent across channels. (Scope = items 1+2+4 from the
scoping discussion; item 3, validate-first, is woven through, not separate.)

## Goal
One continuous effective-charge model Z_eff(v, Z, target-state) that sets the
electronic stopping (∝ (Z_eff/v)²) and is consistent with the ionic/nuclear
channel, from cold matter through warm dense matter to fully-ionized plasma.

## Aperture: the model is NOT decided yet
Brandt–Kitagawa is the leading *cold/WDM* candidate, but keep the choice open —
the plasma regime in particular has more appropriate theory to evaluate first.
Architect Z_eff as **pluggable** (a `--zeff=fixed|ziegler|bk|peter|barriga|auto`
switch) so several models coexist and are compared against data before any becomes
default. Candidate models to survey (Phase 0):
- **Brandt–Kitagawa (1982)** — statistical (Thomas–Fermi) projectile electron
  cloud; cold/WDM; underlies SRIM/Ziegler. Leading cold candidate.
- **Peter & Meyer-ter-Vehn (PRA 43, 1991)** — the canonical charge-state and
  effective-charge treatment for **heavy ions in dense plasma** (ionization balance
  + enhanced stopping). Strong candidate for the plasma branch.
- **Barriga-Carrasco (and co-workers)** — dielectric/RPA-based effective charge and
  charge-state evolution in plasmas; fits our RPA framework. Candidate for the
  plasma branch. (Already cited in the paper's bib.)
- **Ziegler/Betz fit** — the current legacy model; keep as `--zeff=ziegler` for
  back-compatibility and as a cold-matter reference.
Mehlhorn has Barriga-Carrasco, Peter, et al. papers to fold in — review these
before committing to a functional form.

## Current code state (what exists today)
- **Electronic (dedx.f, ~lines 391–449):** velocity-dependent Ziegler/Betz *fit*
  `Z_eff = Z(1 − 1.034·exp[−(137.04·β)/Z^0.69])` = `Z(1 − 1.034·exp[−(v/v₀)/Z^0.69])`,
  applied only for `zzp > 2` (H, He treated as fully stripped). Feeds the stopping
  as `dedx ∝ (zeff/vion)²` and into the loss function `vlhfit(...,zeff,...)`.
  Validated implicitly at Cayzac (carbon, Z=6) but only at high v where Z_eff ≈ Z.
- **Ionic / nuclear (nuclear.py):** `z0p = zp` ("fully stripped"); the pair
  potential screens with `qchar = z0p·zbar`. Full nuclear Z at close approach is
  right for the nucleus–nucleus collision, BUT the projectile electron cloud in the
  Gordon–Kim potential (nuclear_gk.py: D0, rho_B0 built from `z0p`) should carry the
  (Z − Z_eff) *bound* electrons of a partially-stripped projectile — CHECK this.
- **No plasma projectile-ionization model.** The Ziegler form is a cold-matter fit;
  nothing computes the projectile charge state from the plasma T_e/n_e.

## Physics, by regime (the three pieces to unify)
1. **(item 2) Cold / WDM — Brandt–Kitagawa (BK, 1982).** Replace the Ziegler fit
   with the physics-based BK: Thomas–Fermi statistical model of the projectile's own
   electron cloud → fractional ionization `q = 1 − exp(−v_r/(v₀ Z^⅔))` PLUS the
   screening-length / "size effect" term (bound-cloud radius vs target screening).
   Covers Z ≤ 2 naturally; underlies SRIM/Ziegler but from first principles.
2. **(item 1) Hot plasma — ionization balance.** In a plasma the charge state is set
   by collisional-radiative equilibrium with the plasma electrons (T_e, n_e
   dependent), not the velocity-only cold criterion. Plasmas strip more effectively
   than cold matter at equal velocity — the "enhanced stopping" of heavy ions in
   plasma (relevant to the GSI/Volpe proposal). FAC already does the target
   average-atom; the same machinery gives the *projectile* ion's Z_eff(T_e, n_e).
   Dynamic (velocity) stripping may also matter for a fast projectile.
3. **(item 4) Consistency + diagnostics.** Electronic uses Z_eff (net charge).
   Nuclear uses full Z at close approach, with the projectile *bound* cloud
   (Z − Z_eff electrons) in the pair potential (reconcile the GK cloud). Expose
   Z_eff(E) as an output column.

## Proposed architecture
A single `zeff.py`: `zeff(v, Z, target_state)` returning BK for a bound/neutral
target and the plasma-CR value for an ionized one, **blended by the target's own
ionization** (`zbar` from the AA marks where you are on the cold↔plasma axis).
Wire it into BOTH channels; replace the inline Ziegler in dedx.f with a call to it
(or a shared table). Keep a `--zeff=fixed|ziegler|bk|peter|barriga|auto` switch so
the old behavior is reproducible and each candidate model is opt-in and comparable
until one earns default status against data.

## Build order (foundation-first)
- **0. Model survey + pluggable scaffold.** Read the Barriga-Carrasco / Peter &
  Meyer-ter-Vehn / BK papers; decide per-regime functional forms; stand up the
  `zeff.py` module + `--zeff` switch so models are swappable and comparable. Keep the
  aperture open — do NOT hardwire BK.
- **A. Cold/WDM Z_eff (leading candidate BK).** Self-contained upgrade of Ziegler;
  validate against SRIM/Ziegler heavy-ion ranges. Everything downstream blends into
  this.
- **B. Plasma-CR projectile Z_eff (FAC).** Enhanced-stripping regime; validate vs
  heavy-ion-in-plasma data (GSI/Hoffmann-type; ties to the GSI proposal).
- **C. Channel consistency + Z_eff diagnostic output.**

## First target (to confirm at kickoff)
Proposed: **carbon (Z=6) in cold Al** to anchor A (trusted target, C ties to the
Cayzac plasma end later, dense SRIM/PSTAR data), then a heavier ion (Cu or Au) to
stress the low-velocity regime where BK and Ziegler diverge most. Alternative: start
on B (plasma, the novel part for the paper) — deferred to kickoff.

## Validation data
- Cold: SRIM/Ziegler heavy-ion ranges (C, Cu, Au in Al and others); PSTAR for light.
- Plasma: Cayzac 2017 (C in DT, high-v); GSI/Hoffmann heavy-ion-in-plasma
  "enhanced stopping"; see the GSI/Volpe proposal materials.

## Open questions
- GK projectile electron cloud: does it represent Z_eff-consistent bound electrons
  for a partially-stripped projectile? (nuclear_gk.py D0/rho_B0.)
- Cold→plasma blend: continuous in `zbar`, or a physical CR model that reduces to BK
  as T→0? Prefer the latter if tractable with FAC.
- Dynamic vs equilibrium charge state in a fast projectile through a plasma.
