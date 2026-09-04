# Changelog

All notable changes to dedx-erpa are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **Bloch (Z^4) correction no longer over-suppresses sub-thermal-projectile
  stopping in hot plasmas.** For a projectile slower than the thermal electrons
  (v < v_th), the RPA loss number is small, so the additive negative Bloch term
  was a large fraction of it (hitting the -50% cap) and roughly halved the
  stopping -- e.g. a 3.5 MeV alpha in 10 keV DT, giving a ~2.7x too-long range.
  Bloch is a Bethe-regime (fast-ion) correction, so it is now multiplied by a
  thermal switch fth = exp[-(0.4 vth^2/v^2)^2] (vth^2 = 2 te/Hartree) that is ~1
  for v > v_th and ->0 for v << v_th. Verified: cold C/Al/Ag/Au stopping is
  bit-for-bit unchanged; the fast-projectile Bloch effect (e.g. Cayzac, -7%) is
  preserved; the 10 keV alpha range drops 1.247 -> 0.769 g/cm2 (bare-RPA level).
- **Nuclear/ionic stopping in a fully-ionized plasma no longer over-predicts at
  high projectile velocity.** The Gordon-Kim neutral-atom potential (a cold/WDM
  construct) develops a spurious attractive well near the cell radius that breaks
  the Rutherford 1/v^2 limit for a bare-ion plasma, over-predicting the ionic term
  ~30-75x at the 3.5 MeV alpha birth energy. Fully-ionized-plasma runs now use the
  ion-sphere screened-Coulomb potential; the 'total'/range column selects it
  automatically for ionized targets, and a RuntimeWarning fires if `--npot=gk` is
  used on an ionized target. Corrected 3.5 MeV alpha range in DT 100 g/cc, 10 keV
  is consistent with BPS (0.45) and the Zylstra MD fit (0.46).
- **alpha_dt_verify.py: removed a spurious factor of 1/2** in the Chandrasekhar
  friction (`S_field`) that halved the stopping and doubled the range (the 3.5 MeV
  alpha in DT went ~0.89 -> ~0.44 g/cm2, now matching BPS).

### Added
- Output files gain a final **Etot[MeV]** column (total projectile kinetic energy
  = E/A x mass number), so dedx.dat / dedx_nuc.dat give R(E) directly in physical
  units. Appended at the end, so existing 3-column readers are unaffected.
- **`--imass` / `--iwt`**: supply the true ionic-channel isotope masses (e.g.
  D=2.014, T=3.016) for the elastic ion-ion term while the electronic AA, ne, and
  range normalization stay on the electron-density-matched surrogate. The ionic
  term ~ 1/m_target, so a Z=1 hydrogen surrogate over-weights a DT ion channel by
  ~2.2x; `--imass` corrects it (3.5 MeV alpha in DT 100 g/cc: rho_R 0.43 -> 0.51).
- **`--npot=ionsphere`** documented as the correct potential for fully-ionized
  plasmas, with the regime-aware total-column selection and the misuse guard above.
- **examples/**: runnable V&V cases with expected outputs -- cold proton in Al vs
  PSTAR, alpha in DT plasma (ion-sphere + isotope-mass workflow), and a
  potential-choice / Rutherford-limit diagnostic.
- **docs/VV_checklist.md**: a validation / gotchas checklist for new users.

### Changed
- **Range output is now the PHYSICAL (per-ion) range, not per-nucleon.** The
  internal range integral is over energy-per-nucleon (`E/A`), so the range was a
  per-nucleon quantity that needed a x(mass number) correction for multi-nucleon
  projectiles (x4 for an alpha, larger for heavy ions) -- a silent trap for
  everything except protons. `dedx.f` and `nuclear.py` now apply the projectile
  mass factor and write the physical range directly.
- **Output columns are explicitly unit-labeled**: `dedx.dat` header now reads
  `E/A[MeV/u]  dEdx[1e-15eVcm2/atom]  Rphys[mg/cm2]`; `dedx_nuc.dat` gains a units
  comment block. (Verified: 3.5 MeV alpha in cold Al -> 3.43 mg/cm2 = 12.7 um,
  matching the PSTAR-validated proton range at the same E/A to <2%.)

### Added
- `dedx_erpa` package: a dependency-light, `astropy.units`-aware reader for the
  pipeline's results (`read_stopping`, `read_total`, `StoppingResult`,
  `TotalStoppingResult`) for PlasmaPy / scientific-Python interoperability.
- `pyproject.toml` packaging (PEP 621), pip-installable as `dedx-erpa`.
- `tests/` pytest suite covering the IO parser and units API against the
  committed sample result directories.
- GitHub Actions CI (`.github/workflows/ci.yml`): tests on Python 3.10-3.13 plus
  a `ruff` lint job.
- Community docs: `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, this changelog.

### Notes
- The compute pipeline (FAC + Fortran `dedx`) is unchanged and still drives
  result generation via `dedx.py`; the new package wraps the outputs.
