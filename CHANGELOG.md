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
