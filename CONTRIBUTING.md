# Contributing to dedx-erpa

Thanks for your interest in dedx-erpa, an enhanced RPA/LDA model for ion stopping
power from cold matter to high-energy-density plasmas. Contributions, bug
reports, and feature requests are welcome.

## Repository layout

dedx-erpa is part research code, part installable package:

- **Compute pipeline** (`dedx.py`, `dief.py`, `nuclear.py`, `nuclear_gk.py`,
  `rd.py`, `dedx.f`, the `*.py` validation drivers) generates stopping powers.
  It requires `pfac` (the Python interface to [FAC](https://github.com/flexible-atomic-code/fac),
  master branch) and a `gfortran` build of `dedx.f` via `make`.
- **`dedx_erpa/`** is a small, dependency-light package that exposes the results
  through [`astropy.units`](https://docs.astropy.org/en/stable/units/) quantities
  for interoperability with PlasmaPy and the scientific-Python ecosystem. It only
  needs `numpy` and `astropy` and does **not** require FAC.

## Development setup

```bash
python -m pip install -e ".[test]"   # package + pytest
pytest -q                            # runs the FAC-free IO/units tests
```

To work on the compute pipeline as well, install FAC/`pfac` and build the kernel:

```bash
pip install -r requirements.txt      # numpy, scipy, matplotlib
make                                 # builds the dedx Fortran executable
```

## Tests

- The `tests/` suite exercises the units/IO interoperability layer against the
  committed sample result directories (e.g. `ColdAl/`) and runs in CI without
  FAC.
- Pipeline-level validation (cold-NIST nuclear, WDM sweeps, alpha-in-DT) lives in
  the `nuc_test.py`, `nuc_fac_compare.py`, `alpha_dt_verify.py`, and related
  drivers; run these locally where FAC is available.
- Please add a test when you add or change behavior in `dedx_erpa/`.

## Style

- Python 3.10+, [PEP 8](https://peps.python.org/pep-0008/), max line length 88.
- `ruff check dedx_erpa tests` should pass (CI enforces it).

## PlasmaPy interoperability

dedx-erpa aims to be a PlasmaPy **affiliated package**. New public API should
accept and return `astropy.units` quantities and, where natural, interoperate
with `plasmapy.particles`. Please keep that convention when extending
`dedx_erpa/`.

## License

By contributing, you agree that your contributions are licensed under the
project's GNU General Public License v3 or later (see `LICENSE`).
