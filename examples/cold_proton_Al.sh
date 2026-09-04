#!/usr/bin/env bash
# V&V example 1: cold-matter anchor -- proton in solid Al vs NIST PSTAR.
# Electronic + Gordon-Kim nuclear (gk is correct here: cold matter, bound electrons).
# Expected: the range column reproduces PSTAR CSDA to a few percent (paper Fig. 18).
#   e.g. 1 MeV proton in Al -> ~3.9 mg/cm^2 (PSTAR CSDA ~3.9 mg/cm^2).
# Use the fac-env python (has pfac); adjust the path if yours differs.
PY=${PY:-$HOME/miniforge3/envs/fac/bin/python}
$PY dedx.py --zt=13 --aa=2 --d=2.7 --t=0.025 --nuc=1 --npot=gk \
    --emin=1e-3 --emax=10 --mep=60 --od=ColdAl_proton
echo "range column (col 5, mg/cm^2) is CSDA; compare to NIST PSTAR."
