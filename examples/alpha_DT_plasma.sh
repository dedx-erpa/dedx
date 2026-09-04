#!/usr/bin/env bash
# V&V example 2: plasma anchor -- 3.5 MeV alpha in a DT hot spot, done the RIGHT way.
# Z=1 hydrogen surrogate at the DT electron density (40.3 g/cc H = DT 100 g/cc),
# Te=Ti=10 keV, ion-sphere potential (fully ionized) + real D,T ion masses.
# Expected: ion share ~6% at birth; multiply the range column by 2.515/1.008 for the
#   real DT areal density -> rho_R ~ 0.51 g/cm^2 (consistent with BPS 0.45, Zylstra-MD 0.46).
PY=${PY:-$HOME/miniforge3/envs/fac/bin/python}
$PY dedx.py --zp=2 --zt=1 --d=40.3 --t=10000 --ti=10000 --aa=2 --nuc=1 \
    --npot=ionsphere --imass=2.014,3.016 --iwt=1,1 --od=AlphaDT_plasma
echo
echo "WRONG way (for contrast) -- gk on a fully-ionized plasma triggers the guard:"
$PY dedx.py --zp=2 --zt=1 --d=40.3 --t=10000 --ti=10000 --aa=0 --nuc=1 \
    --npot=gk --od=AlphaDT_plasma_gk 2>&1 | grep -i warn || true
