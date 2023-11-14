#!/bin/bash
python dedx.py --zt=5 --aa=2 --mloss=24 --d=2.35 --t=0.025 --od=acB --floss=tst2.dat --zp=2 --emin=5e-3 --emax=5e2
python dedx.py --zt=5 --aa=2 --mloss=24 --d=2.35 --t=1e3 --od=ahB --zp=2 --emin=5e-3 --emax=5e2

python dedx.py --zc='5,7' --wc='1.,1.' --aa=2 --mloss=24 --d=2.35 --t=0.025 --od=acBN --floss=tst2.dat --zp=2 --emin=5e-3 --emax=5e2
python dedx.py --zc='5,7' --wc='1.,1.' --aa=2 --mloss=24 --d=2.35 --t=1e3 --od=ahBN --zp=2 --emin=5e-3 --emax=5e2

