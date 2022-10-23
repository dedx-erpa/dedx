#!/bin/bash
python dedx.py --zt=5 --aa=2 --mloss=24 --d=4.0 --t=2000 --od=tmp/d0B
python dedx.py --zt=5 --aa=2 --mloss=24 --d=40.0 --t=2000 --od=tmp/d1B
python dedx.py --zt=5 --aa=2 --mloss=24 --d=400.0 --t=2000 --od=tmp/d2B
python dedx.py --zt=5 --aa=2 --mloss=24 --d=4000.0 --t=2000 --od=tmp/d3B
python dedx.py --zt=5 --aa=2 --mloss=24 --d=40000.0 --t=2000 --od=tmp/d4B
python dedx.py --zt=5 --aa=2 --mloss=24 --d=400000.0 --t=2000 --od=tmp/d5B
