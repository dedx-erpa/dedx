#!/bin/bash
python bd.py
python bdt.py
python bdt_C.py
python bdt_Al.py
python bdt_Ni.py
python bdt_Au.py
python adt_Al.py
python adt_mylar.py
python cdt_Al.py
python cdt_Ni.py
python ta.py
#to produce radial distribution of the loss function,
#dedx.py needs to be run with option --mout=1.
#this is for cold Al.
python dedx.py --zt=13 --d=2.7 --t=0.025 --floss=tct.dat --od=data/Al --mout=1 --aa=2
