#!/bin/bash

#boron dedx at T=2keV, various densities
python bd.py

#boron dedx at 0.234 and 2.34 g/cc, various temp.
python bdt.py

#carbon dedx, at rho=0.05/0.1/0.5 g/cc, te=0.025/10./20./30. eV
python bdt_C.py

#Al dedx, rho=0.027/0.27/2.7 g/cc, various temp. map out dependence on zbar
python bdt_Al.py

#NI dedx, rho=0.0891/0.891/8.91 g/cc, various temp. map out dependence on zbar
python bdt_Ni.py

#Au dedx, rho=0.193/1.93/19.3 g/cc, various temp. map out dependence on zbar
python bdt_Au.py

#Al dedx, rho=2.7/0.01728/0.00864 g/cc, te=0.025/4.5/15.0 eV.
python adt_Al.py

#Mylar dedx, rho=1.35/0.00864/0.00432 g/cc, te=0.025/3.0/10.0 eV
python adt_mylar.py

#Al dedx, rho=2.7/0.027/0.0027 g/cc, te=0.025/10/20/30/40/48 eV
python cdt_Al.py

#Ni dedx, rho=8.91/0.0891/0.00891 g/cc, te=0.025/10/20/30/42 eV
python cdt_Ni.py

#dedx for single electron temp&dens plasma, 
python ta.py

#to produce radial distribution of the loss function,
#dedx.py needs to be run with option --mout=1.
#this is for cold Al.
python dedx.py --zt=13 --d=2.7 --t=0.025 --floss=tct.dat --od=data/Al --mout=1 --aa=2
