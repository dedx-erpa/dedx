FC=gfortran 
OPTS=-fno-align-commons

dedx: dedx.f
	$(FC) $(OPTS) -o dedx dedx.f

