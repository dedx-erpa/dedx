dedx-erpa  
A package for ion beam stopping power calculations for plasma targets.

Contents:  
1. dief.py, for computing RPA dielectric functions, and stopping powers of  
   ions in a uniform electron gas. Used to tabulate the stopping numbers on  
   a grid of temperature and density grid.  
2. dedx.f, a fortran program to convolve the eRPA stopping powers and the  
   electron density distribution to generate stopping powers of any atomic  
   species.  
3. dedx.py, a driver to compute the electron density distribution in the  
   average atom model, and invoke dedx.f to generate the stopping powers.  
4. t##.dat, a set of tables for proton in uniform electron gas stopping powers  
   calculated with dief.py. used by dedx.py and dedx.f to compute stopping  
   powers for arbitrary electron density distributions.  
5. various utility and example scripts.
6. data/, proton in cold targets for Z=1-92 and select compounds. data  
   for each material is in the sub-directory named after its chemical symbol.  
   dedx.dat contains the dedx and range. dedx.pdf is a plot of the dedx vs E,  
   and range.pdf is a plot of range vs E.  

Instructions for running dedx.py:  

1. Download and install FAC from https://github.com/flexible-atomic-code,  
   which is needed for computing electron density distributions with  
   average atom models   

2. Modify Makefile and compile dedx.f using make  

3. python dedx.py --zp= --zt= ... some necessary options described below:  
--zp= projectile z, default 1  
--zt= target z  
--zc= for compound targets, a comma separated z for individual components.  
--wc= for compund targets, a comma separated weights of components by number.  
--fc= if zc & wc are not given, fc is the chemical formula of the compound.
      e.g., Al2O3 for aluminum oxide.  
--d=  target density in g/cc  
--t=  target temperature in eV
--taa= min temperature used for running average atom model. default 0.5 eV.  
       run aa model with very low temperatures can cause convergence problems.  
       electron density distribution of 0.5 eV aa model is practically the same  
       as the room temperature case of 0.025 eV.
--od= output directory  
--emin= minimum projectile energy in MeV, default 1e-3  
--emax= maximum projectile energy in MeV, default 100.0  
--mep= number of projectile energy points, default 100  
--frho= the file path for the density distribution function to be used in  
        dedx.f. normally, the density is to be computed with average atom  
	model. so no frho needs to be given. but if a density file already  
	exists, it can be used by specifying --aa=0  
--aa= run average atom mode.  
       2, generate electron density distribution by running aa model.  
       1, aa has been run before, just prepare the density distribution using  
          the data from the previous aa run.  
       0, the density distribution file is already present in the output dir.  
--mloss= the mode for computing the stopping power.  
       0, use the fitting formula from the RPA model of Wang et al.  
          PoP, vol. 5, no. 8, pp. 2977, 1998  
       1, use RPA stopping powers without corrections.  
       2, with local field correction (LCF)  
       3, without LFC, but with strong binary collision correction  
       4, with LFC and strong binary collision corrections.  
       11/12/13/14, same as 1/2/3/4, with the addition of Barkas term.  
       21/22/23/24, same as 1/2/3/4, with the addition of Barkas and
       Bloch terms.  
       By default, the bound electron correction term is included.  
       if mloss has a third digit of 1, the bound electron correction
       is omitted  
       The most sensible mode for common calculations would be mloss=24, which  
       is the default if mloss is not specified.

Examples:  
1. Proton in aluminum, solid density at room temperature.  
python dedx.py --zt=13 --aa=2 --d=2.7 --t=0.025 --od=ColdAl  

2. Proton in Mylar (C5H4O2), rho=1.35, te=10.0. average atom model for  
   compounds may take a while to compute.  
python dedx.py --zc='1,6,8' --wc='4,5,2' --aa=2 --d=1.35 --t=10.0 --od=MylarWDM  

3. This is equivalent to example 2.  
python dedx.py --fc=H4C5O2 --aa=2 --d=1.35 --t=10.0 --od=MylarWDM  

After running dedx.py, the output directory contains dedx.dat file.   
the headers starts with '#', and list:  
 nzt = number of constituent atoms  
  zt = target atom z array,  
  wt = target atom weight array,  
  zp = projectile z,  
  rs = unit cell radius,  
  te = electron temperature,  
 rho = material density,  
zbar = mean charge of plasma,  
 mep = number of energy points  
 
the data section has 3 columns,  
      Energy/AMU (MeV)  
      dEdx(10^-15 eV/cm2/atom)  
      Range(mg/cm2)  
