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
   average atom models. The latest version in the master branch is required.  
   Previous releases of FAC is not compatible with dedx-erpa.  

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
--mout= if 1, also ouput radii-dependent contribution to dedx 
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
--v=  verbose level. if 0, the script runs silently. otherwise, progress is printed.

Examples:  
1. Proton in aluminum, solid density at room temperature.  
python dedx.py --zt=13 --aa=2 --d=2.7 --t=0.025 --od=ColdAl  

2. Proton in Mylar (C5H4O2), rho=1.35, te=10.0. average atom model for  
   compounds may take a while to compute.  
python dedx.py --zc='1,6,8' --wc='4,5,2' --aa=2 --d=1.35 --t=10.0 --od=MylarWDM  

3. This is equivalent to example 1.  
python dedx.py --fc=Al --aa=2 --d=2.7 --t=0.025 --od=ColdAl  

4. This is equivalent to example 2.  
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

Nuclear stopping power (cold matter -> WDM)  

dedx.f / dedx.py compute the *electronic* stopping power. The companion module
nuclear.py adds the *nuclear* (elastic ion-ion) stopping power using the
classical two-body model of Faussurier, Blancard & Gauthier, Phys. Plasmas 20,
012705 (2013). The two contributions are returned in the same units
(10^-15 eV cm^2/atom) and add to give the total stopping power, consistent from
cold solids through warm/hot dense matter. See nuclear_model_notes.md for the
equations and validation.  

Enable it by adding options to dedx.py:  
--nuc=  1 to also compute nuclear stopping (default 0).  
--npot= comma-separated pair potentials, subset of gk,ionsphere,yukawa  
        (default gk). gk (full Gordon-Kim) is the recommended model: it builds  
        the projectile-target interaction from the average-atom bound+free  
        electron density with the full electron-gas energy (electrostatic +  
        Thomas-Fermi kinetic + Dirac exchange + PW92 correlation, with the  
        Faussurier/Stein volume-conservation scaling), and reproduces the cold  
        NIST nuclear stopping. The analytic ionsphere/yukawa screen with free  
        electrons only and fall below NIST (Faussurier Fig. 1). The 'total'  
        column uses gk if present.  
--ti=   ion temperature in eV (default = Te). The finite-Ti Maxwellian average  
        (Eq. 6) makes the nuclear stopping negative below a threshold E* in hot  
        plasma (the projectile gains energy from the ion bath). Set --ti=0 for  
        the Ti=0 form (Eq. 10).  
--gkmuffin= GK target density beyond the Wigner-Seitz cell: 1 = muffin-tin  
        interstitial free-electron sea (consistent with the eRPA electronic  
        treatment; default), 0 = isolated neutral atom. The two differ by <0.5%  
        in cold matter and up to ~5% in hot dense matter (see notes).  
--plot= 1 to also save dedx_nuc.pdf: a log-log plot of the electronic,  
        nuclear/ionic (one curve per potential), and total stopping power, plus  
        the total range, with the run conditions in the title.  

This writes dedx_nuc.dat in the output directory, with columns  
      Energy/AMU (MeV)  
      dEdx_e   (electronic, 10^-15 eV cm^2/atom)  
      dEdx_n[<pot>]  (one column per requested potential)  
      dEdx_tot (electronic + chosen potential)  
      Range(mg/cm2) of the total  
and, with --plot=1, the figure dedx_nuc.pdf.  The same plot can be made from an  
existing run with  python -c "import nuclear; nuclear.plot_nuclear('<od>')".  

The combined (electronic + nuclear) model spans cold solids through warm/hot  
dense matter and on to fusion plasmas: at low T the nuclear term is the classical  
elastic ion-ion ("nuclear") stopping; at finite T it becomes the projectile-ion /  
target-ion drag of the Fokker-Planck slowing-down (the same physics as the ion  
term in Mehlhorn, J. Appl. Phys. 52, 6522, 1981).  

Examples:  
1. Proton in cold aluminum, nuclear stopping with all three potentials + plot:  
python dedx.py --zt=13 --d=2.7 --t=0.025 --aa=2 --nuc=1 --npot=gk,ionsphere,yukawa --plot=1 --od=ColdAl  

2. Proton in warm water, Gordon-Kim nuclear stopping only:  
python dedx.py --fc=H2O --d=1.0 --t=10.0 --aa=2 --nuc=1 --npot=gk --od=WarmH2O  

Validation / example drivers:  
python nuc_test.py fig1        # cold Al, three potentials vs NIST/ZBL nuclear  
python nuc_test.py sweep       # Te=1/10/100/1000 eV, nuclear growth + range cut  
python nuc_fac_compare.py      # FAC implementation vs the paper (SCAALP), proton-Al  
python nuc_gap_fig.py          # total = best-kappa electronic + GK nuclear vs PSTAR  
python nuc_ion_verify.py       # finite-Ti ion stopping vs the Fokker-Planck drag  
python alpha_dt_verify.py      # 3.5 MeV alpha in DT: e/ion crossover ~30 keV, range  

Note: the analytic Yukawa potential overestimates the nuclear stopping at high  
projectile energy (its long screened-attractive tail); nuclear stopping there  
is negligible vs electronic, and gk is the recommended potential.  
