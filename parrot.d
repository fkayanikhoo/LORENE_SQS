#################### PHYSICAL PARAMETERS ######################################
1     Relativity parameter: 1 = relativistic computation , 0 = Newtonian
0.45   ent_c : central enthalpy [c^2]
1200.     freq_si : rotation frequency [Hz]
1.0    fact_omega : 1.01 = search for the Keplerian frequency, 1. = otherwise.
2    Requested baryon mass [M_sol] (effective only if mer_mass > mer_max)
1.e32   Requested magnetic moment [A.m^2] (effective only if mer_mass > mer_magmom)
#################### MAGNETIC PARAMETERS ######################################
0.	Requested total charge(conduc=1), charge/baryon (conduc=0)[Lorene unit]
32000.	Requested CFA (current function amplitude) [Lorene unit]
0.	Initial charge (mer =< mer_mag) [Lorene unit]
3500.	Initial CFA    (mer =< mer_mag) [Lorene unit]
10	mer_mag : step at which magnetic quantites are plugged.
15	mer_change_mag : step at which they are increased.
25	mer_fix_mag : step at which they reach their final values.
1	mag_in_eos : use magnetic field value in the EoS (0: false, 1: true)
0	use_magnetisation : include magnetisation terms in equations (0: false,...)
#################### COMPUTATIONAL PARAMETERS #################################
700     mer_max : maximum number of steps
1.e-16   precis : threshold on the enthalpy relative change for ending the computation
30     mer_rot : step at which the rotation is switched on
120.    freq_ini_si : initial rotation frequency [Hz] (switched on at mer = mer_rot)
35      mer_change_omega : step at which f is changed to reach freq_si
45	     mer_fix_omega : step at which f must have reached freq_si
10     delta_mer_kep : number of steps after mer_fix_omega to search for Kepler.
0.3    thres_adapt : threhold on (dH/dr_eq)/dH/dr_pole) for the mapping adaptation
50000   mer_mass : step from which the baryon mass is forced to converge (if negative, variation of Omega)
50000   mer_magmom : step from which the magnetic moment is forced to converge
0.5     aexp_mass : exponent for the increase factor of the central enthalpy
0.5     relax : relaxation factor in the main iteration 
10       mermax_poisson : maximum number of steps in Map_et::poisson
1.5     relax_poisson :  relaxation factor in Map_et::poisson
1.e-14  precis_adapt : precision in Map_et::adapt
0       graph : 1 = graphical outputs during the computation 
#################### MULTI-GRID PARAMETERS ###################################
3	nz : total number of domains
1	nzet : number of domains inside the star
1	nzadapt : number of domains of where the mapping adaptation will be done.
33	nt: number of points in theta (the same in each domain)
1	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
65	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)	
33	1.	<-   nr	  &   min(r)  in domain 1
17	2.	<-   nr   &   min(r)  in domain 2
0.1	enthalpy defining boundary between domains 1 and 2
