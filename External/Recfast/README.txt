=============================================================================================

This directory contains a C++ version of Recfast (v1.0), written by Jens Chluba (Oct 2010). 
This code is based on the C-version by Seager et al. 1999 (this version can be found at 
https://www.cfa.harvard.edu/~sasselov/rec/), but was strongly rewritten and updated to:
 
- perform the normal Recfast computation (Recfast v1.4.2 with normal helium rec)
- the recombination corrections given by Chluba & Thomas 2010. 
- allows running a simple Dark matter annihilation module.

These tasks are achieved using the simple 3-level approach of Recfast. To solve the 
system of ODE's we use a simplified version of the stiff ODE-solver developed by 
Chluba, Vasil & Dursi, 2010. This solver allows us to avoid the switches in the original 
Recfast-code.

When using this code please cite:

Seager et al., 1999, ApJ, 523, L1-L5
Chluba & Thomas, 2010, arXiv:1010.3631v2
Rubino-Martin et al, 2010, MNRAS, 403, 439-452
Chluba, 2010, MNRAS, 402, 1195-1207
Chluba, Vasil & Dursi, 2010, MNRAS, 407, 599-612

as well as

Switzer & Hirata, 2008, Phys.Rev.D. 77, 083006Grin & Hirata, 2010, Phys.Rev.D, 81, 083005Ali-Haimoud & Hirata, 2010, Phys.Rev.D, 82, 063521

=============================================================================================

=============================================================================================

Installation:

To compile the code follow the following steps:

(i) in "Makefile" set "CC" to the C++ compiler that should be used
(ii) type "make" to compile the code. This creates the executable: "Recfast++"

To clean up type "make clean" or "make tidy" (see makefile for difference)

=============================================================================================

=============================================================================================

Running the code:

The code is run by invoking

./Recfast++ runfiles/parameters.dat

The output of the code will be written in "./output/." The filenames will depend on the
chosen runmode.

The parameter file in ./runfiles contains the cosmological parameters and setting for 
the different runmodes.

The entries of "runfiles/parameters.dat" are:

Yp (0.24) 	== helium mass faction
T0 (2.725)	== CMB temperature at z = 0

Om (0.26)	== Omega matter
Ob (0.044)	== Omega Baryon
OL (0.0)	== Omega Lambda ( --> will be determined internally using OL=1-Ok-Om-Orel.
				  Here O_rel== contribution from relativistic species)
Ok (0.0)	== Omega K (curvature parameter Ok=1-O_tot; Ok=0.0 --> spatially flat)

h100 (0.71)	== H0/100
Nnu (3.04)	== effective number of neutrinos
F (1.14)	== Recfast fudge-factor

fDM (0.0)	== DM annihilation efficiency. Typical value would be fDM ~ 2e-24 eV/s 
		   (see Chluba 2010 for details. The corresponding functions are defined 
		    in ./src/DM_annihilation.Recfast.cpp) 
switch (1) 	== 1: includes the recombination correction function given by Chluba & Thomas 2010
		   0: normal Recfast run, with no corrections included. This is equivalent to running
		      Recfastv1.5 with helium-flag=0 and no corrections from Rubino-Martin 2010.
		      Furthermore, we used our on ODE-solver, so that no switching of the ODE system
		      is needed.

=============================================================================================

