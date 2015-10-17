
	This codes implements the multi-point propagator expansion of the matter power spectrum
	up to two loops.

	[1] Installation ::

	MPTbreeze uses the CUBA Library for multidimensional integration to perform two loops
	integrals. The CUBA library can be downloaded from :

        -- http://www.feynarts.de/cuba (and should be simple to install) --

	Note that our integrations have been optimized to version 1.4 of this library (we have no 
	tested it with newer versions).

        We provide a simple compilation script for our code (using ifort) : compile.sh 

	To use it you will have to change the location of libcuba.a to the one in your system, 
	and possibly some flags that might be obsolete in your OS.

	MPTbreeze can also be compiled with other Fortran compilers such as gfortran, but the 
	performance might get worse.

	[2] Input ::

	- Running ./mptbreeze details the running instructions and command line flags

        General inputs are the matter transfer function, relevant cosmological parameters (matter 
	density, spectral tilt, equation of state of DE and amplitude of fluctuations at 8Mpc/h) and 
	the desired redshift output. 

	There are two ways to pass this information to the code : 

	- One is with individual command line inputs. In this case the transfer function must be a 
	two column file (k-Tf)

	- Another is to directly pass the params.ini file used to generate the transfer function 
	in CAMB, which contains much of the relevant information. Additional *needed* information is 
	given as trailing command line options (-sigma8 for sig8, -redshift for output redshift and 
	-filePk for output Pk filename)

	We provide two example scripts, one for each case : run_test.sh and run_test_camb.sh

	Notice that all parameters have default values (declared at the beginning of the code). 
	Make sure you are passing your own.

	[3] Output :: 

	The code evaluates the nonlinear power spectrum from dk = 0.005 h/Mpc up to k_nl = sigv^{-1}
	in steps of size dk (linearly spaced). This election is completely arbitrary and can be
	changed in line 260. 

	Notice however that beyond k ~ k_nl the predictions shows a exponential suppression that
	limits its validity.

	- The output format is given as 

        k P0 P1 P2L 

	corresponding to 0,1 and 2 loops integrals. 

	- The full model prediction would be P0 + P1 + P2 
        
	again one can change this output format trivially

	[4] General Comments and Performance ::

        The speed of the code is mostly controlled by the relative error for the two-loop
        integration, this is set in header-vegas.f (variable -> epsrel). 

	- We recommend epsrel = 0.01 - 0.005 (default value is 0.01) 

	If no (power spectrum) output filename is given the default name is : 
	
	- MPTbreeze_"filename of transfer function"



	** If you find any trouble or you like/dislike something email martincrocce@gmail.com **
