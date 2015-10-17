
	implicit none

	!! Common arguments !! 

	integer   ndim, ncomp, mineval, maxeval, verbose, neval, fail
	real*8    epsrel, epsabs
	parameter (ndim = 5)
	parameter (ncomp = 1)
	parameter (mineval = 0)
	parameter (maxeval = 110000000)
	real*8    integral(ncomp), error(ncomp), prob(ncomp)

	!! VEGAS !!

	parameter (epsrel = 0.01)            !! Relative Error !! 
	parameter (epsabs = 0.0)      
	parameter (verbose = 0*8+0*4+0*2+0*1)
	integer   nstart, nincrease
	parameter (nstart = 4000)    
	parameter (nincrease = 700)

