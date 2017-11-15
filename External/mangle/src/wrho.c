/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Value of summed harmonics at a specified position.

  This is a c interface to fortran function wrho.

   Input: az, el = azimuth and elevation in radians.
	  lmax = maximum harmonic number of harmonics to include.
	  mmax = include azimuthal harmonics -mmax to + mmax
	       = lmax normally, but other values are acceptable.
	  w(i, lm) = spherical transform, dimensioned w(im, nw)
	     w(i,lm), i=1,im, lm=l*(l+1)/2+m+1, l=0,lmax, m=0,l;
	     w(1,lm) is real part, w(2,lm) is imaginary part;
	     NW = ((lmax + 1)(lmax + 2))/ 2 is defined in harmonics.h.
	  lsmooth = smoothing harmonic number
		  = 0. (or < 0.) to skip smoothing.
	  esmooth = exponent of smoothing, ignored if lsmooth <= 0.
		  = 2. for Gaussian smoothing.
  Return value: wrho = sum_lm w_lm Y_lm(az, el)   if lsmooth = 0.
		       or
		       sum_lm w_lm Y_lm exp[- l(l+1) / lsmooth(lsmooth+1) ]
			  if lsmooth > 0. and esmooth = 2. (gaussian smoothing)
		       or
	     sum_lm w_lm Y_lm exp[- [l(l+1) / lsmooth(lsmooth+1)]^(esmooth/2) ]
			  if lsmooth > 0. and general esmooth.
*/
long double wrho(long double az, long double el, int lmax, int mmax, harmonic w[/*NW*/], long double lsmooth, long double esmooth)
{
    int im, nw;
    long double wrho;

    im = IM;
    nw = NW;

    wrho = wrho_(&az, &el, w, &lmax, &mmax, &im, &nw, &lsmooth, &esmooth);

    return(wrho);
}
