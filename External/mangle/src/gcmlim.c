/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Minimum and maximum values of cm = 1-cosl(th) between polygon
  and a unit vector rp.

  This is a c interface to fortran subroutine gcmlim.

   Input: poly is a polygon.
	  *tol = angle within which to merge multiple intersections.
	  rp = unit vector.
  Output: minimum and maximum values of cm = 1-cosl(th).
  Return value:  0 if ok;
		-1 if failed to allocate memory.
*/
int gcmlim(polygon *poly, long double *tol, vec rp, long double *cmmin, long double *cmmax)
{
    /* work arrays */
    int *iord;
    long double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gcmlim: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gcmlim: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }

    /* fortran routine */
    gcmlim_(poly->rp, poly->cm, &poly->np, rp, cmmin, cmmax, tol, phi, iord);

    /* free work arrays */
    free(iord);
    free(phi);

    return(0);
}
