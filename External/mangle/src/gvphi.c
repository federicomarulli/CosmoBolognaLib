/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Point at the centre of that part of a circle which
  (a) lies inside a polygon, and
  (b) contains or is closest to a specified unit vector.

  This is a c interface to fortran subroutine gvphi.

   Input: poly is a polygon.
	  rp, cm define a circle.
	  vi = unit vector.
	  *tol = angle within which to merge multiple intersections.
  Output: angle = length of that part of boundary segment which
		  (a) lies inside the polygon, and
		  (b) contains, or is closest to, unit vector vi.
		= 0. if the boundary lies entirely outside the polygon.
	  v = unit vector at centre of said boundary segment.
  Return value:  0 if ok;
		-1 if failed to allocate memory.
*/
int gvphi(polygon *poly, vec rp, long double cm, vec vi, long double *tol, long double *angle, vec v)
{
    /* work arrays */
    int *iord;
    long double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gvphi: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gvphi: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }

    /* fortran routine */
    gvphi_(angle, v, poly->rp, poly->cm, &poly->np, rp, &cm, vi, tol, phi, iord);

    /* free work arrays */
    free(iord);
    free(phi);

    return(0);
}
