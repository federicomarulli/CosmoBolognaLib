/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Angle along circle centred at unit vector rp,
  with radius th given by cm = 1 - cosl(th),
  lying inside polygon poly.
  The angle is 2 pi if the circle lies entirely inside the polygon.

  This is a c interface to fortran subroutine gphi.

   Input: poly is a polygon.
	  *tol = angle within which to merge multiple intersections.
	  rp = unit vector at the centre of the circle.
	  cm = 1 - cosl(th), where th is the angular radius of the circle.
  Output: *angle = angle along circle, in radians.
  Return value:  0 if ok;
		-1 if failed to allocate memory.
*/
int gphi(polygon *poly, long double *tol, vec rp, long double cm, long double *angle)
{
    /* work arrays */
    int *iord;
    long double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gphi: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gphi: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }

    /* fortran routine */
    gphi_(angle, poly->rp, poly->cm, &poly->np, rp, &cm, tol, phi, iord);

    /* free work arrays */
    free(iord);
    free(phi);

    return(0);
}
