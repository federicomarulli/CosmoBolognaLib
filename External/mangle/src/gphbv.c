/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Correction to boundary and vertex terms per subroutine gspher
  arising from disjoint polygons which abut along a circle.

  This is a c interface to fortran subroutine gphbv.

   Input: poly is a polygon.
	  np = number of caps of 1st polygon.
	  bnd = which cap of the polygon is the abutting boundary.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: bound, vert = as described in comments to gspher.s.f .
  Return value:  0 if ok;
		-1 if could not allocate temporary memory.
*/
int gphbv(polygon *poly, int np, int bnd, long double *tol, long double bound[2], long double vert[2])
{
    int i;
    /* work arrays */
    int *iord;
    long double *phi;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gphbv: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gphbv: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }

    /* translate from c to fortran convention */
    i = bnd + 1;

    /* the fortran routine */
    gphbv_(bound, vert, poly->rp, poly->cm, &poly->np, &np, &poly->np, &i, tol, phi, iord);

    /* free work arrays */
    free(iord);
    free(phi);

    return(0);
}
