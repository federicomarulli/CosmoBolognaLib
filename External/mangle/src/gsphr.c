/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "manglefn.h"
#include "pi.h"

/*------------------------------------------------------------------------------
  Spherical harmonics of polygon.

  This is a simplified c interface to fortran subroutine gspher.
  It returns the spherical harmonics, and does not worry about bound and vert.

   Input: poly is a polygon.
	  lmax = maximum harmonic number.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: w = array containing spherical harmonics of polygon;
	      NW = ((lmax + 1)(lmax + 2))/ 2 is defined in harmonics.h.
  Return value:  0 if ok;
		 1 if fatal error;
		-1 if could not allocate temporary memory.
*/
int gsphr(polygon *poly, int lmax, long double *tol, harmonic w[/*NW*/])
{
    logical ldegen;
    int i, ibv, ier, im, iphi, iw, lmax1, npc, nw, verb;
    long double area, bound[2], darea, vert[2];
    /* work arrays */
    int *iord;
    long double *v, *phw;

    /* determine area without 2 pi ambiguity, and a good value for tol */
    verb = 1;
    ier = garea(poly, tol, verb, &area);
    if (ier) return(ier);

    /* trivial case of zero area */
    if (area == 0.) {
	for (iw = 0; iw < NW; iw++) {
	    for (i = 0; i < IM; i++) w[iw][i] = 0.;
	}

	return(0);
    }

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gsphr: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phw = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phw) {
	fprintf(stderr, "gsphr: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }
    v = (long double *) malloc(sizeof(long double) * (lmax + 1));
    if (!v) {
	fprintf(stderr, "gsphr: failed to allocate memory for %d long doubles\n", lmax + 1);
	return(-1);
    }

    /* parameters */
    lmax1 = lmax + 1;
    im = IM;
    nw = NW;
    npc = 0;
    ibv = 0;
    iphi = 0;

    /* the fortran routine */
    gspher_(&darea, bound, vert, w, &lmax1, &im, &nw, poly->rp, poly->cm, &poly->np, &npc, &ibv, &iphi, tol, phw, iord, v, &ldegen);

    /* monopole harmonic without 2 pi/sqrtl(4 pi) ambiguity */
    w[0][0] = area / sqrtl(4. * PI);

    /* free work arrays */
    free(iord);
    free(phw);
    free(v);

    /* fatal error */
    if (ldegen) return(1);

    return(0);
}

/*------------------------------------------------------------------------------
  Accelerated computation of spherical harmonics of rectangle.

  This is a simplified c interface to fortran subroutine gsphera.
  It returns the spherical harmonics, and does not worry about bound and vert.

  The acceleration involves some overhead, and works only if two or more
  rectangles with the same elmin & elmax are computed in succession.
  The overhead means that the accelerated computation is actually slightly
  slower for just a single rectangle.

   Input: lmax = maximum harmonic number.
  Output: w = array containing spherical harmonics of polygon;
	      NW = ((lmax + 1)(lmax + 2))/ 2 is defined in harmonics.h.
  Return value:  0 if ok;
		-1 if could not allocate temporary memory.
*/
int gsphra(long double azmin, long double azmax, long double elmin, long double elmax, int lmax, harmonic w[/*NW*/])
{
    /* array used for acceleration */
    static long double *dw = 0x0;

    int ibv, im, lmax1, nw;
    long double area, bound[2], vert[2];
    /* work array */
    long double *v;

    /* allocate memory for work arrays */
    v = (long double *) malloc(sizeof(long double) * (lmax + 1));
    if (!v) {
	fprintf(stderr, "gsphra: failed to allocate memory for %d long doubles\n", lmax + 1);
	return(-1);
    }

    /* parameters */
    lmax1 = lmax + 1;
    im = IM;
    nw = NW;
    ibv = 0;

    /* dw contains array that is pre-computed, then used by all rects with same elmin, elmax */
    if (!dw) {
	dw = (long double *) malloc(sizeof(long double) * NW);
	if (!dw) {
	    fprintf(stderr, "gsphra: failed to allocate memory for %d long doubles\n", NW);
	    return(-1);
	}
    }

    /* fortran routine */
    gsphera_(&area, bound, vert, w, &lmax1, &im, &nw, &ibv, &azmin, &azmax, &elmin, &elmax, v, dw);

    /* free work array */
    free(v);

    return(0);
}
