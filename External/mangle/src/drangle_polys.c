/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"
#include "pi.h"

#define TWOPI		(2. * PI)

static int *iord = 0x0;
static long double *cmmin = 0x0, *cmmax = 0x0;

/*------------------------------------------------------------------------------
  Minimum and maximum values of cm = 1-cosl(th) between each of npoly polygons
  and a unit vector rp.

   Input: poly = array of pointers to npoly polygons.
	  npoly = number of polygons in poly array.
	  mtol = initial angular tolerance in radians within which to merge multiple intersections.
	  rp = unit vector.
  Return value: number of polygons done;
		-1 if error.
*/
int cmlim_polys(int npoly, polygon *poly[/*npoly*/], long double mtol, long double rp[3])
{
    int ier, ipoly;
    long double tol;

    /* allocate memory for cmmin, cmmax, iord */
    if (!cmmin) {
	cmmin = (long double *) malloc(sizeof(long double) * npoly);
    } else {
	cmmin = (long double *) realloc(cmmin, sizeof(long double) * npoly);
    }
    if (!cmmin) {
	fprintf(stderr, "cmlim_polys: failed to allocate memory for %d long doubles\n", npoly);
	return(-1);
    }
    if (!cmmax) {
	cmmax = (long double *) malloc(sizeof(long double) * npoly);
    } else {
	cmmax = (long double *) realloc(cmmax, sizeof(long double) * npoly);
    }
    if (!cmmax) {
	fprintf(stderr, "cmlim_polys: failed to allocate memory for %d long doubles\n", npoly);
	return(-1);
    }
    if (!iord) {
	iord = (int *) malloc(sizeof(int) * npoly);
    } else {
	iord = (int *) realloc(iord, sizeof(int) * npoly);
    }
    if (!iord) {
	fprintf(stderr, "cmlim_polys: failed to allocate memory for %d ints\n", npoly);
	return(-1);
    }

    /* min, max distances between rp and each polygon */
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	if (poly[ipoly]->weight == 0.) continue;
	tol = mtol;
	ier = gcmlim(poly[ipoly], &tol, rp, &cmmin[ipoly], &cmmax[ipoly]);
	if (ier) return(-1);
    }

    /* order polygons in increasing order of cmmin */
    findbot(cmmin, npoly, iord, npoly);

    /* number of polygons done */
    return(npoly);
}

/*------------------------------------------------------------------------------
  Angles within mask along circle centred in unit direction rp, with radii th.

   Input: poly = array of pointers to npoly polygons.
	  npoly = number of polygons in poly array.
	  mtol = initial angular tolerance in radians within which to merge multiple intersections.
	  rp = unit vector.
	  nth = number of angular radii.
	  cm = array of 1-cosl(angular radii).
  Output: dr = array containing angles in radians.
  Return value: number of angular radii done;
		-1 if error.
*/
int drangle_polys(int npoly, polygon *poly[/*npoly*/], long double mtol, long double rp[3], int nth, long double cm[/*nth*/], long double dr[/*nth*/])
{
    int ier, ip, ipoly, ith;
    long double angle, tol;

    /* angle within mask at each angular radius */
    for (ith = 0; ith < nth; ith++) {
	/* accumulate angle within each polygon */
	dr[ith] = 0.;
	for (ip = 0; ip < npoly; ip++) {
	    ipoly = iord[ip];
	    /* zero weight polygon contributes nothing */
	    if (poly[ipoly]->weight == 0.) continue;
	    if (cm[ith] <= fabsl(cmmin[ipoly])) {
		/* polygon encloses circle */
		if (cmmin[ipoly] < 0.) {
		  angle = TWOPI;
		/* polygon excludes circle */
		} else if (cmmin[ipoly] >= 0.) {
		  angle = 0.;
		  /* done, given that cmmin are in increasing order */
		  break;
		}
	    } else if (cm[ith] >= fabsl(cmmax[ipoly])) {
		/* circle and polygon enclose each other */
		if (cmmax[ipoly] < 0.) {
		    angle = TWOPI;
		/* circle encloses polygon */
		} else if (cmmax[ipoly] >= 0.) {
		    angle = 0.;
		}
	    } else {
		/* circle intersects boundary of region */
		tol = mtol;
		ier = gphi(poly[ipoly], &tol, rp, cm[ith], &angle);
		if (ier) return(-1);
	    }
	    dr[ith] += poly[ipoly]->weight * angle;
	}
    }

    /* number of angular radii done */
    return(nth);
}
