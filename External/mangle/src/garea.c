/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "logical.h"
#include "manglefn.h"

/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP		4

/*------------------------------------------------------------------------------
  Area of polygon.

  This is a c interface to fortran subroutine garea.

   Input: poly is a polygon.
	  verb = 0 to suppress messages from garea, even fatal ones;
		 1 to allow messages from garea.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *area = area of polygon.
  Return value:  0 if ok;
		 1 if fatal error;
		-1 if failed to allocate memory.
*/
int garea(polygon *poly, long double *tol, int verb, long double *area)
{
    static polygon *dpoly = 0x0;
    logical ldegen;
    int ier, ipmin, ipoly, np;
    long double cmmin, darea;
    /* work arrays */
    int *iord;
    long double *phi;

    /* smallest cap of polygon */
    cmminf(poly, &ipmin, &cmmin);

    /* number of caps in polygon to be passed to garea_ */
    np = (poly->np >= 2 && cmmin > 1.)? poly->np + 1 : poly->np;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * np * 2);
    if (!iord) {
	fprintf(stderr, "garea: failed to allocate memory for %d ints\n", np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * np * 2);
    if (!phi) {
	fprintf(stderr, "garea: failed to allocate memory for %d long doubles\n", np * 2);
	return(-1);
    }

    /* <= 1 caps, or smallest cap has area <= pi */
    if (np == poly->np) {
	/* fortran routine */
  	garea_(area, poly->rp, poly->cm, &poly->np, tol, &verb, phi, iord, &ldegen);

    /* >= 2 caps, and smallest cap has area > pi */
    } else {
	/* make sure dpoly contains enough space */
	ier = room_poly(&dpoly, np, DNP, 0);
	if (ier == -1) {
	    fprintf(stderr, "garea: failed to allocate memory for polygon of %d caps\n", np + DNP);
	    return(-1);
	}

	/* make polygon dpoly with extra cap */
	poly_polyn(poly, poly, ipmin, 1, dpoly);

	/* make the extra cap 1/2 the area of the smallest cap */
	dpoly->cm[poly->np] = cmmin / 2.;

	/* zero area */
	*area = 0.;

	for (ipoly = 0; ipoly < 2; ipoly++) {
	    /* fortran routine */
	    garea_(&darea, dpoly->rp, dpoly->cm, &dpoly->np, tol, &verb, phi, iord, &ldegen);

	    /* accumulate area */
	    *area += darea;

	    /* fatal error */
	    if (ldegen) break;

	    /* the complement of the extra cap */
	    dpoly->cm[poly->np] = - dpoly->cm[poly->np];
	}

    }

    /* free work arrays */
    free(iord);
    free(phi);

    /* fatal error */
    if (ldegen){
      fprintf(stderr,"garea: fatal error in polygon %d, pixel %d\n", poly->id,poly->pixel);
      return(1);
    }
    return(0);
}
