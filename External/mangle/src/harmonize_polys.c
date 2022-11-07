/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"
#include "pi.h"

/* advise how many polygons done if lmax >= this */
#define LMAX_ADVICE		250

/*------------------------------------------------------------------------------
  Spherical harmonics of sum of weighted polygons.

   Input: poly = array of pointers to npoly polygons.
	  npoly = number of polygons in poly array.
	  mtol = initial angular tolerance in radians within which to merge multiple intersections.
	  lmax = maximum harmonic number.
  Output: w = harmonics;
	      NW = ((lmax + 1)(lmax + 2))/ 2 is defined in harmonics.h.
  Return value: number of polygons for which spherical harmonics were computed,
		or -1 if error occurred.
*/
int harmonize_polys(int npoly, polygon *poly[/*npoly*/], long double mtol, int lmax, harmonic w[/*NW*/])
{
    int accelerate, i, ier, ip, ipoly, iq, ir, isrect, iw, naccelerate, ndone, ner, nrect;
    long double azmin, azmax, elmin, elmax, azmn, azmx, elmn, elmx, tol;
    /* work array contains harmonics of single polygon */
    harmonic *dw;
    /* work arrays to deal with possible acceleration */
    int *iord, *ir_to_ip;
    long double *elord;

    /* work arrays */
    dw = (harmonic *) malloc(sizeof(harmonic) * NW);
    if (!dw) {
	fprintf(stderr, "harmonize_polys: failed to allocate memory for %d harmonics\n", NW);
	return(-1);
    }
    iord = (int *) malloc(sizeof(int) * npoly);
    if (!iord) {
	fprintf(stderr, "harmonize_polys: failed to allocate memory for %d ints\n", npoly);
	return(-1);
    }
    ir_to_ip = (int *) malloc(sizeof(int) * npoly);
    if (!ir_to_ip) {
	fprintf(stderr, "harmonize_polys: failed to allocate memory for %d ints\n", npoly);
	return(-1);
    }
    elord = (long double *) malloc(sizeof(long double) * npoly);
    if (!elord) {
	fprintf(stderr, "harmonize_polys: failed to allocate memory for %d long doubles\n", npoly);
	return(-1);
    }

    /* zero harmonics of mask */
    for (iw = 0; iw < NW; iw++) {
	for (i = 0; i < IM; i++) {
	    w[iw][i] = 0.;
	}
    }

    /* determine which polygons are rectangles, for which acceleration may be possible */
    nrect = 0;
    ir = npoly;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	isrect = poly_to_rect(poly[ipoly], &azmin, &azmax, &elmin, &elmax);
	if (isrect && poly[ipoly]->weight != 0.) {
	    ir_to_ip[nrect] = ipoly;
	    elord[nrect] = elmin * 1.e8 + elmax;
	    nrect++;
	} else {
	    ir--;
	    ir_to_ip[ir] = ipoly;
	}
    }
    msg("harmonize_polys: %d polygons are rectangles, for which acceleration may be possible\n",
	nrect);

    /* order rectangles by elmin, elmax */
    findtop(elord, nrect, iord, nrect);

    /* do each polygon */
    ndone = 0;
    naccelerate = 0;
    ner = 0;
    if (lmax >= LMAX_ADVICE) msg("doing polygon number (of %d):\n", npoly);
    for (ip = 0; ip < npoly; ip++) {
	if (lmax >= LMAX_ADVICE) msg(" %d", ip);
	accelerate = 0;
	/* rectangle */
	if (ip < nrect) {
	    ir = iord[ip];
	    ipoly = ir_to_ip[ir];
	    poly_to_rect(poly[ipoly], &azmin, &azmax, &elmin, &elmax);
	    /* does previous rectangle have same elevation limits? */
	    if (ip > 0) {
		iq = iord[ip - 1];
		iq = ir_to_ip[iq];
		poly_to_rect(poly[iq], &azmn, &azmx, &elmn, &elmx);
		/* if so, use acceleration */
		if (elmn == elmin && elmx == elmax) accelerate = 1;
	    }
	    /* if not, does next rectangle have same elevation limits? */
	    if (!accelerate && ip + 1 < nrect) {
		iq = iord[ip + 1];
		iq = ir_to_ip[iq];
		poly_to_rect(poly[iq], &azmn, &azmx, &elmn, &elmx);
		/* if so, worth accelerating */
		if (elmn == elmin && elmx == elmax) accelerate = 1;
	    }
	    /* accelerated computation */
	    if (accelerate) {
		ier = gsphra(azmin, azmax, elmin, elmax, lmax, dw);
		if (ier == -1) return(-1);
	    /* standard computation */
	    } else {
		tol = mtol;
		ier = gsphr(poly[ipoly], lmax, &tol, dw);
		if (ier == -1) return(-1);
	    }
	/* non-rectangle */
	} else {
	    ipoly = ir_to_ip[ip];
	    /* zero weight polygon requires no computation */
	    if (poly[ipoly]->weight == 0.) {
		ndone++;
		continue;
	    } else {
		tol = mtol;
		ier = gsphr(poly[ipoly], lmax, &tol, dw);
		if (ier == -1) return(-1);
	    }
	}
	/* computation failed */
	if (ier) {
	    ner++;
	    if (lmax >= LMAX_ADVICE) msg("\n");
	    fprintf(stderr, "harmonize_polys: computation failed for polygon %d; discard it\n", ipoly);
	/* success */
	} else {
	    naccelerate += accelerate;
	    ndone++;
	    /* increment harmonics of region */
	    for (iw = 0; iw < NW; iw++) {
		for (i = 0; i < IM; i++) {
		    w[iw][i] += dw[iw][i] * poly[ipoly]->weight;
		}
	    }
	}
    }
    if (lmax >= LMAX_ADVICE) msg("\n");

    /* number of computations that were accelerated */
    msg("computation was accelerated for %d rectangles\n", naccelerate);
    /* advise */
    if (ner > 0) {
	fprintf(stderr, "harmonize_polys: discarded %d polygons for which computations failed\n", ner);
    }
    msg("spherical harmonics of %d weighted polygons accumulated\n", ndone);

    /* free work arrays */
    free(dw);
    free(iord);
    free(ir_to_ip);
    free(elord);

    return(ndone);
}
