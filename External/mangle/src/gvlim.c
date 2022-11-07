/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/* number of extra vertices to allocate, to allow for expansion */
#define DNV		4

/*------------------------------------------------------------------------------
  Points on polygon nearest to and farthest from unit direction vi,
  and minimum and maximum values of cm=1-cosl(theta)
  between polygon and unit vector vi.

  This is a wrapper around gvlim,
  that calls gvlim until the arrays are large enough.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  vi = unit vector.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  *vmin_p, vmax_p = pointers to arrays vmin[nv], vmax[nv]
		  giving nearest and farthest points from unit vector vi
		  on each edge; memory for the arrays is allocated.
	  *cmvmin_p, *cmvmax_p = pointer to arrays cmvmin[nv], cmvmax[nv]
		  giving minimum and maximum cm=1-cosl(theta)
		  between each edge and unit vector vi;
		  memory for the arrays is allocated.
	  *cmpmin_p, *cmpmax_p = pointer to arrays cmpmin[np], cmpmax[np]
		  giving minimum and maximum cm=1-cosl(theta)
		  between each circle and unit vector vi;
		  memory for the arrays is allocated.
	  *ipv_p = pointer to array ipv[nv] containing cap number of iv'th edge;
		   that is, ve[iv] lie on cap number ipv;
		   memory for the array is allocated.
	  *gp_p = pointer to array gp[np] giving group circle belongs to;
		  memory for the array is allocated.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  *ev_p = pointer to array ev[nv] of end indices;
		  memory for the array is allocated.
  Return value: 0 if ok;
		1 if fatal degenerate intersection of boundaries;
		-1 if failed to allocate memory.
*/
int gvlims(polygon *poly, int vcirc, long double *tol, vec vi, int *nv, vec **vmin_p, vec **vmax_p, long double **cmvmin_p, long double **cmvmax_p, long double **cmpmin_p, long double **cmpmax_p, int **ipv_p, int **gp_p, int *nev, int *nev0, int **ev_p)
{
    static int nvmax = 0, npmax = 0;
    static int *ipv = 0x0, *gp = 0x0, *ev = 0x0;
    static long double *cmvmin = 0x0, *cmvmax = 0x0, *cmpmin = 0x0, *cmpmax = 0x0;
    static vec *vmin = 0x0, *vmax = 0x0;

    int ier;

    /* putative maximum number of vertices */
    if (poly->np <= 1) {
	*nv = poly->np;
    } else if (poly->np <= 4) {
	*nv = poly->np * (poly->np - 1);
    } else {
	*nv = 6 * (poly->np - 2);
    }

    /* keep trying till the arrays are big enough */
    do {

	/* make sure that allocated arrays contain enough space */
	if (!vmin || !vmax || !cmvmin || !cmvmax || !ipv || !ev || *nv > nvmax) {
	    if (vmin) free(vmin);
	    if (vmax) free(vmax);
	    if (cmvmin) free(cmvmin);
	    if (cmvmax) free(cmvmax);
	    if (ipv) free(ipv);
	    if (ev) free(ev);
	    vmin = (vec *) malloc(sizeof(vec) * (*nv + DNV));
	    if (!vmin) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d vecs\n", *nv + DNV);
		return(-1);
	    }
	    vmax = (vec *) malloc(sizeof(vec) * (*nv + DNV));
	    if (!vmax) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d vecs\n", *nv + DNV);
		return(-1);
	    }
	    cmvmin = (long double *) malloc(sizeof(long double) * (*nv + DNV));
	    if (!cmvmin) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d long doubles\n", *nv + DNV);
		return(-1);
	    }
	    cmvmax = (long double *) malloc(sizeof(long double) * (*nv + DNV));
	    if (!cmvmax) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d long doubles\n", *nv + DNV);
		return(-1);
	    }
	    ipv = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ipv) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d ints\n", *nv + DNV);
		return(-1);
	    }
	    ev = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ev) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d ints\n", *nv + DNV);
		return(-1);
	    }
	    nvmax = *nv + DNV;
	}
	if (!cmpmin || !cmpmax || !gp || poly->np > npmax) { 
	    if (cmpmin) free(cmpmin);
	    if (cmpmax) free(cmpmax);
	    if (gp) free(gp);
	    cmpmin = (long double *) malloc(sizeof(long double) * (poly->np + DNV));
	    if (!cmpmin) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d long doubles\n", poly->np + DNV); 
		return(-1); 
	    }
	    cmpmax = (long double *) malloc(sizeof(long double) * (poly->np + DNV));
	    if (!cmpmax) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d long doubles\n", poly->np + DNV); 
		return(-1); 
	    }
	    gp = (int *) malloc(sizeof(int) * (poly->np + DNV));
	    if (!gp) {
		fprintf(stderr, "gvlims: failed to allocate memory for %d ints\n", poly->np + DNV); 
		return(-1); 
	    }
	    npmax = poly->np + DNV;
	}

	/* compute vertices of polygon */
	ier = gvlim(poly, vcirc, tol, vi, nvmax, nv, vmin, vmax, cmvmin, cmvmax, cmpmin, cmpmax, ipv, gp, nev, nev0, ev);
	if (ier) return(ier);

    } while (*nv > nvmax);

    /* point arguments at allocated arrays */
    *vmin_p = vmin;
    *vmax_p = vmax;
    *cmvmin_p = cmvmin;
    *cmvmax_p = cmvmax;
    *cmpmin_p = cmpmin;
    *cmpmax_p = cmpmax;
    *ev_p = ev;

    return(0);
}

/*------------------------------------------------------------------------------
  Points on polygon nearest to and farthest from unit direction vi,
  and minimum and maximum values of cm=1-cosl(theta)
  between polygon and unit vector vi.

  This is a c interface to fortran subroutine gvlim.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  vi = unit vector.
	  nvmax = dimension of v[nvmax] and ev[nmax].
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  vmin[nv], vmax[nv] = arrays giving nearest and farthest
		    points from unit vector vi on each edge.
	  cmvmin[nv], cmvmax[nv] = minimum and maximum cm=1-cosl(theta)
		    between each edge and unit vector vi.
	  cmpmin[np], cmvmax[np] = minimum and maximum cm=1-cosl(theta)
		    between each circle and unit vector vi.
	  ipv[nv] = cap number of vertices/edges;
		    that is, vertex v[i] and edge points ve[i]
		    lie on cap number ipv.
	  gp[np] = which group of intersecting circles each circle belongs to.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  ev[nev] = end indices of each connected sequence of vertices.
  Return value:  0 if ok;
		 1 if fatal degenerate intersection of boundaries;
		-1 if failed to allocate memory.
*/
int gvlim(polygon *poly, int vcirc, long double *tol, vec vi, int nvmax, int *nv, vec vmin[/*nvmax*/], vec vmax[/*nvmax*/], long double cmvmin[/*nvmax*/], long double cmvmax[/*nvmax*/], long double cmpmin[/*poly->np*/], long double cmpmax[/*poly->np*/], int ipv[/*nvmax*/], int gp[/*poly->np*/], int *nev, int *nev0, int ev[/*nvmax*/])
{
    logical ldegen;
    /* work arrays */
    int *iord, *iwk;
    long double *phi, *wk;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }
    iwk = (int *) malloc(sizeof(int) * nvmax * 4);
    if (!iwk) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d ints\n", nvmax * 4);
	return(-1);
    }
    wk = (long double *) malloc(sizeof(long double) * nvmax);
    if (!wk) {
	fprintf(stderr, "gvlim: failed to allocate memory for %d long doubles\n", nvmax);
	return(-1);
    }

    /* fortran routine */
    gvlim_(vmin, vmax, cmvmin, cmvmax, cmpmin, cmpmax, ipv, gp, ev, &nvmax, nv, nev, nev0, poly->rp, poly->cm, &poly->np, vi, &vcirc, tol, phi, iord, wk, iwk, &ldegen);

    /* number of vertices exceeds putative maximum */
    if (poly->np >= 5 && *nv > 6 * (poly->np - 2)) {
	msg("CONGRATULATIONS!  YOU HAVE DISCOVERED A POLYGON WITH 5 OR MORE CAPS\n");
	msg("(IT HAS %d CAPS) THAT HAS MORE THAN %d VERTICES (IT HAS %d VERTICES).\n", poly->np, 6 * (poly->np - 2), *nv);
	msg("(Either that or you have found a bug.)\n");
	msg("PLEASE EMAIL ME Andrew.Hamilton@colorado.edu THE GOOD NEWS,\n");
	msg("ALONG WITH A POLYGON FILE CONTAINING THE POLYGON THAT DID IT.\n");
	msg("THANKS!\n");
	dump_poly(1, &poly);
	msg("AND THERE'S THE POLYGON FILE I'D LIKE YOU TO SEND.  THANKS!\n");
    }

    /* free work arrays */
    free(iord);
    free(phi);
    free(iwk);
    free(wk);

    /* fatal intersection of boundaries */
    if (ldegen) return(1);

    return(0);
}
