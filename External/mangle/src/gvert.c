/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/* number of extra vertices to allocate, to allow for expansion */
#define DNV		4

/*------------------------------------------------------------------------------
  Points on edges of polygon.

  This is a wrapper around gvert,
  that calls gvert until the arrays are large enough.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  per = 0 or 1 controls meaning of nve.
	  if per = 0:
	     nve = desired number of points on each edge, including vertex;
	  if per = 1:
	     nve = desired number of points per (2 pi) on each edge.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  *ve_p = pointer to array ve[nv][nve] of points on edges;
		  memory for the array is allocated.
	  *angle_p = pointer to array angle[nv] of lengths of edges;
		  memory for the array is allocated.
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
		-1 if could not allocate memory.
*/
int gverts(polygon *poly, int vcirc, long double *tol, int per, int nve, int *nv, vec **ve_p, long double **angle_p, int **ipv_p, int **gp_p, int *nev, int *nev0, int **ev_p)
{
    static int nvmax = 0, nvemax = 0, npmax = 0;
    static int *ipv = 0x0, *gp = 0x0, *ev = 0x0;
    static long double *angle = 0x0;
    static vec *ve = 0x0;

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
	if (!ve || !angle || !ipv || !ev || *nv > nvmax) {
	    if (ve) free(ve);
	    if (angle) free(angle);
	    if (ipv) free(ipv);
	    if (ev) free(ev);
	    if (nve > nvemax) nvemax = nve + DNV;
	    ve = (vec *) malloc(sizeof(vec) * (*nv + DNV) * nvemax);
	    if (!ve) {
		fprintf(stderr, "gverts: failed to allocate memory for %d x %d vecs\n", *nv + DNV, nvemax);
		return(-1);
	    }
	    angle = (long double *) malloc(sizeof(long double) * (*nv + DNV));
	    if (!angle) {
		fprintf(stderr, "gverts: failed to allocate memory for %d long doubles\n", *nv + DNV);
		return(-1);
	    }
	    ipv = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ipv) {
		fprintf(stderr, "gverts: failed to allocate memory for %d ints\n", *nv + DNV);
		return(-1);
	    }
	    ev = (int *) malloc(sizeof(int) * (*nv + DNV));
	    if (!ev) {
		fprintf(stderr, "gverts: failed to allocate memory for %d ints\n", *nv + DNV);
		return(-1);
	    }
	    nvmax = *nv + DNV;
	} else if (nve > nvemax) {
	    if (ve) free(ve);
	    ve = (vec *) malloc(sizeof(vec) * nvmax * (nve + DNV));
	    if (!ve) {
		fprintf(stderr, "gverts: failed to allocate memory for %d x %d vecs\n", nvmax, nve + DNV);
		return(-1);
	    }
	    nvemax = nve + DNV;
	}
	if (!gp || poly->np > npmax) { 
	    if (gp) free(gp);
	    gp = (int *) malloc(sizeof(int) * (poly->np + DNV));
	    if (!gp) {
		fprintf(stderr, "gverts: failed to allocate memory for %d ints\n", poly->np + DNV); 
		return(-1); 
	    }
	    npmax = poly->np + DNV;
	}

	/* compute vertices of polygon */
	ier = gvert(poly, vcirc, tol, nvmax, per, nve, nv, ve, angle, ipv, gp, nev, nev0, ev);
	if (ier) return(1);

    } while (*nv > nvmax);

    /* point ve_p, angle_p, ipv_p, gp_p, and ev_p at ve, angle, ipv, gp, and ev */
    *ve_p = ve;
    *angle_p = angle;
    *ipv_p = ipv;
    *gp_p = gp;
    *ev_p = ev;

    return(0);
}

/*------------------------------------------------------------------------------
  Vertices of polygon.

  This is a c interface to fortran subroutine gvert.

   Input: poly is a polygon.
	  vcirc = 1 to return vertices and midpoints also for bounding circles
		  which have no intersections;
		= 0 not so.
	  nvmax = dimension of ve[nvmax][nve], angle[nvmax], ipv[nvmax], and ev[nvmax].
	  per = 0 or 1 controls meaning of nve.
	  if per = 0:
	     nve = desired number of points on each edge, including vertex;
	  if per = 1:
	     nve = desired number of points per (2 pi) on each edge.
  Input/Output: *tol = angle within which to merge multiple intersections.
  Output: *nv = number of vertices.
	  ve[nv][nve] = points on edges of polygon.
	  angle[nv] = angular lengths of edges of polygon.
	  ipv[nv] = cap number of vertices/edges;
		    that is, vertex v[i] and edge points ve[i]
		    lie on cap number ipv.
	  gp[np] = which group of intersecting circles each circle belongs to.
	  *nev = number of connected sequences of vertices.
	  *nev0 = number of bounding circles which have no intersections.
	  ev[nev] = end indices of each connected sequence of vertices.
  Return value:  0 if ok;
		 1 if fatal intersection of boundaries;
		-1 if failed to allocate memory.
*/
int gvert(polygon *poly, int vcirc, long double *tol, int nvmax, int per, int nve, int *nv, vec ve[/*nvmax * nve*/], long double angle[/*nvmax*/], int ipv[/*nvmax*/], int gp[/*poly->np*/], int *nev, int *nev0, int ev[/*nvmax*/])
{
    int iv;
    logical ldegen;
    /* work arrays */
    int *iord, *iwk;
    long double *phi, *wk;

    /* allocate memory for work arrays */
    iord = (int *) malloc(sizeof(int) * poly->np * 2);
    if (!iord) {
	fprintf(stderr, "gvert: failed to allocate memory for %d ints\n", poly->np * 2);
	return(-1);
    }
    phi = (long double *) malloc(sizeof(long double) * poly->np * 2);
    if (!phi) {
	fprintf(stderr, "gvert: failed to allocate memory for %d long doubles\n", poly->np * 2);
	return(-1);
    }
    iwk = (int *) malloc(sizeof(int) * nvmax * 4);
    if (!iwk) {
	fprintf(stderr, "gvert: failed to allocate memory for %d ints\n", nvmax * 4);
	return(-1);
    }
    wk = (long double *) malloc(sizeof(long double) * nvmax);
    if (!wk) {
	fprintf(stderr, "gvert: failed to allocate memory for %d long doubles\n", nvmax);
	return(-1);
    }

    /* fortran routine */
    gvert_(ve, angle, ipv, gp, ev, &nvmax, nv, &per, &nve, nev, nev0, poly->rp, poly->cm, &poly->np, &vcirc, tol, phi, iord, wk, iwk, &ldegen);

    /* number of vertices exceeds putative maximum */
    if (poly->np >= 5 && *nv > 6 * (poly->np - 2)) {
	msg("CONGRATULATIONS!  YOU HAVE DISCOVERED A POLYGON WITH 5 OR MORE CAPS\n");
	msg("(IT HAS %d CAPS) THAT HAS MORE THAN %d VERTICES (IT HAS %d VERTICES).\n", poly->np, 6 * (poly->np - 2), *nv);
	msg("(Either that or you have found a bug.)\n");
	msg("PLEASE EMAIL ME Andrew.Hamilton@colorado.edu THE GOOD NEWS,\n");
	msg("ALONG WITH A POLYGON FILE CONTAINING THE POLYGON THAT DID IT.\n");
	dump_poly(1, &poly);
	msg("AND THERE'S THE POLYGON FILE I'D LIKE YOU TO SEND.  THANKS!\n");
    }

    /* convert vertex/edge indices from fortran to c convention */
    for (iv = 0; iv < *nv && iv < nvmax; iv++) ipv[iv]--;

    /* free work arrays */
    free(iord);
    free(phi);
    free(iwk);
    free(wk);

    /* fatal intersection of boundaries */
    if (ldegen) return(1);

    return(0);
}
