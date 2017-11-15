/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "manglefn.h"

/* number of extra vertices to allocate, to allow for expansion */
#define DNV		4

/*------------------------------------------------------------------------------
  Point(s) somewhere in the middle of a polygon.
  The number of points is equal to the number of connected sequences of
  vertices, as computed by gverts, except that if there are zero connected
  sequences, then there is a single midpoint.  Zero connected sequences
  happens when the option to return midpoints for non-intersecting circles is
  turned off in gverts, and each boundary of the polygon is a disjoint circle.

  If the connected sequences of vertices delineate holes in the polygon,
  then the number of connected sequences of vertices may exceed
  the number of connected parts of the polygon.  In this case vmid supplies
  more midpoints than there are connected parts of the polygon.

  The intention is to aim for a point (or points) squarely in the middle of
  the polygon, but the algorithm, though mildly paranoid, is not cast iron.

  Designed to be called after gverts().

   Input: poly = pointer to polygon.
          mtol = initial angular tolerance in radians
		 within which to merge multiple intersections.
	  nv = number of vertices of polygon, as output by gverts().
	  nve = number of points per edge, as input to gverts();
		must be even and >= 2.
	  ve = points on edges of polygon, as output by gverts().
	  ipv = cap numbers of edges, as output by gverts().
	  ev = end indices of connected sequences, as output by gverts().
  Output: *nvm = number of middle points in *vm
	       = max(nev, 1).
	  *vm_p = pointer to unit vectors giving middle points;
		  the required memory is allocated.
		  If the algorithm failed to find a point inside the polygon,
		  then the vector for that point is {0.,0.,0.}.
  Return value: 0 if ok (even if no inside point was found);
		-1 if failed to allocate memory.
*/
int vmid(polygon *poly, long double mtol, int nv, int nve, vec ve[/*nv * nve*/], int ipv[/*nv*/], int ev[/*nv*/], int *nvm, vec **vm_p)
{
    const int do_vcirc = 1;
    static int nvmmax = 0;
    static vec *vm = 0x0;

    int i, ier, iev, ip, iv, ivl, ivm, ivu, neva, nev0, nva, scm;
    long double angle, cm, cmbest, s, tol;
    vec  rp, vc, vi;
    int *ipva, *gpa, *eva;
    long double *cmvmin, *cmvmax, *cmpmin, *cmpmax;
    vec *vmin, *vmax;

    /* make sure nve is even and >= 2 */
    if (nve < 2 || nve % 2 != 0) {
	fprintf(stderr, "vmid: nve = %d must be even and >= 2\n", nve);
	fprintf(stderr, "STOP\n");
	exit(1);
    }

    /* number of midpoints */
    for (ivm = 0; ev[ivm] != nv; ivm++);
    *nvm = ivm + 1;

    /* make sure that vm contains enough space */
    if (!vm || *nvm > nvmmax) {
	if (vm) free(vm);
	vm = (vec *) malloc(sizeof(vec) * (*nvm + DNV));
	if (!vm) {
	    fprintf(stderr, "vmid: failed to allocate memory for %d vecs\n", *nvm + DNV);
	    return(-1);
	}
	nvmmax = *nvm + DNV;
    }

    /* polygon is whole sphere: choose point to be north pole */
    if (poly->np == 0) {
	vm[0][0] = 0.;
	vm[0][1] = 0.;
	vm[0][2] = 1.;

    /* polygon has 1 cap: choose point at centre of circle */
    } else if (poly->np == 1) {
	scm = (poly->cm[0] >= 0.)? 1 : -1;
	for (i = 0; i < 3; i++) vm[0][i] = scm * poly->rp[0][i];

    /* polygon has >= 2 caps */
    } else {

	/* index of central point of edge */
	iev = nve / 2;

	/* no vertices implies each boundary is a disjoint circle */
	if (nv == 0) {
	    /* set `edge points' equal to centres of circles */
	    for (ip = 0; ip < poly->np; ip++) {	/* gverts allocated enough memory for this */
		scm = (poly->cm[ip] >= 0.)? 1 : -1;
		for (i = 0; i < 3; i++) ve[ip * nve + iev][i] = scm * poly->rp[ip][i];
	    }
	}

	/* do each connected boundary */
	for (ivm = 0; ivm < *nvm; ivm++) {

	    /* initialize the midpoint to zero */
	    for (i = 0; i < 3; i++) vm[ivm][i] = 0.;

	    /* vertices on boundary */
	    if (nv > 0) {
		ivl = (ivm == 0)? 0 : ev[ivm - 1];
		ivu = ev[ivm];
	    } else {		/* nv = 0 */
		ivl = 0;
		ivu = poly->np;	/* gverts allocated enough memory for this */
	    }

	    /* central vector is summed vector of points at centres of edges of connected boundary */
	    for (i = 0; i < 3; i++) vc[i] = 0.;
	    for (iv = ivl; iv < ivu; iv++) {
		for (i = 0; i < 3; i++) vc[i] += ve[iv * nve + iev][i];
	    }
	    /* length of central vector */
	    s = sqrtl(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);
	    /* s = 0 can happen, though hardly ever */
	    if (s == 0.) {
		ip = ipv[ivl];
		scm = (poly->cm[ip] >= 0.)? 1 : -1;
		for (i = 0; i < 3; i++) vc[i] = scm * poly->rp[ip][i];
	    /* normalize central vector to unit length */
	    } else {
		for (i = 0; i < 3; i++) vc[i] /= s;
	    }

	    /* try several possible central points, and choose the best */
	    cmbest = -1.;
	    for (iv = ivl - 1; iv < ivu; iv++) {
		/* try central vector vc to start with */
		if (iv == ivl - 1) {
		    /* is central vector inside polygon? */
		    if (gptin(poly, vc)) {
			for (i = 0; i < 3; i++) vi[i] = vc[i];
		    } else {
			continue;
		    }
		/* try point on great circle joining central vector vc to centre of each edge */
		} else {
		    /* great circle joining centre of edge of connected boundary to central vector vc */
		    rp_to_gc(ve[iv * nve + iev], vc, rp, &cm);
		    /* central point vi of that part of the great circle inside the polygon, one end of which is the centre of the edge of the connected boundary */
		    tol = mtol;
		    ier = gvphi(poly, rp, cm, ve[iv * nve + iev], &tol, &angle, vi);
		    if (ier == -1) return(-1);
		}
		/* distance of central point to edges */
		tol = mtol;
		ier = gvlims(poly, do_vcirc, &tol, vi, &nva, &vmin, &vmax, &cmvmin, &cmvmax, &cmpmin, &cmpmax, &ipva, &gpa, &neva, &nev0, &eva);
		if (ier == -1) return(-1);
		if (ier) continue;
		/* closest distance of central point to edges */
		cm = 2.;
		for (i = 0; i < nv; i++) {
		    if (cmvmin[i] < cm) cm = cmvmin[i];
		}
		/* best bet so far is the point farthest from any edge */
		if (cm > cmbest) {
		    cmbest = cm;
		    for (i = 0; i < 3; i++) {
			vm[ivm][i] = vi[i];
		    }
		}
	    }

	}

    }

    /* point vm_p at vm array */
    *vm_p = vm;

    return(0);
}

/*------------------------------------------------------------------------------
  Point(s) at the barycentre of the centres of the edges of a polygon.
  The point(s) may lie either inside or outside the polygon.
  The number of points is equal to the number of connected sequences of
  vertices, as computed by gverts, except that if there are zero connected
  sequences, then there is a single point.  Zero connected sequences
  happens when the option to return midpoints for non-intersecting circles is
  turned off in gverts, and each boundary of the polygon is a disjoint circle.

  Designed to be called after gverts().

   Input: poly = pointer to polygon.
	  nv = number of vertices of polygon, as output by gverts().
	  nve = number of points per edge, as input to gverts();
	  ve = points at centres of edges of polygon, as output by gverts().
	  ipv = cap numbers of edges, as output by gverts().
	  ev = end indices of connected sequences, as output by gverts().
  Output: *nvm = number of middle points in *vm
	       = max(nev, 1).
	  *vm = pointer to unit vectors giving barycentres;
		the required memory is allocated.
		If the algorithm failed to find a point inside the polygon,
		then the vector for that point is {0.,0.,0.}.
  Return value: 0 if ok (even if no inside point was found);
		-1 if failed to allocate memory.
*/
int vmidc(polygon *poly, int nv, int nve, vec ve[/*nv * nve*/], int ipv[/*nv*/], int ev[/*nv*/], int *nvm, vec **vm_p)
{
    static int nvmmax = 0;
    static vec *vm = 0x0;

    int i, iev, ip, iv, ivl, ivm, ivu;
    long double s;

    /* number of midpoints */
    for (ivm = 0; ev[ivm] != nv; ivm++);
    *nvm = ivm + 1;

    /* make sure that vm contains enough space */
    if (!vm || *nvm > nvmmax) {
	if (vm) free(vm);
	vm = (vec *) malloc(sizeof(vec) * (*nvm + DNV));
	if (!vm) {
	    fprintf(stderr, "vmidc: failed to allocate memory for %d vecs\n", *nvm + DNV);
	    return(-1);
	}
	nvmmax = *nvm + DNV;
    }

    /* polygon is whole sphere: choose point to be north pole */
    if (poly->np == 0) {
	vm[0][0] = 0.;
	vm[0][1] = 0.;
	vm[0][2] = 1.;

    /* polygon has 1 cap: choose point at centre of circle */
    } else if (poly->np == 1) {
	for (i = 0; i < 3; i++) vm[0][i] = poly->rp[0][i];

    /* no vertices implies each boundary is a disjoint circle */
    } else if (nv == 0) {
	/* choose point at centre of 1st circle */
	for (i = 0; i < 3; i++) vm[0][i] = poly->rp[0][i];

    /* polygon has >= 2 caps */
    } else {

	/* index of central point of edge */
	iev = nve / 2;

	/* do each connected sequence of vertices */
	for (ivm = 0; ivm < *nvm; ivm++) {

	    /* vertices in connected sequence */
	    ivl = (ivm == 0)? 0 : ev[ivm - 1];
	    ivu = ev[ivm];

	    /* central vector is summed vector of points at centres of edges */
	    for (i = 0; i < 3; i++) vm[ivm][i] = 0.;
	    for (iv = ivl; iv < ivu; iv++) {
		for (i = 0; i < 3; i++) {
		    vm[ivm][i] += ve[iv * nve + iev][i];
		}
	    }
	    /* length of central vector */
	    s = sqrtl(vm[ivm][0]*vm[ivm][0] + vm[ivm][1]*vm[ivm][1] + vm[ivm][2]*vm[ivm][2]);
	    /* s = 0 can happen, though hardly ever */
	    if (s == 0.) {
		ip = ipv[ivl];
		for (i = 0; i < 3; i++) vm[ivm][i] = poly->rp[ip][i];
	    /* normalize central vector to unit length */
	    } else {
		for (i = 0; i < 3; i++) vm[ivm][i] /= s;
	    }

	}

    }

    /* point vm_p at vm array */
    *vm_p = vm;

    return(0);
}
