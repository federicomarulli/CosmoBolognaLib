/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include "manglefn.h"

/* initial angular tolerance within which to merge multiple intersections */
extern long double mtol;

/*------------------------------------------------------------------------------
  Convert vertices structure to polygon.
  Assume that vertices are joined by great circles.

  Polygon poly should contain enough room to contain the vertices,
  namely poly->npmax >= vert->nv.
*/
void vert_to_poly(vertices *vert, polygon *poly)
{
    int iv;

    /* number of boundaries of polygon equals number of vertices */
    poly->np = vert->nv;

    /* convert each pair of adjacent vertices to a great circle */
    for (iv = 0; iv < vert->nv; iv++) {
	azel_to_gc(&vert->v[iv], &vert->v[(iv+1) % vert->nv],
	    poly->rp[iv], &poly->cm[iv]);
    }
}

/*------------------------------------------------------------------------------
  Convert edges, stored in vertices structure, to polygon.

   Input: vert = pointer to vertices structure.
	  nve = number of points per edge.
	  ev = cumulative number of points on each connected boundary.
  Output: poly = pointer to polygon structure.
*/
void edge_to_poly(vertices *vert, int nve, int *ev, polygon *poly)
{
    int dev, i, evo, iedge, iv, iv0, iv1, iv2, jv, nedge;

    /* number of edges */
    nedge = vert->nv / nve;

    /* number of boundaries of polygon equals number of edges */
    poly->np = nedge;

    evo = 0;
    iv = 0;
    iedge = 0;
    /* do each connected boundary */
    for (jv = 0; evo < vert->nv; jv++) {
	dev = ev[jv] - evo;
	/* points on connected boundary */
	while (iv < ev[jv]) {
	    i = iv - evo;
	    iv0 = evo + i;
	    iv1 = evo + (i + 1) % dev;
	    iv2 = evo + (i + 2) % dev;
	    if (nve == 1) {
		/* convert pair of adjacent vertices to a great circle */
		azel_to_gc(&vert->v[iv0], &vert->v[iv1],
		    poly->rp[iedge], &poly->cm[iedge]);
	    } else {
		/* convert triple of vertices to circle */
		edge_to_rpcm(&vert->v[iv0], &vert->v[iv1], &vert->v[iv2],
		    poly->rp[iedge], &poly->cm[iedge]);
	    }
	    iv += nve;
	    iedge++;
	}
	evo = ev[jv];
    }
}

/*------------------------------------------------------------------------------
  Convert rectangle to polygon.
  Input: angle [4] = list of angles defining the rectangle, 
         in the order azmin, azmax, elmin, elmax
  Output: poly = pointer to polygon structure
  Polygon poly should contain enough room for the 4 caps of the rectangle
*/
void rect_to_poly(long double angle[4], polygon *poly){
#define ROUND		1.e-5
  long double daz;
  int i;
  int ip =0;
  /*
  angle[0]=azmin;
  angle[1]=azmax;
  angle[2]=elmin;
  angle[3]=elmax;
  */ 
  /* azimuthal extent */
  daz = angle[1] - angle[0];
  /* skip azimuth if it goes full circle */
  if (fabsl(daz - TWOPI) > ROUND) {
    /* phase daz into interval [0, 2 pi) */
    daz = daz - floorl(daz / TWOPI) * TWOPI;
    /* azimuth covers half circle: one cap will do */
    if (fabsl(daz - PI) <= ROUND) {
      i = 0;
      az_to_rpcm(angle[i], i % 2, poly->rp[ip], &poly->cm[ip]);
      ip++;
       /* need two caps */
    } else {
      /* convert azmin, azmax to rp, cm */
      for (i = 0; i < 2; i++) {
	az_to_rpcm(angle[i], i % 2, poly->rp[ip], &poly->cm[ip]);
	ip++;
	
      }
      /* check azimuthal extent is in interval [0, pi] */
      if (daz > PI) {
	fprintf(stderr, " warning:");
	fprintf(stderr, " rectangle has azimuthal extent %.16Lg deg > 180 deg\n",
		places(daz * 180./PI, 14));
      }
    }
  }
  
  /* convert elmin, elmax to rp, cm */
  for (i = 2; i < 4; i++) {    
    /*if one of the el limits is a pole, make a tiny circle around the pole.*/
    /*   
    if(fabsl(angle[i] - ((i==2) ? -PIBYTWO : PIBYTWO)) < ROUND){
      angle[i]+= (i==2)? (5*ROUND) : (-5*ROUND);
    }
    */
    if ((i == 2)? (angle[i] > - PIBYTWO+ROUND) : (angle[i] < PIBYTWO-ROUND)) {
      el_to_rpcm(angle[i], i % 2, poly->rp[ip], &poly->cm[ip]);
       ip++; 
    } 
  } 
  poly->np = ip;
}

/*------------------------------------------------------------------------------
  Convert unit vectors to vertices structure.

   Input: rp = array of nv unit vectors.
	  nv = number of unit vectors.
  Output: pointer to a vertices structure.
*/
void rps_to_vert(int nv, vec rp[/*nv*/], vertices *vert)
{
    int iv;

    /* number of vertices */
    vert->nv = nv;

    /* convert vectors to vertices */
    for (iv = 0; iv < nv; iv++) {
	rp_to_azel(rp[iv], &vert->v[iv]);
    }

    /* phase each vertex to the previous */
    for (iv = 1; iv < nv; iv++) {
	vert->v[iv].az -= rint((vert->v[iv].az - vert->v[iv-1].az) / TWOPI) * TWOPI;
    }
}

/*------------------------------------------------------------------------------
  Convert unit vector to az-el.
*/
void rp_to_azel(vec rp, azel *v)
{
    v->az = atan2l(rp[1], rp[0]);
    v->el = atan2l(rp[2], sqrtl(rp[0] * rp[0] + rp[1] * rp[1]));
}

/*------------------------------------------------------------------------------
  Convert vertex to unit vector.
*/
void azel_to_rp(azel *v, vec rp)
{
    rp[0] = cosl(v->el) * cosl(v->az);
    rp[1] = cosl(v->el) * sinl(v->az);
    rp[2] = sinl(v->el);
}

/*------------------------------------------------------------------------------
  Determine rp, cm for great circle passing through two az-el vertices.
  The great circle goes right-handedly from v0 to v1.
*/
void azel_to_gc(azel *v0, azel *v1, vec rp, long double *cm)
{
    azel *v;
    int iv;
    vec rpv[2];

    for (iv = 0; iv < 2; iv++) {
	v = (iv == 0)? v0: v1;
	rpv[iv][0] = cosl(v->el) * cosl(v->az);
	rpv[iv][1] = cosl(v->el) * sinl(v->az);
	rpv[iv][2] = sinl(v->el);
    }

    rp_to_gc(rpv[0], rpv[1], rp, cm);
}

/*------------------------------------------------------------------------------
  Determine rp, cm for great circle passing through two unit vectors.
  The great circle goes right-handedly from rp0 to rp1.
  If the two unit vectors coincide, suppress the boundary,
  on the assumption that coincident points are redundant
  (for example, a mask-maker may specify a triangle with 4 vertices,
  with 2 vertices being coincident).
*/
void rp_to_gc(vec rp0, vec rp1, vec rp, long double *cm)
{
    int i;
    long double rpa;

    /* cofactors */
    rp[0] = rp0[1]*rp1[2] - rp1[1]*rp0[2];
    rp[1] = rp1[0]*rp0[2] - rp0[0]*rp1[2];
    rp[2] = rp0[0]*rp1[1] - rp1[0]*rp0[1];

    rpa = sqrtl(rp[0]*rp[0] + rp[1]*rp[1] + rp[2]*rp[2]);

    /* indeterminate solution (rp0 and rp1 are same or antipodeal points) */
    if (rpa == 0.) {
    /* suppress boundary, on assumption that coincident points are redundant */
	*cm = 2.;
    /* normalize rp to 1 */
    } else {
	for (i = 0; i < 3; i++) rp[i] /= rpa;
	/* cm = 1 for great circle */
	*cm = 1.;
    }
}

/*------------------------------------------------------------------------------
  Determine rp, cm for circle passing through three az-el points.
  The circle goes right-handedly from v0 to v1 to v2.
*/
void edge_to_rpcm(azel *v0, azel *v1, azel *v2, vec rp, long double *cm)
{
    azel *v=0x0;
    int iv;
    vec rpv[3];

    for (iv = 0; iv < 3; iv++) {
	switch (iv) {
	case 0:	v = v0;	break;
	case 1:	v = v1;	break;
	case 2:	v = v2;	break;
	}
	rpv[iv][0] = cosl(v->el) * cosl(v->az);
	rpv[iv][1] = cosl(v->el) * sinl(v->az);
	rpv[iv][2] = sinl(v->el);
    }

    rp_to_rpcm(rpv[0], rpv[1], rpv[2], rp, cm);
}

/*------------------------------------------------------------------------------
  Determine rp, cm for circle passing through three unit vectors.
  The circle goes right-handedly from rp0 to rp1 to rp2.
  If two of the unit vectors coincide, join with a great circle.
  If three of the unit vectors coincide, suppress the boundary.
*/
void rp_to_rpcm(vec rp0, vec rp1, vec rp2, vec rp, long double *cm)
{
    int coincide, i, j;
    long double det, rpa;
    long double *rpi, *rpj;
    rpi=0x0;
    rpj=0x0;

    /* check whether any two of the three unit vectors coincide */
    coincide = 0;
    for (j = 0; j < 3; j++) {
	switch (j) {
	case 0:	rpj = rp0;	break;
	case 1:	rpj = rp1;	break;
	case 2:	rpj = rp2;	break;
	}
	for (i = 0; i < j; i++) {
	    switch (i) {
	    case 0:	rpi = rp0;	break;
	    case 1:	rpi = rp1;	break;
	    case 2:	rpi = rp2;	break;
	    }
	    /* vector coincide */
	    if (rpi[0] == rpj[0] && rpi[1] == rpj[1] && rpi[2] == rpj[2]) {
		coincide = 1;
		/* set rpj equal to the third vector */
		if (i == 0 && j == 1) {
		    rpj = rp2;
		} else if (i == 0 && j == 2) {
		    rpj = rp1;
		} else if (i == 1 && j == 2) {
		    rpj = rp0;
		}
	    }
	    if (coincide) break;
	}
	if (coincide) break;
    }

    /* if any two vectors coincide, join with great circle */
    if (coincide) {
	rp_to_gc(rpi, rpj, rp, cm);

    /* non-coincident vectors */
    } else {
	/* cofactors, arranged to reduce roundoff */
	rp[0] = rp0[1] * (rp1[2] - rp2[2])
	      + rp1[1] * (rp2[2] - rp0[2])
	      + rp2[1] * (rp0[2] - rp1[2]);
	rp[1] = rp0[0] * (rp2[2] - rp1[2])
	      + rp1[0] * (rp0[2] - rp2[2])
	      + rp2[0] * (rp1[2] - rp0[2]);
	rp[2] = (rp0[0] * rp1[1] - rp1[0] * rp0[1])
	      + (rp1[0] * rp2[1] - rp2[0] * rp1[1])
	      + (rp2[0] * rp0[1] - rp0[0] * rp2[1]);

	/* |rp| */
	rpa = sqrtl(rp[0]*rp[0] + rp[1]*rp[1] + rp[2]*rp[2]);

	/* indeterminate solution (2 of 3 unit vectors coincide: shouldn't happen) */
	if (rpa == 0.) {
	    /* suppress boundary, not knowing what to do with it */
	    *cm = 2.;
	/* normal solution */
	} else {
	    /* determinant */
	    det = (rp0[0] * rp1[1] - rp1[0] * rp0[1]) * rp2[2]
		+ (rp1[0] * rp2[1] - rp2[0] * rp1[1]) * rp0[2]
		+ (rp2[0] * rp0[1] - rp0[0] * rp2[1]) * rp1[2];
	    if (det >= 0.) {
		/* cosl(th) = det / |rp| */
		*cm = 1. - det / rpa;
	    } else {
		/* cosl(th) = - det / |rp| */
		rpa = - rpa;
		*cm = - (1. - det / rpa);
	    }
	    /* normalize rp to 1 */
	    for (i = 0; i < 3; i++) rp[i] /= rpa;
	}
    }
}

/*------------------------------------------------------------------------------
  Convert circle to rp, cm.

   Input: angle = (azimuth, elevation, radius) in radians.
  Output: rp, cm as used by garea, gspher et al.
*/
void circ_to_rpcm(long double angle[3], vec rp, long double *cm)
{
    long double s;

    /* Cartesian coordinates of azimuth, elevation */
    rp[0] = cosl(angle[1]) * cosl(angle[0]);
    rp[1] = cosl(angle[1]) * sinl(angle[0]);
    rp[2] = sinl(angle[1]);
    /* 1 - cosl(radius) = 2 sin^2(radius/2) */
    s = sinl(angle[2] / 2.);
    *cm = s * s * 2.;
    *cm = (angle[2] >= 0.)? *cm : -*cm;
}

/*------------------------------------------------------------------------------
  Convert rp, cm to circle.

   Input: rp, cm as used by garea, gspher et al.
  Output: angle = (azimuth, elevation, radius) in radians.
*/
void rpcm_to_circ(vec rp, long double *cm, long double angle[3])
{
    long double s;

    angle[0] = atan2l(rp[1], rp[0]);
    angle[1] = atan2l(rp[2], sqrtl(rp[0] * rp[0] + rp[1] * rp[1]));
    s = sqrtl(fabsl(*cm) / 2.);
    if (s > 1.) s = 1.;
    angle[2] = 2. * asinl(s);
    angle[2] = (*cm >= 0.)? angle[2] : -angle[2];
}

/*------------------------------------------------------------------------------
  Convert line of constant azimuth to rp, cm.

   Input: az = azimuth in radians.
	  m = 0 for minimum elevation;
	    = 1 for maximum elevation.
  Output: rp, cm as used by garea, gspher et al.
*/
void az_to_rpcm(long double az, int m, vec rp, long double *cm)
{
    /* axis along equator */
    rp[0] = - sinl(az);
    rp[1] = cosl(az);
    rp[2] = 0.;
    /* 1 - cosl(th) = 1 for great circle */
    *cm = 1.;
    /* min, max */
    *cm = (m == 0)? *cm : - *cm;
}

/*------------------------------------------------------------------------------
  Convert line of constant elevation to rp, cm.

   Input: el = elevation in radians.
	  m = 0 for minimum elevation;
	    = 1 for maximum elevation.
  Output: rp, cm as used by garea, gspher et al.
*/
void el_to_rpcm(long double el, int m, vec rp, long double *cm)
{
    /* north pole */
    rp[0] = 0.;
    rp[1] = 0.;
    rp[2] = 1.;
    /* 1 - cosl(th) = 1 - sinl(el) */
    *cm = 1 - sinl(el);
    /* min, max */
    *cm = (m == 0)? *cm : - *cm;
}

/*------------------------------------------------------------------------------
   theta_ij = angle in radians between two unit vectors.
*/
long double thij(vec rpi, vec rpj)
{
    long double cm, th;

    cm = cmij(rpi, rpj);
    th = 2. * asinl(cm / 2.);

    return(th);
}

/*------------------------------------------------------------------------------
   1-cosl(theta_ij) = 2 sin^2(theta_ij/2) between two unit vectors.
*/
long double cmij(vec rpi, vec rpj)
{
    long double cm, dx, dy, dz;

    dx = rpi[0] - rpj[0];
    dy = rpi[1] - rpj[1];
    dz = rpi[2] - rpj[2];
    cm = (dx*dx + dy*dy + dz*dz) / 2.;
    if (cm > 2.) cm = 2.;

    return(cm);
}

/*------------------------------------------------------------------------------
  Determine whether a polygon is a rectangle,
  and if so return <azmin> <azmax> <elmin> <elmax>.

   Input: poly = pointer to polygon.
  Output: *azmin, *azmax, *elmin, *elmax;
	  if *azmin >= *azmax, or *elmin >= *elmax, the rectangle is empty;
	  if *azmin = *azmax + 2*pi, there is no azimuthal constraint.
  Return value: 0 if polygon is not a rectangle,
		1 if polygon is a rectangle,
		2 if polygon is a rectangle with superfluous boundaries.
*/
int poly_to_rect(polygon *poly, long double *azmin, long double *azmax, long double *elmin, long double *elmax)
{
    int iaz, ielmin, ielmax, ip;
    long double az, el;

    *azmin = -TWOPI;
    *azmax = TWOPI;
    *elmin = -PIBYTWO;
    *elmax = PIBYTWO;
    iaz = 0;
    ielmin = 0;
    ielmax = 0;
    for (ip = 0; ip < poly->np; ip++) {
	/* null polygon */
	if (poly->cm[ip] == 0. || poly->cm[ip] <= -2.) return(0);
	/* skip superfluous cap */
	if (poly->cm[ip] >= 2.) continue;
	/* line of constant azimuth */
	if (poly->rp[ip][2] == 0.
	    && (poly->cm[ip] == 1. || poly->cm[ip] == -1.)) {
	    /* az is in (-pi,pi] */
	    az = atan2l(-poly->rp[ip][0], poly->rp[ip][1]);
	    if (poly->cm[ip] == -1.) az += (az <= 0.)? PI : -PI;
	    /* shift azmin, azmax to phase of az, az+pi */
	    if (az >= *azmax) {
		while (az >= *azmax) {
		    *azmin += TWOPI;
		    *azmax += TWOPI;
		}
	    } else if (az + PI <= *azmin) {
		while (az + PI <= *azmin) {
		    *azmin -= TWOPI;
		    *azmax -= TWOPI;
		}
	    }
	    /* revise azmin, azmax */
	    if (az > *azmin) *azmin = az;
	    if (az + PI < *azmax) *azmax = az + PI;
	    iaz++;
	/* line of constant elevation */
	} else if (poly->rp[ip][0] == 0. && poly->rp[ip][1] == 0.) {
	    el = asinl(1. - fabsl(poly->cm[ip]));
	    if (poly->rp[ip][2] < 0.) el = -el;
	    if ((poly->rp[ip][2] > 0. && poly->cm[ip] >= 0.)
		|| (poly->rp[ip][2] < 0. && poly->cm[ip] < 0.)) {
		if (el > *elmin) *elmin = el;
		ielmin++;
	    } else {
		if (el < *elmax) *elmax = el;
		ielmax++;
	    }
	} else {
	    return(0);
	}
    }

    /* azimuthal constraints include nothing */
    if (*azmin > *azmax) {
	*azmin = *azmax = 0.;
    /* azimuthal constraints include everything */
    } else if (*azmin + TWOPI <= *azmax) {
	*azmin = -PI;
	*azmax = PI;
    /* phase azmin to (-pi,pi] */
    } else if (*azmin > PI) {
        *azmin -= TWOPI;
        *azmax -= TWOPI;
    }

    /* should not happen */
    if (iaz > 2 || ielmin > 1 || ielmax > 1) {
	msg("poly_to_rect: polygon contains");
	if (iaz > 2) msg(" %d az boundaries", iaz);
	if (ielmin > 1) msg(" %d elmin boundaries", ielmin);
	if (ielmax > 1) msg(" %d elmax boundaries", ielmax);
	msg("\n");
	return(2);
    }

    return(1);
}

/*------------------------------------------------------------------------------
  Determine whether the first vertex in vert
  is closer to the nearest vertex or to the nearest antivertex of poly.

  Return value: 0 or 1 as first vertex of vert is closer to the nearest
		       vertex or antivertex of poly;
		-1 if error.
*/
int antivert(vertices *vert, polygon *poly)
{
    const int do_vcirc = 0, nve = 1, per = 0;
    int anti, ier, iv, nev, nev0, nv;
    int *ipv, *gp, *ev;
    long double cm, cmmax, cmmin, tol;
    long double *angle;
    vec rp;
    vec *ve;

    /* vert has no vertices */
    if (vert->nv == 0) return(0);

    /* vertices of polygon */
    tol = mtol;
    ier = gverts(poly, do_vcirc, &tol, per, nve, &nv, &ve, &angle, &ipv, &gp, &nev, &nev0, &ev);

    /* error */
    if (ier) return(-1);

    /* poly has less than 3 vertices */
    if (nv < 3) return(0);

    /* convert first vertex to unit vector */
    azel_to_rp(&vert->v[0], rp);

    cmmin = 2.;
    cmmax = 0.;
    for (iv = 0; iv < nv; iv++) {
	cm = cmij(rp, ve[iv]);
	if (cm < cmmin) cmmin = cm;
	if (cm > cmmax) cmmax = cm;
    }
    anti = ((cmmin + cmmax <= 2.)? 0 : 1);

    return(anti);
}
