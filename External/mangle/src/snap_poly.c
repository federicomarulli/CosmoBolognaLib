/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Make almost coincident caps of polygons coincide.

  Input:  fmt = pointer to format structure.
	  npoly = number of polygons to snap.
	  *poly[npoly] = array of npoly pointers to polygon structures.
	  selfsnap = 0 to snap edges of all polygons against each other,
		     1 to snap edges of polygons only against edges of same polygon.
	  axtol = angle in radians [actually 2 sinl(angle/2)]:
		  if angle twixt polar axes of caps <= axtol,
		  then make axis of poly2 cap
		  exactly parallel to axis of poly1 cap.
	  btol = angle in radians:
		 if two axes of caps of poly1 and poly2 are parallel,
		 and if angle between latitudes of caps <= btol,
		 then make latitude of poly2 cap
		 exactly equal to latitude of poly1 cap.
	  thtol = edge tolerance in radians.
	  ytol = edge to length tolerance;
		 if the two vertices and centre point of an edge of poly2 are
		 all closer to a boundary of poly1 than the lesser of
		 (1) thtol, and
		 (2) ytol times the length of the edge,
		 and if in addition at least one of the three points lies
		 inside poly1 (sans said boundary),
		 then make boundary of the poly2 cap equal to that of poly1.
	  mtol = initial tolerance angle for multiple intersections in radians.
	  warnmax = number of times to advise about individual polygon edges being snapped.
  Output: adjusted caps of poly2 (i.e. poly2->rp, poly2->cm).
	  snapped_poly = array of 0 or 1 flagging which polygons were snapped.
			 set to 0x0 on input to ignore.
  Return value: number of caps adjusted.
*/
int snap_polys(format *fmt, int npoly, polygon *poly[/*npoly*/], int selfsnap, long double axtol, long double btol, long double thtol, long double ytol, long double mtol, int warnmax, char snapped_poly[/*npoly*/])
{
    int dnadj, dnadjo, i, j, nadj, pass, snapped, stuck, warn;

    /* initialize snapped polygon flag to zero */
    if (snapped_poly) {
	for (i = 0; i < npoly; i++) snapped_poly[i] = 0;
    }

    nadj = 0;

    /* snap repeatedly, until no more caps snap together */
    pass = 0;
    stuck = 0;
    dnadj = 0;
    do {
	/* snap caps of each pair of polygons in turn, including self-pairs */
	pass++;
	dnadjo = dnadj;
	dnadj = 0;
	warn = 0;
	if (axtol >= 0. || btol >= 0.) {
	    for (i = 0; i < npoly; i++) {
		for (j = i; ((selfsnap)? j == i : j < npoly); j++) {
		    snapped = snap_poly(poly[i], poly[j], axtol, btol);
		    if(snapped==-1){
		      fprintf(stderr, "snap_polys: error in snap_poly for polys %d and %d in pixel %d\n",i,j,poly[i]->pixel);
		      return(-1);
		    }
		    
		    if (snapped) {
			if (warnmax < 0) {
			  if (warn == 0)
			    msg("snap_polys stage 1 pass %d: caps of the following polygons were snapped together:\n", pass);
			  if (warn < warnmax) {
			    if (selfsnap) {
			      msg(" %d", (fmt->newid == 'o')? poly[i]->id : i);
			    } else {
			      msg(" (%d %d)", (fmt->newid == 'o')? poly[i]->id : i, (fmt->newid == 'o')? poly[j]->id : j);
			    }
			  } else if (warn == warnmax) {
			    msg(" ... more\n");
			  }
			}
			if (snapped_poly) {
			    snapped_poly[i] = 1;
			    snapped_poly[j] = 1;
			}
			dnadj += snapped;
			warn++;
		    }
		}
	    }
	}
	if (warnmax > 0 && warn > 0 && warn <= warnmax) msg("\n");
	nadj += dnadj;
	if ((nadj > 0 || !selfsnap) && warnmax) msg("snap_polys stage 1 (axes, latitudes) pass %d: %d caps adjusted\n", pass, dnadj);
	/* avoid infinite loop */
	if (pass > 1 && dnadj >= dnadjo) stuck++;
    } while (dnadj && stuck < 2);
    if (dnadj) {
      if(poly[0]->pixel > 0){
	fprintf(stderr, "snap_polys stage 1: stuck in a loop in pixel %d. continuing ...\n",poly[0]->pixel);
      }
      else{
	fprintf(stderr, "snap_polys stage 1: seem to be stuck in a loop ... exit\n");
      }
    }

    /* trim polygons */
    for (i = 0; i < npoly; i++) {
	trim_poly(poly[i]);
    }

    /* snap repeatedly, until no more caps snap together */
    pass = 0;
    stuck = 0;
    dnadj = 0;
    do {
	/* snap edges of each polygon to caps of each polygon in turn */
	pass++;
	dnadjo = dnadj;
	dnadj = 0;
	warn = 0;
	if (thtol >= 0. && ytol >= 0.) {
	    for (i = 0; i < npoly; i++) {
		for (j = ((selfsnap)? i : 0); ((selfsnap)? j == i : j < npoly); j++) {
		    snapped = snap_polyth(poly[i], poly[j], thtol, ytol, mtol);
		    if(snapped==-1){
		      fprintf(stderr, "snap_polys: error in snap_poly for polys %d and %d in pixel %d\n",i,j,poly[i]->pixel);
		      return(-1);
		    }
		    if (snapped) {
			if (warnmax > 0) {
			    if (warn == 0)
				msg("snap_polys stage 2 pass %d: caps of the following polygons were snapped together:\n", pass);
			    if (warn < warnmax) {
				if (selfsnap) {
				    msg(" %d", (fmt->newid == 'o')? poly[i]->id : i);
				} else {
				    msg(" (%d %d)", (fmt->newid == 'o')? poly[i]->id : i, (fmt->newid == 'o')? poly[j]->id : j);
				}
			    } else if (warn == warnmax) {
				msg(" ... more\n");
			    }
			}
			if (snapped_poly) {
			    snapped_poly[i] = 1;
			    snapped_poly[j] = 1;
			}
			dnadj += snapped;
			warn++;
		    }
		}
	    }
	}
	if (warnmax > 0 && warn > 0 && warn <= warnmax) msg("\n");
	nadj += dnadj;
	if ((nadj > 0 || !selfsnap) && warnmax) msg("snap_polys stage 2 (edges) pass %d: %d caps adjusted\n", pass, dnadj);
	/* avoid infinite loop */
	if (pass > 1 && dnadj >= dnadjo) stuck++;
    } while (dnadj && stuck < 2);
    if (dnadj) {
      if(poly[0]->pixel > 0){
	fprintf(stderr, "snap_polys stage 2: stuck in a loop in pixel %d. continuing ...\n",poly[0]->pixel);
      }
      else{
	fprintf(stderr, "snap_polys stage 2: seem to be stuck in a loop ... exit\n");
      }
    }
    return(nadj);
}

/*------------------------------------------------------------------------------
  Make almost coincident caps of 2 polygons coincide.
  Caps of poly2 are adjusted to equal those of poly1.

   Input: poly1, poly2 = pointers to polygon structures.
	  axtol = angle in radians [actually 2 sinl(angle/2)]:
		  if angle twixt polar axes of caps <= axtol,
		  then make axis of poly2 cap
		  exactly parallel to axis of poly1 cap.
	  btol =  angle in radians:
		  if two axes of caps of poly1 and poly2 are parallel,
		  and if angle between latitudes of caps <= btol,
		  then make latitude of poly2 cap
		  exactly equal to latitude of poly1 cap.
  Output: adjusted caps of poly2 (i.e. poly2->rp, poly2->cm).
  Return value: number of caps adjusted.
*/
int snap_poly(polygon *poly1, polygon *poly2, long double axtol, long double btol)
{
    int adjusted, ip, ip1, ip2, nadj, sp;
    long double cm, dl, drp, dx, dy, dz;

    nadj = 0;
    for (ip1 = 0; ip1 < poly1->np; ip1++) {	/* for each cap of poly1 ... */
	/* superfluous cap */
	if (poly1->cm[ip1] == 0. || fabsl(poly1->cm[ip1]) >= 2.) continue;
	for (ip2 = 0; ip2 < poly2->np; ip2++) {	/* ... and each cap of poly2 */
	    /* superfluous cap */
	    if (poly2->cm[ip2] == 0. || fabsl(poly2->cm[ip2]) >= 2.) continue;
	    for (ip = 0; ip < 2; ip++) {	/* check rp2 = +- rp1 */
		adjusted = 0;
		sp = (ip == 0)? 1 : -1;
		/* [2 sinl(alpha/2)]^2, where alpha is angle twixt axes */
		dx = poly2->rp[ip2][0] - sp * poly1->rp[ip1][0];
		dy = poly2->rp[ip2][1] - sp * poly1->rp[ip1][1];
		dz = poly2->rp[ip2][2] - sp * poly1->rp[ip1][2];
		drp = sqrtl(dx * dx + dy * dy + dz * dz);
		if (drp <= axtol) {		/* axes are nearly parallel */
		    if (!(poly2->rp[ip2][0] == sp * poly1->rp[ip1][0]
		      && poly2->rp[ip2][1] == sp * poly1->rp[ip1][1]
		      && poly2->rp[ip2][2] == sp * poly1->rp[ip1][2])) {
			/* make axis of poly2 cap exactly parallel to poly1
			   (made exactly equal below if caps nearly coincide) */
			poly2->rp[ip2][0] = sp * poly1->rp[ip1][0];
			poly2->rp[ip2][1] = sp * poly1->rp[ip1][1];
			poly2->rp[ip2][2] = sp * poly1->rp[ip1][2];
			adjusted = 1;
		    }
		    /* angle between latitudes of caps */
		    if (sp == 1) {		/* axes are aligned */
			dl = 2. * (asinl(sqrtl(fabsl(poly2->cm[ip2]) / 2.))
			    - asinl(sqrtl(fabsl(poly1->cm[ip1]) / 2.)));
		    } else {			/* axes are anti-aligned */
			dl = 2. * (asinl(sqrtl((2. - fabsl(poly2->cm[ip2])) / 2.))
			    - asinl(sqrtl(fabsl(poly1->cm[ip1]) / 2.)));
		    }
		    if (fabsl(dl) <= btol) {	/* caps nearly coincide */
			if (sp == -1) {
			    /* reflect axis of poly2 cap */
			    poly2->rp[ip2][0] = - poly2->rp[ip2][0];
			    poly2->rp[ip2][1] = - poly2->rp[ip2][1];
			    poly2->rp[ip2][2] = - poly2->rp[ip2][2];
			    adjusted = 1;
			}
			/* set latitude of poly2 cap equal to poly1 */
			cm = (poly2->cm[ip2] >= 0.)?
			    sp * fabsl(poly1->cm[ip1]):
			    - sp * fabsl(poly1->cm[ip1]);
			if (poly2->cm[ip2] != cm) {
			    poly2->cm[ip2] = cm;
			    adjusted = 1;
			}
		    }
		    if (adjusted) {
			nadj++;
			/* no need to test other direction */
			break;
		    }
		}
	    }
	}
    }
    return(nadj);
}
/*------------------------------------------------------------------------------
  Snap edge of poly2 to cap boundary of poly1.
  Caps of poly2 are adjusted to equal those of poly1.

  Input:  poly1, poly2 = pointers to polygon structures.
	  thtol = edge tolerance in radians.
	  ytol = edge to length tolerance;
		 if the two vertices and centre point of an edge of poly2 are
		 all closer to a boundary of poly1 than the lesser of
		 (1) thtol, and
		 (2) ytol times the length of the edge,
		 and if in addition at least one of the three points lies
		 inside poly1 (sans said boundary),
		 then make boundary of the poly2 cap equal to that of poly1.
	  mtol = initial tolerance angle for multiple intersections in radians.
  Output: adjusted caps of poly2 (i.e. poly2->rp, poly2->cm).
  Return value: number of caps adjusted,
		or -1 if error occurred.
*/
int snap_polyth(polygon *poly1, polygon *poly2, long double thtol, long double ytol, long double mtol)
{
    const int per = 0;
    const int nve = 2;

    int adjusted, do_vcirc, i, ier, in, ip1, ip2, iv, ivp, nadj, nev, nev0, nv;
    int *ipv, *gp, *ev;
    long double cm, cm1, dth, dthmax, sp, tol;
    long double *angle;
    vec *v, *ve;

    /* vertices and centres of edges of poly2 */
    do_vcirc = 0;
    tol = mtol;
    ier = gverts(poly2, do_vcirc, &tol, per, nve, &nv, &ve, &angle, &ipv, &gp, &nev, &nev0, &ev);
    if (ier != 0) return(-1);

    /* convert angle of each edge to scalar length angle * sinl(theta) */
    for (iv = 0; iv < nv; iv++) {
	ip2 = ipv[iv];
	cm = fabsl(poly2->cm[ip2]);
	angle[iv] = angle[iv] * sqrtl(cm * (2. - cm));
    }

    nadj = 0;
    /* for each edge of poly2 ... */
    for (iv = 0; iv < nv; iv++) {
	ivp = (iv + 1) % nv;
	ip2 = ipv[iv];

	/* ... and each axis of poly1 */
	for (ip1 = 0; ip1 < poly1->np; ip1++) {
	    adjusted = 0;

	    /* distance from edge of poly2 to cap of poly1 */
	    cm1 = poly1->cm[ip1];
            poly1->cm[ip1] = 2.;        /* suppress cap of poly1 */
	    in = 0;
	    dthmax = 0.;
	    for (i = 0; i < 3; i++) {
		/* vertex, centre point, vertex of edge of poly2 */
		v = &ve[(iv * nve + i) % (nv * nve)];
		in |= gptin(poly1, *v);	/* in if any one point is in */
		cm = cmij(*v, poly1->rp[ip1]);
		dth = 2. * (sqrtl(cm/2.) - sqrtl(fabsl(cm1/2.)));
		dth = fabsl(dth);	/* angle from point to cap of poly1 */
		if (dth > dthmax) dthmax = dth;
	    }
            poly1->cm[ip1] = cm1;       /* restore cap of poly1 */

	    /* three points of poly2 edge are all close to boundary of poly1 */
	    if (in && dthmax <= thtol && dthmax <= ytol * angle[iv]) {
		sp = poly1->rp[ip1][0] * poly2->rp[ip2][0] + poly1->rp[ip1][1] * poly2->rp[ip2][1] + poly1->rp[ip1][2] * poly2->rp[ip2][2];
		sp = (sp >= 0.)? 1. : -1.;
		if (!(poly2->rp[ip2][0] == poly1->rp[ip1][0]
		  && poly2->rp[ip2][1] == poly1->rp[ip1][1]
		  && poly2->rp[ip2][2] == poly1->rp[ip1][2])) {
		    /* make axis of poly2 cap exactly equal to that of poly1 */
		    poly2->rp[ip2][0] = poly1->rp[ip1][0];
		    poly2->rp[ip2][1] = poly1->rp[ip1][1];
		    poly2->rp[ip2][2] = poly1->rp[ip1][2];
		    adjusted = 1;
		}
		/* set latitude of poly2 cap equal to that of poly1 */
		cm = (poly2->cm[ip2] >= 0.)?
		    sp * fabsl(poly1->cm[ip1]):
		    - sp * fabsl(poly1->cm[ip1]);
		if (poly2->cm[ip2] != cm) {
		    poly2->cm[ip2] = cm;
		    adjusted = 1;
		}
		if (adjusted) nadj++;
	    }

	}

    }

    /* trim adjusted polygon */
    if (nadj > 0) trim_poly(poly2);

    return(nadj);
}
