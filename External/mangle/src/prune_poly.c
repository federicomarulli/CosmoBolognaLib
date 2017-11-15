/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Remove all superfluous caps from polygon.

  If all you want to do is to stop garea, gspher et al from complaining,
  then use trim_poly(), not prune_poly().

  After suppressing obviously superfluous caps with trim_poly(),
  which makes garea et al happy, prune_poly applies garea to detect whether
  there are any redundant caps enclosing the entire polygon, in which case it
  removes those caps, in addition to removing caps suppressed by trim_poly().

  Null polygons are replaced with a single null cap.
  Note that a polygon with no caps is the whole sphere, not a null polygon.

   Input: poly is a pointer to a polygon.
	  mtol = initial angular tolerance in radians
		 within which to merge multiple intersections.
  Output: poly with all superfluous caps removed;
	  the number of caps is changed.
  Return value: -1 if error;
		0 if nothing changed;
		1 if something changed;
		2 if nothing changed, and polygon is null;
		3 if polygon was changed to null polygon.
*/
int prune_poly(polygon *poly, long double mtol)
{
    int i, ier, ip, iret, jp, verb;
    long double area, area_tot, cm, tol;

    /* first cut */
    iret = trim_poly(poly);

    /* trim_poly detected null polygon */
    if (iret >= 2) return(iret);

    /* area of intersection */
    tol = mtol;
    //if(tol > 0.01) printf("prune_poly: tol = %Lf\n", tol);
    verb = 1;
    ier = garea(poly, &tol, verb, &area_tot);
    if (ier) {
      printf("prune_poly: received error message from garea, mtol = %Lf\n", tol);
      return(-1);
    }

    /* null polygon */
    if (area_tot == 0.) {
	poly->rp[0][0] = 0.;
	poly->rp[0][1] = 0.;
	poly->rp[0][2] = 1.;
	poly->cm[0] = 0.;
	poly->np = 1;
	return(3);
    }

    /* test whether suppressing cap changes area or not */
    verb = 0;
    for (ip = 0; ip < poly->np; ip++) {
	if (poly->cm[ip] >= 2.) continue;	/* cap is already superfluous */
	cm = poly->cm[ip];			/* save latitude */
	poly->cm[ip] = 2.;			/* suppress cap */
	tol = mtol;
	//if(tol > 0.01) printf("prune_poly: tol = %Lf\n", tol);
	ier = garea(poly, &tol, verb, &area);	/* area sans cap */
	if (ier == -1) {
	  printf("prune_poly: received error message from garea\n");
	  return(-1);
	}
	if (ier || area != area_tot) {		/* cap affects area */
	    poly->cm[ip] = cm;			/* so restore cap */
	    if(ier) printf("prune_poly: omission of a cap caused garea failure, so cap is being restored\n");
	}
    }

    /* remove superfluous caps */
    ip = 0;
    for (jp = 0; jp < poly->np; jp++) {
	/* copy down cap */
	if (poly->cm[jp] < 2.) {
	    for (i = 0; i < 3; i++) {
		poly->rp[ip][i] = poly->rp[jp][i];
	    }
	    poly->cm[ip] = poly->cm[jp];
	    ip++;
	/* skip superfluous cap */
	} else {
	    iret = 1;
	}
    }
    poly->np = ip;

    return(iret);
}

/*------------------------------------------------------------------------------
  Suppress obviously superfluous caps from polygon,
  by setting cm = 2 for second of two coincident caps.
  In addition, detect null caps (those which contain nothing),
  and complementary caps (those which exclude each other,
  in which case replace the polygon with a single null cap.

  Two caps are considered coincident if their axes (rp) and latitudes (cm)
  are EXACTLY equal.  Caps with axes pointing in opposite directions are
  not detected.

  Two caps are considered complementary if their axes (rp) are EXACTLY equal,
  and their latitudes (cm) are EXACTLY opposing (cm of one is -cm of the other).
  Axes pointing in opposite directions are not detected.

  All this makes garea, gspher et al happy, provided that near coincident caps,
  including those with axes pointing in opposite directions, have been
  modified by `snap' to coincide exactly, with coaligned axes.

   Input: poly is a pointer to a polygon.
  Output: poly with obviously superfluous caps suppressed;
	  the number and order of caps remains unchanged
	  UNLESS polygon is replaced by null polygon.
  Return value: 0 if nothing changed;
		1 if one or more caps were suppressed;
		2 if nothing changed, and polygon is null;
		3 if polygon was changed to null polygon.
*/
int trim_poly(polygon *poly)
{
    int ip, iret, jp;

    /* initialize return value to no change */
    iret = 0;

    /* check for cap which excludes everything */
    for (jp = 0; jp < poly->np; jp++) {
	if (poly->cm[jp] == 0. || poly->cm[jp] <= -2.) {
	    if (poly->np == 1		/* polygon is already single null cap */
		&& poly->rp[0][0] == 0.
		&& poly->rp[0][1] == 0.
		&& poly->rp[0][2] == 1.
		&& poly->cm[0] == 0.) {
		return(2);
	    } else {			/* change polygon to single null cap */
		poly->rp[0][0] = 0.;
		poly->rp[0][1] = 0.;
		poly->rp[0][2] = 1.;
		poly->cm[0] = 0.;
		poly->np = 1;
		return(3);
	    }
	}
    }

    /* for each cap jp, check for coincident caps */
    for (jp = 0; jp < poly->np; jp++) {
	/* don't check superfluous cap */
	if (poly->cm[jp] >= 2.) continue;
	for (ip = jp+1; ip < poly->np; ip++) {
	    /* don't check superfluous cap */
	    if (poly->cm[ip] >= 2.) continue;
	    /* cap axes coincide */
	    if (poly->rp[ip][0] == poly->rp[jp][0]
		&& poly->rp[ip][1] == poly->rp[jp][1]
		&& poly->rp[ip][2] == poly->rp[jp][2]) {
		/* suppress coincident cap ip */
		if (poly->cm[ip] == poly->cm[jp]) {
		    poly->cm[ip] = 2.;
		    iret = 1;
		} else if (poly->cm[ip] == - poly->cm[jp]) {
		/* complementary cap means polygon is null */
		    poly->rp[0][0] = 0.;
		    poly->rp[0][1] = 0.;
		    poly->rp[0][2] = 1.;
		    poly->cm[0] = 0.;
		    poly->np = 1;
		    return(3);
		}
	    }
	}
    }

    return(iret);
}

/*------------------------------------------------------------------------------
  Suppress obviously superfluous caps from polygon,
  by setting cm = 2 for second of two coincident caps.

  Similar to trim_poly(), but does not attempt to detect null or
  complementary caps.

  Two caps are considered coincident if their axes (rp) and latitudes (cm)
  are EXACTLY equal.  Caps with axes pointing in opposite directions are
  not detected.

  In general this is not enough to make garea, gspher et al happy.

   Input: poly is a pointer to a polygon.
  Output: poly with obviously superfluous caps suppressed;
	  the number and order of caps remains unchanged.
  Return value: 0 if nothing changed;
		1 if one or more caps were suppressed.
*/
int touch_poly(polygon *poly)
{
    int ip, iret, jp;

    /* initialize return value to no change */
    iret = 0;

    /* for each cap jp, check for coincident caps */
    for (jp = 0; jp < poly->np; jp++) {
	/* don't check superfluous cap */
	if (poly->cm[jp] >= 2.) continue;
	for (ip = jp+1; ip < poly->np; ip++) {
	    /* don't check superfluous cap */
	    if (poly->cm[ip] >= 2.) continue;
	    /* cap axes coincide */
	    if (poly->rp[ip][0] == poly->rp[jp][0]
		&& poly->rp[ip][1] == poly->rp[jp][1]
		&& poly->rp[ip][2] == poly->rp[jp][2]) {
		/* suppress coincident cap ip */
		if (poly->cm[ip] == poly->cm[jp]) {
		    poly->cm[ip] = 2.;
		    iret = 1;
		}
	    }
	}
    }

    return(iret);
}
