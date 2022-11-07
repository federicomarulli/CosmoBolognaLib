/*------------------------------------------------------------------------------
  © A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include "manglefn.h"

/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP		4

/*------------------------------------------------------------------------------
  Partition disconnected polygon into connected polygons.
  Identifies the groups of boundaries of the polygon
  (two circles are friends, belonging to the same group, if they intersect,
  and friends of friends are friends),
  calls partition_gpoly to partition each group polygon into its parts,
  and combines parts from each group.

  The value of the flag overwrite_original determines whether the
  original polygon in poly is overwritten or not.

   Input: *poly is a polygon.
	  npolys = maximum number of polygons available in polys array.
	  mtol = initial angular tolerance within which to merge multiple intersections.
	  all_oneboundary = 2 to lasso all one-boundary polygons,
			    1 to lasso only those one-boundary polygons
			      with more caps than vertices,
			    0 never to lasso one-boundary polygons;
			    in all cases, all multi-boundary polygons are lassoed.
	  adjust_lasso = how to tighten lasso:
		       = 0 as tight as possible,
			 1 for balkanize,
			 2 for ransack.
	  force_split = 1 to force a polygon to be split
			  even if no boundary can be lassoed;
			0 otherwise.
	  overwrite_original = 2 to overwrite original polygon poly in all cases,
				 whether or not partitioning succeeds;
			       1 to overwrite poly only if partitioning succeeds;
			       0 never to overwrite original poly.
  Output: (*poly and) polys[i], i = 0 to *npoly-1, are the parts of *poly.
	  *npoly = number of polygons in polys;
		   if *npoly > npolys, then there was not enough space.
  Return value: -1 if error occurred;
		 0 if *poly was fully partitioned into its parts;
		 1 if *poly was partially partioned;
		   the parts of *poly constitute a valid set of parts,
		   whose union equals the input *poly,
		   though some parts are disconnected.
*/
int partition_poly(polygon **poly, int npolys, polygon *polys[/*npolys*/], long double mtol, int all_oneboundary, int adjust_lasso, int force_split, int overwrite_original, int *npoly)
{
    const int do_vcirc = 1;
    const int per = 0;
    const int nve = 2;
    int dnp, dnpoly, ier, iev, igp, ip, ipoly, jpoly, kpoly, nev, nev0, ngp, np, nret, nv, verb;
    long double area, areag, atol, tol;
    int *ipv, *gp_tmp, *ev;
    long double *angle;
    vec *ve;
    /* work arrays */
    int *gp = 0x0, *gpg = 0x0;
    polygon *gpoly = 0x0;

    /* initialize return value to normal */
    nret = 0;

    /* initialize number of polygon parts in polys to zero */
    *npoly = 0;

    /* call gverts to determine which groups connected boundaries belong to */
    tol = mtol;
    ier = gverts(*poly, do_vcirc, &tol, per, nve, &nv, &ve, &angle, &ipv, &gp_tmp, &nev, &nev0, &ev);
    if (ier) goto error;

    /* no boundaries */
    if (nev == 0) return(0);

    /* polygon has 1 connected boundary, and not too many caps */
    if (!all_oneboundary && nev == 1 && (*poly)->np <= nv + 1) return(0);

    /* allocate memory for gp and gpg */
    gp = (int *) malloc(sizeof(int) * (*poly)->np);
    if (!gp) {
	fprintf(stderr, "partition_poly: failed to allocate memory for %d ints\n", (*poly)->np);
	return(-1);
    }
    gpg = (int *) malloc(sizeof(int) * (*poly)->np);
    if (!gpg) {
	fprintf(stderr, "partition_poly: failed to allocate memory for %d ints\n", (*poly)->np);
	return(-1);
    }

    /* copy group numbers, because they will be destroyed by subsequent call to gverts() */
    for (ip = 0; ip < (*poly)->np; ip++) {
	gp[ip] = gp_tmp[ip];
    }

    /* distinct groups */
    ngp = 0;
    for (iev = 0; iev < nev; iev++) {
	if (ngp == 0 || gp[ipv[ev[iev] - 1]] != gpg[ngp - 1]) {
	    gpg[ngp] = gp[ipv[ev[iev] - 1]];
	    ngp++;
	}
    }

    /* accept error messages from garea */
    verb = 1;

    /* each group of connected boundaries */
    for (igp = 0; igp < ngp; igp++) {

	/* polygon contains only one group of connected boundaries */
	if (ngp == 1) {
	    /* point group polygon gpoly at original polygon */
	    gpoly = *poly;

	/* polygon contains more than one group of connected boundaries */
	} else {
	    /* make sure group polygon gpoly contains enough space */
	    np = (*poly)->np;
	    dnp = 0;
	    ier = room_poly(&gpoly, np, dnp, 0);
	    if (ier == -1) goto out_of_memory;

	    /* make group polygon */
	    group_poly(*poly, gp, gpg[igp], gpoly);

	}

	/* partition group polygon */
	ier = partition_gpoly(gpoly, npolys - *npoly, &polys[*npoly], tol, all_oneboundary, adjust_lasso, force_split, &dnpoly);
	if (ier == -1) goto error;

	/* flag that partitioning was only partly successful */
	if (ier == 1) nret = 1;

	/* just one group */
	if (ngp == 1) {
	    /* number of polygon parts */
	    *npoly = dnpoly;

	    /* not enough polygons */
	    if (*npoly > npolys) return(0);

	/* two or more groups */
	} else {
	    /* first group */
	    if (igp == 0) {
		if (dnpoly == 0) {
		    /* not enough polygons */
		    if (*npoly >= npolys) {
			(*npoly)++;
			return(0);
		    }

		    /* move group polygon into cumulative polygon polys[*npoly] */
		    free_poly(polys[*npoly]);
		    polys[*npoly] = gpoly;
		    gpoly = 0x0;
		    (*npoly)++;
		} else {
		    *npoly = dnpoly;
		}

	    /* subsequent groups */
	    } else {
		dnp = DNP;

		/* index of next polygon after cumulative and group polygons */
		kpoly = *npoly + dnpoly;

		/* accumulate polygons */
		for (ipoly = 0; ipoly < *npoly; ipoly++) {

		    /* area of cumulative polygon polys[ipoly] */
		    atol = tol;
		    ier = garea(polys[ipoly], &atol, verb, &area);
		    if (ier) goto error;

		    /* not enough polygons */
		    if (kpoly >= npolys) {
			*npoly = kpoly + 1;
			return(0);
		    }

		    /* intersection of cumulative polygon polys[ipoly] with group polygon gpoly */
		    np = polys[ipoly]->np + gpoly->np;
		    ier = room_poly(&polys[kpoly], np, dnp, 0);
		    if (ier == -1) goto out_of_memory;
		    poly_poly(polys[ipoly], gpoly, polys[kpoly]);

		    /* area of intersection of cumulative polygon polys[ipoly] with group polygon gpoly */
		    ier = garea(polys[kpoly], &atol, verb, &areag);
		    if (ier) goto error;

		    /* paranoid check: group polygon gpoly does not intersect cumulative polygon polys[ipoly]: should not happen */
		    if (areag == 0.) {
			fprintf(stderr, "partition_poly: group %d does not intersect part %d of polygon with input id %d; should not happen; continuing ...\n", igp, ipoly, polys[ipoly]->id);
			/* goto error; */
			/* group polygon gpoly encloses cumulative polygon polys[ipoly], so offers no further constraint */
		    } else if (areag == area) {
			continue;
		    }

		    /* group polygon was not partitioned */
		    if (dnpoly == 0) {
			/* replace cumulative polygon polys[ipoly] with intersection polygon polys[kpoly] */
			free_poly(polys[ipoly]);
			polys[ipoly] = polys[kpoly];
			polys[kpoly] = 0x0;

		    /* group polygon was partitioned */
		    } else {
			/* intersect cumulative polygon polys[ipoly] with each part of group polygon */
			for (jpoly = *npoly; jpoly < *npoly + dnpoly; jpoly++) {
			    /* not enough polygons */
			    if (kpoly >= npolys) {
				*npoly = kpoly + 1;
				return(0);
			    }
			    np = polys[ipoly]->np + gpoly->np;
			    ier = room_poly(&polys[kpoly], np, dnp, 0);
			    if (ier == -1) goto out_of_memory;
			    poly_poly(polys[ipoly], polys[jpoly], polys[kpoly]);
			    /* increment number of polygons */
			    kpoly++;
			}

			/* replace cumulative polygon polys[ipoly] with last intersected part polys[kpoly - 1] */
			free_poly(polys[ipoly]);
			kpoly--;
			polys[ipoly] = polys[kpoly];
			polys[kpoly] = 0x0;

		    }

		}

		if (dnpoly > 0) {
		    /* free parts of group polygon */
		    for (jpoly = *npoly; jpoly < *npoly + dnpoly; jpoly++) {
			free_poly(polys[jpoly]);
		    }
		    /* decrement number of polygons */
		    kpoly -= dnpoly;
		    /* move down cumulative polygons */
		    for (ipoly = *npoly; ipoly < kpoly; ipoly++) {
			polys[ipoly] = polys[ipoly + dnpoly];
		    }
		    /* nullify vacated polygons */
		    for (ipoly = kpoly; ipoly < kpoly + dnpoly; ipoly++) {
			polys[ipoly] = 0x0;
		    }
		}

		/* revised number of polygons */
		*npoly = kpoly;

	    }

	}

    }

    /* polygon contains just one group of connected boundaries */
    if (ngp == 1) {
	/* point original polygon at group polygon */
	*poly = gpoly;

    /* polygon contains more than one group of connected boundaries */
    } else {
	/* free group polygon */
	free_poly(gpoly);

    }

    /* move final polygon part to poly */
    if (*npoly > 0
	&& (overwrite_original == 2
	|| (overwrite_original == 1 && ier == 0))) {

	/* free input poly */
	free_poly(*poly);

	/* point poly at last polygon part */
	*poly = polys[*npoly - 1];

	/* nullify last polygon part */
	polys[*npoly - 1] = 0x0;

	/* decrement number of polygons */
	(*npoly)--;
    }

    if (gp) free(gp);
    if (gpg) free(gpg);

    return(nret);

    /* ---------------- error returns ---------------- */
    error:
    if (gp) free(gp);
    if (gpg) free(gpg);
    return(-1);

    out_of_memory:
    fprintf(stderr, "partition_poly: failed to allocate memory for polygon of %d caps\n", np + dnp);
    if (gp) free(gp);
    if (gpg) free(gpg);
    return(-1);
}

/*------------------------------------------------------------------------------
  Partition group polygon into connected polygons
  by calling part_poly repeatedly until the group polygon is fully partitioned,
  or until partitioning fails.

  If the group polygon has two or more connected boundaries none of which
  can be lassoed, then the force_split option controls whether part_poly
  should be forced to split the polygon in two.
  If a split is forced, then each of the resulting two polygons is subjected
  to further partitioning, or is in turn forcibly split if necessary.
  If too many forcible splits occur, then it is assumed that the procedure
  is not converging, and the routine bails out.

   Input: gpoly is a polygon all of whose circles belong to a single group.
	  npolys = maximum number of polygons available in polys array.
	  mtol = parameter passed to part_poly.
	  all_oneboundary = parameter passed to part_poly.
	  adjust_lasso = parameter passed to part_poly.
	  force_split = 1 to force part_poly to split a polygon
			  even if no boundary can be lassoed;
			0 otherwise.
  Output: polys[i], i = 0 to *npoly-1, are the parts of gpoly;
	      if return value = 0, then:
		  if *npoly = 0, then:
		      input gpoly already consists of a single connected part
		      that needs no partitioning;
		  if *npoly > 0, then:
		      gpoly was fully partitioned into its parts;
	      if return value = 1, then:
		  if force_split = 0, then:
		      gpoly was partially partitioned,
		      and polys[*npoly-1] contains those parts of the
		      input gpoly that were not partitioned successfully;
		      in spite of the failure,
		      polys[i], i = 0 to *npoly-1
		      constitute a valid set of parts,
		      whose union equals the input gpoly;
		  if force_split = 1, then:
		      the attempt to partition gpoly was abandoned after
		      gpoly was forcibly split too many times;
		      polys[i], i = 0 to *npoly-1
		      contain the parts of gpoly obtained so far;
		      in spite of the failure,
		      polys[i], i = 0 to *npoly-1
		      constitute a valid set of parts,
		      whose union equals the input gpoly,
		      though some parts are disconnected.
	  *npoly = number of polygons in polys;
		   if *npoly > npolys, then there was not enough space.
  Return value: -1 if error occurred;
		 0 if gpoly was successfully partitioned;
		 1 if gpoly was not fully partitioned.
*/
int partition_gpoly(polygon *gpoly, int npolys, polygon *polys[/*npolys*/], long double mtol, int all_oneboundary, int adjust_lasso, int force_split, int *npoly)
{
/* bail out if number of forcibly split polygons to partition exceeds this maximum */
#define	NFORCEMAX	200
    int dnpoly, dnpoly_try, iforce, ier, ipoly, iprune, jpoly, nforce;
    int do_poly[NFORCEMAX];
    polygon *poly;

    /* initialize number of polygon parts in polys to zero */
    *npoly = 0;

    /* number of forcibly split polygons to partition */
    nforce = 0;

    /* partition these polygons */
    for (iforce = -1; iforce < nforce; iforce++) {

	/* start with input group polygon */
	if (iforce == -1) {
	    poly = gpoly;
	/* subsequent polygons that have been split forcibly */
	} else {
	    poly = polys[do_poly[iforce]];
	}

	/* partition polygon repeatedly, as long as progress is made */
	while (1) {
	    /* partition polygon by lassoing its boundaries with circles */
	    ier = part_poly(poly, npolys - *npoly, &polys[*npoly], mtol, all_oneboundary, adjust_lasso, force_split, &dnpoly, &dnpoly_try);
	    if (ier == -1) return(-1);

	    /* increment number of polygons made */
	    *npoly += dnpoly;

	    /* not enough polygons */
	    if (*npoly > npolys) return(0);

	    /* partitioning was partially successful */
	    if (dnpoly > 0 && dnpoly < dnpoly_try && ier == 0) {
		/* move last polygon into poly and repeat */
		(*npoly)--;
		poly = polys[*npoly];
		polys[*npoly] = 0x0;

	    /* partitioning succeeded, or failed completely */
	    } else {
		/* poly was left unpartitioned, and is not the original group polygon */
		if (dnpoly == 0 && *npoly > 0) {
		    /* move poly into next polygon */
		    polys[*npoly] = poly;
		    /* nullify source of poly, if it is a forcibly split polygon */
		    if (iforce >= 0 && poly == polys[do_poly[iforce]]) polys[do_poly[iforce]] = 0x0;
		    (*npoly)++;
		}
		/* break out of loop */
		break;
	    }
	}

	/* part_poly forcibly split polygon into two parts, which need further partitioning */
	if (ier == 1) {
	    if (nforce + 2 <= NFORCEMAX) {
		for (ipoly = *npoly - 2; ipoly < *npoly; ipoly++) {
		    /* prune polygon that needs further partitioning */
		    iprune = prune_poly(polys[ipoly], mtol);
		    if (iprune == -1) {
			fprintf(stderr, "partition_gpoly: failed to prune forcibly split part %d of polygon with input id %d\n", nforce, polys[ipoly]->id);
			goto error;
		    }
		    if (iprune >= 2) {
			fprintf(stderr, "partition_poly: forcibly split part %d of polygon with input id %d has zero area; should not happen; continuing ...\n", nforce, polys[ipoly]->id);
			dump_poly(2, &polys[*npoly - 2]);
			continue;
		    }
		    /* flag polygon for further partitioning */
		    do_poly[nforce] = ipoly;
		    nforce++;
		}

	    /* too many forcible splits: bail out */
	    } else {
		/* advise */
		msg("partition_gpoly: unable to lasso parts of a polygon even after it has been\n");
		msg("split forcibly many times; bailing out with %d polygon parts.\n", *npoly);
		msg("CONGRATULATIONS!  YOU HAVE FOUND A POLYGON THAT BEATS MANGLE.\n");
		msg("PLEASE EMAIL ME Andrew.Hamilton@colorado.edu THE GOOD NEWS\n");
		msg("ALONG WITH A POLYGON FILE CONTAINING THE POLYGON THAT DID IT.\n");
		/* dump the polygon to a polygon file */
		dump_poly(1, &gpoly);
		msg("AND THERE'S THE POLYGON FILE I'D LIKE YOU TO SEND.  THANKS!\n");
		/* number of forcibly split polygons that have been partitioned */
		nforce = iforce + 1;
		/* break out of partitioning loop */
		break;

	    }
	}
    }

    /* flag that partitioning was only partially successful */
    if (dnpoly < dnpoly_try) ier = 1;

    /* remove forcibly split polygons that were partitioned */
    if (nforce > 0) {
	iforce = 0;
	jpoly = 0;
	for (ipoly = 0; ipoly < *npoly; ipoly++) {
	    /* free forcibly split polygons */
	    if (iforce < nforce && ipoly == do_poly[iforce]) {
		if (polys[ipoly]) free_poly(polys[ipoly]);
		iforce++;
	    /* move down polygons */
	    } else {
		polys[jpoly] = polys[ipoly];
		jpoly++;
	    }
	}
	/* nullify vacated polygons */
	for (ipoly = jpoly; ipoly < *npoly; ipoly++) {
	    polys[ipoly] = 0x0;
	}
	/* revise number of polygons */
	*npoly = jpoly;
    }

    return(ier);

    /* ---------------- error returns ---------------- */
    error:
    return(-1);
}

/*------------------------------------------------------------------------------
  Partition a polygon by lassoing its connected boundaries with circles.

  Normally the input polygon poly would be a group polygon,
  all of whose connected boundaries belong to a single group.
  If the boundaries of the input polygon belong to more than one group,
  then the polygon is not simply-connected, and each non-simply-connected part
  of the polygon will contain two or more boundaries.
  The algorithm will attempt to lasso all these boundaries, even though
  lassoing boundaries of a non-simply-connected part of the polygon
  cannot partition the polygon.

  The routine attempts to lasso all the connected boundaries of a polygon,
  unless the polygon has only a single connected boundary.
  If the polygon has a single connected boundary, then
  the all_oneboundary option controls whether or not this boundary is lassoed.

  An attempted lasso is discarded if it lies fully inside or fully outside
  all the caps of the polygon.  It would be incorrect to retain a lasso
  that lies fully inside the polygon, and it would be superfluous to retain
  a lasso that fully encloses the polygon.
  In the normal case where the input polygon is a group polygon,
  a lasso can lie fully inside or outside the caps of the group polygon
  if the group polygon has a single connected boundary.

  If the input polygon has two or more connected boundaries none of which
  can be lassoed, then the force_split option controls whether the
  routine gives up, or else forcibly splits the input polygon into two parts,
  each of which will require further splitting.

   Input: poly is a polygon.
	  npolys = maximum number of polygons available in polys array.
	  mtol = initial angular tolerance within which to merge multiple intersections.
	  all_oneboundary = 2 to lasso all one-boundary polygons,
			    1 to lasso only those one-boundary polygons
			      with more caps than vertices,
			    0 never to lasso one-boundary polygons;
			    in all cases, all multi-boundary polygons are lassoed.
	  adjust_lasso = how to tighten lasso:
		       = 0 as tight as possible,
			 1 for balkanize,
			 2 for ransack.
	  force_split = 1 to force a polygon to be split
			  even if no boundary can be lassoed;
			0 otherwise.
  Output: polys[i], i = 0 to npoly-1, are the connected parts of poly;
              if return value = 0, then:
		  if *npoly = *npoly_try = 0, then:
		      input poly already consists of a single connected part
		      that needs no partitioning;
		  if *npoly = *npoly_try > 0, then:
		      poly was fully partitioned into its parts;
		  if *npoly < *npoly_try then:
		      poly was partially partitioned,
		      and poly[npoly - 1] contains the parts of poly
		      that were not partitioned successfully;
	      if return value = 1,
		    which can occur only if force_split = 1, then:
		    poly contains at least 2 connected boundaries,
		    none of which could be lassoed successfully,
		    and instead poly was split forcibly into two,
		    with the two parts in polys[0] and polys[1].
	  *npoly = number of polygons in polys;
		   if *npoly > npolys, then there was not enough space.
	  *npoly_try = attempted number of polygons.
  Return value: -1 if error occurred;
		 0 for a normal return;
		 1 if poly was split forcibly when no boundary could be lassoed;
		   can only occur if force_split = 1.
*/
int part_poly(polygon *poly, int npolys, polygon *polys[/*npolys*/], long double mtol, int all_oneboundary, int adjust_lasso, int force_split, int *npoly, int *npoly_try)
{
/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP		4
    static polygon *extracap = 0x0;
    const int do_vcirc = 1;
    const int per = 0;
    const int nve = 2;
    const int itmax = 30;
    int dnp, found, i, ier, iev, ip, it, iv, ivm, ivmax_that, ivmax_this, ivmin_that, ivmin_this, nev, nev0, np, nret, nv, nvm;
    int *ipv, *gp, *ev;
    long double *angle;
    vec *ve, *vm;
    long double *cmvmin, *cmvmax, *cmpmin, *cmpmax;
    vec *vmax, *vmin;
    long double cmbest, cme, cmforce, dth, dthbest, dthforce, dthm, dthp, s, th, thm, tol;
    long double cmpmax_all, cmpmin_all, cmvmax_that, cmvmax_this, cmvmin_that, cmvmin_this, thmax_that, thmax_this, thmin_that, thmin_this;
    vec v, vmbest, vmforce;

    /* initialize return value to normal */
    nret = 0;

    /* initialize number of polygon parts in polys to zero */
    *npoly = 0;

    /* initialize attempted number of polygon parts to zero */
    *npoly_try = 0;

    /* vertices and centres of edges of polygon */
    tol = mtol;
    ier = gverts(poly, do_vcirc, &tol, per, nve, &nv, &ve, &angle, &ipv, &gp, &nev, &nev0, &ev);
    if (ier) return(-1);

    /* no boundaries */
    if (nev == 0) return(nret);

    /* polygon has just 1 connected boundary */
    if (nev == 1) {
	if (all_oneboundary == 0) {
	    return(nret);
	} else if (all_oneboundary == 1) {
	    /* polygon has not too many caps */
	    if (poly->np <= nv + 1) return(nret);
	}
    }

    /* barycentres of connected boundaries of polygon */
    ier = vmidc(poly, nv, nve, ve, ipv, ev, &nvm, &vm);
    if (ier == -1) return(-1);

    /* number of polygons to try to split into */
    *npoly_try = nvm;

    /* initialize lasso to be used to force split */
    if (force_split) dthforce = -4.;

    /* attempt to partition polygon around each barycentre vm[ivm] */
    for (ivm = 0; ivm < nvm; ivm++) {

	/* initialize best lasso so far */
	if (force_split) dthbest = -4.;

	/* repeat until find isolating lasso */
	it = 0;
	do {

	    /* points on each edge nearest to and farthest from vm[ivm] */
	    ier = gvlims(poly, do_vcirc, &tol, vm[ivm], &nv, &vmin, &vmax, &cmvmin, &cmvmax, &cmpmin, &cmpmax, &ipv, &gp, &nev, &nev0, &ev);
	    if (ier == -1) return(-1);
	    if (ier) break;

	    /* distances that exclude/enclose all circles */
	    cmpmin_all = 2.;
	    cmpmax_all = 0.;
	    for (ip = 0; ip < poly->np; ip++) {
		if (cmpmin[ip] <= cmpmin_all) {
		    cmpmin_all = cmpmin[ip];
		}
		if (cmpmax[ip] >= cmpmax_all) {
		    cmpmax_all = cmpmax[ip];
		}
	    }

	    /* distances that exclude/enclose this connected boundary */
	    ivmax_this = -1;
	    ivmin_this = -1;
	    cmvmin_this = 2.;
	    cmvmax_this = 0.;
	    for (iv = (ivm == 0)? 0 : ev[ivm - 1]; iv < ev[ivm]; iv++) {
		if (cmvmin[iv] <= cmvmin_this) {
		    ivmin_this = iv;
		    cmvmin_this = cmvmin[iv];
		}
		if (cmvmax[iv] >= cmvmax_this) {
		    ivmax_this = iv;
		    cmvmax_this = cmvmax[iv];
		}
	    }

	    /* distances that exclude/enclose other connected boundaries */
	    ivmin_that = -1;
	    ivmax_that = -1;
	    cmvmin_that = 2.;
	    cmvmax_that = 0.;
	    for (iev = 0; iev < nev; iev++) {
		if (iev == ivm) continue;
		for (iv = (iev == 0)? 0 : ev[iev - 1]; iv < ev[iev]; iv++) {
		    if (cmvmin[iv] <= cmvmin_that) {
			ivmin_that = iv;
			cmvmin_that = cmvmin[iv];
		    }
		    if (cmvmax[iv] >= cmvmax_that) {
			ivmax_that = iv;
			cmvmax_that = cmvmax[iv];
		    }
		}
	    }

	    /* angles corresponding to cmmin_this, cmvmax_this, cmvmin_that, and cmvmax_that */
	    thmin_this = 2. * asinl(sqrtl(cmvmin_this / 2.));
	    thmax_this = 2. * asinl(sqrtl(cmvmax_this / 2.));
	    thmin_that = 2. * asinl(sqrtl(cmvmin_that / 2.));
	    thmax_that = 2. * asinl(sqrtl(cmvmax_that / 2.));

	    dthp = thmin_that - thmax_this;
	    dthm = thmin_this - thmax_that;

	    dth = (dthp >= dthm)? dthp : dthm;

	    /* found lasso that isolates this connected boundary of polygon */
	    if (dth >= - tol) {

		/* if (it >= 10) printf("%21.15Lg %21.15Lg %21.15Lg %4d %2d%21.15Lg %21.15Lg %21.15Lg %21.15Lg %21.15Lg %21.15Lg\n", vm[ivm][0], vm[ivm][1], vm[ivm][2], ivm, it, thmin_this, thmax_this, thmin_that, thmax_that, dthp, dthm); */

		/* thmin_that >= thmax_this - tol */
		if (dthp >= dthm) {
		    /* isolating lasso */
		    th = (thmax_this + thmin_that) / 2.;
		    switch (adjust_lasso) {
		    /* as tight as possible */
		    case 0:	thm = thmax_this;		break;
		    /* for balkanize: tiny angles give garea problems */
		    case 1:	thm = thmax_this + .001;	break;
		    /* for ransack: want tight lasso */
		    case 2:	thm = thmax_this * 1.05;	break;
		    }
		    /* tighten lasso */
		    if (th > thm) th = thm;
		    if (th < PI) {
			s = sinl(th / 2.);
			cme = 2. * s * s;
		    } else {
			cme = 2.;
		    }

		    /* discard lasso that completely encloses all circles */
		    if (cme >= cmpmax_all) {
			/* printf("%21.15Lg %21.15Lg %21.15Lg %21.15Lg\n", vm[ivm][0], vm[ivm][1], vm[ivm][2], cme); */
			/* decrement number of polygons to try for */
			(*npoly_try)--;
			/* break out of search loop */
			break;
		    }

		/* thmin_this >= thmax_that - tol */
		} else {
		    /* isolating lasso */
		    th = (thmax_that + thmin_this) / 2.;
		    switch (adjust_lasso) {
		    /* as tight as possible */
		    case 0:	thm = thmin_this;		break;
		    /* for balkanize: tiny angles give garea problems */
		    case 1:	thm = thmin_this - .001;	break;
		    /* for ransack: want tight lasso */
		    case 2:	thm = thmin_this / 1.05;	break;
		    }
		    /* tighten lasso */
		    if (th < thm) th = thm;
		    if (th > 0.) {
			s = sinl(th / 2.);
			cme = 2. * s * s;
		    } else {
			cme = 0.;
		    }

		    /* discard lasso that completely encloses all circles */
		    if (cme <= cmpmin_all) {
			/* printf("%21.15Lg %21.15Lg %21.15Lg %21.15Lg\n", vm[ivm][0], vm[ivm][1], vm[ivm][2], cme); */
			/* decrement number of polygons to try for */
			(*npoly_try)--;
			/* break out of search loop */
			break;
		    }

		    cme = - cme;

		}

		/* not enough polygons for a new one */
		if (*npoly >= npolys) {
		    (*npoly)++;
		    return(0);
		}

		/* put isolating lasso into new polygon */
		np = 1;
		dnp = DNP;
		ier = room_poly(&extracap, np, dnp, 0);
		if (ier == -1) goto out_of_memory;
		for (i = 0; i < 3; i++) {
		    extracap->rp[0][i] = vm[ivm][i];
		}

		extracap->cm[0] = cme;
		extracap->np = 1;

		/* make sure new polygon contains enough space */
		np = poly->np + 1;
		dnp = 0;
		ier = room_poly(&polys[*npoly], np, dnp, 0);
		if (ier == -1) goto out_of_memory;

		/* combination of poly with new circle */
		poly_poly(poly, extracap, polys[*npoly]);

		/* increment number of polygons */
		(*npoly)++;

		/* flag found isolating boundary */
		found = 1;

	    /* failed to find isolating lasso */
	    } else {

		dthp = thmin_that - thmin_this;
		dthm = thmax_this - thmax_that;

		if (dthp >= dthm) {
		    /* paranoid check that ivmax_this and ivmin_that were initialized */
		    if (ivmax_this == -1 || ivmin_that == -1) {
			/* error should never happen */
			fprintf(stderr, "partition_poly: ivmax_this = %d ivmin_that = %d should be in interval [0, %d]\n", ivmax_this, ivmin_that, nv);
			return(-1);
		    }

		    /* record the best lasso so far */
		    if (force_split) {
			if (dthp > dthbest) {
			    dthbest = dthp;
			    /* circle that does not enclose connected boundary,
			       but does exclude other connected boundaries */
			    th = (thmin_this + thmin_that) / 2.;
			    thm = thmin_that * .999;
			    if (th < thm) th = thm;
			    s = sinl(th / 2.);
			    cmbest = 2. * s * s;
			    for (i = 0; i < 3; i++) {
				vmbest[i] = vm[ivm][i];
			    }
			}
		    }

		    /* vector to this from that */
		    for (i = 0; i < 3; i++) {
			v[i] = vmax[ivmax_this][i] - vmin[ivmin_that][i];
		    }
		    s = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		    /* translate centre point along vector */
		    /* the 0.02 puts centre just beyond equal distance */
		    for (i = 0; i < 3; i++) {
			v[i] = vm[ivm][i] + ((cmvmax_this - cmvmin_that) / s + 0.02 * (it + 1)) * v[i];
		    }
		    s = sqrtl(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		    for (i = 0; i < 3; i++) {
			vm[ivm][i] = v[i] / s;
		    }

		} else {
		    /* paranoid check that ivmin_this and ivmax_that were initialized */
		    if (ivmin_this == -1 || ivmax_that == -1) {
			/* error should never happen */
			fprintf(stderr, "partition_poly: ivmin_this = %d ivmax_that = %d should be in interval [0, %d]\n", ivmin_this, ivmax_that, nv);
			return(-1);
		    }

		    /* record the best lasso so far */
		    if (force_split) {
			if (dthm > dthbest) {
			    dthbest = dthm;
			    /* circle that does not enclose connected boundary,
			       but does exclude other connected boundaries */
			    th = PI - (thmax_this + thmax_that) / 2.;
			    thm = (PI - thmax_that) * .999;
			    if (th < thm) th = thm;
			    s = cosl(th / 2.);
			    cmbest = - 2. * s * s;
			    for (i = 0; i < 3; i++) {
				vmbest[i] = vm[ivm][i];
			    }
			}
		    }

		    /* vector to that from this */
		    for (i = 0; i < 3; i++) {
			v[i] = vmax[ivmax_that][i] - vmin[ivmin_this][i];
		    }
		    s = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		    /* translate centre point along vector */
		    /* the 0.02 puts centre just beyond equal distance */
		    for (i = 0; i < 3; i++) {
			v[i] = vm[ivm][i] + ((cmvmax_that - cmvmin_this) / s + 0.02 * (it + 1)) * v[i];
		    }
		    s = sqrtl(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		    for (i = 0; i < 3; i++) {
			vm[ivm][i] = v[i] / s;
		    }

		}

		/* flag failed to find isolating boundary */
		found = 0;

	    }

	} while (!found && it++ < itmax);

	/* record the best lasso so far */
	if (force_split && !found) {
	    if (dthbest > dthforce) {
		dthforce = dthbest;
		cmforce = cmbest;
		for (i = 0; i < 3; i++) {
		    vmforce[i] = vmbest[i];
		}
	    }
	}

    }

    /* no polygons were lassoed */
    if (*npoly == 0) {

	/* go with original polygon */
	if (*npoly_try <= 1) {
	    *npoly_try = 0;
	    return(nret);
	}

	/* split polygon along best circle so far, even though it is not a lasso */
	if (force_split) {
	    /* not enough polygons for a new one */
	    if (*npoly >= npolys) {
		(*npoly)++;
		return(0);
	    }

	    /* put circle into new polygon */
	    np = 1;
	    dnp = DNP;
	    ier = room_poly(&extracap, np, dnp, 0);
	    if (ier == -1) goto out_of_memory;
	    for (i = 0; i < 3; i++) {
		extracap->rp[0][i] = vmforce[i];
	    }
	    extracap->cm[0] = cmforce;
	    extracap->np = 1;

	    /* make sure new polygon contains enough space */
	    np = poly->np + 1;
	    dnp = 0;
	    ier = room_poly(&polys[*npoly], np, dnp, 0);
	    if (ier == -1) goto out_of_memory;

	    /* combination of polygon with new circle */
	    poly_poly(poly, extracap, polys[*npoly]);

	    /* increment number of polygons */
	    (*npoly)++;

	    /* flag that split was forced */
	    nret = 1;
	}

    }

    /* if some polygons were lassoed and others not, add complement of all new caps to polygon */
    if (*npoly > 0 && *npoly < *npoly_try) {
	    /* not enough polygons for a new one */
	if (*npoly >= npolys) {
	    (*npoly)++;
	    return(0);
	}

	/* put complement of all new caps into new polygon */
	np = *npoly;
	dnp = DNP;
	ier = room_poly(&extracap, np, dnp, 0);
	if (ier == -1) goto out_of_memory;
	for (np = 0; np < *npoly; np++) {
	    ip = polys[np]->np - 1;
	    for (i = 0; i < 3; i++) {
		extracap->rp[np][i] = polys[np]->rp[ip][i];
	    }
	    extracap->cm[np] = - polys[np]->cm[ip];
	}
	extracap->np = *npoly;

	/* make sure new polygon contains enough space */
	np = poly->np + *npoly;
	dnp = 0;
	ier = room_poly(&polys[*npoly], np, dnp, 0);
	if (ier == -1) goto out_of_memory;

	/* poly with complement of new caps from other polygons */
	poly_poly(poly, extracap, polys[*npoly]);

	/* increment number of polygons */
	(*npoly)++;
    }

    /* trim new polygons to suppress obviously coincident caps */
    for (ip = 0; ip < *npoly; ip++) {
	trim_poly(polys[ip]);
    }

    return(nret);

    /* ---------------- error returns ---------------- */
    out_of_memory:
    fprintf(stderr, "part_poly: failed to allocate memory for polygon of %d caps\n", np + dnp);
    return(-1);
}
