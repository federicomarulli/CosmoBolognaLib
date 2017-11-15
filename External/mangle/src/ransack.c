/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqm:c:r:s:e:u:p:";

/* allocate polygons as a global array */
polygon *poly_global[NPOLYSMAX];

/* local functions */
void	usage(void);
#ifdef	GCC
int	ransack(char *, format *, int, int npolysmax, polygon *[npolysmax]);
int	lasso_poly(polygon **, int npolys, polygon *[npolys], long double, int *);
#else
int	ransack(char *, format *, int, int npolysmax, polygon *[/*npolysmax*/]);
int	lasso_poly(polygon **, int npolys, polygon *[/*npolys*/], long double, int *);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, np, npoly, npolys,i;
    polygon **poly;
    poly=poly_global;

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind >= 1) {
            fprintf(stderr, "%s requires at least 2 arguments: polygon_infile, and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- ransack ----------------\n");

    /* advise data format */
    advise_fmt(&fmt);

    /* tolerance angle for multiple intersections */
    if (mtol != 0.) {
	scale(&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %Lg%c will be treated as coincident\n", mtol, munit);
	scale(&mtol, munit, 'r');
	munit = 'r';
    }

    /* warn about seed not being set */
    if (seed_read == 0) {
	msg("warning: seed was not set on command line: using default seed %d\n", seed);
    }

    /* read polygons */
    npoly = 0;
    nfiles = argc - 1 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &poly[npoly]);
	if (npolys == -1) exit(1);
	npoly += npolys;
    }
    if (nfiles >= 2) {
	msg("total of %d polygons read\n", npoly);
    }
    if (npoly == 0) {
	msg("STOP\n");
	exit(0);
    }
 
   if (snapped==0 || balkanized==0) {
     msg("WARNING: 'snapped' and 'balkanized' keywords not found in all input files.\n");
     msg("Running ransack on polygons that are not snapped and balkanized may give misleading results.\n");
   }

    /* random points in polygons */
    ifile = argc - 1;
    np = ransack(argv[ifile], &fmt, npoly, NPOLYSMAX, poly);
    if (np == -1) exit(1);
 
    for(i=0;i<npoly;i++){
      free_poly(poly[i]);
    }

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("ransack [-d] [-q] [-c<seed>] [-r<n>] [-m<a>[u]] [-s<n>] [-e<n>] [-u<inunit>[,<outunit>]] [-p[+|-][<n>]] [-i<f>[<n>][u]] polygon_infile1 [polygon_infile2 ...] outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Generate random az, el positions within mask defined by poly.
  The results are written to out_filename.

   Input: out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  npoly = number of polygons in poly array.
	  npolysmax = maximum number of polygons in poly array.
	  poly = array of pointers to polygons.
	  mtol = initial tolerance angle for multiple intersections.
  Return value: number of random points generated,
		or -1 if error occurred.
*/
int ransack(char *out_filename, format *fmt, int npoly, int npolysmax, polygon *poly[/*npolysmax*/])
{
/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP			4
/* length of state vector for random number generator */
#define STATELEN		256
    static char state[STATELEN], stateo[STATELEN];
#define AZEL_STR_LEN		32
    char output[] = "output";
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int dnp, dnwl, i, idmin, idmax, idwidth, ier, in, inull, ip, ipmin, ipoly, iprune, irandom, lassoed, np, nwl, tries, verb, width, k;
    int *dlasso=0x0, *lasso=0x0;
    long double area, cmmin, cmi, phi, rpoly, si, tol, w, wcum, x, y, z;
    long double *wpoly;
    vec rp, xi, yi;
    azel v;
    char *out_fn;
    FILE *outfile;

    /* open out_filename for writing */
    if (!out_filename || strcmp(out_filename, "-") == 0) {
	outfile = stdout;
	out_fn = output;
    } else {
	outfile = fopen(out_filename, "w");
	if (!outfile) {
	    fprintf(stderr, "ransack: cannot open %s for writing\n", out_filename);
	    goto error;
	}
	out_fn = out_filename;
    }

    /* advise angular units */
    if (fmt->outunit != fmt->inunit) {
	msg("units of output az, el angles will be ");
	switch (fmt->outunit) {
#include "angunit.h"
	}
	msg("\n");
    }

    /* initialize random number generator used by ransack() */
    initstate(seed, state, STATELEN);
    /* initialize random number generator used by ikrand() */
    initstate(seed, stateo, STATELEN);

    /* prune polygons, discarding those with zero weight * area */
    msg("pruning %d polygons ...\n", npoly);
    ier = 0;
    inull = 0;
    np = 0;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	/* zero weight polygon */
	if (poly[ipoly]->weight == 0.) {
            inull++;
	    free_poly(poly[ipoly]);
	    poly[ipoly] = 0x0;
	} else {
	/* prune polygon */
	    iprune = prune_poly(poly[ipoly], mtol);
	    /* error */
	    if (iprune == -1) {
		ier++;
		free_poly(poly[ipoly]);
		poly[ipoly] = 0x0;
		fprintf(stderr, "ransack: failed to prune polygon %d; discard it\n", ipoly);
		/* goto error; */
	    /* zero area polygon */
	    } else if (iprune >= 2) {
		inull++;
		free_poly(poly[ipoly]);
		poly[ipoly] = 0x0;
	    } else {
		np++;
	    }
	}
    }
   /*copy down non-null polygons*/
    k=0;
    for(ipoly = 0; ipoly < npoly; ipoly++){
      if(poly[ipoly]){
	poly[k++]=poly[ipoly];
      }
    }
    /*after copying non-null polygons, k should be equal to np */
    if(k!=np){
      fprintf(stderr, "ransack: should be left with %d non-null polygons, but actually have %d\n",np,k);
    }

    /*nullify the rest of the array, but don't free, since pointers have been copied above*/
    for(ipoly=np; ipoly < npoly; ipoly++){
      poly[ipoly]=0x0;
    }

    if (ier > 0) {
	msg("discarding %d unprunable polygons\n", ier);
    }
    if (inull > 0) {
	msg("discarding %d polygons with zero weight * area\n", inull);
    }
    /* number of polygons with finite weight * area */
    npoly = np;

    /* no polygons */
    if (npoly == 0) {
	fprintf(stderr, "ransack: no polygons to generate random points inside!\n");
	goto error;
    }

    /* pre-lasso polygons if there are many random points */
    if (nrandom >= npoly) {
	msg("lassoing %d polygons ...\n", npoly);

	/* lasso each polygon */
	np = npoly;
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    ier = lasso_poly(&poly[ipoly], npolysmax - np, &poly[np], mtol, &dnp);
	    if (ier == -1) {
		fprintf(stderr, "ransack: UHOH at polygon %d; continuing ...\n", poly[ipoly]->id);
	    }

	    /* lassoed polygons are an improvement over original polygon */
	    if (dnp > 0) {
		/* check whether exceeded maximum number of polygons */
		if (np + dnp > npolysmax) {
		    fprintf(stderr, "ransack: total number of polygons exceeded maximum %d\n", npolysmax);
		    fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
		    goto error;
		}

		/* decrement dnp by 1 */
		dnp--;

		/* increment number of polygons */
		np += dnp;

		/* move last polygon part into poly[ipoly] */
		free_poly(poly[ipoly]);
		poly[ipoly] = poly[np];
		poly[np] = 0x0;
	    }
	}

	/* revised number of polygons */
	npoly = np;

	/* flag that all polygons have been lassoed */
	lassoed = 1;

    /* two few random points to make it worth pre-lassoing */
    } else {
	/* flag that all polygons have not been lassoed */
	lassoed = 0;

    }

    /* allocate memory for wpoly array */
    nwl = npoly;
    wpoly = (long double *) malloc(sizeof(long double) * nwl);
    if (!wpoly) {
        fprintf(stderr, "ransack: failed to allocate memory for %d long doubles\n", nwl);
        goto error;
    }
    if (!lassoed) {
	/* allocate memory for lasso and dlasso arrays */
	lasso = (int *) malloc(sizeof(int) * nwl);
	if (!lasso) {
	    fprintf(stderr, "ransack: failed to allocate memory for %d ints\n", nwl);
	    goto error;
	}
	dlasso = (int *) malloc(sizeof(int) * nwl);
	if (!dlasso) {
	    fprintf(stderr, "ransack: failed to allocate memory for %d ints\n", nwl);
	    goto error;
	}

	/* initialize dlasso array to zero */
	for (ipoly = 0; ipoly < nwl; ipoly++) dlasso[ipoly] = 0;
    }

    /* largest width of polygon id number */
    idmin = 0;
    idmax = 0;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	if (poly[ipoly]->id < idmin) idmin = poly[ipoly]->id;
	if (poly[ipoly]->id > idmax) idmax = poly[ipoly]->id;
    }
    idmin = ((idmin < 0)? floorl(log10l((long double)-idmin)) + 2 : 1);
    idmax = ((idmax > 0)? floorl(log10l((long double)idmax)) + 1 : 1);
    idwidth = ((idmin > idmax)? idmin : idmax);

    /* write header */
    wrangle(0., fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
    width = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    }
    fprintf(outfile, "%*s\t%*s\t%*s\n", width, az_str, width, el_str, idwidth, "id");

    /* accept error messages from garea */
    /* unprunable polygons were already discarded, so garea should give no errors */
    verb = 1;

    /* cumulative area times weight of polygons */
    w = 0.;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	/* skip null polygons */
	if (poly[ipoly]) {
	    /* area of polygon */
	    tol = mtol;
	    ier = garea(poly[ipoly], &tol, verb, &area);
	    if (ier) goto error;
	    /* accumulate weight times area */
	    w += poly[ipoly]->weight * area;
	}
	wpoly[ipoly] = w;
    }
    wcum = w;

    /* random points */
    if (strcmp(out_fn, output) != 0) {
	msg("generating %d random points from seed %u in %d polygons ...\n", nrandom, seed, npoly);
    }
    for (irandom = 0; irandom < nrandom; irandom++) {

	/* random number in interval [0, 1) wcum */
	setstate(state);
	rpoly = drandom() * wcum;
	setstate(stateo);

	/* which polygon to put random point in */
	ipoly = search(npoly, wpoly, rpoly);

	/* guard against roundoff */
	if (ipoly >= npoly) {
	    fprintf(stderr, "ransack: %d should be < %d (i.e. %.15Lg < %.15Lg)\n", ipoly, npoly, rpoly, wpoly[npoly - 1]);
	    ipoly = npoly - 1;
	}

	/* all polygons have not been lassoed */
	if (!lassoed) {

	    /* polygon has not yet been lassoed */
	    if  (dlasso[ipoly] == 0) {

		/* lasso polygon */
		ier = lasso_poly(&poly[ipoly], npolysmax - np, &poly[np], mtol, &dnp);
		if (ier == -1) {
		    fprintf(stderr, "ransack: UHOH at polygon %d; continuing ...\n", poly[ipoly]->id);
		}

		/* go with original polygon */
		if (dnp == 0) {
		    /* lasso, dlasso */
		    lasso[ipoly] = ipoly;
		    dlasso[ipoly] = 1;

		/* lassoed polygons are an improvement over original */
		} else {
		    /* check whether exceeded maximum number of polygons */
		    if (np + dnp > npolysmax) {
			fprintf(stderr, "ransack: total number of polygons exceeded maximum %d\n", npolysmax);
			fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
			goto error;
		    }

		    /* just one lassoed polygon */
		    if (dnp == 1) {
			/* move last polygon part into poly[ipoly] */
			free_poly(poly[ipoly]);
			poly[ipoly] = poly[np];
			poly[np] = 0x0;

			/* lasso, dlasso */
			lasso[ipoly] = ipoly;
			dlasso[ipoly] = 1;

		    /* more than one lassoed polygon */
		    } else {
			/* enlarge memory for wpoly, lasso, and dlasso arrays */
			if (np + dnp > nwl) {
			    dnwl = dnp + 1024;
			    wpoly = (long double *) realloc(wpoly, sizeof(long double) * (nwl + dnwl));
			    if (!wpoly) {
				fprintf(stderr, "ransack: failed to reallocate memory for %d long doubles\n", nwl + dnwl);
				goto error;
			    }
			    lasso = (int *) realloc(lasso, sizeof(int) * (nwl + dnwl));
			    if (!lasso) {
				fprintf(stderr, "ransack: failed to reallocate memory for %d ints\n", nwl + dnwl);
				goto error;
			    }
			    dlasso = (int *) realloc(dlasso, sizeof(int) * (nwl + dnwl));
			    if (!dlasso) {
				fprintf(stderr, "ransack: failed to reallocate memory for %d ints\n", nwl + dnwl);
				goto error;
			    }

			    /* initialize new part of lasso and dlasso arrays to inconsistent values */
			    for (ipoly = nwl; ipoly < nwl + dnwl; ipoly++) lasso[ipoly] = 1;
			    for (ipoly = nwl; ipoly < nwl + dnwl; ipoly++) dlasso[ipoly] = 0;

			    /* revised size of wpoly, lasso, and dlasso arrays */
			    nwl += dnwl;
			}

			/* lasso, dlasso */
			lasso[ipoly] = np;
			dlasso[ipoly] = dnp;

			/* cumulative weight times area of lassoed polygons */
			w = (ipoly == 0)? 0. : wpoly[ipoly-1];
			for (ip = np; ip < np + dnp; ip++) {
			    /* area of polygon */
			    tol = mtol;
			    ier = garea(poly[ip], &tol, verb, &area);
			    if (ier) goto error;
			    /* accumulate area times weight */
			    w += poly[ip]->weight * area;
			    wpoly[ip] = w;
			}

			/* increment number of polygons */
			np += dnp;
		    }

		}

	    }

	    /* polygon was partitioned into at least two */
	    if (dlasso[ipoly] >= 2) {
		/* which polygon to put random point in */
		ip = search(dlasso[ipoly], &wpoly[lasso[ipoly]], rpoly);

		/* guard against roundoff */
		if (ip >= lasso[ipoly] + dlasso[ipoly]) {
		    fprintf(stderr, "ransack: %d should be < %d (i.e. %.15Lg < %.15Lg)\n", ip, lasso[ipoly] + dlasso[ipoly], rpoly, wpoly[lasso[ipoly] + dlasso[ipoly] - 1]);
		    ip = lasso[ipoly] + dlasso[ipoly] - 1;
		}

		/* revised polygon number to put random point in */
		ipoly = ip;
	    }
	}

	/* smallest cap of polygon */
	cmminf(poly[ipoly], &ipmin, &cmmin);

	/* random point within polygon */
	tries = 0;
	do {
	    tries++;
	    /* random point within smallest cap */
	    setstate(state);
	    phi = TWOPI * drandom();
	    cmi = cmmin * drandom();
	    setstate(stateo);
	    /* coordinates of random point in cap frame */
	    si=sqrtl(cmi * (2. - cmi));
	    x = si * cosl(phi);
	    y = si * sinl(phi);
	    z = 1. - cmi;
	    /* polygon has caps */
	    if (poly[ipoly]->np > 0) {
		if (poly[ipoly]->cm[ipmin] < 0.) z = -z;
		/* Cartesian axes with z-axis along cap axis */
		gaxisi_(poly[ipoly]->rp[ipmin], xi, yi);
		/* coordinates of random point */
		for (i = 0; i < 3; i++) rp[i] = x * xi[i] + y * yi[i] + z * poly[ipoly]->rp[ipmin][i];
		/* whether random point is inside polygon */
		in = gptin(poly[ipoly], rp);
	    /* polygon has no caps, so is the whole sphere */
	    } else {
		rp[0] = x;
		rp[1] = y;
		rp[2] = z;
		in = 1;
	    }
	} while (!in);

	/* convert unit vector to az, el */
	rp_to_azel(rp, &v);
	v.az -= floorl(v.az / TWOPI) * TWOPI;

	/* convert az and el from radians to output units */
	scale_azel(&v, 'r', fmt->outunit);

	/* write result */
	wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
	wrangle(v.el, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, el_str);
	fprintf(outfile, "%s\t%s\t%*d\n", az_str, el_str, idwidth, poly[ipoly]->id);
	/* fprintf(outfile, "%s %s %d %d %d %Lg %Lg %Lg %Lg %d %d\n", az_str, el_str, irandom, ipoly, tries, wcum, rpoly / wcum, area, TWOPI * cmmin / area, ipmin, poly[ipoly]->np); */

    }

    /* advise */
    if (outfile != stdout) {
	msg("ransack: %d random positions written to %s\n", nrandom, out_fn);
    }

    return(nrandom);

    /* error returns */
    error:
    return(-1);
}

/*------------------------------------------------------------------------------
  Lasso polygon,
  keeping the lassoed parts only if the sum of the areas of lassos is
  sufficiently less than the area of the tightest cap of the original polygon.

  Output: *np = number of lassoed parts;
              = 0 to retain original polygon.
  Return value: same as partition_poly:
		-1 if error occurred;
		 0 ok;
		 1 if *poly was only partially partioned.
*/
int lasso_poly(polygon **poly, int npolys, polygon *polys[/*npolys*/], long double mtol, int *np)
{
/* part_poly should lasso all one-boundary polygons */
#define ALL_ONEBOUNDARY		2
/* how part_poly should tighten lasso */
#define ADJUST_LASSO		2
/* part_poly should not force polygon to be split even if no part can be lassoed */
#define FORCE_SPLIT		0
/* partition_poly should never overwrite original polygons */
#define	OVERWRITE_ORIGINAL	0
    int ier, ip, ipmin;
    long double cmmin, cmmino, cmmint;

    /* area/(2 pi) of smallest cap of polygon */
    cmminf(*poly, &ipmin, &cmmino);

    /* lasso polygon */
    ier = partition_poly(poly, npolys, polys, mtol, ALL_ONEBOUNDARY, ADJUST_LASSO, FORCE_SPLIT, OVERWRITE_ORIGINAL, np);

    /* polygon was successfully lassoed */
    if (ier == 0) {
	if (*np > 0) {
	    /* not enough polygons */
	    if (*np > npolys) return(0);

	    /* area/(2 pi) of combined smallest caps of lassoed polygon */
	    cmmint = 0.;
	    for (ip = 0; ip < *np; ip++) {
		cmminf(polys[ip], &ipmin, &cmmin);
		cmmint += cmmin;
	    }

	    /* lassoed polygons are a genuine improvement */
	    if ((*np == 1 && cmmint < cmmino) || cmmint <= .9 * cmmino) {

	    /* lassoed polygons are too large to bother with */
	    } else {
		*np = 0;

	    }
	}

    } else {
	*np = 0;
    }

    return(ier);
}
