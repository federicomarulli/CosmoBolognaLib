/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqm:s:e:i:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);
#ifdef	GCC
int	rrcoeffs(int npoly, polygon *[npoly], long double *, long double [2], long double[2]);
#else
int	rrcoeffs(int npoly, polygon *[/*npoly*/], long double *, long double [2], long double[2]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly, npolys, nws,i;
    long double area, bound[2], vert[2];
    polygon **polys;
    polys=polys_global;

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile, and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- rrcoeffs ----------------\n");

    msg("WARNING: coefficients 2 and 3 may be incorrect because the contribution from point abutments between polygons is not yet implemented.\n");
    msg("However, coefficients 0 and 1 are good.\n");

    /* tolerance angle for multiple intersections */
    if (mtol != 0.) {
	scale(&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %Lg%c will be treated as coincident\n", mtol, munit);
	scale(&mtol, munit, 'r');
	munit = 'r';
    }

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 1 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &polys[npoly]);
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

    /* area, bound, and vert of region */
    npoly = rrcoeffs(npoly, polys, &area, bound, vert);
    if (npoly == -1) exit(1);

    /* advise area */
    msg("area of (weighted^2) region is %.15Lg str\n", area);

    /* write polygons */
    ifile = argc - 1;
    nws = wrrrcoeffs(argv[ifile], area, bound, vert);
    if (nws == -1) exit(1);

    for(i=0;i<npoly;i++){
      free_poly(polys[i]);
    }

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("rrcoeffs [-d] [-q] [-l<n>] [-m<a>[u]] [-s<n>] [-e<n>] [-i<f>[<n>][u]] polygon_infile1 [polygon_infile2 ...] outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Coefficients of series expansion of correlation <WW> at angular separation th
  <WW> = 2 pi area
	 - 4 bound[0] sinl(th/2)
	 + 2 vert[0] sin^2(th/2)
	 + (2/3 bound[1] + 8/9 vert[1]) sin^3(th/2) + ...

  Return value: number of polygons for which area, bound, and vert were computed.
		    or -1 if error occurred.
*/
int rrcoeffs(int npoly, polygon *poly[/*npoly*/], long double *area, long double bound[2], long double vert[2])
{
    int ier, ip, ipoly, jp, jpoly, ndone, ner, np;
    long double darea, dbound[2], dvert[2], toli, tolj, ww;
    harmonic dw[1];
    polygon *polyij = 0x0;

    /* initialize area, bound, and vert to zero */
    *area = 0.;
    bound[0] = 0.;
    bound[1] = 0.;
    vert[0] = 0.;
    vert[1] = 0.;

    ndone = 0;
    ner = 0;
    lmax = 0;

    /* each polygon */
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	ww = poly[ipoly]->weight * poly[ipoly]->weight;

	/* zero weight polygon requires no computation */
	if (ww == 0.) {
	    ndone++;
	    continue;

	/* contribution to correlation from self-correlation of polygons */
	} else {
	    /* compute area, bound, and vert */
	    toli = mtol;
	    ier = gspher(poly[ipoly], lmax, &toli, &darea, dbound, dvert, dw);
	    if (ier == -1) return(-1);

	}

	/* computation failed */
	if (ier) {
	    ner++;
	    fprintf(stderr, "rrcoeffs: computation failed for polygon %d; discard it\n", ipoly);

	/* success */
	} else {
	    ndone++;
	    /* increment area, bound, and vert */
	    *area += darea * ww;
	    bound[0] += dbound[0] * ww;
	    bound[1] += dbound[1] * ww;
	    vert[0] += dvert[0] * ww;
	    vert[1] += dvert[1] * ww;

	}

	/* contribution to correlation from abutting polygons */
	for (jpoly = 0; jpoly < ipoly; jpoly++) {
	    ww = - 2. * poly[ipoly]->weight * poly[jpoly]->weight;
	    if (ww == 0.) continue;

	    /* look for boundary dividing poly[ipoly] and poly[jpoly] */
	    for (ip = 0; ip < poly[ipoly]->np; ip++) {

		for (jp = 0; jp < poly[jpoly]->np; jp++) {

		    /* poly[ipoly] and poly[jpoly] abut along a common boundary */
		    if (poly[ipoly]->cm[ip] == - poly[jpoly]->cm[jp]
			&& poly[ipoly]->rp[ip][0] == poly[jpoly]->rp[jp][0]
			&& poly[ipoly]->rp[ip][1] == poly[jpoly]->rp[jp][1]
			&& poly[ipoly]->rp[ip][2] == poly[jpoly]->rp[jp][2]) {

			/* make sure polyij contains enough space for intersection */
			np = poly[ipoly]->np + poly[jpoly]->np;
			ier = room_poly(&polyij, np, DNP, 0);
			if (ier == -1) {
			    fprintf(stderr, "rrcoeffs: failed to allocate memory for polygon of %d caps\n", np + DNP);
			    return(-1);
			}

			/* make polygon which is the intersection of the 2 polygons */
			poly_poly(poly[ipoly], poly[jpoly], polyij);

			/* suppress abutting boundary from poly[jpoly] */
			polyij->cm[poly[ipoly]->np + jp] = 2.;

			/* compute bound and vert */
			tolj = toli;
			ier = gphbv(polyij, poly[ipoly]->np, ip, &tolj, dbound, dvert);
			if (ier == -1) return(-1);

			/* increment bound and vert */
			bound[0] += dbound[0] * ww;
			bound[1] += dbound[1] * ww;
			vert[0] += dvert[0] * ww;
			vert[1] += dvert[1] * ww;

		    }

		}

	    }

	}

    }

    /* advise */
    if (ner > 0) {
	msg("discarded %d polygons for which computations failed\n", ner);
    }
    msg("area, bound, and vert accumulated for %d weighted polygons\n", ndone);

    return(ndone);
}
