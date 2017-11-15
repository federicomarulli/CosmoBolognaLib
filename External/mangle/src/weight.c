/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqz:m:s:e:v:p:i:o:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);
#ifdef	GCC
int	weight(int npoly, polygon *[npoly], char *);
#else
int	weight(int npoly, polygon *[/*npoly*/], char *);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly, npolys,i;
    polygon **polys;
    polys=polys_global;

    /* default output format */
    fmt.out = keywords[POLYGON];

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile and polygon_outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    /* survey must have been specified */
    if (!survey) {
	fprintf(stderr, "%s requires -z<survey> option to specify the name of a survey,\n", argv[0]);
	fprintf(stderr, "%*s or the name of a file containing a list of weights.\n", (int)strlen(argv[0]), "");
	fprintf(stderr, "Please look in weight_fn.c for the names of known surveys.\n");
	exit(1);
    }

    msg("---------------- weight ----------------\n");

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

    /* weight polygons */
    npoly = weight(npoly, polys, survey);

    ifile = argc - 1;
    npoly = wrmask(argv[ifile], &fmt, npoly, polys);
    if (npoly == -1) exit(1);

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
    printf("weight [-d] [-q] -z<survey> [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn|-vp] [-p[+|-][<n>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Weight polygons.

   Input: poly = array of pointers to polygons.
	  npoly = pointer to number of polygons.
	  survey = name of survey, or of filename containing list of weights.
  Output: polys = array of pointers to polygons;
  Return value: number of polygons weighted,
		or -1 if error occurred.
*/
int weight(int npoly, polygon *poly[/*npoly*/], char *survey)
{
    const int per = 0;
    const int nve = 2;

    int do_vcirc, i, imid, ipoly, iverts, ivm, nev, nev0, nomid, nv, nvm, nzero;
    int *ipv, *gp, *ev;
    long double tol;
    long double *angle;
    vec *ve, *vm;
    azel v;

    nomid = 0;
    nzero = 0;
    for (ipoly = 0; ipoly < npoly; ipoly++) {
	/* vertices of polygon */
	do_vcirc = 0;
	tol = mtol;
	iverts = gverts(poly[ipoly], do_vcirc, &tol, per, nve, &nv, &ve, &angle, &ipv, &gp, &nev, &nev0, &ev);
	if (iverts != 0) return(-1);
	/* point somewhere in the middle of the polygon */
	imid = vmid(poly[ipoly], tol, nv, nve, ve, ipv, ev, &nvm, &vm);
	if (imid == -1) return(-1);
	/* check found a point inside the polygon */
	imid = 0;
	for (ivm = 0; ivm < nvm; ivm++) {
	    if (vm[ivm][0] != 0. || vm[ivm][1] != 0. || vm[ivm][2] != 0.) {
		imid = 1;
		if (ivm > 0) for (i = 0; i < 3; i++) vm[0][i] = vm[ivm][i];
		break;
	    }
	}
	/* found a point */
        if (imid == 1) {
	    /* convert unit vector to az, el */
	    rp_to_azel(*vm, &v);
	    /* scale angles from radians to degrees */
	    scale_azel(&v, 'r', 'd');
	    /* weight at that point */
	    poly[ipoly]->weight = weight_fn(v.az, v.el, survey);
	    if (poly[ipoly]->weight == 0.) nzero++;
	/* failed to find a point */
	} else {
	  //call weight_fn to stay in right place in weight file if reading from weights from a file
	    weight_fn(v.az, v.el, survey);
	    if (nomid == 0) msg("weight: failed to find interior point for the following polygons:\n");
	    msg(" %d", (fmt.newid == 'n')? ipoly : poly[ipoly]->id);
	    nomid++;
	}
    }
    if (nomid > 0) msg("\n");
    if (nomid > 0) {
	msg("weight: failed to find interior point for %d polygons\n", nomid);
        msg("FAILURE TO FIND INTERIOR POINT PROBABLY MEANS YOU HAVE A WEIRD-SHAPED POLYGON.\n");
	msg("PLEASE FILL IN THE CORRECT WEIGHT BY HAND.\n");
    }

    /* assign new polygon id numbers in place of inherited ids */
    if (fmt.newid == 'n') {
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    poly[ipoly]->id = ipoly;
	}
    }

    if (fmt.newid == 'p') {
      for (ipoly = 0; ipoly < npoly; ipoly++) {
	poly[ipoly]->id = poly[ipoly]->pixel;
      }
    }

    /* warn about zero weights */
    if (nzero > 0) msg("weight: %d polygons have zero weight\n", nzero);

    return(npoly);
}
