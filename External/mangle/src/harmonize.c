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
const char *optstr = "dql:m:s:e:i:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly, npolys, nws,i;
    long double area;
    harmonic *w;
    polygon **polys;
    polys=polys_global;

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile, and Wlm_outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- harmonize ----------------\n");

    /* advise harmonic number */
    msg("maximum harmonic number %d\n", lmax);

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
    
    if (snapped==0 || balkanized==0) {
      msg("WARNING: 'snapped' and 'balkanized' keywords not found in all input files.\n");
      msg("Running harmonize on polygons that are not snapped and balkanized may give misleading results.\n");
    }

    /* allocate array containing spherical harmonics of complete mask */
    w = (harmonic *) malloc(sizeof(harmonic) * NW);
    if (!w) {
        fprintf(stderr, "harmonize: failed to allocate memory for %d harmonics\n", NW);
        exit(1);
    }

    /* spherical harmonics of region */
    npoly = harmonize_polys(npoly, polys, mtol, lmax, w);
    if (npoly == -1) exit(1);

    /* advise area */
    area = w[0][0] * 2. * sqrtl(PI);
    msg("area of (weighted) region is %.15Lg str\n", area);

    /* write polygons */
    ifile = argc - 1;
    nws = wrspher(argv[ifile], lmax, w);
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
    printf("harmonize [-d] [-q] [-l<lmax>] [-m<a>[u]] [-s<n>] [-e<n>] [-i<f>[<n>][u]] polygon_infile1 [polygon_infile2 ...] Wlm_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"
