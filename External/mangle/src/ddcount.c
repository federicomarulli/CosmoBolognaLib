/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef TIME
#include <time.h>
#endif
#include "inputfile.h"
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqs:e:u:p:i:";

/* allocate polygons as a global array */
polygon *poly_global[NPOLYSMAX];

/* declared in rdmask */
extern inputfile file;

/* local functions */
void	usage(void);
#ifdef	GCC
long	ddcount(char *, char *, char *, format *, int npoly, polygon *[npoly]);
#else
long	ddcount(char *, char *, char *, format *, int npoly, polygon *[/*npoly*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly,i;
    long np;
    polygon **poly;
    poly=poly_global;

    /* parse arguments */
    parse_args(argc, argv);

    /* three input and one output filename required as arguments */
    nfiles = argc - optind;
    if (nfiles != 4) {
	if (optind > 1 || nfiles >= 1) {
	    fprintf(stderr, "%s requires 4 arguments: polygon_infile, azel_infile, th_infile, and dd_outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- ddcount ----------------\n");

    /* advise data format */
    advise_fmt(&fmt);


    /* read polygons */
    ifile = optind;
    npoly = rdmask(argv[optind], &fmt, NPOLYSMAX, poly);
    if (npoly == -1) exit(1);
    if (npoly == 0) {
	msg("STOP\n");
	exit(0);
    }

    /* pair counts */
    np = ddcount(argv[optind + 1], argv[optind + 2], argv[optind + 3], &fmt, npoly, poly);
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
    printf("ddcount [-d] [-q] [-s<n>] [-e<n>] [-u<inunit>] [-p[+|-][<n>]] [-i<f>[<n>][u]] polygon_infile azel_infile th_infile dd_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Counts of pairs in bins bounded by radii th, centred at az el.

  Anglar positions az, el are read from azel_in_filename,
  angular radii th are read from th_in_filename,
  or from azel_in_filename if th_in_filename is null,
  and the results are written to out_filename.

  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: azel_in_filename = name of file to read az, el from;
			"" or "-" means read from standard input.
	  th_in_filename = name of file to read th from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  npoly = number of polygons in poly array.
	  poly = array of pointers to polygons.
  Return value: number of distinct pairs counted,
		or -1 if error occurred.
*/
long ddcount(char *azel_in_filename, char *th_in_filename, char *out_filename, format *fmt, int npoly, polygon *poly[/*npoly*/])
{
#define AZEL_STR_LEN	32
    char input[] = "input", output[] = "output";
    /* maximum number of angular angular radii: will expand as necessary */
    static int nthmax = 0;
    /* maximum number of az-el points: will expand as necessary */
    static int nazelmax = 0;
    static int *dd = 0x0, *id = 0x0, *iord = 0x0;
    static long double *cm = 0x0, *th = 0x0;
    static azel *v = 0x0;
    static vec *rp = 0x0;

#ifdef TIME
    clock_t time;
#endif
    char inunit;
    char *word, *next;
    char th_str[AZEL_STR_LEN];
    int i, iazel, idi, ird, ith, j, jazel, manyid, nazel, nid, noid, nth;
    int *id_p;
    long np;
    long double az, cmm, el, s, t;
    char *out_fn;
    FILE *outfile;

    /* open th_in_filename for reading */
    if (strcmp(th_in_filename, "-") == 0) {
	file.file = stdin;
	file.name = input;
    } else {
	file.file = fopen(th_in_filename, "r");
	if (!file.file) {
	    fprintf(stderr, "cannot open %s for reading\n", th_in_filename);
	    return(-1);
	}
	file.name = th_in_filename;
    }
    file.line_number = 0;

    inunit = (fmt->inunit == 'h')? 'd' : fmt->inunit;
    msg("will take units of input th angles in %s to be ", file.name);
    switch (inunit) {
#include "angunit.h"
    }
    msg("\n");

    /* read angular radii th from th_in_filename */
    nth = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;
	/* read angular radius from line */
	ird = rdangle(file.line, &next, inunit, &t);
	/* error */
	if (ird < 1) {
	    /* retry if nothing read, otherwise break */
	    if (nth > 0) break;
	/* ok */
	} else if (ird == 1) {
	    if (nth >= nthmax) {
		if (nthmax == 0) {
		    nthmax = 64;
		} else {
		    nthmax *= 2;
		}
		/* (re)allocate memory for th array */
		th = (long double *) realloc(th, sizeof(long double) * nthmax);
		if (!th) {
		    fprintf(stderr, "ddcount: failed to allocate memory for %d long doubles\n", nthmax);
		    return(-1);
		}
	    }
	    /* store th */
	    th[nth] = t;
	    nth++;
	}
    }

    if (file.file != stdin) {
	/* close th_in_filename */
	fclose(file.file);
	/* advise */
	msg("%d angular radii read from %s\n", nth, file.name);
    }

    if (nth == 0) return(nth);

    /* open azel_in_filename for reading */
    if (!azel_in_filename || strcmp(azel_in_filename, "-") == 0) {
	file.file = stdin;
	file.name = input;
    } else {
	file.file = fopen(azel_in_filename, "r");
	if (!file.file) {
	    fprintf(stderr, "cannot open %s for reading\n", azel_in_filename);
	    return(-1);
	}
	file.name = azel_in_filename;
    }
    file.line_number = 0;

    /* advise input angular units */
    msg("will take units of input az, el angles in %s to be ", file.name);
    switch (fmt->inunit) {
#include "angunit.h"
    }
    msg("\n");

    /* read angular positions az, el from azel_in_filename */
    nazel = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &az);
	/* skip header */
	if (ird != 1 && nazel == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &el);
	/* skip header */
	if (ird != 1 && nazel == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* (re)allocate memory for array of az-el points */
	if (nazel >= nazelmax) {
	    if (nazelmax == 0) {
		nazelmax = 64;
	    } else {
		nazelmax *= 2;
	    }
	    v = (azel *) realloc(v, sizeof(azel) * nazelmax);
	    if (!v) {
		fprintf(stderr, "ddcount: failed to allocate memory for %d az-el points\n", nazelmax);
		return(-1);
	    }
	}

	/* record az-el */
	v[nazel].az = az;
	v[nazel].el = el;

        /* increment number of az-el points */
	nazel++;
    }

    if (file.file != stdin) {
	/* close azel_in_filename */
	fclose(file.file);
	/* advise */
	msg("%d angular positions az, el read from %s\n", nazel, file.name);
    }

    if (nazel == 0) return(nazel);

    /* open out_filename for writing */
    if (!out_filename || strcmp(out_filename, "-") == 0) {
	outfile = stdout;
	out_fn = output;
    } else {
	outfile = fopen(out_filename, "w");
	if (!outfile) {
	    fprintf(stderr, "cannot open %s for writing\n", out_filename);
	    return(-1);
	}
	out_fn = out_filename;
    }

    /* (re)allocate memory */
    cm = (long double *) realloc(dd, sizeof(long double) * nth);
    if (!cm) {
	fprintf(stderr, "ddcount: failed to allocate memory for %d long doubles\n", nth);
	return(-1);
    }
    dd = (int *) realloc(dd, sizeof(int) * nth);
    if (!dd) {
	fprintf(stderr, "ddcount: failed to allocate memory for %d ints\n", nth);
	return(-1);
    }
    id = (int *) realloc(id, sizeof(int) * nazel);
    if (!id) {
	fprintf(stderr, "ddcount: failed to allocate memory for %d ints\n", nazel);
	return(-1);
    }
    rp = (vec *) realloc(rp, sizeof(vec) * nazel);
    if (!rp) {
	fprintf(stderr, "ddcount: failed to allocate memory for %d unit vectors\n", nazel);
	return(-1);
    }
    iord = (int *) realloc(iord, sizeof(int) * nazel);
    if (!iord) {
	fprintf(stderr, "ddcount: failed to allocate memory for %d ints\n", nazel);
	return(-1);
    }

    /* store 1 - cosl(th) in cm array */
    for (ith = 0; ith < nth; ith++) {
	t = th[ith];
	scale(&t, inunit, 'r');
	s = sinl(t / 2.);
	cm[ith] = 2. * s * s;
    }

    /* convert az-el angles to radians */
    for (iazel = 0; iazel < nazel; iazel++) {
	scale_azel(&v[iazel], fmt->inunit, 'r');
    }

    /* convert az-el points to unit vectors */
    for (iazel = 0; iazel < nazel; iazel++) {
	azel_to_rp(&v[iazel], rp[iazel]);
    }

    /* polygon id number(s) of az-el points */
    msg("figuring polygon id number(s) of each az-el point ...");
    noid = 0;
    manyid = 0;
    for (iazel = 0; iazel < nazel; iazel++) {
	nid = poly_id(npoly, poly, v[iazel].az, v[iazel].el, &id_p);
	if (nid == 0) {
	    noid++;
	} else if (nid > 1) {
	    manyid++;
	}
	/* store first polygon id of point */
	if (nid == 0) {
	    id[iazel] = -1;
	} else {
	    id[iazel] = id_p[0];
	}
    }
    msg(" done\n");
    if (noid > 0) {
	msg("%d az-el points lie outside the angular mask: discard them\n", noid);
    }
    if (manyid > 0) {
	msg("%d az-el points lie inside more than one polygon: use only first polygon\n", manyid);
    }

    /* order az-el points in increasing order of polygon id */
    finibot(id, nazel, iord, nazel);

    /* write header */
    fprintf(outfile, "th(%c):", inunit);
    for (ith = 0; ith < nth; ith++) {
	wrangle(th[ith], inunit, fmt->outprecision, AZEL_STR_LEN, th_str);
	fprintf(outfile, "\t%s", th_str);
    }
    fprintf(outfile, "\n");

    msg("counting pairs ...");

#ifdef TIME
    printf(" timing ... ");
    fflush(stdout);
    time = clock();
#endif

    /* loop over az-el points */
    idi = -1;
    nid = 0;
    np = 0;
    for (iazel = 0; iazel < nazel; iazel++) {
	i = iord[iazel];
	/* skip points outside mask */
	if (id[i] == -1) continue;
	/* new polygon */
	if (id[i] != idi) {
	    idi = id[i];
	    /* reset pair counts to zero */
	    for (ith = 0; ith < nth; ith++) dd[ith] = 0;
	}
//printf(" %d", idi);
	/* az-el neighbours within same polygon */
	for (jazel = iazel + 1; jazel < nazel; jazel++) {
	    j = iord[jazel];
	    /* exit at new polygon */
	    if (id[j] != idi) break;
	    /* 1 - cosl(th_ij) */
	    cmm = cmij(rp[i], rp[j]);
	    /* ith such that cm[ith-1] <= cmm < cm[ith] */
	    ith = search(nth, cm, cmm);
	    /* increment count in this bin */
	    if (ith < nth) dd[ith]++;
	    /* increment total pair count */
	    np++;
	}
	/* write counts for this polygon */
	if (iazel + 1 == nazel || id[iord[iazel + 1]] != idi) {
	    fprintf(outfile, "%d", idi);
	    for (ith = 0; ith < nth; ith++) {
		fprintf(outfile, "\t%d", dd[ith]);
	    }
	    fprintf(outfile, "\n");
	    fflush(outfile);
	    nid++;
//printf("\n");
	}
    }

#ifdef TIME
    time = clock() - time;
    printf("done in %Lg sec\n", (float)time / (float)CLOCKS_PER_SEC);
#else
    msg("\n");
#endif

    /* advise */
    if (outfile != stdout) {
	fclose(outfile);
	msg("%d distinct pairs in %d th-bins x %d polygons written to %s\n", np, nth, nid, out_fn);
    }

    return(np);
}
