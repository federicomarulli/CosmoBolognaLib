/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef TIME
#include <time.h>
#endif
#include "inputfile.h"
#include "manglefn.h"
#include "defaults.h"

/* redefine default angular unit for output DR angles to radians */
#undef OUTUNIT
#define OUTUNIT		'r'

/* getopt options */
const char *optstr = "dqm:hs:e:u:p:i:";

/* allocate polygons as a global array */
polygon *poly_global[NPOLYSMAX];

/* declared in rdmask */
extern inputfile file;

/* local functions */
void	usage(void);
#ifdef	GCC
int	drangle(char *, char *, char *, format *, int npoly, polygon *[npoly]);
#else
int	drangle(char *, char *, char *, format *, int npoly, polygon *[/*npoly*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    char *th_in_filename;
    int nfiles, np, npoly,i;
    polygon **poly;
    poly=poly_global;

    /* set angular unit for output DR angles to default */
    fmt.outunit = OUTUNIT;

    /* parse arguments */
    parse_args(argc, argv);

    /* two or three input and one output filename required as arguments */
    nfiles = argc - optind;
    if (!(nfiles == 3 || nfiles == 4)) {
	if (optind > 1 || nfiles >= 1) {
	    fprintf(stderr, "%s requires 3 arguments: polygon_infile, azel+th_infile, and dr_outfile\n", argv[0]);
	    fprintf(stderr, "%*s       or 4 arguments: polygon_infile, azel_infile, th_infile, and dr_outfile\n", (int)strlen(argv[0]), " ");
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    /* summary option only possible with 3 input and 1 output files */
    if (summary && nfiles != 4) {
	fprintf(stderr, "%s with summary option -h requires 3 input files and 1 output file:\n", argv[0]);
	fprintf(stderr, "%*s polygon_infile azel_infile th_infile dr_outfile\n", (int)strlen(argv[0]), " ");
	usage();
	exit(1);
    }

    msg("---------------- drangle ----------------\n");

    /* advise data format */
    advise_fmt(&fmt);


    /* read polygons */
    npoly = rdmask(argv[optind], &fmt, NPOLYSMAX, poly);
    if (npoly == -1) exit(1);
    if (npoly == 0) {
	msg("STOP\n");
	exit(0);
    }
    if (snapped==0 || balkanized==0) {
      fprintf(stderr, "Error: input polygons must be snapped and balkanized before using drangle.\n");
      fprintf(stderr, "If your polygons are already snapped and balkanized, add the 'snapped' and\n'balkanized' keywords at the beginning of each of your input polygon files.\n");
      exit(1);
    }

    /* name of file containing angular radii th, if present */
    th_in_filename = (nfiles == 3)? 0x0 : argv[optind + 2];

    /* angles */
    np = drangle(argv[optind + 1], th_in_filename, argv[argc - 1], &fmt, npoly, poly);
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
    printf("drangle [-d] [-q] [-h] [-m<a>[u]] [-s<n>] [-e<n>] [-u<inunit>[,<outunit>]] [-p[+|-][<n>]] [-i<f>[<n>][u]] polygon_infile azel[th]_infile [th_infile] dr_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Angles within mask along circles centred at az el, radii th.

  Anglar positions az, el are read from azel_in_filename,
  angular radii th are read from th_in_filename,
  or from azel_in_filename if th_in_filename is null,
  and the results are written to out_filename.

  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: azel_in_filename = name of file to read az, el from;
			"" or "-" means read from standard input.
	  th_in_filename = name of file to read th from;
			null means read from azel_in_filename;
			"-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  npoly = number of polygons in poly array.
	  poly = array of pointers to polygons.
  Return value: number of lines written,
		or -1 if error occurred.
*/
int drangle(char *azel_in_filename, char *th_in_filename, char *out_filename, format *fmt, int npoly, polygon *poly[/*npoly*/])
{
#define AZEL_STR_LEN	32
    char input[] = "input", output[] = "output";
    /* maximum number of angular angular radii: will expand as necessary */
    static int nthmax = 0;
    static int ndrmax = 0;
    static long double *th = 0x0, *cm = 0x0, *drsum = 0x0;
    static long double *dr = 0x0;

#ifdef TIME
    clock_t time;
#endif
    char inunit, outunit;
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN], th_str[AZEL_STR_LEN], dr_str[AZEL_STR_LEN];
    int ier, ird, ith, len, lenth, np, nt, nth;
    long double rp[3], s, t;
    azel v;
    char *out_fn;
    FILE *outfile;

    inunit = (fmt->inunit == 'h')? 'd' : fmt->inunit;

    /* read angular radii from th_in_filename */
    if (th_in_filename) {

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

	msg("will take units of input th angles in %s to be ", file.name);
	switch (inunit) {
#include "angunit.h"
	}
	msg("\n");

	/* read angular radii from th_in_filename */
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
			fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", nthmax);
			return(-1);
		    }
		}
		/* copy angular radius into th array */
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

	/* (re)allocate memory for th array */
	th = (long double *) realloc(th, sizeof(long double) * nth);
	if (!th) {
	    fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", nth);
	    return(-1);
	}
	/* (re)allocate memory for cm array */
	cm = (long double *) realloc(cm, sizeof(long double) * nth);
	if (!cm) {
	    fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", nth);
	    return(-1);
	}
	/* (re)allocate memory for drsum array */
	drsum = (long double *) realloc(drsum, sizeof(long double) * nth);
	if (!drsum) {
	    fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", nth);
	    return(-1);
	}
	nthmax = nth;
    }

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

    /* advise input angular units */
    if (th_in_filename) {
	msg("will take units of input az, el angles in %s to be ", file.name);
    } else {
	msg("will take units of input az, el, and th angles in %s to be ", file.name);
    }
    switch (fmt->inunit) {
#include "angunit.h"
    }
    if (!th_in_filename && fmt->inunit == 'h') {
	msg(", th in deg");
    }
    msg("\n");

    /* advise output angular unit */
    outunit = (fmt->outunit == 'h')? 'd' : fmt->outunit;
    msg("units of output DR angles will be ");
    switch (outunit) {
#include "angunit.h"
    }
    msg("\n");

#ifdef TIME
    printf("timing ... ");
    fflush(stdout);
    time = clock();
#endif

    /* length of az, el and th, dr strings to be written to output */
    t = 0.;
    wrangle(t, fmt->inunit, fmt->outprecision, AZEL_STR_LEN, az_str);
    len = strlen(az_str);
    t = 0.;
    wrangle(t, inunit, fmt->outprecision, AZEL_STR_LEN, th_str);
    lenth = strlen(th_str);

    /* write header */
    if (!summary) {
	if (fmt->inunit == 'h') {
	    sprintf(az_str, "az(hms)");
	    sprintf(el_str, "el(dms)");
	    sprintf(th_str, "th(d):");
	} else {
	    sprintf(az_str, "az(%c)", fmt->inunit);
	    sprintf(el_str, "el(%c)", fmt->inunit);
	    sprintf(th_str, "th(%c):", fmt->inunit);
	}
	if (th_in_filename) {
	    fprintf(outfile, "%*s %*s %*s\n", len, " ", len, " ", lenth, th_str);
	}
	fprintf(outfile, "%*s %*s", len, az_str, len, el_str);
	if (th_in_filename) {
	    for (ith = 0; ith < nth; ith++) {
		wrangle(th[ith], inunit, fmt->outprecision, AZEL_STR_LEN, th_str);
		fprintf(outfile, " %s", th_str);
	    }
	    fprintf(outfile, "\n");
	} else {
	    fprintf(outfile, " %*s\n", lenth, th_str);
	}
    }

    /* initialize th and drsum */
    if (th_in_filename) {
	/* convert th from input units to radians */
	for (ith = 0; ith < nth; ith++) {
	    scale(&th[ith], inunit, 'r');
	}

	/* initialize sum of dr to zero */
	for (ith = 0; ith < nth; ith++) {
	    drsum[ith] = 0.;
	}
    }

    /* interpretive read/write loop */
    np = 0;
    nt = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &v.az);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &v.el);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to radians */
	scale_azel(&v, fmt->inunit, 'r');

	/* read th */
	if (!th_in_filename) {
	    nth = 0;
	    while (1) {
		word = next;
		ird = rdangle(word, &next, inunit, &t);
		/* done */
		if (ird < 1) break;
		/* ok */
		if (nth >= nthmax) {
		    if (nthmax == 0) {
			nthmax = 64;
		    } else {
			nthmax *= 2;
		    }
		    th = (long double *) realloc(th, sizeof(long double) * nthmax);
		    if (!th) {
			fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", nthmax);
			return(-1);
		    }
		    /* (re)allocate memory for cm array */
		    cm = (long double *) realloc(cm, sizeof(long double) * nthmax);
		    if (!cm) {
			fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", nthmax);
			return(-1);
		    }
		}
		th[nth] = t;
		/* increment count of angular radii */
		nth++;
	    }
	    /* convert th from input units to radians */
	    for (ith = 0; ith < nth; ith++) {
		scale(&th[nth], inunit, 'r');
	    }
	}

	/* allocate memory for dr */
	if (nth > ndrmax) {
	    ndrmax = nth;
	    dr = (long double *) realloc(dr, sizeof(long double) * ndrmax);
	    if (!dr) {
		fprintf(stderr, "drangle: failed to allocate memory for %d long doubles\n", ndrmax);
		return(-1);
	    }
	}

	/* unit vector corresponding to angular position az, el */
	azel_to_rp(&v, rp);

	/* limiting cm = 1-cosl(th) values to each polygon */
	ier = cmlim_polys(npoly, poly, mtol, rp);
	if (ier == -1) return(-1);

	/* cm = 1 - cosl(th) */
	for (ith = 0; ith < nth; ith++) {
	    s = sinl(th[ith] / 2.);
	    cm[ith] = 2. * s * s;
	}

	/* angles about rp direction at radii th */
	ier = drangle_polys(npoly, poly, mtol, rp, nth, cm, dr);
	if (ier == -1) return(-1);

	/* sum of dr at radii th */
	if (th_in_filename) {
	    for (ith = 0; ith < nth; ith++) {
		drsum[ith] += dr[ith];
	    }
	}

	/* convert az and el from radians to original input units */
	scale_azel(&v, 'r', fmt->inunit);

	/* write result */
	if (!summary) {
	    wrangle(v.az, fmt->inunit, fmt->outprecision, AZEL_STR_LEN, az_str);
	    wrangle(v.el, fmt->inunit, fmt->outprecision, AZEL_STR_LEN, el_str);
	    fprintf(outfile, "%s %s", az_str, el_str);
	    for (ith = 0; ith < nth; ith++) {
		scale(&dr[ith], 'r', outunit);
		wrangle(dr[ith], outunit, fmt->outprecision, AZEL_STR_LEN, dr_str);
		fprintf(outfile, " %s", dr_str);
	    }
	    fprintf(outfile, "\n");
	    fflush(outfile);
	}

        /* increment counters of results */
	np++;
	nt += nth;

	/* warn about a potentially huge output file */
	if (np == 100 && !summary && th_in_filename && outfile != stdout) {
	    msg("hmm, looks like %s could grow pretty large ...\n", out_fn);
	    msg("try using the -h switch if you only want a summary in the output file\n");
	}
    }

    /* write sum of dr */
    if (summary) {
	sprintf(th_str, "th(%c)", inunit);
	sprintf(dr_str, "DR(%c)", outunit);
	fprintf(outfile, "%*s %*s\n", lenth, th_str, lenth, dr_str);
	for (ith = 0; ith < nth; ith++) {
	    scale(&th[ith], 'r', inunit);
	    wrangle(th[ith], inunit, fmt->outprecision, AZEL_STR_LEN, th_str);
	    scale(&drsum[ith], 'r', outunit);
	    wrangle(drsum[ith] / (long double)np, outunit, fmt->outprecision, AZEL_STR_LEN, dr_str);
	    fprintf(outfile, "%s %s\n", th_str, dr_str);
	}
	fflush(outfile);
    } else if (th_in_filename) {
	fprintf(outfile, "%*s", 2 * len + 1, "Average:");
	for (ith = 0; ith < nth; ith++) {
	    scale(&drsum[ith], 'r', outunit);
	    wrangle(drsum[ith] / (long double)np, outunit, fmt->outprecision, AZEL_STR_LEN, dr_str);
	    fprintf(outfile, " %s", dr_str);
	}
	fprintf(outfile, "\n");
	fflush(outfile);
    }

    /* close azel_in_filename */
    if (file.file != stdin) {
	fclose(file.file);
    }

#ifdef TIME
    time = clock() - time;
    printf("done in %Lg sec\n", (float)time / (float)CLOCKS_PER_SEC);
#endif

    /* advise */
    if (outfile != stdout) {
	fclose(outfile);
	if (summary) {
	    msg("drangle: header + %d lines written to %s\n", nth, out_fn);
	} else {
	    if (th_in_filename) {
		msg("drangle: %d x %d = ", nth, np);
	    } else {
		msg("drangle: total of ");
	    }
	    msg("%d angles at %d positions written to %s\n", nt, np, out_fn);
	}
    }

    return(nt);
}
