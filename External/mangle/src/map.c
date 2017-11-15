/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "inputfile.h"
#include "manglefn.h"
#include "defaults.h"

/* redefine default maximum harmonic */
#undef LMAX
#define LMAX		MAXINT

/* getopt options */
const char *optstr = "dqw:l:g:x:u:p:";

/* local functions */
void	usage(void);
#ifdef	GCC
int	map(char *, char *, format *, int lmax, harmonic w[NW], long double, long double);
#else
int	map(char *, char *, format *, int lmax, harmonic w[/*NW*/], long double, long double);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int nmap, nws;
    harmonic *w_p;

    /* lmax will be read from Wlm_filename, or set by command line */
    lmax = LMAX;

    /* parse arguments */
    parse_args(argc, argv);

    /* one input and output filename required as arguments */
    if (argc - optind != 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires 2 arguments: azel_infile, and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    /* file containing harmonics must have been specified */
    if (!Wlm_filename) {
	fprintf(stderr, "%s requires -w<Wlm_filename> option to specify file containing harmonics\n", argv[0]);
	exit(1);
    }

    msg("---------------- map ----------------\n");

    /* advise */
    if (lsmooth == 0.) {
	msg("no smoothing\n");
    } else {
	msg("smoothing harmonic number lsmooth = %Lg\n", lsmooth);
	if (esmooth != 2.) {
	    msg("smoothing exponent = %Lg\n", esmooth);
	}
    }

    /* read harmonics */
    nws = rdspher(Wlm_filename, &lmax, &w_p);
    if (nws == -1) exit(1);

    /* map */
    nmap = map(argv[argc - 2], argv[argc - 1], &fmt, lmax, w_p, lsmooth, esmooth);
    if (nmap == -1) exit(1);

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("map [-d] [-q] -w<Wlmfile> [-l<lmax>] [-g<lsmooth>] [-u<inunit>[,<outunit>]] [-p[+|-][<n>]] azel_infile outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Map.  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: in_filename = name of file to read from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  w = array containing spherical harmonics.
	  lmax = maximum harmonic number.
	  lsmooth = smoothing harmonic number (0. = no smooth).
	  esmooth = smoothing exponent (2. = gaussian).
  Return value: number of items written,
		or -1 if error occurred.
*/
int map(char *in_filename, char *out_filename, format *fmt, int lmax, harmonic w[/*NW*/], long double lsmooth, long double esmooth)
{
/* precision of map values written to file */
#define PRECISION	8
#define AZEL_STR_LEN	32
    inputfile file = {
	'\0',	/* input filename */
	0x0,	/* input file stream */
	'\0',	/* line buffer */
	64,		/* size of line buffer (will expand as necessary) */
	0,		/* line number */
	0		/* maximum number of characters to read (0 = no limit) */
    };
    char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int ird, len, mmax, nmap, width;
    long double rho;
    azel v;
    char *out_fn;
    FILE *outfile;

    /* open in_filename for reading */
    if (!in_filename || strcmp(in_filename, "-") == 0) {
	file.file = stdin;
	file.name = input;
    } else {
	file.file = fopen(in_filename, "r");
	if (!file.file) {
	    fprintf(stderr, "cannot open %s for reading\n", in_filename);
	    return(-1);
	}
	file.name = in_filename;
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

    /* advise angular units */
    msg("will take units of input az, el angles in %s to be ", file.name);
    switch (fmt->inunit) {
#include "angunit.h"
    }
    msg("\n");
    if (fmt->outunit != fmt->inunit) {
	msg("units of output az, el angles will be ");
	switch (fmt->outunit) {
#include "angunit.h"
	}
	msg("\n");
    }

    /* setting mmax = lmax ensures complete set of azimuthal harmonics */
    mmax = lmax;

    /* width of map value */
    width = PRECISION + 6;

    /* write header */
    v.az = 0.;
    wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
    len = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    }
    fprintf(outfile, "%*s %*s %*s\n", len, az_str, len, el_str, width - 4, "wrho");

    /* interpretive read/write loop */
    nmap = 0;
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
	if (ird != 1 && nmap == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &v.el);
	/* skip header */
	if (ird != 1 && nmap == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to radians */
	scale_azel(&v, fmt->inunit, 'r');

	/*
	  The entire of map.c is an interface to the next line of code.
	  Bizarre, huh?
	*/
        /* compute the value of the window function at this point */
	rho = wrho(v.az, v.el, lmax, mmax, w, lsmooth, esmooth);

	/* convert az and el from radians to output units */
	scale_azel(&v, 'r', fmt->outunit);

	/* write result */
	wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
	wrangle(v.el, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, el_str);
	fprintf(outfile, "%s %s %- #*.*Lg\n", az_str, el_str, width, PRECISION, rho);
	fflush(outfile);

        /* increment counter of results */
	nmap++;
    }

    if (outfile != stdout) {
	msg("map: %d values written to %s\n", nmap, out_fn);
    }

    return(nmap);
}
