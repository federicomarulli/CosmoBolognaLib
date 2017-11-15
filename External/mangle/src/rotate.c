/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "inputfile.h"
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqf:u:p:";

/* local functions */
void	usage(void);
int	rotate(char *, char *, format *);
void	rotate_azel(format *, azel *, azel *);

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int itr, np;

    /* parse arguments */
    parse_args(argc, argv);
    /* parse option <fopt> to -f<fopt> */
    if (fopt) itr = parse_fopt();

    /* one input and one output filename required as arguments */
    if (argc - optind != 2) {
	if (optind > 1 || argc - optind >= 1) {
	    fprintf(stderr, "%s requires 2 arguments: polygon_infile and polygon_outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- rotate ----------------\n");

    /* rotate */
    np = rotate(argv[argc - 2], argv[argc - 1], &fmt);
    if (np == -1) exit(1);

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("rotate [-d] [-q] [-f[<inframe>[,<outframe>]|<azn>,<eln>,<azp>[u]]] [-u<inunit>[,<outunit>]] [-p[+|-][<n>]] azel_infile azel_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
*/
#include "parse_fopt.c"

/*------------------------------------------------------------------------------
  Rotate az, el positions from one frame to another.
  The az, el positions are read from in_filename,
  and rotated az, el positions are written to out_filename.
  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: in_filename = name of file to read from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
  Return value: number of lines written,
		or -1 if error occurred.
*/
int rotate(char *in_filename, char *out_filename, format *fmt)
{
#ifndef BUFSIZE
#  define	BUFSIZE		64
#endif
#define AZEL_STR_LEN	32
    static inputfile file = {
	'\0',		/* input filename */
	0x0,		/* input file stream */
	'\0',		/* line buffer */
	BUFSIZE,	/* size of line buffer (will expand as necessary) */
	0,		/* line number */
	0		/* maximum number of characters to read (0 = no limit) */
    };
    char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int ird, len, np;
    long double circle;
    azel vi, vf;
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

    /* advise custom transformation */
    if (fmt->outframe == -1) {
	/* multiple transformation */
	if (strchr(fopt, ':')) {
	    msg(" -f%s is equivalent to\n", fopt);
	    msg(" -f%.16Lg,%.16Lg,%.16Lg%c\n",
		places(fmt->azn, 14), places(fmt->eln, 14), places(fmt->azp, 14), fmt->trunit);
	/* single transformation */
	} else {
	    msg("rotate -f%.16Lg,%.16Lg,%.16Lg%c\n",
		fmt->azn, fmt->eln, fmt->azp, fmt->trunit);
	}
    /* advise standard transformation */
    } else if (fmt->inframe != fmt->outframe) {
	msg("rotate from %s to %s\n",
	    frames[fmt->inframe], frames[fmt->outframe]);
	msg(" -f%s,%s is equivalent to\n",
	    frames[fmt->inframe], frames[fmt->outframe]);
        msg(" -f%.16Lg,%.16Lg,%.16Lg%c\n",
	    places(fmt->azn, 14), places(fmt->eln, 14), places(fmt->azp, 14), fmt->trunit);
    }

    /* angular units */
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

    /* write header */
    vf.az = 0.;
    wrangle(vf.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
    len = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    }
    fprintf(outfile, "%*s %*s\n", len, az_str, len, el_str);

    /* interpretive read/write loop */
    np = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &vi.az);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &vi.el);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* identity: treat specially to avoid loss of precision in scaling */
	if (fmt->inframe == fmt->outframe) {
	    /* output angles = input angles */
	    vf.az = vi.az;
	    vf.el = vi.el;

	    circle = 360.;
	    scale(&circle, 'd', fmt->inunit);

	    /* phase az */
	    switch (fmt->outphase) {
	    case '+':	if (vf.az < 0.) vf.az += circle;		break;
	    case '-':	if (vf.az > circle / 2.) vf.az -= circle;	break;
	    }

	    /* convert az and el from input to output units */
	    scale_azel(&vf, fmt->inunit, fmt->outunit);

	/* normal rotation */
	} else {
	    /* convert az and el from input units to degrees */
	    scale_azel(&vi, fmt->inunit, 'd');

	    /* rotate az and el */
	    rotate_azel(fmt, &vi, &vf);

	    /* phase az */
	    switch (fmt->outphase) {
	    case '+':	if (vf.az < 0.) vf.az += 360.;		break;
	    case '-':	if (vf.az > 180.) vf.az -= 360.;	break;
	    }

	    /* convert az and el from degrees to output units */
	    scale_azel(&vf, 'd', fmt->outunit);
	}

	/* write result */
	wrangle(vf.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
	wrangle(vf.el, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, el_str);
	fprintf(outfile, "%s %s\n", az_str, el_str);
	fflush(outfile);

        /* increment counters of results */
	np++;
    }

    if (outfile != stdout) {
	msg("rotate: %d positions written to %s\n", np, out_fn);
    }

    return(np);
}

/*------------------------------------------------------------------------------
   Interface to fortran routines that actually do the rotation.
*/
void rotate_azel(format *fmt, azel *vi, azel *vf)
{
    /* custom */
    if (fmt->outframe == -1) {
	azel_(&vi->az, &vi->el, &fmt->azn, &fmt->eln, &fmt->azp, &vf->az, &vf->el);

    /* built-ins */
    } else {
	fframe_(&fmt->inframe, &vi->az, &vi->el, &fmt->outframe, &vf->az, &vf->el);

    }
}
