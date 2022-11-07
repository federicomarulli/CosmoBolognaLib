/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputfile.h"
#include "manglefn.h"

/* verbosity */
extern int verbose;

/*------------------------------------------------------------------------------
  Call a user-supplied weight function,
  or read a user-supplied file,
  to determine the weight of the region at angular position az, el in deg.
================================================================================
In c, the user-supplied function should take the form:

long double weight(long double *az, long double *el)

verbose is available as a global variable that can be declared with
extern int verbose;
or you can use the msg() routine for printing messages.
================================================================================
In fortran, the same user-supplied function should take the form
(note the trailing _ added to the function name, as is standard when
calling fortran routines from c):

      real*8 function weight_(az, el)
      real*8 az, el
================================================================================

   Input: az, el = azimuth, elevation in degrees.
	  survey = name of survey, or name of file containing list of weights.
  Output: weight at angular position az, el;
	  if survey is null, then weight is set to 1.
*/
long double weight_fn(long double az, long double el, char *survey)
{
    long double weight;

    /* null survey */
    if (!survey) {
	weight = 1.;

    /* 2dF 10k quasar survey */
    } else if (strcmpl(survey, "2QZ10k") == 0) {
	weight = twoqz_(&az, &el, &verbose);

    /* 2dF 100k galaxy survey */
    } else if (strcmpl(survey, "2dF100k") == 0) {
	weight = twodf100k_(&az, &el);

    /* 2dF 230k galaxy survey */
    } else if (strcmpl(survey, "2dF230k") == 0) {
	weight = twodf230k_(&az, &el);

    /* try reading weights from file */
    } else {
	weight = rdweight(survey);

    }

    return(weight);
}

/*------------------------------------------------------------------------------
  Try to read weight from file.
  After arbitrary header lines, the file should contain a list of weights,
  one weight in each line.

  If the file contains only one weight, then that weight will be returned
  every time.  Otherwise a new weight will be returned each time, and it is
  an error if the file does not contain enough weights.

  Input: survey = name of file containing weights.
  Return value: weight.
*/
long double rdweight(char *survey)
{
#ifndef BUFSIZE
#  define	BUFSIZE		64
#endif
    static int init = 0, nweight = 0;
    static long double weight;
    static inputfile file = {
	'\0',		/* input filename */
	0x0,		/* input file stream */
	'\0',		/* line buffer */
	BUFSIZE,	/* size of line buffer (will expand as necessary) */
	0,		/* line number */
	0		/* maximum number of characters to read (0 = no limit) */
    };

    int ird, iscan;

    /* first call */
    if (!init) {
	/* test if survey is name of file */
	file.file = fopen(survey, "r");
	/* unknown survey */
	if (!file.file) {
	    fprintf(stderr, "weight_fn: unknown file or survey %s\n", survey);
	    fprintf(stderr, "please look in weight_fn.c for the names of known surveys\n");
	    exit(1);
	}
	file.name = survey;
	/* survey appears to be a file */
	msg("will read weights for each polygon from %s\n", survey);
	init = 1;
    }

    /* every call */
    if (nweight >= 0) do {
	/* read line of data */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) exit(1);
	/* EOF */
	if (ird == 0) {
	    /* flag there's only one weight */
	    if (nweight == 1) {
		nweight = -1;
		break;
	    /* premature EOF */
	    } else {
		fprintf(stderr, "weight_fn: unexpected EOF at line %d of %s\n", file.line_number, survey);
		exit(1);
	    }
	}
	/* read contents of line into weight */
	iscan = sscanf(file.line, "%Lg", &weight);
	nweight++;
    } while (!iscan);

    return(weight);
}
