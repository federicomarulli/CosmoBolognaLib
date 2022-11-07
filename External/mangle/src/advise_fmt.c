/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <string.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Advise data format.

   Input: fmt = pointer to format structure.
*/
void advise_fmt(format *fmt)
{
    int n;
    char plural[] = "s";

    /* format of input data */
    if (fmt->in) {
	if (fmt->in[strlen(fmt->in) - 1] == 's') {
	    strcpy(plural, "");
	} else {
	    strcpy(plural, "s");
	}
	if (fmt->single == 1) {
	    msg("input data format will be defined by keywords in infiles\n");
	} else {		/* fmt->single == 0 */
	    msg("(initial) input data format is %s with", fmt->in);
	    if (strcmp(fmt->in, "edges") == 0) {
		msg(" %d points/edge and", fmt->innve);
		n = fmt->nn * fmt->innve * fmt->n;
	    } else {
		n = fmt->nn * fmt->n;
	    }
	    if (fmt->n == 0) {
		msg(" variable no. of %s%s per line\n", fmt->in, plural);
	    } else {
		msg(" %d %s%s (%d numbers) per line\n", fmt->n, fmt->in, plural, n);
	    }
	}

	/* angular unit of input data */
	if (fmt->inunitp && !(strcmp(fmt->in, "polygon") == 0 || strcmp(fmt->in, "Region") == 0)) {
	    msg("will take units of any angles in input polygon files to be %c (", fmt->inunitp);
	    switch (fmt->inunitp) {
#include "angunit.h"
	    }
	    msg(")\n");
	}

    }

    /* format of output data */
    if (fmt->out) {
	msg("output data format will be %s", fmt->out);
	if (strcmp(fmt->out, "edges") == 0) msg(" with %d points/edge", fmt->outnve);
	msg("\n");
	if (strcmp(fmt->out, "vertices") == 0 || strcmp(fmt->out, "edges") == 0)
	     msg("WARNING: %s output format loses information\n", fmt->out);

	/* angular unit of output data */
	if (fmt->outunitp && !(strcmp(fmt->out, "polygon") == 0 || strcmp(fmt->out, "Region") == 0)) {
	    msg("units of angles in output polygon files will be %c (", fmt->outunitp);
	    switch (fmt->outunitp) {
#include "angunit.h"
	    }
	    msg(")\n");
	}

    }

}
