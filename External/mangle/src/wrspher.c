/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Write spherical harmonics.

   Input: filename = name of file to write to;
		     "" or "-" means write to standard output.
	  lmax = maximum harmonic number.
	  w = array containing harmonics;
	      NW = ((lmax + 1)(lmax + 2))/ 2 is defined in harmonics.h.
  Return value: number of (complex) harmonics written,
		or -1 if error occurred.
*/
int wrspher(char *filename, int lmax, harmonic w[/*NW*/])
{
/* precision with which harmonics are written */
#define PRECISION	16
    int i, iw, width;
    FILE *file;

    /* open filename for writing */
    if (!filename || strcmp(filename, "-") == 0) {
	file = stdout;
    } else {
	file = fopen(filename, "w");
	if (!file) {
	    fprintf(stderr, "cannot open %s for writing\n", filename);
	    return(-1);
	}
    }

    /* width of each number */
    width = PRECISION + 7;

    /* write */
    fprintf(file, "%12d %12d %12d\n", lmax, IM, NW);
    for (iw = 0; iw < NW; iw++) {
	for (i = 0; i < IM; i++) {
	    fprintf(file, " %- #*.*Lg", width, PRECISION, w[iw][i]);
	}
	fprintf(file, "\n");
    }

    /* advise */
    msg("%d x %d harmonics up to lmax = %d written to %s\n",
	IM, NW, lmax, (file == stdout)? "output": filename);

    /* close file */
    if (file != stdout) fclose(file);

    return(NW);
}
