/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "manglefn.h"

#define WHERE		fprintf(stderr, "rdspher: at line %d of %s\n", line_number, fn)

/*------------------------------------------------------------------------------
  Read spherical harmonics.

   Input: filename = name of file to read from;
		     "" or "-" means read from standard input.
	  lmax = maximum harmonic number.
	  w = array containing harmonics.
  Return value: number of (complex) harmonics read,
		or -1 if error occurred.
*/
int rdspher(char *filename, int *lmax_p, harmonic **w_p)
{
    char input[] = "input";
    char *fn;
    static harmonic *w;
    int i, im, iscan, iw, lmax, nw;
    FILE *file;

    /* open filename for reading */
    if (!filename || strcmp(filename, "-") == 0) {
	file = stdin;
	fn = input;
    } else {
	file = fopen(filename, "r");
	if (!file) {
	    fprintf(stderr, "rdspher: cannot open %s for reading\n", filename);
	    return(-1);
	}
	fn = filename;
    }

    /* read number of harmonics */
    iscan = fscanf(file, "%d %d %d", &lmax, &im, &nw);
    if (iscan != 3) {
	fprintf(stderr, "rdspher: at line 1 of %s\n", fn);
	fprintf(stderr, " expecting 3 integers\n");
	return(-1);
    }
    msg("lmax = %d in %s\n", lmax, fn);

    if (lmax > *lmax_p) lmax = *lmax_p;

    /* allocate memory for array containing spherical harmonics */
    w = (harmonic *) malloc(sizeof(harmonic) * NW);
    if (!w) {
	fprintf(stderr, "rdspher: failed to allocate memory for %d harmonics\n", NW);
        return(-1);
    }

    /* zero harmonics */
    for (iw = 0; iw < NW; iw++) {
	for (i = 0; i < IM; i++) {
	    w[iw][i] = 0.;
	}
    }

    /* read harmonics */
    for (iw = 0; iw < NW; iw++) {
	for (i = 0; i < im; i++) {
	    iscan = fscanf(file, "%Lg", &w[iw][i]);
	    if (iscan != 1) {
		fprintf(stderr, "rdspher: error reading line %d of %s\n", iw + 2, fn);
		return(-1);
	    }
	}
    }

    /* advise */
    msg("harmonics up to lmax = %d read from %s\n", lmax, filename);

    /* close file */
    if (file != stdin) fclose(file);

    /* point lmax_p at maximum harmonic */
    *lmax_p = lmax;
    /* point w_p at spherical harmonics */
    *w_p = w;

    return(NW);
}
