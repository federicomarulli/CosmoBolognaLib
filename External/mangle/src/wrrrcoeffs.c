/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Write first few coefficients rrcoeff[i]
  of the power series expansion in sinl(th/2)
  of the self-correlation of mask W at angular separation th

              //
     <W W> =  || W(n) W(n') delta(n.n' - cos th) do do'
              //
				       i
	   = sum rrcoeff[i] [sinl(th/2)]
	      i

   Input: filename = name of file to write to;
		     "" or "-" means write to standard output.
	  area = area of mask in steradians.
	  bound, vert = coefficients computed by subroutine gspher.s.f for each polygon,
			with corrections from gphbv.s.f for polygons that abut along a finite boundary,
			and further corrections NOT YET IMPLEMENTED for polygons that touch at isolated points.
  Return value: number of coefficients written,
		or -1 if error occurred.
*/
int wrrrcoeffs(char *filename, long double area, long double bound[2], long double vert[2])
{
/* precision with which harmonics are written */
#define PRECISION	16
/* number of coefficients */
#define NCOEFF		4
    int icoeff, width;
    long double rrcoeff[NCOEFF];
    FILE *file;

    /* first 4 coefficients of power series expansion of <RR> */
    rrcoeff[0] = area * 2. * PI;
    rrcoeff[1] = - bound[0] * 4.;
    rrcoeff[2] = vert[0] * 2.;
    rrcoeff[3] = bound[1] * 2./3. + vert[1] * 8./9.;

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
    fprintf(file, "%*s\n", PRECISION, "WWcoeffs");
    for (icoeff = 0; icoeff < NCOEFF; icoeff++) {
	fprintf(file, "%- #*.*Lg\n", width, PRECISION, rrcoeff[icoeff]);
    }

    /* advise */
    msg("%d coefficients of series expansion of <WW> written to %s\n",
	NCOEFF, (file == stdout)? "output": filename);

    /* close file */
    if (file != stdout) fclose(file);

    return(NCOEFF);
}
