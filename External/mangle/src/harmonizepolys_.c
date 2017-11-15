/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdlib.h>
#include "manglefn.h"

/* polygons declared in rdmask_() */
extern int npolys;
extern polygon *polys[];

/*------------------------------------------------------------------------------
  Simplified fortran interface to harmonize_polys routine.
      real *8 mtol
      integer lmax
      real *8 w(IM,NW)
      call harmonizepolys(mtol, lmax, w)
*/
void harmonizepolys_(long double *mtol, int *lmax, harmonic w[])
{
    int ndone;

    ndone = harmonize_polys(npolys, polys, *mtol, *lmax, w);
    if (ndone == -1) exit(1);
}
