/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdlib.h>
#include "manglefn.h"

/* polygons declared in rdmask_() */
extern int npolys;
extern polygon *polys[];

/*------------------------------------------------------------------------------
  Simplified fortran interface to cmlim_polys routine.
      real*8 mtol
      real*8 rp(3)
      call cmlimpolys(mtol, rp)
*/
void cmlimpolys_(long double *mtol, vec rp)
{
    int ndone;

    ndone = cmlim_polys(npolys, polys, *mtol, rp);
    if (ndone == -1) exit(1);
}

/*------------------------------------------------------------------------------
  Simplified fortran interface to drangle_polys routine.
      real*8 mtol
      real*8 rp(3)
      integer nth
      real*8 cm(nth),dr(nth)
      call dranglepolys(mtol, rp, nth, cm, dr)
*/
void dranglepolys_(long double *mtol, vec rp, int *nth, long double cm[/**nth*/], long double dr[/**nth*/])
{
    int ndone;

    ndone = drangle_polys(npolys, polys, *mtol, rp, *nth, cm, dr);
    if (ndone == -1) exit(1);
}
