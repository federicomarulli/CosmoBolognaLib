/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Determine whether unit vector lies inside polygon.

  This is a c interface to fortran logical function gptin.

   Input: poly is a polygon.
	  rp = unit vector.
  Return value: 1 if in;
		0 if not in.
*/
int gptin(polygon *poly, vec rp)
{
    /* fortran routine */
    if (gptin_(poly->rp, poly->cm, &poly->np, rp)) return(1);
    return(0);
}
