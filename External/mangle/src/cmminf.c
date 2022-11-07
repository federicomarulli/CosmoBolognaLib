/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Find smallest cap of polygon.
*/
void cmminf(polygon *poly, int *ipmin, long double *cmmin)
{
    int ip;
    long double cmi;

    *cmmin = 2.;
    for (ip = 0; ip < poly->np; ip++) {
	cmi = (poly->cm[ip] >= 0.)? poly->cm[ip] : 2. + poly->cm[ip];
	if (cmi <= *cmmin) {
		*ipmin = ip;
		*cmmin = cmi;
	}
    }
}
