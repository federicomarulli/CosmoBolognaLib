/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Find position within ordered array by binary chop.
  Assumes points are uncorrelated between calls
  (otherwise it would be better to start search from previous point).

  Return value: i such that array[i-1] <= point < array[i]
		0 if point < array[0]
		n if point >= array[n-1]
*/
int search(int n, long double array[/*n*/], long double point)
{
    int i, im, ip;

    /* point below minimum */
    if (point < array[0]) {
	return(0);

    /* point above maximum */
    } else if (point >= array[n-1]) {
	return(n);

    /* point between limits */
    } else {
	im = 0;
	ip = n-1;
	/* binary chop */
	while (im + 1 < ip) {
	    i = (im + ip) / 2;
	    if (point < array[i]) {
		ip = i;
	    } else if (point >= array[i]) {
		im = i;
	    }
	};
	return(ip);
    }
}
