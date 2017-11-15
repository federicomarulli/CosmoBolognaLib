/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Id numbers of polygons containing position az, el.

   Input: poly = array of pointers to npoly polygons.
	  npoly = number of polygons in poly array.
	  az, el = angular position in radians.
  Output: id_p = pointer to array containing id numbers of polygons;
		 the required memory is allocated.
  Return value: number of polygons that contain az, el position.
*/
int poly_id(int npoly, polygon *poly[/*npoly*/], long double az, long double el, int **id_p, long double **weight_p)
{
/* number of extra polygon id numbers to allocate, to allow for expansion */
#define DNID		16
    static int nidmax = 0;
    static int *id = 0x0;
    static long double *weight = 0x0;

    int ipoly, nid;
    long double rp[3];

    /* unit vector corresponding to angular position az, el */
    rp[0] = cosl(el) * cosl(az);
    rp[1] = cosl(el) * sinl(az);
    rp[2] = sinl(el);

    nid = 0;

    /* keep trying till the id array is big enough */
    do {
        /* make sure that allocated id array contain enough space */
        if (!id || nid > nidmax) {
            if (id) free(id);
            if (weight) free(weight);
            id = (int *) malloc(sizeof(int) * (nid + DNID));
	    weight = (long double *) malloc(sizeof(long double) * (nid + DNID));
            if (!id) {
                fprintf(stderr, "poly_id: failed to allocate memory for %d ints\n", nid + DNID);
                return(-1);
            }
            if (!weight) {
                fprintf(stderr, "poly_id: failed to allocate memory for %d long doubles\n", nid + DNID);
                return(-1);
            }
	    nidmax = nid + DNID;
	}

	nid = 0;
	/* do each polygon in turn */
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	  /* id number of each polygon that contains az, el position */
	  if (gptin(poly[ipoly], rp)) {
	    if (nid < nidmax){ 
	      id[nid] = poly[ipoly]->id;
	      weight[nid] = poly[ipoly]->weight;
	    }
	    nid++;
	  }
	}
	
    } while (nid > nidmax);
    
    /* point id_p at id array */
    *id_p = id;
    /* point weight_p at weight array */
    *weight_p = weight;

    /* number of polygons containing az, el position */
    return(nid);
}
