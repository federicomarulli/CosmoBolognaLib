/*-----------------------------------------------------------------------------
© A J S Hamilton 2001
-----------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Make almost coincident caps of polygons coincide.

  Input: npoly = number of polygons to snap.
	 *poly[npoly] = array of npoly pointers to polygon structures.
	 mtol = tolerance angle for multiple intersections.
	 fmt = pointer to format structure.
	 axtol, btol, thtol, ytol = tolerance angles (see documentation).
	 selfsnap = determines whether or not to snap edges only against edges of
	            the same polygon.
  Return value: number of caps adjusted.
*/
int snap(int npoly, polygon *poly[/*npoly*/], long double mtol, format *fmt, long double axtol, long double btol, long double thtol, long double ytol, int selfsnap)
{
#define WARNMAX		8
  int i, inull, iprune, nadj, dnadj, warnmax;
  int *start;
  int *total;
  int p, max_pixel, ier;
  
  /* start by sorting polygons by pixel number*/
  poly_sort(npoly,poly,'p');

  /* allocate memory for pixel info arrays start and total */ 
  /* if only self-snapping, don't use pixelization */
  max_pixel=(selfsnap)? 1 : poly[npoly-1]->pixel+1; 
  start = (int *) malloc(sizeof(int) * max_pixel);
  if (!start) {
    fprintf(stderr, "snap_polys: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }
  total = (int *) malloc(sizeof(int) * max_pixel);
  if (!total) {
    fprintf(stderr, "snap_polys: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }
  
  /* if we're only doing self-snapping, don't use the pixelization info */ 
  if(selfsnap){
    start[0]=0;
    total[0]=npoly;
  }
  else{
    /* build lists of starting indices of each pixel and total number of polygons in each pixel*/
    ier=pixel_list(npoly, poly, max_pixel, start, total);
    if (ier == -1) {
      fprintf(stderr, "snap: error building pixel index lists\n");
      return(-1);
    }
  }
  
  /*turn off warning messages if using more than one pixel*/
  warnmax= (max_pixel<=1) ? WARNMAX : 0; 

  /* snap edges of polygons to each other */
  nadj=0;
  for(p=0;p<max_pixel;p++){
    if(total[p]==0) continue;
    dnadj=snap_polys(fmt, total[p], &poly[start[p]], selfsnap, axtol, btol, thtol, ytol, mtol,((selfsnap)? warnmax : warnmax/2),0x0);
    if(dnadj==-1) return(-1);
    nadj+=dnadj;
  }
  
  /* prune polygons */
  inull = 0;
  for (i = 0; i < npoly; i++) {
    iprune = prune_poly(poly[i], mtol);
    if (iprune >= 2) {
      if (WARNMAX > 0 && inull == 0)
	msg("warning from snap: the following polygons have zero area:\n");
      if (inull < WARNMAX) {
	msg(" %d", (fmt->newid == 'o')? poly[i]->id : i);
      } else if (inull == WARNMAX) {
	msg(" ... more\n");
      }
      inull++;
    }
  }
  if (WARNMAX > 0 && inull > 0 && inull <= WARNMAX) msg("\n");
  if (inull > 0) msg("snap: %d snapped polygons have zero area (but are being retained)\n", inull);
  
  /* assign new polygon id numbers */
  if (fmt->newid == 'n') {
    for (i = 0; i < npoly; i++) {
      poly[i]->id = i;
    }
  }

  if (fmt->newid == 'p') {
    for (i = 0; i < npoly; i++) {
      poly[i]->id = poly[i]->pixel;
    }
  }

  
  msg("snap: total of %d caps adjusted\n", nadj);
  
  return(nadj);
}
