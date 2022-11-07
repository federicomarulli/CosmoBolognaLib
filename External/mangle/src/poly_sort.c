/*------------------------------------------------------------------------------
© M E C Swanson 2005
------------------------------------------------------------------------------*/
#include <stdlib.h>
#include "manglefn.h"

/*polygon comparison functions*/
int poly_cmp_pixel(polygon **poly1, polygon **poly2)
{
  int pixel=(*poly1)->pixel - (*poly2)->pixel;
  return(pixel);
}
int poly_cmp_id(polygon **poly1, polygon **poly2)
{
  int id=(*poly1)->id - (*poly2)->id;
  return(id);
}
int poly_cmp_weight(polygon **poly1, polygon **poly2)
{
  long double weight=(*poly1)->weight - (*poly2)->weight;
  return(weight);
}

/*polygon sorting function
  sorts an array polys of n polygon pointers by either pixel number ('p'), 
  id number ('i') or weight ('w').
*/
void poly_sort(int npoly, polygon *poly[],char key){
  if(key=='p'){
    mysort(poly, npoly, sizeof(polygon *),poly_cmp_pixel);
  }
  else if(key=='i'){
    mysort(poly, npoly, sizeof(polygon *),poly_cmp_id);
  }
  else if(key=='w'){
    mysort(poly, npoly, sizeof(polygon *),poly_cmp_weight);
  }
  else{
    fprintf(stderr,"sort key %c not recognized.  Array can't be sorted\n", key);
  }

}

/* Function pixel_list generates arrays start[] and total[] from a polygon array which has been
   sorted by pixel number.
   inputs: 
   npoly: number of polygons in polys
   poly[]: array of pointers to polygons (must be sorted by pixel number)
   max_pixel: highest allowed pixel number in the polygon array
   outputs:
   start[k] contains the starting index for the polygons in pixel k
   total[k] contains the total number of polygons in pixel k
   returns 0 if successful, -1 if there's an error

  */ 
int pixel_list(int npoly, polygon *poly[], int max_pixel, int start[], int total[]){
  int i,j,k,k_old;

  /* initialize output arrays */

  for(i=0;i<max_pixel;i++){
    start[i]=0;
    total[i]=0;
  }
  k_old=-1;
  for(j=0;j<npoly;j++){
    k=poly[j]->pixel;
    if(k<k_old){
      fprintf(stderr, "Error in pixel_list: polygon array not sorted.  Please use poly_sort first.\n");
      return(-1);
    }
    if(k>max_pixel){
      fprintf(stderr, "Error in pixel_list: polygon %d is in pixel %d > max_pixel %d\n", j, k, max_pixel);
      return(-1);
    }
    
    if(k==k_old){
      total[k]++;
    }
    else if(k>k_old){
      start[k]=j;
      total[k]++;
    }
    k_old=k;
  }
  return(0);
}

