/*--------------------------------------------------------------------
Imported from SDSSPix software by J C Hill 08/14/06.
Get the original at http://lahmu.phyast.pitt.edu/~scranton/SDSSPix/
--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "manglefn.h"

unsigned long nx0, ny0;
long double pi, deg2Rad, rad2Deg, strad2Deg, etaOffSet;
long double surveyCenterRA, surveyCenterDEC, node, etaPole;

void assign_parameters()
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg, etaOffSet;
  extern long double surveyCenterRA, surveyCenterDEC, node, etaPole;

  /* With these base resolutions, we can achieve nearly equal area pixels, with
     nearly square pixels at lambda = 30 degrees.  A resolution factor of 
     4 gives us pixels as wide as a stripe.  The seeing/reddening/sky maps
     are made with resolution=256. */

  nx0 = 36;
  ny0 = 13;

  pi = 2.0*asinl(1.0);
  deg2Rad = pi/180.0;
  rad2Deg = 180.0/pi;
  strad2Deg = 360.0*360.0/(4.0*pi*pi);

  /* These parameters are necessary for translation of LAMBDA-ETA coordinates
     into x-y-z vectors in a way consistent with translations from RA-DEC.  We
     also need etaOffSet to compensate for the fact that the stripes stradle
     the ETA = 0 meridian, rather than being bounded by it. */

  etaOffSet = 1.25;
  surveyCenterRA = 185.0;
  surveyCenterDEC = 32.5;
  node = deg2Rad*(surveyCenterRA - 90.0);
  etaPole = deg2Rad*surveyCenterDEC;
}

void pix2ang(int resolution, unsigned long pixnum, long double *lambda, long double *eta)
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  long nx, ny, i, j;

  /* This module takes a pixel index and converts it into LAMBDA-ETA
     coordinates.  Notice that we require that -180 < ETA < 180 and use the
     offset in ETA. 

     Also in IDL: pix2ang.pro
  */
  
  nx = nx0*resolution;
  ny = ny0*resolution;

  j = pixnum/nx;
  i = pixnum - nx*j;

  *eta = rad2Deg*(2.0*pi*(i+0.5))/nx + etaOffSet;
  if (*eta >= 180.0) *eta -= 360.0;
  *lambda = 90.0 - rad2Deg*acosl(1.0-2.0*(j+0.5)/ny);
}

void ang2pix(int resolution, long double lambda, long double eta, unsigned long *pixnum)
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  unsigned long nx, ny, i, j;
  long double eta2;

  /* The complement to pix2ang, this module converts LAMBDA-ETA to pixel
     index.  This single number uniquely identifies a pixel on the sky. 
     
     Again, we need to subtract off the etaOffSet to make the coordinate 
     system work. 

     Also in IDL: ang2pix.pro
*/

  nx = nx0*resolution;
  ny = ny0*resolution;

  eta -= etaOffSet;

  eta *= deg2Rad;

  if (eta >= 0.0) {
    eta2 = eta;
  } else {
    eta2 = eta + 2.0*pi;
  }

  i = nx*eta2/(2.0*pi);

  lambda = (90.0 - lambda)*deg2Rad;

  if (lambda >= pi) {
    j = ny - 1;
  } else {
    j = ny*((1.0 - cosl(lambda))/2.0);
  }

  *pixnum = nx*j + i;
}

void pix2ang_radec(int resolution, unsigned long pixnum, long double *ra, long double *dec) {

  void csurvey2eq(long double lambda, long double eta, long double *ra, long double *dec);
  long double lam, eta, ra_tmp, dec_tmp;

  /* Same as pix2ang, but returning RA-DEC coordinates instead of LAMBDA-ETA */

  pix2ang(resolution,pixnum,&lam,&eta);

  csurvey2eq(lam,eta,&ra_tmp,&dec_tmp);

  *ra = ra_tmp;
  *dec = dec_tmp;
}

void ang2pix_radec(int resolution, long double ra, long double dec, unsigned long *pixnum) {

  void eq2csurvey(long double ra, long double dec, long double *lambda, long double *eta);
  unsigned long tmp_pixnum;
  long double lambda, eta;

  /* Same as ang2pix, but taking RA-DEC coordinates instead of LAMBDA-ETA */

  eq2csurvey(ra,dec,&lambda,&eta);
 
  ang2pix(resolution,lambda,eta,&tmp_pixnum);

  *pixnum = tmp_pixnum;
}

void csurvey2eq(long double lambda, long double eta, long double *ra, long double *dec) {

  long double x, y, z;

  /* Conversion from LAMBDA-ETA to RA-DEC coordinates */

  x = -1.0*sinl(lambda*deg2Rad);
  y = cosl(lambda*deg2Rad)*cosl(eta*deg2Rad+etaPole);
  z = cosl(lambda*deg2Rad)*sinl(eta*deg2Rad+etaPole);

  *ra = (atan2l(y,x) + node)/deg2Rad;
  *dec = asinl(z)/deg2Rad;
}

void eq2csurvey(long double ra, long double dec, long double *lambda, long double *eta) {

  long double x, y, z;

  /* Conversion from RA-DEC to LAMBDA-ETA coordinates */

  x = cosl(deg2Rad*ra-node)*cosl(deg2Rad*dec);
  y = sinl(deg2Rad*ra-node)*cosl(deg2Rad*dec);
  z = sinl(deg2Rad*dec);

  *lambda = -1.0*asinl(x)/deg2Rad;
  *eta = (atan2l(z,y) - etaPole)/deg2Rad;

  if (*eta < -180.0) *eta += 360.0;
  if (*eta > 180.0) *eta -= 360.0;
}

/* void downsample(int resolution, gsl_matrix *inmap, gsl_matrix *outmap)
{
  extern unsigned long nx0, ny0;
  unsigned long i, j, nx, ny;

  nx = nx0*resolution;
  ny = ny0*resolution;

  for (i=0;i<nx/2;i++) {
    for (j=0;j<ny/2;j++) {
      outmap->data[i*outmap->tda+j] = 
	0.25*(inmap->data[2*i*inmap->tda+2*j] + 
	      inmap->data[2*i*inmap->tda+2*j+1] + 
	      inmap->data[(2*i+1)*inmap->tda+2*j] + 
	      inmap->data[(2*i+1)*inmap->tda+2*j+1]);
    }
  }
} */
  
/* void upsample(int resolution, gsl_matrix *inmap, gsl_matrix *outmap)
{
  extern unsigned long nx0, ny0;
  unsigned long i, j, nx, ny;

  nx = nx0*resolution;
  ny = ny0*resolution;

  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      outmap->data[2*i*outmap->tda+2*j] = inmap->data[i*inmap->tda+j];
      outmap->data[2*i*outmap->tda+2*j+1] = inmap->data[i*inmap->tda+j];
      outmap->data[(2*i+1)*outmap->tda+2*j] = inmap->data[i*inmap->tda+j];
      outmap->data[(2*i+1)*outmap->tda+2*j+1] = inmap->data[i*inmap->tda+j];
    }
  }
} */

void superpix(int hi_resolution, unsigned long hi_pixnum, 
	      int lo_resolution, unsigned long *lo_pixnum) 
{
  extern unsigned long nx0, ny0;
  unsigned long nx_hi, ny_hi, nx_lo, ny_lo, i, j, ratio;

  /* Takes pixel index at high resolution and returns the pixel which contains
     it at a lower resolution.

     Also in IDL: superpix.pro
  */

  if (hi_resolution < lo_resolution) {
    printf("Can't go from low resolution to higher resolution.\n");
    exit(1);
  }

  if (lo_resolution == 0) *lo_pixnum = 0;
  else {
    nx_hi = nx0*hi_resolution;
    ny_hi = ny0*hi_resolution;
    nx_lo = nx0*lo_resolution;
    ny_lo = ny0*lo_resolution;
  
    ratio = hi_resolution/lo_resolution;

    j = hi_pixnum/nx_hi;
    i = hi_pixnum - nx_hi*j;

    i /= ratio;
    j /= ratio;
   
    *lo_pixnum = nx_lo*j + i;
  }
}

void subpix(int resolution, unsigned long pixnum, unsigned long *sub_pixnum1, 
	    unsigned long *sub_pixnum2, unsigned long *sub_pixnum3, 
	    unsigned long *sub_pixnum4)
{
  extern unsigned long nx0, ny0;
  unsigned long nx_hi, ny_hi, nx_lo, ny_lo, i, j;

  /* Reverse of superpix.  

  **Only works for long doubled resolution**

  */

  
  nx_hi = 2*nx0*resolution;
  ny_hi = 2*ny0*resolution;
  nx_lo = nx0*resolution;
  ny_lo = ny0*resolution;
  
  j = pixnum/nx_lo;
  i = pixnum - nx_lo*j; 
   
  *sub_pixnum1 = nx_hi*(2*j) + 2*i;
  *sub_pixnum2 = nx_hi*(2*j) + 2*i + 1;
  *sub_pixnum3 = nx_hi*(2*j + 1) + 2*i;
  *sub_pixnum4 = nx_hi*(2*j + 1) + 2*i + 1;

}

void pix_bound(int resolution, unsigned long pixnum, 
	       long double *lammin, long double *lammax, long double *etamin, long double *etamax)
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  unsigned long nx, ny, i, j;

  /* Returns ETA-LAMBDA boundaries for a given pixel index.

  Also in IDL: pix_bound.pro

  */

  nx = nx0*resolution;
  ny = ny0*resolution;

  j = pixnum/nx;
  i = pixnum - nx*j;
  
  *etamin = rad2Deg*(2.0*pi*i)/nx + etaOffSet;
  if (*etamin >= 180.0) *etamin -= 360.0;
  *etamax = rad2Deg*(2.0*pi*(i+1))/nx + etaOffSet;
  if (*etamax >= 180.0) *etamax -= 360.0;
  *lammin = 90.0 - rad2Deg*acosl(1.0 - 2.0*(j+1)/ny);
  *lammax = 90.0 - rad2Deg*acosl(1.0 - 2.0*j/ny);

}

long double pix_area(int resolution, unsigned long pixnum)
{
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  long double lammin, lammax, etamin, etamax;

  /* Returns the area of a given pixel index and resolution */

  pix_bound(resolution,pixnum,&lammin,&lammax,&etamin,&etamax);

  return strad2Deg*(deg2Rad*(etamax-etamin))*
    (sinl(deg2Rad*lammax)-sinl(deg2Rad*lammin));
}

void pix2xyz(int resolution, unsigned long pixnum, 
	     long double *x, long double *y, long double *z)
{
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  extern long double surveyCenterRA, surveyCenterDEC, node, etaPole;
  long double lam, eta;

  pix2ang(resolution,pixnum,&lam,&eta);

  *x = -1.0*sinl(lam*deg2Rad);
  *y = cosl(lam*deg2Rad)*cosl(eta*deg2Rad+etaPole);
  *z = cosl(lam*deg2Rad)*sinl(eta*deg2Rad+etaPole);
}

void area_index(int resolution, long double lammin, long double lammax, long double etamin, 
		long double etamax, unsigned long *x_min, unsigned long *x_max, 
		unsigned long *y_min, unsigned long *y_max)
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  void ang2pix(int resolution, long double lambda, long double eta, 
	       unsigned long *pixnum);
  unsigned long nx, ny, pixnum;

  /* Given a range in LAMBDA and ETA, returns the x-y boundaries.  Pixel
     indices in this range can be found by taking 

     nx = nx0*resolution;
     ny = ny0*resolution;
  
     n_pixel = (x_max - x_min + 1)*(y_max - y_min + 1);

     index_array = gsl_vector_long_alloc(n_pixel);

     k = 0;
     for (j=y_min;j<=y_max;j++) {
       for (i=x_min;i<=x_max;i++) {
         index_array->data[k] = nx*j + i;
         k++;
       }
     }

     Also in IDL: area_index.pro
  */


  nx = nx0*resolution;
  ny = ny0*resolution;

  ang2pix(resolution,lammax,etamin,&pixnum);
  *y_min = pixnum/nx;
  *x_min = pixnum - nx*(*y_min);
  
  ang2pix(resolution,lammin,etamax,&pixnum);
  *y_max = pixnum/nx;
  *x_max = pixnum - nx*(*y_max);

  *y_min += 1;
  *y_max -= 1;
  
}

void area_index_stripe(int resolution, int stripe, 
		       unsigned long *x_min, unsigned long *x_max, 
		       unsigned long *y_min, unsigned long *y_max)
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  void ang2pix(int resolution, long double lambda, long double eta, 
	       unsigned long *pixnum);
  void primary_bound(int stripe, long double *lammin, long double *lammax, 
		     long double *etamin, long double *etamax);
  unsigned long nx, ny, pixnum;
  long double lammin, lammax, etamin, etamax;


  /* Similar to area_index, but this returns the x-y range of the pixels in 
     the primary region of a given stripe.

     Also in IDL: area_index_stripe.pro
  */

  nx = nx0*resolution;
  ny = ny0*resolution;

  primary_bound(stripe,&lammin,&lammax,&etamin,&etamax);

  /*printf("Stripe %i: %Lf %Lf %Lf %Lf\n",stripe,
    lammin,lammax,etamin,etamax);*/

  ang2pix(resolution,lammax,etamin,&pixnum);
  *y_min = pixnum/nx;
  *x_min = pixnum - nx*(*y_min);

  ang2pix(resolution,lammin,etamax,&pixnum);
  *y_max = pixnum/nx;
  *x_max = pixnum - nx*(*y_max);

  *y_min += 1;
  *y_max -= 1;

  /* printf("Found %u (%u x %u) pixels within stripe boundary...\n",
     (*x_max - *x_min + 1)*(*y_max - *y_min + 1),
     (*x_max - *x_min + 1),(*y_max - *y_min + 1));*/

}
  
/* void sort_mask_resolution(gsl_vector_int *resolution_array, 
			  gsl_vector_ulong *pixnum_array, unsigned long n_mask)
{
  gsl_permutation *pixel_index;
  gsl_vector_ulong *tmp_pixnum_array; 
  gsl_vector_int *tmp_resolution_array;
  unsigned long i, j;

  tmp_pixnum_array = gsl_vector_ulong_alloc(n_mask);
  tmp_resolution_array = gsl_vector_int_alloc(n_mask);
  pixel_index = gsl_permutation_alloc(n_mask);

  gsl_sort_vector_int_index(pixel_index,resolution_array);

  for (i=0;i<n_mask;i++) {
    j = pixel_index->data[i];
    tmp_pixnum_array->data[i] = pixnum_array->data[j];
    tmp_resolution_array->data[i] = resolution_array->data[j];
  }

  for (i=0;i<n_mask;i++) {
    pixnum_array->data[i] = tmp_pixnum_array->data[i];
    resolution_array->data[i] = tmp_resolution_array->data[i];
  }

  gsl_vector_int_free(tmp_resolution_array);
  gsl_vector_ulong_free(tmp_pixnum_array);
  gsl_permutation_free(pixel_index);

}

void sort_mask_pixnum(gsl_vector_ulong *pixnum_array, 
		      gsl_vector_int *resolution_array, unsigned long n_mask,
		      gsl_vector_int *resolution_region_array, 
		      gsl_vector_ulong *resolution_start_array,
		      gsl_vector_ulong *resolution_finish_array, int n_res)
{
  unsigned long i, j, k, n, n_sub_mask;
  gsl_vector_ulong *tmp_pixnum_array, *total_pixnum_array;
  gsl_permutation *pixel_index;

  total_pixnum_array = gsl_vector_ulong_alloc(n_mask);

  for (k=0;k<n_res;k++) {
    if (resolution_start_array->data[k] == resolution_finish_array->data[k]) {
      j = resolution_start_array->data[k];
      total_pixnum_array->data[j] = pixnum_array->data[j];
    } else {
      n_sub_mask = resolution_finish_array->data[k] - 
	resolution_start_array->data[k] + 1;
      
      tmp_pixnum_array = gsl_vector_ulong_alloc(n_sub_mask);
      pixel_index = gsl_permutation_alloc(n_sub_mask);
      
      j = resolution_start_array->data[k];
      for (i=0;i<n_sub_mask;i++) {
	tmp_pixnum_array->data[i] = pixnum_array->data[j];
	j++;
      }
  
      gsl_sort_vector_ulong_index(pixel_index,tmp_pixnum_array);

      j = resolution_start_array->data[k];
      for (i=0;i<n_sub_mask;i++) {
	n = pixel_index->data[i];
	total_pixnum_array->data[j] = tmp_pixnum_array->data[n];
	j++;
      }
      
      gsl_vector_ulong_free(tmp_pixnum_array);
      gsl_permutation_free(pixel_index);
    }
  }
   
  for (i=0;i<n_mask;i++) 
    pixnum_array->data[i] = total_pixnum_array->data[i];

  gsl_vector_ulong_free(total_pixnum_array);

}

int find_n_res(gsl_vector_int *resolution_array, unsigned long n_mask)
{
  int n_res, resolution;
  unsigned long i, not_sorted;

  n_res = 0;
  not_sorted = 0;

  for (i=1;i<n_mask;i++) {
    if (resolution_array->data[i-1] > resolution_array->data[i]) {
      not_sorted = 1;
    }
  }
  
  if (not_sorted == 1) {
    n_res = -1;
  } else {
    resolution = resolution_array->data[0];
    n_res = 1;
    for (i=1;i<n_mask;i++) {
      if (resolution_array->data[i] != resolution) {
	n_res++;
	resolution = resolution_array->data[i];
      }
    }
  }
  
  return n_res;
}

long find_n_superpix(int superpix_resolution, gsl_vector_ulong *pixnum_array,
		    gsl_vector_int *resolution_array, unsigned long n_mask)
{
  gsl_vector_ulong *superpix_array;
  unsigned long superpixnum;
  unsigned long n_superpix;
  unsigned long i,j,k;

  n_superpix = 0;

  superpix_array = gsl_vector_ulong_alloc(n_mask);

  for (i=0;i<n_mask;i++) 
    superpix(resolution_array->data[i],pixnum_array->data[i],
	     superpix_resolution,&superpix_array->data[i]);

  gsl_sort_vector_ulong(superpix_array);

  n_superpix = 1;
  for (i=1;i<n_mask;i++) {
    if (superpix_array->data[i] != superpix_array->data[i-1]) {
      n_superpix++;
    }
  }
  
  return n_superpix;
}


void find_resolution_bounds(gsl_vector_int *resolution_array, 
			    unsigned long n_mask,
			    gsl_vector_int *resolution_region_array, 
			    gsl_vector_ulong *resolution_start_array,
			    gsl_vector_ulong *resolution_finish_array)
{
  unsigned long i, n_res;
  
  resolution_region_array->data[0] = resolution_array->data[0];
  resolution_start_array->data[0] = 0;
  resolution_finish_array->data[0] = 0;
  n_res = 0;

  for (i=1;i<n_mask;i++) {
    if (resolution_array->data[i] == resolution_array->data[i-1]) {
      resolution_finish_array->data[n_res] = i;
    } else {
      n_res++;
      resolution_region_array->data[n_res] = resolution_array->data[i];
      resolution_start_array->data[n_res] = i;
      resolution_finish_array->data[n_res] = i;
    }
  }
}

void find_superpix_bounds(gsl_vector_ulong *superpix_array, 
			  unsigned long n_mask,
			  gsl_vector_ulong *superpix_region_array, 
			  gsl_vector_ulong *superpix_start_array,
			  gsl_vector_ulong *superpix_finish_array)
{
  unsigned long i, n_superpix;
  
  superpix_region_array->data[0] = superpix_array->data[0];
  superpix_start_array->data[0] = 0;
  superpix_finish_array->data[0] = 0;
  n_superpix = 0;

  for (i=1;i<n_mask;i++) {
    if (superpix_array->data[i] == superpix_array->data[i-1]) {
      superpix_finish_array->data[n_superpix] = i;
    } else {
      n_superpix++;
      superpix_region_array->data[n_superpix] = superpix_array->data[i];
      superpix_start_array->data[n_superpix] = i;
      superpix_finish_array->data[n_superpix] = i;
    }
  }
}


void make_resolution_struct(gsl_vector_ulong *pixnum_array, 
			    gsl_vector_int *resolution_array, 
			    unsigned long n_pixel, 
			    resolution_struct *res_struct, int n_res)
{
  unsigned long i,j;
  gsl_vector_ulong *resolution_start_array, *resolution_finish_array;
  gsl_vector_int *resolution_region_array;
  int k; 

  resolution_region_array = gsl_vector_int_alloc(n_res);
  resolution_start_array = gsl_vector_ulong_alloc(n_res);
  resolution_finish_array = gsl_vector_ulong_alloc(n_res);

  find_resolution_bounds(resolution_array,n_pixel,resolution_region_array,
			 resolution_start_array,resolution_finish_array);

  sort_mask_pixnum(pixnum_array,resolution_array,n_pixel,
		   resolution_region_array,resolution_start_array,
		   resolution_finish_array,n_res);
  
  for (k=0;k<n_res;k++) {
    res_struct[k].start = resolution_start_array->data[k];
    res_struct[k].finish = resolution_finish_array->data[k];
    res_struct[k].n_pixel = 
      res_struct[k].finish - res_struct[k].start + 1;
    res_struct[k].resolution = resolution_region_array->data[k];
    printf("%u pixels with resolution of %i\n",res_struct[k].n_pixel,
       res_struct[k].resolution);
    res_struct[k].pixnum = gsl_vector_ulong_alloc(res_struct[k].n_pixel);
    j = res_struct[k].start;
    for (i=0;i<res_struct[k].n_pixel;i++) {
      res_struct[k].pixnum->data[i] = pixnum_array->data[j];
      j++;
    }
  }

  gsl_vector_int_free(resolution_region_array);
  gsl_vector_ulong_free(resolution_start_array);
  gsl_vector_ulong_free(resolution_finish_array);

}

void make_superpix_struct(int superpix_resolution,
			  gsl_vector_ulong *pixnum_array, 
			  gsl_vector_int *resolution_array, 
			  unsigned long n_pixel, 
			  superpixnum_struct *superpix_struct, 
			  unsigned long n_superpix)
{
  unsigned long i,j, k;
  gsl_vector_ulong *superpix_start_array, *superpix_finish_array;
  gsl_vector_ulong *superpix_region_array, *tmp_pixnum_array, *superpix_array;
  gsl_vector_ulong *tmp_superpix_array;
  gsl_vector_int *tmp_resolution_array;
  gsl_permutation *pixel_index;

  superpix_region_array = gsl_vector_ulong_alloc(n_superpix);
  superpix_start_array = gsl_vector_ulong_alloc(n_superpix);
  superpix_finish_array = gsl_vector_ulong_alloc(n_superpix);
  superpix_array = gsl_vector_ulong_alloc(n_pixel);
  tmp_superpix_array = gsl_vector_ulong_alloc(n_pixel);
  tmp_pixnum_array = gsl_vector_ulong_alloc(n_pixel);
  tmp_resolution_array = gsl_vector_int_alloc(n_pixel);
  pixel_index = gsl_permutation_alloc(n_pixel);


  for (i=0;i<n_pixel;i++) {
    superpix(resolution_array->data[i],pixnum_array->data[i],
	     superpix_resolution,&superpix_array->data[i]);
    printf("%d\n",superpix_array->data[i]);
  }

  gsl_sort_vector_ulong_index(pixel_index,superpix_array);
  
  for (i=0;i<n_pixel;i++) {
    k = pixel_index->data[i];
    tmp_resolution_array->data[i] = resolution_array->data[k];
    tmp_pixnum_array->data[i] = pixnum_array->data[k];
    tmp_superpix_array->data[i] = superpix_array->data[k];
  }
    
  for (i=0;i<n_pixel;i++) {
    resolution_array->data[i] = tmp_resolution_array->data[i];
    pixnum_array->data[i] = tmp_pixnum_array->data[i];
    superpix_array->data[i] = tmp_superpix_array->data[i];
  }

  gsl_vector_ulong_free(tmp_superpix_array);
  gsl_vector_ulong_free(tmp_pixnum_array);
  gsl_vector_int_free(tmp_resolution_array);
  gsl_permutation_free(pixel_index);

  printf("Finding superpix bounds...\n");
  
  find_superpix_bounds(superpix_array,n_pixel,superpix_region_array,
		       superpix_start_array,superpix_finish_array);

  for (k=0;k<n_superpix;k++) {
    printf("%d %d %d\n",superpix_start_array->data[k],
       superpix_finish_array->data[k],superpix_region_array->data[k]);
    superpix_struct[k].n_pixel = 
      superpix_finish_array->data[k] - superpix_start_array->data[k] + 1;
    superpix_struct[k].resolution = superpix_resolution;
    superpix_struct[k].superpixnum = superpix_region_array->data[k];

    tmp_pixnum_array = gsl_vector_ulong_alloc(superpix_struct[k].n_pixel);
    tmp_resolution_array = gsl_vector_int_alloc(superpix_struct[k].n_pixel);
    pixel_index = gsl_permutation_alloc(superpix_struct[k].n_pixel);

    j = 0;
    for (i=superpix_start_array->data[k];
	 i<=superpix_finish_array->data[k];i++) {
      tmp_resolution_array->data[j] = resolution_array->data[i];
      j++;
    }

    gsl_sort_vector_int_index(pixel_index,tmp_resolution_array);
    
    for (i=0;i<superpix_struct[k].n_pixel;i++) {
      j = pixel_index->data[i] + superpix_start_array->data[k];
      tmp_pixnum_array->data[i] = pixnum_array->data[j];
      tmp_resolution_array->data[i] = resolution_array->data[j];
      printf("%i %u\n",tmp_resolution_array->data[i],
	 tmp_pixnum_array->data[i]);
    }

    superpix_struct[k].n_res = 
      find_n_res(tmp_resolution_array, superpix_struct[k].n_pixel);
    
    printf("Found %i resolutions in superpixel %u\n",
       superpix_struct[k].n_res,k);

    if (!(superpix_struct[k].res_struct=
	  malloc(superpix_struct[k].n_res*sizeof(resolution_struct)))) {
      printf("Couldn't allocate memory...\n");
      exit(1);
    }

    make_resolution_struct(tmp_pixnum_array, tmp_resolution_array, 
			   superpix_struct[k].n_pixel, 
			   superpix_struct[k].res_struct,
			   superpix_struct[k].n_res);

    gsl_vector_ulong_free(tmp_pixnum_array);
    gsl_vector_int_free(tmp_resolution_array);
    gsl_permutation_free(pixel_index);
  }

  gsl_vector_ulong_free(superpix_region_array);
  gsl_vector_ulong_free(superpix_start_array);
  gsl_vector_ulong_free(superpix_finish_array);
  gsl_vector_ulong_free(superpix_array);


}

void rand_pixel_position(int resolution, unsigned long pixnum, 
			 long double *lambda, long double *eta)
{
  extern unsigned long nx0, ny0;
  extern long double pi, deg2Rad, rad2Deg, strad2Deg; 
  void pix_bound(int resolution, unsigned long pixnum, 
		 long double *lammin, long double *lammax, 
		 long double *etamin, long double *etamax);
  long double z_min, z, z_length, lammin, lammax, etamin, etamax, eta_length;

  pix_bound(resolution,pixnum,&lammin,&lammax,&etamin,&etamax);

  eta_length = etamax - etamin;
  z_min = sinl(deg2Rad*lammin);
  z_length = 2.0/(ny0*resolution);

  *eta = eta_length*gsl_rng_uniform(mt19937_rand) + etamin;
  z = z_length*gsl_rng_uniform(mt19937_rand) + z_min;
  *lambda = asinl(z)/deg2Rad;
  
}
*/

long double stripe_inclination(int stripe) 
{
  return 2.5*(stripe-10);
}

void primary_bound(int stripe, long double *lammin, long double *lammax, 
		   long double *etamin, long double *etamax)
{
  long double inc;
  long double stripe_inclination(int stripe);
  
  inc = stripe_inclination(stripe);
  
  *etamin = inc - 32.5 - 1.25 + 0.0000001;
  *etamax = inc - 32.5 + 1.25 - 0.0000001;

  *lammin = 10000.0;
  *lammax = -10000.0;

  if (stripe == 9) {
    *lammin = -58.8;
    *lammax = 53.5;
  } 
  if (stripe == 10) {
    *lammin = -63.0;
    *lammax = 64.95;
  }
  if (stripe == 11) {
    *lammin = -60.4;
    *lammax = 55.2;
  }
  if (stripe == 12) {
    *lammin = -64.1;
    *lammax = 56.6;
  }
  if (stripe == 13) {
    *lammin = -62.15;
    *lammax = 57.8;
  } 
  if (stripe == 14) {
    *lammin = -62.4;
    *lammax = 58.9;
  }
  if (stripe == 15) {
    *lammin = -64.95;
    *lammax = 59.8;
  } 
  if (stripe == 16) {
    *lammin = -63.1;
    *lammax = 60.6;
  } 
  if (stripe == 17) {
    *lammin = -63.4;
    *lammax = 61.2;
  } 
  if (stripe == 18) {
    *lammin = -63.6;
    *lammax = 61.8;
  } 
  if (stripe == 19) {
    *lammin = -63.7;
    *lammax = 62.3;
  }
  if (stripe == 20) {
    *lammin = -63.8;
    *lammax = 62.7;
  }
  if (stripe == 21) {
    *lammin = -63.7;
    *lammax = 63.1;
  }
  if (stripe == 22) {
    *lammin = -63.7;
    *lammax = 63.3;
  }
  if (stripe == 23) {
    *lammin = -63.5;
    *lammax = 63.5;
  }
  if (stripe == 24) {
    *lammin = -63.3;
    *lammax = 63.7;
  }
  if (stripe == 25) {
    *lammin = -63.1;
    *lammax = 63.7;
  }
  if (stripe == 26) {
    *lammin = -62.7;
    *lammax = 63.8;
  }
  if (stripe == 27) {
    *lammin = -64.75;
    *lammax = 63.7;
  } 
  if (stripe == 28) {
    *lammin = -65.55;
    *lammax = 63.6;
  }
  if (stripe == 29) {
    *lammin = -62.8;
    *lammax = 63.4;
  }
  if (stripe == 30) {
    *lammin = -63.0;
    *lammax = 63.1;
  }
  if (stripe == 31) {
    *lammin = -60.4;
    *lammax = 62.8;
  }
  if (stripe == 32) {
    *lammin = -60.0;
    *lammax = 63.35;
  }
  if (stripe == 33) {
    *lammin = -60.0;
    *lammax = 61.9;
  }
  if (stripe == 34) {
    *lammin = -59.25;
    *lammax = 61.8;
  } 
  if (stripe == 35) {
    *lammin = -55.2;
    *lammax = 60.95;
  }
  if (stripe == 36) {
    *lammin = -53.6;
    *lammax = 61.65;
  }
  if (stripe == 37) {
    *lammin = -52.3;
    *lammax = 58.8;
  }
  if (stripe == 76) {
    *lammin = -27.95;
    *lammax = 48.5;
  }
  if (stripe == 82) {
    *lammin = -60.0;
    *lammax = 60.0;
  }
  if (stripe == 86) {
    *lammin = -61.8;
    *lammax = 55.7;
  }

  if (*lammin > 1000.0) *lammin = -60.0;
  if (*lammax < -1000.0) *lammax = 60.0;
}
/*
void hunt(gsl_vector *xx, long double x, long *jlo)
{
  long up, down, mid, n;
  int ascnd;
  
  n = xx->size;

  if ((x < xx->data[0]) || (x > xx->data[n-1])) {
    if (x < xx->data[0]) {
      *jlo = n+1;
    } else {
      *jlo = n-1;
    } 
  } else {
    
    down = -1;
    up = n;
    
    while (up-down > 1) {
      mid = down + (up-down)/2;
      if (x >= xx->data[mid]) {
        down = mid;
      } else {
        up = mid;
      }
    }
    *jlo = down;
  }
  
}


void ihunt(gsl_vector_int *xx, int x, long *jlo)
{
  long up, down, mid, n;
  int ascnd;
  
  n = xx->size;

  if ((x < xx->data[0]) || (x > xx->data[n-1])) {
    *jlo = n+1;
  } else {
    
    down = -1;
    up = n;
    
    while (up-down > 1) {
      mid = down + (up-down)/2;
      if (x >= xx->data[mid]) {
        down = mid;
      } else {
        up = mid;
      }
    }
    *jlo = down;
  }
  
}


void lhunt(gsl_vector_ulong *xx, unsigned long x, unsigned long *jlo)
{
  long up, down, mid, n;
  int ascnd;
  
  n = xx->size;

  if ((x < xx->data[0]) || (x > xx->data[n-1])) {
    *jlo = n+1;
  } else {
    
    down = -1;
    up = n;
    
    while (up-down > 1) {
      mid = down + (up-down)/2;
      if (x >= xx->data[mid]) {
        down = mid;
      } else {
        up = mid;
      }
    }
    *jlo = down;
  }
  
}

*/
