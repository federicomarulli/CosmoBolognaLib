/*------------------------------------------------------------------------------
© M E C Swanson 2005
------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <math.h>
#include "pi.h"
#include "manglefn.h"

/* Function which_pixel returns the pixel number for a given azimuth and 
   elevation angle.
   inputs:
   az: azimuth angle (in radians)
   el: elevation angle (in radians)
   res: desired resolution of the pixel to be returned
   scheme: pixelization scheme
   returns the number of the pixel containing the point, or -1 if error occurs
*/

int which_pixel(long double az, long double el, int res, char scheme)
{
  int n,m,pix,base_pix,i;
  unsigned long pixnum;
  long double az_check, el_check;
  int *parent_pixels;

  if(az<0){
    az+=TWOPI;
  }
  
  if(az>TWOPI || az<0){
    fprintf(stderr, "error in which_pixel: az must lie between 0 and 2*PI.\n");
    return(-1);
  }
  if(el>PIBYTWO || el<-PIBYTWO){
    fprintf(stderr, "error in which_pixel: el must lie between -PI/2 and PI/2.\n");
    return(-1);
  }
  if(res<0){
    fprintf(stderr, "error in which_pixel: resolution must be an integer >=0.\n");
    return(-1);
  }
  

  if(scheme=='s'){
    // this scheme divides up the sphere by rectangles in az and el, and is numbered 
    // such that the resolution is encoded in each pixel number.  The whole sky is pixel 0,
    // pixels 1, 2, 3, and 4 are each 1/4 of the sky (resolution 1), pixels 5-20 are each 
    // 1/16 of the sky (resolution 2), etc.


    if(az==TWOPI) az=0;
    n=(sinl(el)==1) ? 0 : ceill((1-sinl(el))/2*powl(2,res))-1;
    
    m=floorl(az/(TWOPI)*powl(2,res));
    base_pix=powl(2,res)*n+m;
    pix=pixel_start(res,scheme)+base_pix; 
    return(pix);
  }
  else if(scheme=='d'){
    assign_parameters();
    if(az==TWOPI) az=0;
    /* ang2pix_radec takes az and el in degrees */
    az *= (180.0/PI);
    el *= (180.0/PI);

    if(res==0){
      return(0);
    }

    else if(res==1){
      parent_pixels = (int *) malloc(sizeof(int) * (3));
      if(!parent_pixels){
	fprintf(stderr, "error in which_pixel: failed to allocate memory for 3 integers\n");
	return(-1);
      }

      ang2pix_radec(1, az, el, &pixnum);
      pix = (int)pixnum;

      i = get_parent_pixels(pix+pixel_start(2,scheme), parent_pixels, scheme);
      if(i==1){
	fprintf(stderr, "error in which_pixel: get_parent_pixels failed\n");
	return(-1);
      }

      return(parent_pixels[1]);
    }

    else{
      ang2pix_radec((int)powl(2,res-2), az, el, &pixnum);
      pix = (int)pixnum;

       /* check
       pix2ang_radec((int)powl(2,res-1), pixnum, &az_check, &el_check);
       printf("pix2ang_radec(pixnum = %d) = %Lf, %Lf\n", (int)pixnum, az_check, el_check); */

      return(pix+pixel_start(res, scheme));
    }
  }
  else{
    fprintf(stderr, "error in which_pixel: pixel scheme %c not recognized.\n", scheme);
    return(-1);  
  }
}

/* Function get_parent_pixels generates a list of the parent pixels for a given child pixel.
   inputs:
   pix_c: number of child pixel 
   scheme: pixelization scheme
   output:
   pix_p[]: array containing parent pixels.  pix_p[r] is the parent pixel of resolution r.
   returns 0 on success, 1 on error
*/

int get_parent_pixels(int pix_c, int pix_p[], char scheme){
  int m,n,res,base_pix,i,j;
  unsigned long pixp;
  //long double res_d;

  if(pix_c<0){
    fprintf(stderr, "error in get_parent_pixels: %d is not a valid pixel number\n",pix_c);
    return(1);
  }

  res=get_res(pix_c, scheme);
  if(res==-1) return (1);

  if(scheme=='s'){
    // this scheme divides up the sphere by rectangles in az and el, and is numbered 
    // such that the resolution is encoded in each pixel number.  The whole sky is pixel 0,
    // pixels 1, 2, 3, and 4 are each 1/4 of the sky (resolution 1), pixels 5-20 are each 
    // 1/16 of the sky (resolution 2), etc.
    
    base_pix=pix_c-pixel_start(res,scheme);
    m=base_pix % (int)(powl(2,res));
    n=(base_pix-m)/powl(2,res);

    for(i=res;i>=0;i--){
      //put pixel number into array
      pix_p[i]=pixel_start(i,scheme)+(int)(powl(2,i))*n+m;
      //make child pixel into next parent pixel
      n=n/2;
      m=m/2;
    }
    return(0);
  }
  else if(scheme=='d'){
    assign_parameters();
    //printf("res = %d\n", res);
    //printf("res1 (1) = %d\n", res1);
    //if (res >= 1) { 
    //  res_d = (long double)(logl((long double)res)/logl(2.0))+1;
    //  res1 = (int)(res_d + 0.1);
    //}
    //printf("res1 (2) = %d\n", res1);
    pix_p[res]=pix_c;
    for(i=res;i>2;i--){
      //printf("args to superpix: %d, %d, %d\n", (int)powl(2,i-1), pix_p[i]-pixel_start(i, scheme), (int)powl(2,i-2));
      superpix((int)powl(2,i-2), (unsigned long)pix_p[i]-(unsigned long)pixel_start(i, scheme), (int)powl(2,i-3), &pixp);
      //printf("pixp = %d\n", (int)pixp);
      pix_p[i-1] = (int)pixp + pixel_start(i-1, scheme);
    }
    for(j=0;j<=5;j++){
      for(i=118+j*72;i<=152+j*72;i+=2){
	if(pix_p[2]==i || pix_p[2]==i+1 || pix_p[2]==i+36 || pix_p[2]==i+37) pix_p[1]=(i-118-j*36)/2+1;
      }
    }
    for(i=550;i<=585;i++){
      if(pix_p[2]==i) {
	pix_p[1]=(i-114)/4;
      }
    }
    pix_p[0]=0;
    return(0);
  }
  else{
    fprintf(stderr, "error in get_parent_pixels: pixel scheme %c not recognized.\n", scheme);
    return(1);   
  }
}

/* Function pixel_start returns the starting pixel number for the pixels of a given resolution
   inputs:
   res: resolution 
   scheme: pixelization scheme
   returns the starting pixel number for the set of pixels of resolution res, 
   or -1 if error occurs
*/

int pixel_start(int res, char scheme){

  //int res1;
  //long double res_d;
  
  if(res<0){
    fprintf(stderr, "error in pixel_start: %d not a valid resolution.\n", res);
    return(-1);
  }
  
  if(scheme=='s'){
    return ((int)((powl(4,res)-1)/3));
  }
  else if(scheme=='d'){
    //res_d = (long double)(logl((long double)res)/logl(2.0))+1;
    //res1 = (int)(res_d + 0.1);
    //printf("pixel_start: res = %d\n", res);
    if(res==0) return(0);
    else if(res==1) return(1);
    else return (468*(((int)powl(4,res-1)-4)/12)+118);
  }
  else{
    fprintf(stderr, "error in pixel_start: pixel scheme %c not recognized.\n", scheme);
    return(-1);  
  }
  
}  
