/*------------------------------------------------------------------------------
© M E C Swanson 2005
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include "pi.h"
#include "manglefn.h"


/* Function get_pixel takes a pixel number and returns a pointer to a polygon 
   containing that pixel.
   inputs:
   pix: pixel number
   scheme: pixelization scheme
   returns pointer to polygon containing pixel, or 0x0 if an error occurs
*/

polygon *get_pixel(int pix, char scheme){
  int m,n,res,base_pix,i,ier,pix_c[4];
  long double azmax, azmin, elmax, elmin;
  long double lammin, lammax, etamin, etamax;
  long double angle[4], lammin_c[4], lammax_c[4], etamin_c[4], etamax_c[4];
  azel v[4],v_r[4];
  polygon *pixel;
  
  if(pix<0){
    fprintf(stderr, "error in get_pixel: %d not a valid pixel number.\n", pix);
    return(0x0);
  }

  res=get_res(pix,scheme);
  if(res==-1) return (0x0);
  
  if(scheme=='s'){
    // this scheme divides up the sphere by rectangles in az and el, and is numbered 
    // such that the resolution is encoded in each pixel number.  The whole sky is pixel 0,
    // pixels 1, 2, 3, and 4 are each 1/4 of the sky (resolution 1), pixels 5-20 are each 
    // 1/16 of the sky (resolution 2), etc.
  
    base_pix=pix-pixel_start(res,scheme);
   
    m=base_pix % (int)(powl(2,res));
    n=(base_pix-m)/powl(2,res);
    azmin=TWOPI/powl(2,res)*m;
    azmax=TWOPI/powl(2,res)*(m+1);
    elmin=asinl(1-2.0/powl(2,res)*(n+1));
    elmax=asinl(1-2.0/powl(2,res)*n);
    
    angle[0]=azmin;
    angle[1]=azmax;
    angle[2]=elmin;
    angle[3]=elmax;

    pixel=new_poly(4);
    if (!pixel) {
      fprintf(stderr, "error in get_pixel: failed to allocate memory for polygon of 4 caps\n");
      return(0x0);
    }
    
    //    printf("az range: %Lf - %Lf, el range: %Lf - %Lf\n", azmin, azmax, elmin,elmax);
    rect_to_poly(angle,pixel);
    pixel->pixel=pix;
 
    if(!pixel){
      fprintf(stderr, "error in get_pixel: polygon is NULL.\n");
    }
    
    return(pixel);
    
  }
  else if(scheme=='d'){
    //this is the SDSSPix pixelization scheme; see http://lahmu.phyast.pitt.edu/~scranton/SDSSPix/
    //for more details

    assign_parameters();

    if(res==0){
      pixel=new_poly(0);
      pixel->weight=1.;
      pixel->pixel=0;

      return(pixel);
    }

    else if(res==1){
      ier = get_child_pixels(pix, pix_c, scheme);
      if(ier==1){
	fprintf(stderr, "error in get_pixel: get_child_pixels failed\n");
	return(0x0);
      }

      for(i=0;i<=3;i++){
        pix_bound(1, (unsigned long)pix_c[i] - (unsigned long)pixel_start(2, scheme), &lammin_c[i], &lammax_c[i], &etamin_c[i], &etamax_c[i]);
      }

      for(i=0;i<=3;i++){
        if(lammin_c[i]<=lammin_c[(i+1)%4] && lammin_c[i]<=lammin_c[(i+2)%4] && lammin_c[i]<=lammin_c[(i+3)%4]) lammin=lammin_c[i];
	if(lammax_c[i]>=lammax_c[(i+1)%4] && lammax_c[i]>=lammax_c[(i+2)%4] && lammax_c[i]>=lammax_c[(i+3)%4]) lammax=lammax_c[i];
	if(etamin_c[i]<=etamin_c[(i+1)%4] && etamin_c[i]<=etamin_c[(i+2)%4] && etamin_c[i]<=etamin_c[(i+3)%4]) etamin=etamin_c[i];
	if(etamax_c[i]>=etamax_c[(i+1)%4] && etamax_c[i]>=etamax_c[(i+2)%4] && etamax_c[i]>=etamax_c[(i+3)%4]) etamax=etamax_c[i];
      }

      /* fix pixels which span eta=0 */
      for(i=0;i<=5;i++){
	if(pix==9+18*i || pix==113){
	  etamin=etamin_c[0];
	  etamax=etamax_c[3];
	}
      }

      if(lammin < -90.) lammin = -90.;
      if(lammax > 90.) lammax = 90.;
      if(etamin < 0.) etamin += 360.;
      if(etamax > 360.) etamax -= 360.;
      lammin *= (PI/180.0);
      lammax *= (PI/180.0);
      etamin *= (PI/180.0);
      etamax *= (PI/180.0);

      angle[0]=etamin;
      angle[1]=etamax;
      angle[2]=lammin;
      angle[3]=lammax;

      pixel=new_poly(4);
      if (!pixel) {
        fprintf(stderr, "error in get_pixel: failed to allocate memory for polygon of 4 caps\n");
        return(0x0);
      }

      rect_to_poly(angle,pixel);
      pixel->pixel=pix;

      if(!pixel){
        fprintf(stderr, "error in get_pixel: polygon is NULL.\n");
        return(0x0);
      }

      for(i=0;i<pixel->np;i++){
	rp_to_azel(pixel->rp[i], &(v[i]));
	v[i].az *= (180./PI);
	v[i].el *= (180./PI);
	if(v[i].az < -180.) v[i].az += 360.;
	if(v[i].az > 180.) v[i].az -= 360.;
	if(v[i].el < -90.) v[i].el = -90.;
	if(v[i].el > 90.) v[i].el = 90.;
	csurvey2eq(v[i].el, v[i].az, &(v_r[i].az), &(v_r[i].el));
	if(v_r[i].az < 0.) v_r[i].az += 360.;
	if(v_r[i].az > 360.) v_r[i].az -= 360.;
	if(v_r[i].el < -90.) v_r[i].el = -90.;
	if(v_r[i].el > 90.) v_r[i].el = 90.;
	v_r[i].az *= (PI/180.);
	v_r[i].el *= (PI/180.);
	azel_to_rp(&(v_r[i]), pixel->rp[i]);
      }
      
      return(pixel);
    }

    else{
      pix_bound((int)powl(2,res-2), (unsigned long)pix - (unsigned long)pixel_start(res, scheme), &lammin, &lammax, &etamin, &etamax);

      if(lammin < -90.) lammin = -90.;
      if(lammax > 90.) lammax = 90.;
      if(etamin < 0.) etamin += 360.;
      if(etamax > 360.) etamax -= 360.;
      lammin *= (PI/180.0);
      lammax *= (PI/180.0);
      etamin *= (PI/180.0);
      etamax *= (PI/180.0);
  
      angle[0]=etamin;
      angle[1]=etamax;
      angle[2]=lammin;
      angle[3]=lammax;

      pixel=new_poly(4);
      if (!pixel) {
        fprintf(stderr, "error in get_pixel: failed to allocate memory for polygon of 4 caps\n");
        return(0x0);
      }

      rect_to_poly(angle,pixel);
      pixel->pixel=pix;
 
      if(!pixel){
        fprintf(stderr, "error in get_pixel: polygon is NULL.\n");
        return(0x0);
      }

      for(i=0;i<pixel->np;i++){
	rp_to_azel(pixel->rp[i], &(v[i]));
	v[i].az *= (180./PI);
	v[i].el *= (180./PI);
	if(v[i].az < -180.) v[i].az += 360.;
	if(v[i].az > 180.) v[i].az -= 360.;
	if(v[i].el < -90.) v[i].el = -90.;
	if(v[i].el > 90.) v[i].el = 90.;
	csurvey2eq(v[i].el, v[i].az, &(v_r[i].az), &(v_r[i].el));
	if(v_r[i].az < 0.) v_r[i].az += 360.;
	if(v_r[i].az > 360.) v_r[i].az -= 360.;
	if(v_r[i].el < -90.) v_r[i].el = -90.;
	if(v_r[i].el > 90.) v_r[i].el = 90.;
	v_r[i].az *= (PI/180.);
	v_r[i].el *= (PI/180.);
	azel_to_rp(&(v_r[i]), pixel->rp[i]);
      }
      
      return(pixel);
    }
  }
  else{
    fprintf(stderr, "error in get_pixel: pixel scheme %c not recognized.\n", scheme);
    return(0x0);   
  }
}

/* Function get_child_pixels takes a pixel number and calculates the numbers of the 
   child pixels of that pixel
   inputs: 
   pix_p: parent pixel number
   scheme: pixelization scheme
   outputs:
   pix_c[4]: array containing the pixel numbers of the 4 child pixels
   returns 0 on success, 1 if an error occurs
*/

int get_child_pixels(int pix_p, int pix_c[], char scheme){
  int mp,np,res,base_pix,i;
  unsigned long pix_c0, pix_c1, pix_c2, pix_c3;

  if(pix_p<0){
    fprintf(stderr, "error in get_child_pixels: %d is not a valid pixel number\n",pix_p);
    return(1);
  }

  res=get_res(pix_p, scheme);
  //printf("get_res(pix_p = %d) = %d\n", pix_p, res);
  if(res==-1) return (1);

  if(scheme=='s'){
    // this scheme divides up the sphere by rectangles in az and el, and is numbered 
    // such that the resolution is encoded in each pixel number.  The whole sky is pixel 0,
    // pixels 1, 2, 3, and 4 are each 1/4 of the sky (resolution 1), pixels 5-20 are each 
    // 1/16 of the sky (resolution 2), etc.
    
    base_pix=pix_p-pixel_start(res,scheme);
    mp=base_pix % (int)(powl(2,res));
    np=(base_pix-mp)/powl(2,res);
    
    //child pixels will have nc=2*np or 2*np+1, mc=2*mp or 2*mp+1, res_c=res+1
    //for first child pixel (nc=2*np, mc=2*mp), the base pixel number is given by
    //base_pix_c = 2^res_c*nc+mc = 2^(res+1)*2*np+2*mp=2^res*4*np+2*mp
    //combine this with base_pix_p=2^res*np+mp and extra resolution term 4^res 
    //to get formula for the number for the first child pixel number below
    
    pix_c[0]=pix_p+powl(4,res)+powl(2,res)*3*np+mp;
    pix_c[1]=pix_c[0]+1;
    pix_c[2]=pix_c[0]+powl(2,res+1);
    pix_c[3]=pix_c[2]+1;
    return(0);
  }
  else if(scheme=='d'){
    assign_parameters();
    if (pix_p==0){
      for(i=0;i<=116;i++) pix_c[i]=i+1;
      return(0);
    }
    else if (pix_p>=1 && pix_p<=108){
      for(i=1;i<=91;i+=18){
	if(pix_p>=i && pix_p<=i+17) {
	  pix_c[0]=4*(i-1)+2*(pix_p-i)+pixel_start(res+1,scheme);
	  break;
	}
      }
      pix_c[1]=pix_c[0]+1;
      pix_c[2]=pix_c[0]+36;
      pix_c[3]=pix_c[0]+37;
      return(0);
    }
    else if (pix_p>=109 && pix_p<=117){
      pix_c[0]=5*pix_p+(114-pix_p);
      for(i=1;i<=3;i++) pix_c[i]=pix_c[0]+i;
      return(0);
    }
    else {
      subpix((int)powl(2,res-2), (unsigned long)(pix_p-pixel_start(res, scheme)), &pix_c0, &pix_c1, &pix_c2, &pix_c3);
      pix_c[0] = (int)pix_c0 + pixel_start(res+1, scheme);
      pix_c[1] = (int)pix_c1 + pixel_start(res+1, scheme);
      pix_c[2] = (int)pix_c2 + pixel_start(res+1, scheme);
      pix_c[3] = (int)pix_c3 + pixel_start(res+1, scheme);
      return(0);
    }
  }
  else{
    fprintf(stderr, "error in get_child_pixels: pixel scheme %c not recognized.\n", scheme);
    return(1);   
  }
}

/* Function get_res takes a pixel number and returns the resolution 
   implied by that pixel number.
   inputs:
   pix: pixel number
   scheme: pixelization scheme
   returns the resolution of the pixel, or -1 if an error occurs
*/

int get_res(int pix, char scheme){
  int res;
  
  if(pix<0){
    fprintf(stderr, "error in get_res: %d not a valid pixel number.\n", pix);
    return(-1);
  }

  if(scheme=='s'){
    for(res=0;pix>=powl(4,res);res++){
      pix-=(int)powl(4,res);
    }
    return(res);
  }
  /* else if(scheme=='h'){
    if(pix==0){
      res=0;
      return(res);
    }
    else if(pix>0){
      for(res=2;pix>(int)(12*(powl(4,res)-4)/12);res++){
      }
    return((int)(res-1));  
    }
    } */
  else if(scheme=='d'){
    if(pix==0) return(0);
    else if(pix>=1 && pix <=117) return(1);
    else pix-=117;
    for(res=2;pix>(int)powl(4,res-2)*468;res++){
      pix-=(int)powl(4,res-2)*468;
    }
    // sdss_res increases by factors of 2 instead of increments of 1
    // sdss_res = (int)powl(2,res-1);
    return(res);
  }
  else{
    fprintf(stderr, "error in get_res: pixel scheme %c not recognized.\n", scheme);
    return(-1);  
  } 
  return(-1);
}
