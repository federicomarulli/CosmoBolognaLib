#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "manglefn.h"
#include "defaults.h"

#define ARGLEN 10
/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

int main(int argc, char *argv[])
{
  int ifile, ipoly, nfiles, npoly;
   
  int i, pixel, res, n, m, pixel_num;
  long double ra, dec;
  char scheme;
  int *child_pix;
  int children;
  int *parent_pix;
   polygon **polys;
   polys=polys_global;

  /* default output format */
  fmt.out = keywords[POLYGON];
  
  
  /*    
  if(argc<5){
    msg("enter the arguments for which_pixel as command line arguments:\n ra, dec, resolution, and pixelization scheme.\n");
    exit(1);
  } 
  else{
    ra=atof(argv[1]);
    dec=atof(argv[2]);
    res=atoi(argv[3]);
    scheme=argv[4][0];
    
    scale(&ra, 'd', 'r');
    scale(&dec, 'd','r');
    
    pixel = which_pixel(ra,dec,res,scheme);
    printf("pixel=%i\n",pixel);  
    return(0);
  }
  */
  
  
  /*  
  if(argc<3){
    msg("enter the arguments for get_child_pixels as command line arguments:\n pixel number and pixelization scheme.\n");
    exit(1);
  }
  else{
    pixel_num=atoi(argv[1]);
    scheme=argv[2][0];

    //allocate memory for child_pix array
    if(scheme=='d' && pixel_num==0){
      child_pix=(int *) malloc(sizeof(int) * 117);
      children=117;
      if(!child_pix){
	fprintf(stderr, "get_child_pixels: failed to allocate memory for %d integers\n", 117);
	return(-1);
      }
    }
    else{
      child_pix=(int *) malloc(sizeof(int) * 4);
      children=4;
      if(!child_pix){
	fprintf(stderr, "get_child_pixels: failed to allocate memory for %d integers\n", 4);
	return(-1);
      }
    }

    get_child_pixels(pixel_num,child_pix,scheme);
    printf("parent pixel = %d\n", pixel_num);

    for(i=0;i<children;i++){
      printf("child pixel %d = %d\n", i+1, child_pix[i]);
    }
    return(0);
  }
  */  

  /*  
    if(argc<3){
      msg("enter the arguments for get_parent_pixels as command line arguments:\n pixel number and pixelization scheme.\n");
      exit(1);
    }
    else{
      pixel_num=atoi(argv[1]);
      scheme=argv[2][0];
      res=get_res(pixel_num, scheme);
      printf("res=%d\n",res);
      
      if(pixel_num==0 && scheme=='d'){
        parent_pix = (int *) malloc(sizeof(int) * (168));
	if (!parent_pix){
	  fprintf(stderr, "test: failed to allocate memory for 168 integers\n");
	  exit(1);
        }
      }
      else{
        parent_pix = (int *) malloc(sizeof(int) * (res+1));
        if (!parent_pix) {
	  fprintf(stderr, "test: failed to allocate memory for %d integers\n", res);
	  exit(1);
        }
      }
      get_parent_pixels(pixel_num,parent_pix,scheme);
      printf("child pixel = %d\nparent pixels =", pixel_num);
      for(i=0;i<res;i++){
	printf(" %d, ",parent_pix[i]);
      }
      printf("and %d\n",parent_pix[res]);
      free(parent_pix);
      return(0);
    }
  
  */

  
   
   if(argc<4){
    msg("enter as command line arguments:\n resolution, pixelization scheme, and name of output file\n");
    return(1);
  }
  else{
    res=atoi(argv[1]);
    scheme=argv[2][0];
    
    if(scheme=='s'){
      npoly=powl(4,res);
    }

    if(scheme=='d'){
      if(res==0) npoly=1;
      else if(res==1) npoly=117;
      else{
	npoly=468*powl(4,res-2);
      }
    }

    npoly=1;
    polys[0]=get_pixel(400, scheme);

    /*    for(ipoly=0;ipoly<npoly;ipoly++){
      //pixel_num=ipoly+(int)((powl(4,res)+1)/3);
      pixel_num=ipoly+pixel_start(res, scheme);
      polys[ipoly]=get_pixel(pixel_num,scheme);
      
      m=ipoly % (int)(powl(2,res));
      n=(ipoly-m)/powl(2,res);
      
      // polys[ipoly]->weight=(n+m) % 2;
      polys[ipoly]->id=ipoly;
      polys[ipoly]->weight=(long double)ipoly/(long double)npoly;
    }
    */
    
    ifile = argc - 1;
    
    advise_fmt(&fmt);
    
    npoly = wrmask(argv[ifile], &fmt, npoly, polys);
    
    if (npoly == -1) exit(1);

    for(ipoly=0;ipoly<npoly;ipoly++){
    free_poly(polys[ipoly]);
    }

      
      return(0);
  }
   
    

}
