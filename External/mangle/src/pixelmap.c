/*------------------------------------------------------------------------------
© M E C Swanson 2006
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqm:s:e:v:p:P:i:o:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);
#ifdef	GCC
int	pixelmap(int *npoly, polygon *[*npoly]);
#else
int	pixelmap(int *npoly, polygon *[/**npoly*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nadj, nfiles, npoly, npolys, res_max_temp,i;
    char scheme_temp;
    polygon **polys;
    polys=polys_global;

    /* default output format */
    fmt.out = keywords[POLYGON];
    /* default is to renumber output polygons with pixel numbers as id numbers */
    fmt.newid = 'p';

    /* parse arguments */
    parse_args(argc, argv);

    /* at least one input and output filename required as arguments */
    if (argc - optind < 2) {
	if (optind > 1 || argc - optind == 1) {
	    fprintf(stderr, "%s requires at least 2 arguments: polygon_infile and polygon_outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- pixelmap ----------------\n");

    /* tolerance angle for multiple intersections */
    if (mtol != 0.) {
	scale(&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %Lg%c will be treated as coincident\n", mtol, munit);
	scale(&mtol, munit, 'r');
	munit = 'r';
    }

    /* save res_max as defined on command line rather than using value from file*/
    /* value of scheme in file will override scheme defined on command line */
    res_max_temp=res_max;
    scheme_temp=scheme;

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 1 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &polys[npoly]);
	if (npolys == -1) exit(1);
	npoly += npolys;
    }
    if (nfiles >= 2) {
        msg("total of %d polygons read\n", npoly);
    }
    if (npoly == 0) {
	msg("STOP\n");
	exit(0);
    }
    
    res_max=res_max_temp;
    msg("pixelization scheme %c, making map at resolution %d\n", scheme, res_max);
    
    if (snapped==0 || balkanized==0) {
      fprintf(stderr, "Error: input polygons must be snapped and balkanized before using pixelmap.\n");
      fprintf(stderr, "If your polygons are already snapped and balkanized, add the 'snapped' and\n'balkanized' keywords at the beginning of each of your input polygon files.\n");
      exit(1);
    }

    /* pixelmap polygons */
    nadj = pixelmap(&npoly, polys);
    if (nadj == -1) exit(1);

    ifile = argc - 1;
    npoly = wrmask(argv[ifile], &fmt, npoly, polys);
    if (npoly == -1) exit(1);
 
    for(i=0;i<npoly;i++){
      free_poly(polys[i]);
    }

    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("pixelmap [-d] [-q] [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn|-vp] [-p[+|-][<n>]] [-P[scheme][<p>][,<r>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Take pixelized polygons, find the average weight within each pixel, and return a set of polygons consisting of the pixels weighted with the average weight.

   Input: poly = array of pointers to polygons.
	  npoly = pointer to number of polygons.
  Output: polys = array of pointers to polygons;
  Return value: number of polygons discarded by pixelmapping,
		or -1 if error occurred.
*/
int pixelmap(int *npoly, polygon *poly[/**npoly*/])
{
  int i, j, nadj, k, kstart,kend,numpix;
  int *start;
  int *total;
  int *parent_pixels;
  int begin, end, p,max_pixel,min_pixel, ier, verb,res1,res2;
  long double tol,area, tot_area;
  long double *av_weight;
  long double *av_weight0;

  poly_sort(*npoly,poly,'p');
  min_pixel = poly[0]->pixel; 
  max_pixel = poly[*npoly-1]->pixel+1; 
  res1=get_res(min_pixel,scheme);
  res2=get_res(max_pixel,scheme);

  if(res1<res_max){
    fprintf(stderr,"pixelmap: there are pixels in the mask with a lower resolution than the desired pixelmap resolution %d.  The desired pixelmap resolution can be set with the -P option.\n",res_max);
    fprintf(stderr,"Before using pixelmap, use pixelize with the -P0,r option to pixelize the entire mask to the desired resolution r.\n");
    return(-1);
  }

 
  /* allocate memory for pixel info arrays start and total */ 
  start = (int *) malloc(sizeof(int) * max_pixel);
  if (!start) {
    fprintf(stderr, "pixelmap: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }
  total = (int *) malloc(sizeof(int) * max_pixel);
  if (!total) {
    fprintf(stderr, "pixelmap: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }

  /* build lists of starting indices of each pixel and total number of polygons in each pixel*/
  ier=pixel_list(*npoly, poly, max_pixel, start, total);
  if (ier == -1) {
    fprintf(stderr, "pixelmap: error building pixel index lists\n");
    return(-1);
  } 

  //allocate memory for parent pixels array
  parent_pixels = (int *) malloc(sizeof(int) * (res2+1));
  if (!parent_pixels) {
    fprintf(stderr, "pixelmap: failed to allocate memory for %d integers\n", res2+1);
    return(-1);
  }

  //kstart=number of first pixel at desired output resolution
  //kend=number of last pixel at desired output resolution
  if(res_max==-1){
    kstart=pixel_start(res1,scheme);
    kend=pixel_start(res2+1,scheme)-1;
  }
  else{
    kstart=pixel_start(res_max,scheme);
    kend=pixel_start(res_max+1,scheme)-1;
  }
  av_weight0= (long double *) malloc(sizeof(long double) * (kend-kstart+1) );
  if (!av_weight0) {
    fprintf(stderr, "pixelmap: failed to allocate memory for %d integers\n", kend-kstart+1 );
    return(-1);
  }
  //make av_weight an array indexed by the pixel number 
  av_weight=av_weight0-kstart;

  //set av_weight array to 0 initially
  for(k=kstart;k<=kend;k++){
    av_weight[k]=0;
  }

  nadj = 0;
  verb=1;
   
  /*find average weight of polygons within each pixel*/
  for(p=min_pixel;p<max_pixel;p++){
    begin=start[p];
    end=start[p]+total[p];
    ier=get_parent_pixels(p,parent_pixels,scheme);
    if(ier) return(-1);
    
    //set k to the pixel at the desired output resolution, or to the pixel number if using
    //existing resolution

    k=(res_max==-1) ? p : parent_pixels[res_max];
    
    for (i = begin; i < end; i++) {
      if (!poly[i]) continue;
      tol=mtol;
      ier = garea(poly[i], &tol, verb, &area);
      if(ier==1 || ier == -1){
	fprintf(stderr, "error %d in garea in polygon %d\n", ier, poly[i]->id);
	continue;
      }
      
      av_weight[k]+=poly[i]->weight * area;
    }
  }
  
  //replace polygons in input array with non-zero weight pixels
  j=0;
  for(k=kstart;k<=kend;k++){
    if(av_weight[k]==0) continue;
    free_poly(poly[j]);
    poly[j]=get_pixel(k,scheme);
    tol=mtol;
    ier = garea(poly[j], &tol, verb, &tot_area);
    if(ier==1 || ier == -1){
      fprintf(stderr, "pixelmap: error in garea in pixel %d\n",p);
      continue;
    }
    poly[j]->weight=av_weight[k]/tot_area;
    j++;
    if(j> *npoly ){
      fprintf(stderr,"pixelmap: number of pixels with non-zero weight exceeds number of polygons.\n");
      fprintf(stderr, "Try running unify on your mask to remove zero-weight polygons before using pixelmap.\n");
    }
  }

  numpix=j;
  
  for(j=numpix; j< *npoly; j++){
    free_poly(poly[j]);
    poly[j] = 0x0;
    nadj++;
  }

  *npoly=numpix;

  free(start);
  free(total);
  free(parent_pixels);
  free(av_weight0);
  
  /* assign new polygon id numbers */
  if (fmt.newid == 'n') {
    for (i = 0; i < *npoly; i++) {
      poly[i]->id = i;
    }
  }
  
  if (fmt.newid == 'p') {
    for (i = 0; i < *npoly; i++) {
      poly[i]->id = poly[i]->pixel;
    }
  }
  
  /* advise */
  msg("pixelmap: %d pixels in map\n", numpix);

  return(nadj);
}
