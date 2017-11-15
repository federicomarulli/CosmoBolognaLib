/*------------------------------------------------------------------------------
  © M E C Swanson 2005
  ------------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"
//#include <mcheck.h>

/* getopt options */
const char *optstr = "dqm:s:e:v:p:P:i:o:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);
#ifdef	GCC
int	pixelize(int npoly, polygon *[npoly], int npolys, polygon *[npolys]);
int	pixel_loop(int pix, int n, polygon *[n], int out_max, polygon *[out_max]);

#else
int	pixelize(int npoly, polygon *[/*npoly*/], int npolys, polygon *[/*npolys*/]);
int	pixel_loop(int pix, int n, polygon *[/*n*/], int out_max, polygon *[/*out_max*/]);

#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
  int ifile, nfiles, npoly, npolys, i, res_max_temp;
  char scheme_temp;
  polygon **polys;
  polys=polys_global;

  //   mtrace();

  /* default output format */
  fmt.out = keywords[POLYGON];
  /* default is to renumber output polygons with old id numbers */
  fmt.newid = 'o';

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
  
  msg("---------------- pixelize ----------------\n");
  
  // snap angles 
  scale(&axtol, axunit, 's');
  scale(&btol, bunit, 's');
  scale(&thtol, thunit, 's');
  axunit = 's';
  bunit = 's';
  thunit = 's';
  //  msg("snap angles: axis %Lg%c latitude %Lg%c edge %Lg%c\n", axtol, axunit, btol, bunit, thtol, thunit);
  scale(&axtol, axunit, 'r');
  scale(&btol, bunit, 'r');
  scale(&thtol, thunit, 'r');
  axunit = 'r';
  bunit = 'r';
  thunit = 'r';
  
  /* tolerance angle for multiple intersections */
  if (mtol != 0.) {
    scale(&mtol, munit, 's');
    munit = 's';
    msg("multiple intersections closer than %Lg%c will be treated as coincident\n", mtol, munit);
    scale(&mtol, munit, 'r');
    munit = 'r';
  }
  
  msg("pixelization scheme %c, maximum resolution %d\n", scheme, res_max);
  msg("maximum number of polygons allowed in each pixel: %d\n", polys_per_pixel);
  scheme_temp=scheme;
  res_max_temp=res_max;
  

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

  if(scheme!=scheme_temp || res_max!=res_max_temp){
    msg("warning: pixelization information in input file is being discarded\n");
    scheme=scheme_temp;
    res_max=res_max_temp;
  }

  /* pixelize polygons */
  npolys = pixelize(npoly, polys, NPOLYSMAX - npoly, &polys[npoly]);
  if (npolys == -1) exit(1);

  pixelized=1;
  if(polys_per_pixel>0) res_max=-1;

  /* write polygons */
  ifile = argc - 1;
  npolys = wrmask(argv[ifile], &fmt, npolys, &polys[npoly]);
  if (npolys == -1) exit(1);
  /* memmsg(); */

  for(i=0;i<npoly+npolys;i++){
    free_poly(polys[i]);
  }

  return(0);
}

/*------------------------------------------------------------------------------
 */
void usage(void)
{
  printf("usage:\n");
  //  printf("pixelize [-d] [-q] [-a<a>[u]] [-b<a>[u]] [-t<a>[u]] [-y<r>] [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn|-vp] [-p[+|-][<n>]] [-P[scheme][<r>][,<p>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
 printf("pixelize [-d] [-q] [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn|-vp] [-p[+|-][<n>]] [-P[scheme][<p>][,<r>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
 */
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Pixelize: split polygons against a pre-defined pixel map such that each polygon is in only one pixel
  
  Input: npoly = number of polygons.
  poly = array of pointers to polygons.
  npolys = maximum number of output polygons.
  Output: polys = array of pointers to polygons.
  Return value: number of disjoint connected polygons,
  or -1 if error occurred.
*/
int pixelize(int npoly, polygon *poly[/*npoly*/], int npolys, polygon *polys[/*npolys*/])
{
  /* part_poly should lasso one-boundary polygons only if they have too many caps */
#define ALL_ONEBOUNDARY		1
  /* how part_poly should tighten lasso */
#define ADJUST_LASSO		1
  /* part_poly should force polygon to be split even if no part can be lassoed */
#define	FORCE_SPLIT		1
  /* partition_poly should overwrite all original polygons */
#define OVERWRITE_ORIGINAL	2
#define WARNMAX			8
   char *snapped_polys = 0x0;
   int isnap,j, nadj;
  int dn, dnp, failed, i, ier, inull, ip, iprune, m, n, np;

  msg("pruning input polygons\n");

  /* start by pruning all input polygons */
  np = 0;
  inull = 0;
  for (i = 0; i < npoly; i++) {
    iprune = prune_poly(poly[i], mtol);
    /* error */
    if (iprune == -1) {
      fprintf(stderr, "pixelize: initial prune failed at polygon %d\n", poly[i]->id);
      return(-1);
    }
    /* zero area polygon */
    if (iprune >= 2) {
      if (WARNMAX > 0 && inull == 0) msg("warning from pixelize: the following polygons have zero area & are being discarded:\n");
      if (inull < WARNMAX) {
	msg(" %d", (fmt.newid == 'o')? poly[i]->id : i);
      } else if (inull == WARNMAX) {
	msg(" ... more\n");
      }
      inull++;
    } else {
      np++;
    }
  }
  if (WARNMAX > 0 && inull > 0 && inull <= WARNMAX) msg("\n");
  if (inull > 0) {
    msg("pixelize: %d polygons with zero area are being discarded.\n", inull);
  }

  /* number of polygons */
  msg("pixelizing %d polygons ...\n", np);

  /* set all input polygons to be in pixel 0 (whole sky)*/
  inull=0;
  for (i = 0; i < npoly; i++) {
    if(poly[i]->pixel!=0){
      poly[i]->pixel = 0;
      if(WARNMAX>0 && inull ==0)  msg("warning from pixelize: following polygons are being re-set to be in pixel 0:\n");
      if (inull < WARNMAX) {
	msg(" %d", (fmt.newid == 'o')? poly[i]->id : i);
      } else if (inull == WARNMAX) {
	msg(" ... more\n");
      }
      inull++;     		       
    }
  }
  if (WARNMAX > 0 && inull > 0 && inull <= WARNMAX) msg("\n");
  if (inull > 0) {
    msg("pixelize: %d polygons have been re-set to be in pixel 0.\n", inull);
  }

  /* nullify all output polygons */
  for (i = 0; i < npolys; i++) {
    polys[i] = 0x0;
  }

  msg("pixelize stage 1 (fragment each polygon so it is in only one pixel):\n");
  
  /*call recursive pixel_loop to split polygons into pixels*/
  n=pixel_loop(0,npoly,poly,npolys,polys);
  if(n==-1) return(-1);

  dnp=n-np;
  np=n;
  
  msg("added %d polygons to make %d\n", dnp, np);
  
  /* partition disconnected polygons into connected parts  */
  msg("pixelize stage 2 (partition disconnected polygons into connected parts):\n");
  m = n;
  dnp = 0;
  ip = 0;
  failed = 0;
  for (i = 0; i < m; i++) {
    /* skip null polygons */
    if (!polys[i] || (polys[i]->np > 0 && polys[i]->cm[0] == 0.)) continue;
    /* partition disconnected polygons */
    ier = partition_poly(&polys[i], npolys - n, &polys[n], mtol, ALL_ONEBOUNDARY, ADJUST_LASSO, FORCE_SPLIT, OVERWRITE_ORIGINAL, &dn);
    /* error */
    if (ier == -1) {
      fprintf(stderr, "pixelize: UHOH at polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : ip);
      continue;
      /* return(-1); */
      /* failed to partition polygon into desired number of parts */
    } else if (ier == 1) {
      fprintf(stderr, "pixelize: failed to partition polygon %d fully; partitioned it into %d parts\n", (fmt.newid == 'o')? polys[i]->id : ip, dn + 1);
      failed++;
    }
    /* increment index of next subset of fragments */
    n += dn;
    /* increment polygon count */
    np += dn;
    dnp += dn;
    /* check whether exceeded maximum number of polygons */
    if (n > npolys) {
      fprintf(stderr, "pixelize: total number of polygons exceeded maximum %d\n", npoly + npolys);
      fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
      return(-1);
    }
    ip++;
  }
  msg("added %d polygons to make %d\n", dnp, np);
  
  if (failed > 0) {
    msg("pixelize: failed to split %d polygons into desired number of connected parts\n", failed);
    msg(".............................................................................\n");
    msg("Failure to split polygon probably means:\n");
    msg("either (1) you forgot to run snap on all your input polygon files;\n");
    msg("    or (2) the polygon is too small for the numerics to cope with;\n");
    msg("    or (3) you have a weird-shaped polygon.\n");
    msg("You may ignore this warning message if the weights of polygons in the input\n");
    msg("polygon file(s) are already correct, and you do not want to reweight them.\n");
    msg("Similarly, you may ignore this warning message if you do want to reweight the\n");
    msg("polygons, but the weights of the different parts of each unsplit polygon are\n");
    msg("the same.  If you want to reweight the different parts of an unsplit polygon\n");
    msg("with different weights, then you will need to split that polygon by hand.\n");
    msg("Whatever the case, the output file of pixelized polygons constitutes\n");
    msg("a valid mask with each polygon in only one pixel, and is safe to use.\n");
    msg(".............................................................................\n");
  }
    
  /*
  
  // prune 
  msg("pruning ... \n");
  j = 0;
  inull = 0;
  for (i = 0; i < n; i++) {
    iprune = prune_poly(polys[i], mtol);
    if (iprune == -1) {
      fprintf(stderr, "pixelize: failed to prune polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : j);
       return(-1); 
    }
    if (iprune >= 2) {
      free_poly(polys[i]);
      polys[i] = 0x0;
      inull++;
    } else {
      polys[j] = polys[i];
      j++;
    }
  }
  if (inull > 0) msg("pixelize: %d pixelized polygons have zero area, and are being discarded\n", inull);
  n = j;

   
      
  //allocate snapped_polys array 
  snapped_polys = (char *) malloc(sizeof(char) * n);
  if (!snapped_polys) {
    fprintf(stderr, "pixelize: failed to allocate memory for %d characters\n", n);
    return(-1);
  }
  
  //snap edges of each polygon
  selfsnap = 1;
  nadj = snap_polys(&fmt, n, polys, selfsnap, axtol, btol, thtol, ytol, mtol, WARNMAX, snapped_polys);
  if(nadj==-1){
    msg("pixelize: error snapping pixelized polygons\n");
    return(-1);
  }

  //number of polygons whose edges were snapped
  isnap = 0;
  for (i = 0; i < n; i++) if (snapped_polys[i]) isnap++;
  if (isnap > 0) msg("pixelize: edges of %d pixelized polygons were snapped\n", isnap);
  
  //prune snapped polygons
  j = 0;
  inull = 0;
  for (i = 0; i < n; i++) {
    if (snapped_polys[i]) {
      iprune = prune_poly(polys[i], mtol);
      if (iprune == -1) {
	fprintf(stderr, "pixelize: failed to prune polygon %d; continuing ...\n", (fmt.newid == 'o')? polys[i]->id : j);
	// return(-1);
      }
      if (iprune >= 2) {
	free_poly(polys[i]);
	polys[i] = 0x0;
	inull++;
      } else {
	polys[j] = polys[i];
	j++;
      }
    } else {
      polys[j] = polys[i];
      j++;
    }
  }
  if (inull > 0) msg("pixelize: %d snapped polygons have zero area, and are being discarded\n", inull);
  n = j;
  
  //free snapped_polys array
  free(snapped_polys);
  */

  if(n!=-1){
    /* sort polygons by pixel number */
    poly_sort(n, polys, 'p');
  }
  
  /* assign new polygon id numbers in place of inherited ids */
  if (fmt.newid == 'n') {
    for (i = 0; i < n; i++) {
      polys[i]->id = i;
    }
  }

  if (fmt.newid == 'p') {
    for (i = 0; i < n; i++) {
      polys[i]->id = polys[i]->pixel;
    }
  }


  return(n);
}

/*
  Function pixel_loop takes a list of all of the polygons in the input pixel and then splits 
  them into the four child pixels of the input pixel.  It then recursively calls itself on
  each of the child pixels, until the desired level of pixelization is reached.
  Inputs: 
  pix: input pixel number
  n = number of polygons.
  input = array of pointers to polygons.
  out_max = maximum number of output polygons.
  Output:
  output = array of pointers to polygons.
  Return value: number of polygons written to output array,
  or -1 if error occurred.
*/

int pixel_loop(int pix, int n, polygon *input[/*n*/], int out_max, polygon *output[/*out_max*/]){
  int *child_pix,children;
  int i,j,k,m,out,nout;
  int ier, iprune, np;
  polygon *pixel;
  polygon **poly;
 
  //allocate memory for work array of polygon pointers
  poly=(polygon **) malloc(sizeof(polygon *) * n);
  if(!poly){
    fprintf(stderr, "pixel_loop: failed to allocate memory for %d polygon pointers\n",n);
    return(-1);
  }
  
  // allocate memory for child_pix array
  if(pix==0 && scheme=='d'){
    child_pix=(int *) malloc(sizeof(int) * 117);
    children=117;
    if(!child_pix){
      fprintf(stderr, "pixel_loop: failed to allocate memory for 117 integers\n");
      return(-1);
    }
  }
  else{
    child_pix=(int *) malloc(sizeof(int) * 4);
    children=4;
    if(!child_pix){
      fprintf(stderr, "pixel_loop: failed to allocate memory for %d integers\n", 4);
      return(-1);
    }
  }

  get_child_pixels(pix, child_pix, scheme);
  out=0;
  
  for(i=0;i<children;i++){    
    /*get the current child pixel*/
    pixel=get_pixel(child_pix[i], scheme);
 
    if(!pixel){
      fprintf(stderr, "error in pixel_loop: could not get pixel %d\n", child_pix[i]); 
      return(-1);
    }
    
    /*loop through input polygons to find the ones that overlap with current child pixel*/
    for(j=0;j<n;j++){
      /* skip null polygons */
      if (input[j]->np > 0 && input[j]->cm[0] == 0.){
	poly[j] = 0x0;
	continue;
      }

      np=input[j]->np+pixel->np;
      poly[j]=new_poly(np);
      if(!poly[j]){
	fprintf(stderr, "error in pixel_loop: failed to allocate memory for polygon of %d caps\n", np);
	return(-1);
      }
      /*set poly[j] to the intersection of input[j] and current child pixel*/
      poly_poly(input[j],pixel,poly[j]);
      poly[j]->pixel=pixel->pixel;
    
      iprune = prune_poly(poly[j], mtol);
      if (iprune == -1) {
	fprintf(stderr, "pixelize: failed to prune polygon for pixel %d; continuing ...\n", poly[j]->pixel);
	//return(-1);
      }
      /*if polygon is null, get rid of it*/
      if (iprune >= 2) {
	free_poly(poly[j]);
	poly[j] = 0x0;
      }       
    }
    
    /*copy down non-null polygons*/
    k=0;
    for(j=0;j<n;j++){
      if(poly[j]){
	poly[k++]=poly[j];
      }
    }
    m=k;
    /*nullify the rest of the array, but don't free, since pointers have been copied above*/
    for(j=m;j<n;j++){
      poly[j]=0x0;
    }
    
    /*if we're below the max resolution, recursively call pixel_loop on the current child pixel */ 
    if(m>polys_per_pixel && get_res(child_pix[i],scheme)<res_max){
      //printf("calling pixel loop for pixel %d with %d polygons\n",child_pix[i],m);
      nout=pixel_loop(child_pix[i],m,poly,out_max-out,&output[out]);
      if(nout==-1) return(-1);
      out+=nout;
    }
    else{
      for(k=0;k<m;k++){
	/* check whether exceeded maximum number of polygons */
	if (out >= out_max) {
	  fprintf(stderr, "pixel_loop: total number of polygons exceeded maximum %d\n", NPOLYSMAX);
	  fprintf(stderr, "if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
	  return(-1);
	}
	/*make sure output polygon has enough room */
	ier = room_poly(&output[out], poly[k]->np, DNP, 0);
	if (ier == -1) {
	  fprintf(stderr, "error in pixel_loop: failed to allocate memory for polygon of %d caps\n", poly[i]->np + DNP);
	  return(-1);
	}

	/*copy polygon to output array*/
	copy_poly(poly[k],output[out]);
	out++;    
      }
    }
    /*free up memory for next child pixel*/
    
    free_poly(pixel);
    for(j=0;j<n;j++){
      free_poly(poly[j]);
    }
  }
  free(child_pix);
  free(poly);
  return out;
}
