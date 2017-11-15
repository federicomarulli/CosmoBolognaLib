/*--------------------------------------------------------------------
(C) J C Hill 2006
--------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include "pi.h"
#include "manglefn.h"
#include "defaults.h"

/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP             4

/* getopt options */
const char *optstr = "dqa:b:t:y:m:s:e:p:i:o:H";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void     usage(void);
#ifdef  GCC
int     rasterize(int nhealpix_poly, int npoly, polygon *[npoly], int nweights, long double [nweights]);
#else
int     rasterize(int nhealpix_poly, int npoly, polygon *[/*npoly*/], int nweights, long double [/*nweights*/]);
#endif

/*--------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
  int ifile, nfiles, npoly, npolys, nhealpix_poly, nhealpix_polys, j, k, nweights, nweight;
  long double *weights;
  
  polygon **polys;
  polys=polys_global;

  /* default output format */
  //fmt.out = keywords[HEALPIX_WEIGHT];
  fmt.out = keywords[POLYGON];

  /* parse arguments */
  parse_args(argc, argv);

  /* at least two input and one output filenames required as arguments */
  if (argc - optind < 3) {
      if (optind > 1 || argc - optind == 1 || argc - optind == 2) {
         fprintf(stderr, "%s requires at least 3 arguments: polygon_infile1, polygon_infile2, and polygon_outfile\n", argv[0]);
         usage();
         exit(1);
     } else {
         usage();
         exit(0);
     }
  }

  msg("---------------- rasterize ----------------\n");

  /* snap angles */
  scale(&axtol, axunit, 's');
  scale(&btol, bunit, 's');
  scale(&thtol, thunit, 's');
  axunit = 's';
  bunit = 's';
  thunit = 's';
  msg("snap angles: axis %Lg%c latitude %Lg%c edge %Lg%c\n", axtol, axunit, btol, bunit, thtol, thunit);
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

  /* advise data format */
  advise_fmt(&fmt);

  /* read polygons from polygon_infile1 (healpix pixels, or some other 'rasterizer' pixels) */
  /* the id numbers of these polygons should match the pixel numbers of this pixelization scheme;
     for example, if you are using HEALPix, the id numbers should match the HEALPix pixel numbers
     in the NESTED scheme */
  nhealpix_poly = 0;
  ifile = optind;
  nhealpix_polys = rdmask(argv[ifile], &fmt, NPOLYSMAX - nhealpix_poly, &polys[nhealpix_poly]);
  if (nhealpix_polys == -1) exit(1);
  nhealpix_poly += nhealpix_polys;

  if (nhealpix_poly == 0) {
     msg("STOP\n");
     exit(0);
  }

  /* Input rasterizer polygons need not be balkanized if they are non-overlapping by construction,
     which is the case for the HEALPix polygons.  This is a special case - all other mangle functions
     that require balkanization require all input files to be balkanized.  To avoid getting an error
     here, increment the 'balkanized' counter here if the rasterizer polygons are not balkanized. */
  if (balkanized == 0) {
    balkanized++;
  }

  /* set nweights equal to maximum id number in rasterizer file */
  nweights = 0;
  for (k = 0; k < nhealpix_poly; k++) {
    if (polys[k]->id > nweights) nweights = polys[k]->id;
  }

  /* read polygons from polygon_infile2, polygon_infile3, etc. */
  npoly = nhealpix_poly;
  nfiles = argc - 2 - optind;
  for (ifile = optind + 1; ifile < optind + 1 + nfiles; ifile++) {
      npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &polys[npoly]);
      if (npolys == -1) exit(1);
      npoly += npolys;
  }
  if (nfiles >= 2) {
    msg("total of %d polygons read from mask files\n", npoly-nhealpix_poly);
  }
  if (npoly-nhealpix_poly == 0) {
    msg("STOP\n");
    exit(0);
  }

    if (snapped==0 || balkanized==0) {
      fprintf(stderr, "Error: input polygons must be snapped and balkanized before rasterization.\n");
      fprintf(stderr, "If your polygons are already snapped and balkanized, add the 'snapped' and\n'balkanized' keywords at the beginning of each of your input polygon files.\n");
      exit(1);
    }


  /* allocate memory for weights array */
  weights = (long double *) malloc(sizeof(long double) * (nweights));
  if (!weights) {
     fprintf(stderr, "rasterize: failed to allocate memory for %d long doubles\n", nweights);
     exit(1);
  }

  /* initialize weights array to 0 */
  for (k = 0; k < nweights; k++) weights[k] = 0.;

  /* rasterize */
  nweight = rasterize(nhealpix_poly, npoly, polys, nweights, weights);
  if (nweight == -1) exit(1);

  /* copy new weights to original rasterizer polygons */
  for (k = 0; k < nhealpix_poly; k++) {
    for (j = 0; j < nweights; j++) {
      if (polys[k]->id == j+1) {
	polys[k]->weight = weights[j];
	break;
      }
    }
  }


  ifile = argc - 1;
  if (strcmp(fmt.out, "healpix_weight") == 0) {
    nweight = wr_healpix_weight(argv[ifile], &fmt, nweights, weights);
    if (nweight == -1) exit(1);
  }
  else {
    nweight = wrmask(argv[ifile], &fmt, nhealpix_poly, polys);
    if (nweight == -1) exit(1);
  }

  /* free array */
  for(k = 0; k < npoly; k++){
    free_poly(polys[k]);
  }

  return(0);

}

/*-------------------------------------------------------------------------
*/
void usage(void)
{
     printf("usage:\n");
     printf("rasterize [-d] [-q] [-a<a>[u]] [-b<a>[u]] [-t<a>[u]] [-y<r>] [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn] [-p[+|-][<n>]] [-i<f>[<n>][u]] [-o<f>[u]] [-H] polygon_infile1 polygon_infile2 [polygon_infile3 ...] polygon_outfile\n");
#include "usage.h"
}

/*-------------------------------------------------------------------------
*/
#include "parse_args.c"

/*-------------------------------------------------------------------------
  Rasterize a mask of input polygons against a mask of rasterizer polygons.

  Input: nhealpix_poly = number of rasterizer polygons.
         npoly = total number of polygons in input array.
	 polys = array of pointers to polygons.
	 nweights = number of weights in output array.
  Output: weights = array of rasterizer weights.
  Return value: number of weights in array,
                or -1 if error occurred.
*/

int rasterize(int nhealpix_poly, int npoly, polygon *polys[/*npoly*/], int nweights, long double weights[/*nweights*/])
{
  int min_pixel, max_pixel, ier, ier_h, ier_i, i, j, ipix, ipoly, begin_r, end_r, begin_m, end_m, verb, np, iprune;
  int *start_r, *start_m, *total_r, *total_m;
  long double *areas, area_h, area_i, tol;

  static polygon *polyint = 0x0;

  /* make sure weights are all zero for rasterizer pixels */
  for (i = 0; i < nhealpix_poly; i++) {
      polys[i]->weight = 0.;
  }

  /* allocate memory for rasterizer areas array */
  areas = (long double *) malloc(sizeof(long double) * (nweights));
  if (!areas) {
    fprintf(stderr, "rasterize: failed to allocate memory for %d long doubles\n", nweights);
    exit(1);
  }

  /* initialize rasterizer areas array to 0 */
  for (i = 0; i < nweights; i++) areas[i] = 0.;

  /* allow error messages from garea */
  verb = 1;

  /* find areas of rasterizer pixels for later use */
  for (i = 1; i <= nweights; i++) {
    for (j = 0; j < nhealpix_poly; j++) {
      if (polys[j]->id == i) {
        tol = mtol;
        ier_h = garea(polys[j], &tol, verb, &area_h);
        if (ier_h == 1) {
          fprintf(stderr, "fatal error in garea\n");
          exit(1);
        }
        if (ier_h == -1) {
          fprintf(stderr, "failed to allocate memory in garea\n");
          exit(1);
        }
        areas[i-1] += area_h;
      }
    }
  }

  /* sort arrays by pixel number */
  poly_sort(nhealpix_poly, polys, 'p');
  poly_sort(npoly-nhealpix_poly, &(polys[nhealpix_poly]), 'p');

  /* allocate memory for pixel info arrays start_r, start_m, total_r, and total_m */
  min_pixel = polys[0]->pixel;
  max_pixel = (polys[nhealpix_poly-1]->pixel+1>polys[npoly-1]->pixel+1)?(polys[nhealpix_poly-1]->pixel+1):(polys[npoly-1]->pixel+1);
  start_r = (int *) malloc(sizeof(int) * max_pixel);
  if (!start_r) {
    fprintf(stderr, "rasterize: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }
  start_m = (int *) malloc(sizeof(int) * max_pixel);
  if (!start_m) {
    fprintf(stderr, "rasterize: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }
  total_r = (int *) malloc(sizeof(int) * max_pixel);
  if (!total_r) {
    fprintf(stderr, "rasterize: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }
  total_m = (int *) malloc(sizeof(int) * max_pixel);
  if (!total_m) {
    fprintf(stderr, "rasterize: failed to allocate memory for %d integers\n", max_pixel);
    return(-1);
  }

  /* build lists of starting indices of each pixel and total number of polygons in each pixel */
  ier = pixel_list(nhealpix_poly, polys, max_pixel, start_r, total_r);
  if (ier == -1) {
    fprintf(stderr, "rasterize: error building pixel index lists for rasterizer polygons\n");
    return(-1);
  }

  ier = pixel_list(npoly-nhealpix_poly, &(polys[nhealpix_poly]), max_pixel, start_m, total_m);
  if (ier == -1) {
    fprintf(stderr, "rasterize: error building pixel index lists for input mask polygons\n");
    return(-1);
  }

  /* correction due to the start_m array's offset */
  for (i = min_pixel; i < max_pixel; i++) {
    start_m[i] += nhealpix_poly;
  }

  /* compute intersection of each input mask polygon with each rasterizer polygon */
  for (ipix = min_pixel; ipix < max_pixel; ipix++) {
    begin_r = start_r[ipix];
    end_r = start_r[ipix] + total_r[ipix];
    begin_m = start_m[ipix];
    end_m = start_m[ipix] + total_m[ipix];

    for (ipoly = begin_m; ipoly < end_m; ipoly++) {
      /* disregard any null polygons */
      if (!polys[ipoly]) continue;

      for (i = begin_r; i < end_r; i++) {

	/* make sure polyint contains enough space for intersection */
	np = polys[ipoly]->np + polys[i]->np;
	ier = room_poly(&polyint, np, DNP, 0);
	if (ier == -1) goto out_of_memory;

	poly_poly(polys[ipoly], polys[i], polyint);

	/* suppress coincident boundaries, to make garea happy */
	iprune = trim_poly(polyint);

	/* intersection of polys[ipoly] and polys[i] is null polygon */
	if (iprune >= 2) area_i = 0.;

	else {
	  tol = mtol;
	  ier_i = garea(polyint, &tol, verb, &area_i);
	  if (ier_i == 1) {
	    fprintf(stderr, "fatal error in garea\n");
	    return(-1);
	  }
	  if (ier_i == -1) {
	    fprintf(stderr, "failed to allocate memory in garea\n");
	    return(-1);
	  }
	}
	    
	weights[(polys[i]->id)-1] += (area_i)*(polys[ipoly]->weight);
      }
    }
  }

  for (i=0; i<nweights; i++) {
    if(areas[i]!=0){
      weights[i] = weights[i]/areas[i];
    }
    else{
      weights[i]=0;
      fprintf(stderr,"WARNING: rasterize: area of rasterizer polygon %d is zero.  Assigning zero weight.\n",i);
    }
  }

  return(i+1);

  /* ----- error return ----- */
  out_of_memory:
  fprintf(stderr, "rasterize: failed to allocate memory for polygon of %d caps\n", np + DNP);
  return(-1);

}
