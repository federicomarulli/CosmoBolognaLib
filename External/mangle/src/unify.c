/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqm:s:e:v:p:Ui:o:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);
int	unify_poly(polygon **, polygon *);
#ifdef	GCC
int	unify(int *npoly, polygon *[*npoly]);
#else
int	unify(int *npoly, polygon *[/**npoly*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nadj, nfiles, npoly, npolys,i;
    polygon **polys;
    polys=polys_global;

    /* default output format */
    fmt.out = keywords[POLYGON];
    /* default is to renumber output polygons with new id numbers */
    fmt.newid = 'n';

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

    msg("---------------- unify ----------------\n");

    /* tolerance angle for multiple intersections */
    if (mtol != 0.) {
	scale(&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %Lg%c will be treated as coincident\n", mtol, munit);
	scale(&mtol, munit, 'r');
	munit = 'r';
    }

    if(unpixelize){
      msg("removing pixelization info by unifying across the whole mask\n");
    }
    else{
      msg("only unifying within each pixel - to unify across the whole mask, use -U\n");
    }

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

    if (snapped==0 || balkanized==0) {
      fprintf(stderr, "Error: input polygons must be snapped and balkanized before unification.\n");
      fprintf(stderr, "If your polygons are already snapped and balkanized, add the 'snapped' and\n'balkanized' keywords at the beginning of each of your input polygon files.\n");
      exit(1);
    }


    /* unify polygons */
    nadj = unify(&npoly, polys);
    if (nadj == -1) exit(1);

    if(unpixelize) pixelized=0;
    
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
    printf("unify [-d] [-q] [-m<a>[u]] [-s<n>] [-e<n>] [-vo|-vn|-vp] [-p[+|-][<n>]] [-U] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Unify two polygons poly1 and poly2, if possible.

   Input: *poly1, poly2 = pointers to polygon structures.
  Output: if poly1 and poly2 were unified, then *poly1 contains unified polygon;
	  if poly1 and poly2 were not unified, then *poly1 is unchanged.
  Return value: -1 = error occurred;
	        0 = poly1 and poly2 were not unified;
	        1 = poly1 and poly2 were unified.
*/
int unify_poly(polygon **poly1, polygon *poly2)
{
/* number of extra caps to allocate to polygon, to allow for expansion */
#define DNP		4
    static polygon *poly = 0x0;

    int bnd, bndin, bndout, bnd1, bnd2, i, ier, i1, i2, np, verb;
    int np1, np2;
    long double area, areain, tol;
    polygon *polyin, *polyout;

    bnd = 0;
    bnd1 = 0;
    bnd2 = 0;
    /* look for single boundary dividing poly1 and poly2 */
    for (i1 = 0; i1 < (*poly1)->np; i1++) {
	if (bnd >= 2) break;
	for (i2 = 0; i2 < poly2->np; i2++) {
	    if ((*poly1)->cm[i1] == - poly2->cm[i2]
		&& (*poly1)->rp[i1][0] == poly2->rp[i2][0]
		&& (*poly1)->rp[i1][1] == poly2->rp[i2][1]
		&& (*poly1)->rp[i1][2] == poly2->rp[i2][2]) {
		bnd++;
		if (bnd >= 2) break;
		bnd1 = i1;
		bnd2 = i2;
	    }
	}
    }

    /* poly1 and poly2 are not separated by a single boundary */
    if (bnd != 1) return(0);
    
    np1=(*poly1)->np;
    np2=poly2->np;
    np=np1+np2;

    /* make sure poly contains enough space for intersection */
    
    ier = room_poly(&poly, np, DNP, 0);
    if (ier == -1) goto out_of_memory;

    /* check whether areas of poly1 and poly2 equal subareas of unified poly */
    for (i = 0; i < 2; i++) {
	if (i == 0) {
	    polyin = poly2;
	    polyout = *poly1;
	    bndin = bnd2;
	    bndout = poly2->np + bnd1;
	} else {
	    polyin = *poly1;
	    polyout = poly2;
	    bndin = bnd1;
	    bndout = (*poly1)->np + bnd2;
	}
	/* area of polyin */
	tol = mtol;
	verb = 1;
	ier = garea(polyin, &tol, verb, &areain);
	if (ier) goto error;
	/* intersection of polyin and polyout (in that order!) */
	poly_poly(polyin, polyout, poly);
	/* suppress coincident boundaries */
	touch_poly(poly);
	/* suppress excluding boundary of unified poly */
	poly->cm[bndout] = 2.;
	/* subarea of unified poly */
	verb = 0;
	ier = garea(poly, &tol, verb, &area);
	if (ier == -1) goto error;
	if (ier || area != areain) break;
    }

    /* do not unify poly1 and poly2 */
    if (ier || area != areain) return(0);

    /* suppress dividing boundary */
    poly->cm[bndin] = 2.;

    /* prune unified polygon */
    if (prune_poly(poly, mtol) == -1) return(-1);

    /* make sure poly1 contains enough space */
    np = poly->np;
    ier = room_poly(poly1, np, DNP, 0);
    if (ier == -1) goto out_of_memory;

    /* copy unified polygon into poly1 */
    copy_poly(poly, *poly1);

    return(1);

    /* ---------------- error returns ---------------- */
    error:
    return(-1);

    out_of_memory:
    fprintf(stderr, "unify_poly: failed to allocate memory for polygon of %d caps\n", np + DNP);
    return(-1);
}

/*------------------------------------------------------------------------------
  Unify polygons.

   Input: poly = array of pointers to polygons.
	  npoly = pointer to number of polygons.
  Output: polys = array of pointers to polygons;
  Return value: number of polygons unified,
		or -1 if error occurred.
*/
int unify(int *npoly, polygon *poly[/**npoly*/])
{
#define WARNMAX		8
    int dnadj, i, j, unified, nadj, pass, warnmax;
    int *start;
    int *total;
    int begin, end, p, max_pixel, ier;

    /* allocate memory for pixel info arrays start and total */ 
    /* if unpixelizing, don't use pixelization info */
    max_pixel = (unpixelize)? 1: poly[*npoly-1]->pixel+1; 
    start = (int *) malloc(sizeof(int) * max_pixel);
    if (!start) {
      fprintf(stderr, "unify: failed to allocate memory for %d integers\n", max_pixel);
      return(-1);
    }
    total = (int *) malloc(sizeof(int) * max_pixel);
    if (!total) {
      fprintf(stderr, "unify: failed to allocate memory for %d integers\n", max_pixel);
      return(-1);
    }

    /* if we're unpixelizing, don't use the pixelization info */ 
    if(unpixelize){
      start[0]=0;
      total[0]=*npoly;
    }
    else{
      /* build lists of starting indices of each pixel and total number of polygons in each pixel*/
      ier=pixel_list(*npoly, poly, max_pixel, start, total);
      if (ier == -1) {
	fprintf(stderr, "unify: error building pixel index lists\n");
	return(-1);
      }
    }
    
    /*turn off warning messages if using more than one pixel*/
    warnmax= (max_pixel<=1) ? WARNMAX : 0; 

    nadj = 0;

    /* discard polygons that have zero weight */
    dnadj = 0;
    for (i = 0; i < *npoly; i++) {
	if (!poly[i]) continue;
	if (poly[i]->weight == 0.) {
	  if (WARNMAX > 0 && dnadj == 0)
	    msg("unify: the following polygons have zero weight and are being discarded:\n");
	  if (dnadj < WARNMAX) {
	    msg(" %d", poly[i]->id);
	  } else if (dnadj == WARNMAX) {
	    msg(" ... more\n");
	  }
	  dnadj++;
	  free_poly(poly[i]);
	  poly[i] = 0x0;	  
	}
    }
    if (WARNMAX > 0 && dnadj > 0 && dnadj <= WARNMAX) msg("\n");
    if (dnadj > 0) msg("unify: %d polygons with zero weight were discarded\n", dnadj);
    nadj += dnadj;

    /* discard polygons that are null */
    dnadj = 0;
    for (i = 0; i < *npoly; i++) {
	if (!poly[i]) continue;
	/* prune polygon (should really already have been done by balkanize) */
	if (prune_poly(poly[i], mtol) == -1) {
	    fprintf(stderr, "unify_poly: failed to prune polygon %d; continuing\n", poly[i]->id);
	    continue;
	}
	if (poly[i]->np > 0 && poly[i]->cm[0] == 0.) {
	  if (WARNMAX > 0 && dnadj == 0)
	    msg("unify: the following polygons have zero area and are being discarded:\n");
	  if (dnadj < WARNMAX) {
	    msg(" %d", poly[i]->id);
	  } else if (dnadj == WARNMAX) {
	    msg(" ... more\n");
	  }
	  dnadj++;
	  free_poly(poly[i]);
	  poly[i] = 0x0;	  
        }
    }
    if (WARNMAX > 0 && dnadj > 0 && dnadj <= WARNMAX) msg("\n");
    if (dnadj > 0) msg("unify: %d polygons with zero area were discarded\n", dnadj);
    nadj += dnadj;

    /*unify polygons within each pixel*/
    /* unify repeatedly, until no more unification occurs */
    for(p=0;p<max_pixel;p++){
      begin=start[p];
      end=start[p]+total[p];
      
      pass = 0;
      do {
	pass++;
	dnadj = 0;
	/* try unifying each polygon in turn ... */
	
	for (i = begin; i < end; i++) {
	  if (!poly[i]) continue;
	  
	  /* ... with another polygon */
	  for (j = i+1; j < end; j++) {
	    if (!poly[j]) continue;
	    /* only unify polygons with the same weight */
	    if (poly[i]->weight != poly[j]->weight) continue;
	    
	    /* if applying old ids, then only unify polygons with the same id */
	    if (fmt.newid == 'o' && poly[i]->id != poly[j]->id) continue;
	    
	    /* try unifying polygons */
	    unified = unify_poly(&poly[i], poly[j]);
	    if (unified == -1) {
	      if (warnmax/2 > 0 && dnadj > 0 && dnadj <= warnmax/2) msg("\n");
	      fprintf(stderr, "unify_poly: failed to unify polygons %d & %d; continuing\n", poly[i]->id, poly[j]->id);
	      continue;
	    }
	    /* polygons were unified */
	    if (unified) {
	      if(warnmax > 0){
		if (warnmax/2 > 0 && dnadj == 0)
		  msg("unify pass %d: the following polygons are being unified:\n", pass);
		if (dnadj < warnmax/2) {
		  msg(" (%d %d)", poly[i]->id, poly[j]->id);
		} else if (dnadj == warnmax/2) {
		  msg(" ... more\n");
		}
	      }
	      free_poly(poly[j]);
	      poly[j] = 0x0;
	      dnadj++;
	    }
	  }
	}
	if (warnmax/2 > 0 && dnadj > 0 && dnadj <= warnmax/2) msg("\n");
	if(warnmax) msg("unify pass %d: %d polygons unified\n", pass, dnadj);
	nadj += dnadj;
	
      } while (dnadj);
    }
    
    free(start);
    free(total);

    /* copy down polygons */
    j = 0;
    for (i = 0; i < *npoly; i++) {
	if (poly[i]) {
	    poly[j] = poly[i];
	    j++;
	}
    }
    *npoly = j;

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

    /* if unpixelizing, set all pixel numbers to zero */
    if(unpixelize) {
      for (i = 0; i < *npoly; i++) {
	poly[i]->pixel = 0;
      }
    }

    /* advise */
    msg("unify: %d polygons discarded or unified\n", nadj);
    
    return(nadj);
}
