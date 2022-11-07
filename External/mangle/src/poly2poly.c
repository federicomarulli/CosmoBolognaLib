/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqm:j:J:k:K:ns:e:v:p:i:o:";

/* allocate polygons as a global array */
polygon *polys_global[NPOLYSMAX];

/* local functions */
void	usage(void);
#ifdef	GCC
int	intersect_poly(int npoly1, polygon *[npoly1], int npoly2, polygon *[npoly2], long double);
#else
int	intersect_poly(int npoly1, polygon *[/*npoly1*/], int npoly2, polygon *[/*npoly2*/], long double);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, ipoly, nfiles, npoly, npolys, i;
    polygon **polys;
    polys=polys_global;

    /* default output format */
    fmt.out = keywords[POLYGON];

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

    msg("---------------- poly2poly ----------------\n");

    /* tolerance angle for multiple intersections */
    if (mtol != 0.) {
	scale (&mtol, munit, 's');
	munit = 's';
	msg("multiple intersections closer than %Lg%c will be treated as coincident\n", mtol, munit);
	scale (&mtol, munit, 'r');
	munit = 'r';
    }

    /* weight limits */
    if (is_weight_min && is_weight_max) {
	/* min <= max */
	if (weight_min <= weight_max) {
	    msg("will keep only polygons with weights inside [%Lg, %Lg]\n", weight_min, weight_max);
	/* min > max */
	} else {
	    msg("will keep only polygons with weights >= %Lg or <= %Lg\n", weight_min, weight_max);
	    msg("         (only polygons with weights outside (%Lg, %Lg))\n", weight_max, weight_min);
	}
    } else if (is_weight_min) {
	msg("will keep only polygons with weights >= %Lg\n", weight_min);
    } else if (is_weight_max) {
	msg("will keep only polygons with weights <= %Lg\n", weight_max);
    }
    /* area limits */
    if (is_area_min && is_area_max) {
	/* min <= max */
	if (area_min < area_max) {
	    msg("will keep only polygons with areas inside [%Lg, %Lg]\n", area_min, area_max);
	/* min > max */
	} else {
	    msg("will keep only polygons with areas >= %Lg or <= %Lg\n", area_min, area_max);
	    msg("         (only polygons with areas outside (%Lg, %Lg))\n", area_max, area_min);
	}
    } else if (is_area_min) {
	msg("will keep only polygons with areas >= %Lg\n", area_min);
    } else if (is_area_max) {
	msg("will keep only polygons with areas <= %Lg\n", area_max);
    }
    /* id limits */
    if (is_id_min && is_id_max) {
	/* min <= max */
	if (id_min < id_max) {
	    msg("will keep only polygons with ids inside [%d, %d]\n", id_min, id_max);
	/* min > max */
	} else {
	    msg("will keep only polygons with ids >= %d or <= %d\n", id_min, id_max);
	    msg("         (only polygons with ids outside (%d, %d))\n", id_max, id_min);
	}
    } else if (is_id_min) {
	msg("will keep only polygons with areas >= %d\n", id_min);
    } else if (is_id_max) {
	msg("will keep only polygons with areas <= %d\n", id_max);
    }
    /* pixel limits */
    if (is_pixel_min && is_pixel_max) {
	/* min <= max */
	if (pixel_min < pixel_max) {
	    msg("will keep only polygons with pixel numbers inside [%d, %d]\n", pixel_min, pixel_max);
	/* min > max */
	} else {
	    msg("will keep only polygons with pixel numbers >= %d or <= %d\n", pixel_min, pixel_max);
	    msg("         (only polygons with pixel numbers outside (%d, %d))\n", pixel_max, pixel_min);
	}
    } else if (is_pixel_min) {
	msg("will keep only polygons with pixel numbers >= %d\n", pixel_min);
    } else if (is_pixel_max) {
	msg("will keep only polygons with pixel numbers <= %d\n", pixel_max);
    }

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 1 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &polys[npoly]);
	if (npolys == -1) exit(1);
	/* intersect polygons of infile1 with those of subsequent infiles */
	if (ifile > optind && intersect) {
	    npoly = intersect_poly(npoly, polys, npolys, &polys[npoly], mtol);
	    if (npoly == -1) exit(1);
	/* increment number of polygons */
	} else {
	    npoly += npolys;
	}
    }
    if (nfiles >= 2 && !intersect) {
        msg("total of %d polygons read\n", npoly);
    }
    if (npoly == 0) {
	msg("STOP\n");
	exit(0);
    }

    /* apply new id numbers to output polygons */
    if (fmt.newid == 'n') {
	for (ipoly = 0; ipoly < npoly; ipoly++) {
	    polys[ipoly]->id = ipoly;
	}
    }

    if (fmt.newid == 'p') {
      for (ipoly = 0; ipoly < npoly; ipoly++) {
	polys[ipoly]->id = polys[ipoly]->pixel;
      }
    }

    /* write polygons */
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
    printf("poly2poly [-d] [-q] [-m<a>[u]] [-j[<min>][,<max>]] [-J[<min>][,<max>]] [-k[min][,<max>]] [-K[min][,<max>]] [-n] [-s<n>] [-e<n>] [-vo|-vn|-vp] [-p[+|-][<n>]] [-i<f>[<n>][u]] [-o<f>[u]] polygon_infile1 [polygon_infile2 ...] polygon_outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Intersect polygons of poly1 with any polygon(s) of poly2 having the
  same id number.

  This subroutine implements the -n option of poly2poly.
*/
int intersect_poly(int npoly1, polygon *poly1[/*npoly1*/], int npoly2, polygon *poly2[/*npoly2*/], long double mtol)
{
    int ier, inull, iprune, i, j, k, np;

    /* intersect each poly1 with any poly2 having same id number */
    for (i = 0; i < npoly1; i++) {
	for (j = 0; j < npoly2; j++) {
	    if (poly1[i]->id == poly2[j]->id && poly1[i]->pixel == poly2[j]->pixel ) {
		/* make sure poly1 contains enough space for intersection */
		np = poly1[i]->np + poly2[j]->np;
		ier = room_poly(&poly1[i], np, 0, 1);
		if (ier == -1) goto out_of_memory;

		/* intersection of poly1 and poly2 */
		poly_poly(poly1[i], poly2[j], poly1[i]);
	    }
	}
    }

    /* free poly2 polygons */
    for (j = 0; j < npoly2; j++) {
        free_poly(poly2[j]);
	poly2[j] = 0x0;
    }

    /* prune poly1 polygons */
    j = 0;
    inull = 0;
    for (i = 0; i < npoly1; i++) {
          iprune = prune_poly(poly1[i], mtol);
        if (iprune == -1) {
	    fprintf(stderr, "intersect_poly: failed to prune polygon %d; continuing ...\n", (fmt.newid == 'o')? poly1[i]->id : j);
	}
	if (iprune >= 2) {
	    free_poly(poly1[i]);
	    poly1[i] = 0x0;
	    inull++;
	} else {
	    j++;
	}
   }

   /*copy down non-null polygons*/
    k=0;
    for(i = 0; i < npoly1; i++){
      if(poly1[i]){
	poly1[k++]=poly1[i];
      }
    }
    /*after copying non-null polygons, k should be equal to j */
    if(k!=j){
      fprintf(stderr, "intersect_poly: should be left with %d non-null polygons, but actually have %d\n",j,k);
    }

    /*nullify the rest of the array, but don't free, since pointers have been copied above*/
    for(i=j; i < npoly1; i++){
      poly1[i]=0x0;
    }

    if (inull > 0) msg("%d intersected polygons have zero area, and are being discarded\n", inull);
    npoly1 = j;

    return(npoly1);

    /* ---------------- error returns ---------------- */
    out_of_memory:
    fprintf(stderr, "intersect_poly: failed to allocate memory for polygon of %d caps\n", np);
    return(-1);
}
