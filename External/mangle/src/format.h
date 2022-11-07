/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#ifndef FORMAT_H
#define FORMAT_H

#include <stdio.h>

/*
  Structure defining format of data.
  If your format needs more descriptives, shove 'em in here.
  DON'T FORGET TO SET THE DEFAULT VALUES IN defaults.h,
  AND TO UPDATE copy_format().
*/
typedef struct {
    char *in;		/* keyword defining the input data format */
    char *out;		/* keyword defining the output data format */
    size_t skip;	/* skip first skip characters of each line */
    size_t end;		/* read only up to end'th character of each line */
    char single;	/* keyword defines precisely one polygon */
    int n;		/* the number of thingys defined by keyword */
    int nn;		/* the number of thingys per thingy */
    int innve;		/* the input number of points per edge */
    int outper;		/* controls interpretation of outnve */
    int outnve;		/* the output number of points per edge */
    int id;		/* id number of current polygon */
    char newid;		/* whether to use old or new id number */
    int pixel;          /* pixel that current polygon is in */ 
    long double weight;	/* weight of current polygon */
    char inunitp;	/* angular units of input polygon data */
    char outunitp;	/* angular units of output polygon data */
    int inframe;	/* angular frame of input az, el data */
    int outframe;	/* angular frame of output az, el data */
    char inunit;	/* angular units of input az, el data */
    char outunit;	/* angular units of output az, el data */
    int outprecision;	/* digits after decimal point in output angles */
    char outphase;	/* '-' or '+' to make output azimuth in interval (-pi, pi] or [0, 2 pi) */
    long double azn;		/* azimuth of new pole wrt original frame */
    long double eln;		/* elevation of new pole wrt original frame
			 = elevation of original pole wrt new frame */
    long double azp;		/* azimuth of original pole wrt new frame */
    char trunit;	/* angular units of transformation angles */
    int nweights;       /* the total number of weights/polygons, for use with healpix_weight input files and rasterize */
} format;

#endif	/* FORMAT_H */
