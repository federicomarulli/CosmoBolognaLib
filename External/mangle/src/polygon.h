/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#ifndef POLYGON_H
#define POLYGON_H

typedef long double vec[3];

typedef struct {		/* polygon structure */
  int np;			/* number of caps of polygon */
  int npmax;			/* dimension of allocated rp and cm arrays */
  vec *rp;			/* pointer to array rp[np][3] of axis coords */
  long double *cm;			/* pointer to array cm[np] of 1 - cosl(theta) */
  int id;			/* id number of polygon */
  int pixel;                    /* pixel that polygon is in */
  long double weight;		/* weight of polygon */
} polygon;

#endif	/* POLYGON_H */
