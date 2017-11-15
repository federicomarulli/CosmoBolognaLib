/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#ifndef VERTICES_H
#define VERTICES_H

typedef struct {		/* azel structure */
    long double az;			/* azimuth */
    long double el;			/* elevation */
} azel;

typedef struct {		/* vertices structure */
    int nv;			/* number of azel vertices */
    int nvmax;			/* dimension of allocated v array */
    azel *v;			/* array v[nv] of azel structures */
} vertices;

#endif	/* VERTICES_H */
