/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/#include <limits.h>
#ifndef LOGICAL_H
#define LOGICAL_H

/* define a c type logical equivalent to fortran's logical, which (on most
systems) is a 32 bit integer
 */
#if LONG_MAX==2147483647
typedef long int logical;	/* type of fortran logical, according to f2c */
#elif  INT_MAX==2147483647
typedef int logical;            /* if long is too big, make logical an int */
#else
#error Could not define a 32-bit integer type for "logical"
#endif

#endif	/* LOGICAL_H */
