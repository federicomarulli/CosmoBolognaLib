/*------------------------------------------------------------------------------
© M E C Swanson 2008
------------------------------------------------------------------------------*/
/* determines which sorting function to use depending on OS */
#ifndef POLYSORT_H
#define POLYSORT_H

#ifdef LINUX
#define mysort qsort
#endif
#ifdef SUN
#define mysort qsort
#endif
#ifdef MACOSX
#define mysort mergesort
#endif

#ifndef mysort
#define mysort qsort
#endif

#endif	/* POLYSORT_H */
