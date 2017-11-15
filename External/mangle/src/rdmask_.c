/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "manglefn.h"
#include "defaults.h"

/* global declaration of polygons here */
int npolys;
polygon *polys[NPOLYSMAX];

/*------------------------------------------------------------------------------
  Simplified fortran interface to rdmask routine.
      call rdmask()
*/
void rdmask_(void)
{
    char name[256];

    printf(" enter INPUT polygon file:\n");
    scanf("%256s", name);
    npolys = rdmask(name, &fmt, NPOLYSMAX, polys);
    if (npolys == -1) exit(1);
}
