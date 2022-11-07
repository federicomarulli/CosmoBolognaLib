/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  c interface to fortran subroutines in findtop.s.f
*/
void findtop(long double a[], int na, int iord[], int nb)
{
    int i;

    findtop_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void findbot(long double a[], int na, int iord[], int nb)
{
    int i;

    findbot_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void findtpa(long double a[], int na, int iord[], int nb)
{
    int i;

    findtpa_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void findbta(long double a[], int na, int iord[], int nb)
{
    int i;

    findbta_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void finitop(int a[], int na, int iord[], int nb)
{
    int i;

    finitop_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void finibot(int a[], int na, int iord[], int nb)
{
    int i;

    finibot_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void finitpa(int a[], int na, int iord[], int nb)
{
    int i;

    finitpa_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}

void finibta(int a[], int na, int iord[], int nb)
{
    int i;

    finibta_(a, &na, iord, &nb);
    /* convert from fortran to c convention */
    for (i = 0; i < nb; i++) iord[i]--;
}
