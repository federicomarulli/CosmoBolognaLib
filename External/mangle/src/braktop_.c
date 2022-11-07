/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  c interface to fortran subroutines in braktop.s.f
*/
void braktop(long double aa, int *ia, long double a[], int n, int l)
{
    braktop_(&aa, ia, a, &n, &l);
}

void brakbot(long double aa, int *ia, long double a[], int n, int l)
{
    brakbot_(&aa, ia, a, &n, &l);
}

void braktpa(long double aa, int *ia, long double a[], int n, int l)
{
    braktpa_(&aa, ia, a, &n, &l);
}

void brakbta(long double aa, int *ia, long double a[], int n, int l)
{
    brakbta_(&aa, ia, a, &n, &l);
}
