/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Compare strings converted to lower case.
*/
int strcmpl(const char *s1, const char *s2)
{
    int iret;
    size_t l, len1, len2;
    char *l1, *l2;

    len1 = strlen(s1) + 1;
    len2 = strlen(s2) + 1;

    l1 = (char *) malloc(sizeof(char) * len1);
    l2 = (char *) malloc(sizeof(char) * len2);

    for (l = 0; l < len1; l++) *(l1 + l) = tolower(*(s1 + l));
    for (l = 0; l < len2; l++) *(l2 + l) = tolower(*(s2 + l));

    iret = strcmp(l1, l2);

    free(l1);
    free(l2);

    return(iret);
}

/*------------------------------------------------------------------------------
  Compare first n, at most, characters of strings converted to lower case.
*/
int strncmpl(const char *s1, const char *s2, size_t n)
{
    int iret;
    size_t l, len1, len2;
    char *l1, *l2;

    len1 = strlen(s1) + 1;
    len2 = strlen(s2) + 1;

    if (len1 > n) len1 = n;
    if (len2 > n) len2 = n;

    l1 = (char *) malloc(sizeof(char) * len1);
    l2 = (char *) malloc(sizeof(char) * len2);

    for (l = 0; l < len1; l++) *(l1 + l) = tolower(*(s1 + l));
    for (l = 0; l < len2; l++) *(l2 + l) = tolower(*(s2 + l));

    iret = strncmp(l1, l2, n);

    free(l1);
    free(l2);

    return(iret);
}
