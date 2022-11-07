/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Read angle from word.

   Input: word = pointer to string.
	  unit = units of angle;
		 this is needed to detect 'h' unit [hms(RA) & dms(Dec)],
		 not to change the units of the angle, which is unchanged.
  Output: *next = pointer to character immediately after word;
		  *next may be the same pointer as word.
	  angle = angle.
  Return value: 1 = ok,
		otherwise error.
*/
int rdangle(char *word, char **next, char unit, long double *angle)
{
    const char *blank = " \t\n\r";
    const char *number = "+-.0123456789eE";

    int i, ird, sgn;
    unsigned int deg, min;
    long double sec;
    char *ch;

    /* find first non-blank character */
    ch = word;
    while (*ch && strchr(blank, *ch)) ch++;
    if (!*ch) return(0);

    /* units are hours, minutes, seconds (Right Ascension)
       and degrees, minutes, seconds (Declination) */
    if (unit == 'h') {
	sgn = 1;
	if (*ch == '-') {
	    sgn = -1;
	    ch++;
	}
	ird = sscanf(ch, "%d %d %Lf", &deg, &min, &sec);
	*angle = sgn * ((sec/60. + min)/60. + deg);
	for (i = 0; i < ird; i++) {
	    while (*ch && strchr(blank, *ch)) ch++;
	    while (*ch && strchr(number, *ch)) ch++;
	}
	if (ird == 3) {
	   ird = 1;
	} else {
	   ird = 0;
	}
    } else {
	ird = sscanf(ch, "%Lf", angle);
	while (*ch && strchr(number, *ch)) ch++;
    }

    /* point next at character after word */
    *next = ch;

    return(ird);
}
