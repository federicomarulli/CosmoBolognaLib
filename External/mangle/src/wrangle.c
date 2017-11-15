/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Write angle into string in appropriate format.
  In particular, this routine supports the hms dms format for RA & Dec.

  All angles written by mangle go through this routine.

   Input: angle = angle.
          unit = format in which to write angle;
		 this is only used to decide the format, not to rescale the
		 angle, which remains unchanged;
                 the 'h' unit signifies hms(RA) and dms(Dec);
		 otherwise angle is written to str as a long double.
	  precision = number of digits after decimal point in output angles.
  Output: str = pointer to string containing the angle.
	  str_len = length of string.
*/
void wrangle(long double angle, char unit, int precision, size_t str_len, char str[/*str_len*/])
{
/* default number of significant digits */
#define DIGITS		9
    char sign;
    int hour, min;
    int width;
    long double a, sec;

    if (unit == 'h') {
	sign = (angle < 0.)? '-' : ' ';
	a = fabsl(angle);
	hour = floorl(a);
	a = (a - hour) * 60.;
	min = floorl(a);
	a = (a - min) * 60.;
	sec = a;
	if (precision < 0) precision = DIGITS - 7;
	width = precision + 3;
	if (precision == 0) width--;
	snprintf(str, str_len, "%c%02d %02d %0*.*Lf",
	    sign, hour, min, width, precision, sec);
    } else {
	switch (unit) {
	case 'r':	/* radians */
	default:
	    if (precision < 0) precision = DIGITS - 1;
	    width = precision + 3;
	    break;
	case 'd':	/* degrees */
	case '°':
	    if (precision < 0) precision = DIGITS - 3;
	    width = precision + 5;
	    break;
	case 'm':	/* arcminutes */
	case '\'':
	case '´':
	    if (precision < 0) precision = DIGITS - 5;
	    width = precision + 7;
	    break;
	case 's':	/* arcseconds */
	case '"':
	case '¨':
	    if (precision < 0) precision = DIGITS - 7;
	    width = precision + 9;
	    break;
	}
	if (precision == 0) width--;
	snprintf(str, str_len, "%*.*Lf",
	    width, precision, angle);
    }
}
