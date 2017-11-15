/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include "manglefn.h"

#ifdef GCC
# include <stdarg.h>
#endif
#ifdef LINUX
# include <stdarg.h>
#endif
#ifdef SUN
# include <sys/varargs.h>
#endif

extern int verbose;

/*------------------------------------------------------------------------------
  Print messages.
*/
void msg(char *fmt, ...)
{
    va_list args;

    if (verbose) {
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	fflush(stdout);
    }
}
