/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputfile.h"

/*------------------------------------------------------------------------------
  Read line into buffer, expanding memory as needed.

   Input: file = pointer to inputfile structure.
  Return value:  1 = ok;
		 0 = EOF;
		-1 = error.
*/
int rdline(inputfile *file)
{
#define BUFSIZE0	64
#define WHERE	fprintf(stderr, "rdline: at line %d of %s:", file->line_number, file->name)
    char *line_new;

    /* adjust bufsize to desired number of characters to read */
    if (file->end > 0 && file->bufsize != file->end + 1) {
	if (file->line) free(file->line);
	file->line = 0x0;
	file->bufsize = file->end + 1;
    }

    /* allocate memory for line buffer */
    if (!file->line) {
	if (file->bufsize <= 0) file->bufsize = BUFSIZE0;
	file->line = (char *) malloc(sizeof(char) * file->bufsize);
	if (!file->line) {
	    WHERE;
	    fprintf(stderr, " line too long (%zd characters)\n", file->bufsize);
	    return(-1);
	}
    }

    /* read line of data */
    if (!fgets(file->line, file->bufsize, file->file)) return(0);

    /* increment line number */
    file->line_number++;

    /* expand memory for line buffer if necessary */
    if (file->end == 0) {
	while (file->line[strlen(file->line) - 1] != '\n') {
	    file->bufsize *= 2;
	    line_new = (char *) malloc(sizeof(char) * file->bufsize);
	    if (!line_new) {
		WHERE;
		fprintf(stderr, " line too long (%zd characters)\n", file->bufsize);
		return(-1);
	    }
	    strncpy(line_new, file->line, file->bufsize/2 - 1);
	    free(file->line);
	    file->line = line_new;
	    if (!fgets(file->line + file->bufsize/2 - 1, file->bufsize/2 + 1, file->file)) {
		WHERE;
		fprintf(stderr, " missing newline at EOF\n");
		return(0);
	    }
	}
    }

    /* success */
    return(1);
}
