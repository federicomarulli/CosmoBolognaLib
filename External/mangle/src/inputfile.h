/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#ifndef INPUTFILE_H
#define INPUTFILE_H

#include <stdio.h>

typedef struct {
    char *name;			/* filename */
    FILE *file;			/* file stream */
    char *line;			/* line buffer */
    size_t bufsize;		/* size of line buffer (will expand as necessary) */
    unsigned int line_number;	/* line number */
    unsigned int end;		/* maximum number of characters to read (0 = no limit) */
} inputfile;

int	rdline(inputfile *);

#endif	/* INPUTFILE_H */
