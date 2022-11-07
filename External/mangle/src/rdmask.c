/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------
This subroutine recognizes the following formats:

<keyword> <id> ( <int> caps, <long double> weight ...

polygon
-------
polygon <id> ( <n> caps, <long double> weight):
<rp_00> <rp_10> <rp_20> <cm_0>
<rp_01> <rp_11> <rp_21> <cm_1>
...
<rp_0n> <rp_1n> <rp_2n> <cm_n>
Each rp_i is a unit vector defining the north polar axis of the cap,
while cm = 1 - cosl(theta) defines the cap latitude, theta being the polar angle.
Positive cm designates the region north of the latitude,
while negative cm designates the region south of the latitude.

circle
------
<az_0> <el_0> <th_0> ... <az_n> <el_n> <th_n>
Each circle is defined by the azimuth and elevation of its axis,
and by the angular radius, the polar angle theta, of the circle about that axis.

vertices
--------
<az_0> <el_0> <az_1> <el_1> ... <az_n> <el_n>
Each vertex is defined by an azimuth and an elevation,
and they are joined by great circles.
The vertices wind right-handedly about the enclosed region, as in
   3  0
   2  1
 <- az <-
with elevation increasing upward, and azimuth increasing to the LEFT
(as it does when you look at the sky with the north pole at zenith).
A rectangle would join the vertices 0 to 3 above as
   0 azmin 1 elmin 2 azmax 3 elmax 0
There MUST be at least 3 vertices: a polygon bounded by 1 vertex
is ill-defined, and a polygon bounded by 2 vertices is null.
The interior angle at each vertex MUST be <= 180 deg (convex).

edges
-----
<az_0> <el_0> <az_1> <el_1> ... <az_n> <el_n>
This is like vertices, but with the addition of one or more points on each edge.
The vertices and edge-points wind right-handedly about the polygon, as in
      0
    5   1
  4   3   2
  <- az <-
with elevation increasing upward, and azimuth increasing to the LEFT
(as it does when you look at the sky with the north pole at zenith).

rectangle
---------
<azmin_0> <azmax_0> <elmin_0> <elmax_0> ... <azmin_n> <azmax_n> <elmin_n> <elmax_n>
Each rectangle is bounded by a minimum and maximum azimuth and elevation.
The azimuthal extent of the rectangle MUST NOT exceed 180 deg.

healpix_weight
--------------
<weight>
Simply a list of the weight of each HEALPix pixel at some resolution (i.e., nside); we read in
these weights and tack them onto the correct HEALPix polygon.
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputfile.h"
#include "manglefn.h"

#define WHERE		fprintf(stderr, "rdmask: at line %d of %s:", file.line_number, file.name)

inputfile file = {
    '\0',	/* input name */
    0x0,	/* input file stream */
    '\0',	/* line buffer */
    512,	/* size of line buffer (will expand as necessary) */
    0,		/* line number */
    0		/* maximum number of characters to read (0 = no limit) */
};

static int blankline = 0;	/* blank lines encountered where numbers expected */

/* dictionary of keywords */
extern char *keywords[];

/* min, max weights to keep (only used here in healpix_weight; other formats
   discard these unwanted polygons in wrmask.c) */
extern int is_weight_min, is_weight_max;

/* min, max areas to keep (only used here in healpix_weight; other formats
   discard these unwanted polygons in wrmask.c) */
extern int is_area_min, is_area_max;

/* global pixelization variables */
extern int res_max;
extern char scheme;
extern int pixelized;
extern int snapped;
extern int balkanized;

/* counter for input files*/
extern int infiles;

/* local functions */
char	*get_keyword(char *, char **, format *);
int	new_fmt(char *, char **, format *);
char	*get_word(char *, const char *, int, size_t *);
int	get_n(format *, char *);
int	get_nang(format *fmt, char *, int);
polygon	*get_poly(format *);
polygon	*rd_poly(format *);
polygon *rd_hpix(format *);
polygon	*rd_circ(format *);
polygon	*rd_edge(format *);
polygon	*rd_rect(format *);

/*------------------------------------------------------------------------------
  Read mask of polygons from file.

  The format is determined by a keyword and associated parameters.
  The format may or may not be initialized by the calling program.
  Reading continues in the current format until a line is read
  whose first word is a keyword.
  If a keyword is encountered, the format is changed accordingly.

   Input: name = name of file to read from;
		     "" or "-" means read from standard input.
	  fmt = pointer to format structure.
	  npolys = maximum number of polygons to read.
  Output: polys = polygons read.
  Return value: number of polygons read,
		or -1 if error occurred.
*/
int rdmask(char *name, format *fmt, int npolys, polygon *polys[/*npolys*/])
{
    char *input = "input";
    char *line_rest, *word;
    int ird;
    int npoly = 0;
    format in_fmt;
    polygon *poly = 0x0;

    /* store original contents of fmt */
    copy_format(fmt, &in_fmt);

    /* open name for reading */
    if (!name || strcmp(name, "-") == 0) {
	file.file = stdin;
	file.name = input;
    } else {
	file.file = fopen(name, "r");
	if (!file.file) {
	    fprintf(stderr, "rdmask: cannot open %s for reading\n", name);
	    return(-1);
	}
	file.name = name;
    }
    file.line_number = 0;
    file.end = fmt->end;

    /* read data until hit EOF */
    while (1) {
	/* read line of data */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) goto error;
	/* EOF */
	if (ird == 0) break;

	/* look for keyword as first word in line */
	word = get_keyword(file.line, &line_rest, fmt);
	/* initialize to new format */
	if (word) ird = new_fmt(word, &line_rest, fmt);
	if (ird == -1) goto error;
	
	/* read polygon */
	if ((word && fmt->single) || (!word && !fmt->single)) {
	  
	  poly = get_poly(fmt);
  
	    if (poly) {
		if (npoly >= npolys) {
		    fprintf(stderr, "rdmask: number of polygons exceeds maximum %d\n", NPOLYSMAX);
		    fprintf(stderr, " if you need more space, enlarge NPOLYSMAX in defines.h, and recompile\n");
		    goto error;
		}
		polys[npoly] = poly;
		npoly++;
	    }
	}
    }
    infiles++;

    /* for HEALPix weights, only 1 input file allowed (if you have multiple files,
       just combine them)*/
    if(fmt->in == keywords[HEALPIX_WEIGHT] && infiles >= 2) {
       fprintf(stderr, "error: only 1 input file allowed for HEALPix_weight format!\n");
       fprintf(stderr, "combine your HEALPix_weight files into one file before running any commands\n");
       return(-1);
    }

    /* check whether too few HEALPix weights were read */
    if(fmt->in == keywords[HEALPIX_WEIGHT] && fmt->id < fmt->nweights) {
        fprintf(stderr, "error: the number of HEALPix weights in the input file is less than %d (or the file does not end in a newline)\n", fmt->nweights);
        goto error;
    }

    /* if there are pixelized files, make sure all files are pixelized */ 
    if(pixelized>0 && infiles!=pixelized){
       fprintf(stderr, "error: some input files are pixelized and some are not.\n");
       fprintf(stderr, "all input files must have consistent pixelization.\n");
      goto error;
    }
    /* if there are snapped files, make sure all files are snapped */ 
    if(snapped>0 && infiles!=snapped){
       msg("WARNING: some input files are snapped and some are not.\n");
       snapped=0;
    }
    /* if there are balkanized files, make sure all files are balkanized */ 
    if(balkanized>0 && infiles!=balkanized){
      msg("WARNING: some input files are balkanized and some are not.\n");
      balkanized=0;
    }

    /* advise */
    msg("%d polygons read from %s\n", npoly, file.name);

    /* close file */
    if (file.file != stdin) fclose(file.file);

    /* warn about all blank file */
    if (npoly == 0 && blankline > 0) {
	msg("Is format correct?  Maybe check -s<n> -e<n> -i<f>[<n>][u] command line options?\n");
    }

    /* restore format as it was on input */
    copy_format(&in_fmt, fmt);

    return(npoly);

    /* ---------------- error returns ---------------- */
    error:

    /* restore format as it was on input */
    copy_format(&in_fmt, fmt);

    return(-1);
}

/*------------------------------------------------------------------------------
  Determine whether first word in string is a keyword.

   Input: str = string.
  Output: *str_rest = string following keyword.
  Return value: pointer to string containing keyword,
		or null if no match.
*/
char *get_keyword(char *str, char **str_rest, format *fmt)
{
    const char *blank = " \t\n\r";
    char c;
    int id, ird, n;

    char *word;
    int ikey;
    size_t word_len;

    /* first word in str */
    word = get_word(str, blank, 0, &word_len);
    if (!word) return(0x0);

    /* compare first word in str against keywords */
    for (ikey = 0; keywords[ikey]; ikey++) {
	if (strncmp(word, keywords[ikey], strlen(keywords[ikey])) == 0) {
	   *str_rest = word + word_len;
	   return(keywords[ikey]);
	}
    }

    /* format not yet specified, or spolygon */
    if (!fmt->in || strcmp(fmt->in, "spolygon") == 0) {
	/* check for line starting with 2 integers, indicating spolygon format */
	ird = sscanf(str, "%d %d%[ \t\n]", &id, &n, &c);
	if (ird == 3 && n >= 0) {
	    *str_rest = word;
	    return(keywords[SPOLYGON]);
	}
    }

    /* check for possibly truncated keyword */
    for (ikey = 0; keywords[ikey]; ikey++) {
	if (strlen(word) >= 4 && strncmp(word, keywords[ikey], 4) == 0) {
	    msg("at line %d of %s: is  %.*s  supposed to be the keyword  %s ?\n", file.line_number, file.name, word_len, word, keywords[ikey]);
	}
    }

    /* no match */
    return(0x0);
}

/*------------------------------------------------------------------------------
  Change format according to contents of line_rest.

   Input: keyword = keyword of new format.
	  *line_rest = string containing possible arguments of format.
	  fmt = pointer to format structure.
  Output: contents of format structure fmt revised.
  Return value: 0 = ok;
		-1 = error.
*/
int new_fmt(char *keyword, char **line_rest, format *fmt)
{
    /* const char *Region_fmt = "%d ( %d caps, %d holes):"; */
    const char *Region_fmt = "%d%*[^0-9]%d%*[^0-9]%d";
    /* const char *generic_fmt = "%d ( %d caps, %Lf weight ...:"; */
    //const char *generic_fmt = "%d%*[^0-9]%d%*[^0-9-.]%Lf";
    //const char *edges_fmt = "%d%*[^0-9]%d%*[^0-9]%d%*[^0-9-.]%Lf";
    const char *healpix_weight_fmt = "%d";
    const char *skip_fmt = "%d";
    const char *end_fmt = "%d";
    const char *unit_fmt = " %c";
    const char *pix_fmt = "%d%c";
    char *word;
    char *blank = " \t\n\r";
    int ird, iscan, nholes, i, flag;
    size_t word_len;
    long double temp_pixel;
    int res_max_temp;
    char scheme_temp;

    /* defaults */
    fmt->id = 0;
    fmt->pixel = 0;
    fmt->weight = 1.;
    
    iscan=0;
    
    /* Regions */
    if (strcmp(keyword, "Region") == 0) {
	fmt->in = keyword;
	ird = 0;
	word = *line_rest;
	do {
	    word = get_word(word, blank, 0, &word_len);
	    if (!word) break;
	    switch (ird) {
	    case 0:	iscan = sscanf(word, "%d", &fmt->id);	break;
	    case 1:	iscan = sscanf(word, "%d", &fmt->n);	break;
	    case 2:	iscan = sscanf(word, "%d", &nholes);	break;
	    }
	    if (iscan == 1) ird++;
	    word += word_len;
	} while (word && ird < 3);
	if (ird != 3) {
	    WHERE;
	    fprintf(stderr, " expecting line of format:\n");
	    fprintf(stderr, " %s %s\n", keyword, Region_fmt);
            return(-1);
	}
	/* format defines only one polygon */
	fmt->single = 1;
    
    /* healpix_weight */
    } else if (strcmp(keyword, "healpix_weight") == 0) {
        fmt->in = keyword;
	if(is_weight_min || is_weight_max) {
	  fprintf(stderr, "error: -j option invalid for HEALPix_weight format!\n");
	  fprintf(stderr, "if you want to discard certain weights, just set them to 0 in the input file\n");
	  return(-1);
	}
	if(is_area_min || is_area_max) {
	  fprintf(stderr, "error: -k option invalid for HEALPix_weight format!\n");
	  fprintf(stderr, "HEALPix pixels are equal-area, so area limits serve no purpose\n");
	  return(-1);
	}
	fmt->n = 1;
	fmt->single = 0;
	ird = 0;
	word = *line_rest;
	do {
	   word = get_word(word, blank, 0, &word_len);
	   if (!word) break;
	   switch (ird) {
	   case 0: iscan = sscanf(word, "%d", &fmt->nweights);           break;
	   }
	   if (iscan == 1) ird++;
	   word += word_len;
	} while (word && ird < 1);
	if (ird != 1) {
	   WHERE;
	   fprintf(stderr, " expecting line of format:\n");
	   fprintf(stderr, "%s %s\n", keyword, healpix_weight_fmt);
	   return(-1);
	}
	/* catch wrongly-formatted healpix_weight input file */
	if(fmt->nweights == 0) {
	  fprintf(stderr, "error: input file contains 0 weights???\n");
	  fprintf(stderr, "need line of format:\n");
	  fprintf(stderr, "%s %s\n", keyword, healpix_weight_fmt);
	  return(-1);
	}
	/* check whether nweights is a valid number based on the HEALPix pixelization scheme */
	for (i=1; i<=8192; i=2*i) {
       	    if (12*powl(i,2) == fmt->nweights) {
	       flag = 1;
	       break;
	    }
	    else {
	       flag = 0;
	       continue;
	    }
	}
	if (fmt->nweights == 1) flag = 1;
	if (flag == 0) {
  	   fprintf(stderr, "error: input file does not contain a valid number of weights according\n");
	   fprintf(stderr, "to the HEALPix pixelization scheme (see http://healpix.jpl.nasa.gov)\n");
	   return(-1);
	   }

    /* polygons, spolygons, circles, rectangles, vertices */
    } else if (strcmp(keyword, "polygon") == 0
	|| strcmp(keyword, "spolygon") == 0
	|| strcmp(keyword, "circle") == 0
	|| strcmp(keyword, "edges") == 0
	|| strcmp(keyword, "rectangle") == 0
	|| strcmp(keyword, "vertices") == 0) {
	fmt->in = keyword;
	if (strcmp(keyword, "edges") == 0) {
	    ird = 0;
	    word = *line_rest;
	    do {
		word = get_word(word, blank, 0, &word_len);
		if (!word) break;
		switch (ird) {
		case 0:	iscan = sscanf(word, "%d", &fmt->id);		break;
		case 1:	iscan = sscanf(word, "%d", &fmt->innve);	break;
		case 2:	iscan = sscanf(word, "%d", &fmt->n);		break;
		case 3:	iscan = sscanf(word, "%Lf", &fmt->weight);	break;
		}
		if (iscan == 1) ird++;
		word += word_len;
	    } while (word && ird < 4);
	    /* default number of edges per line is variable (fmt->n = 0) */
	    if (ird < 3) fmt->n = 0;
	    /* default number of points/edge is 2 */
	    if (ird < 2) fmt->innve = 2;
	} else {
	    ird = 0;
	    word = *line_rest;
	    do {
		word = get_word(word, blank, 0, &word_len);
		if (!word) break;
		switch (ird) {
		case 0:	iscan = sscanf(word, "%d", &fmt->id);          break;
		case 1:	iscan = sscanf(word, "%d", &fmt->n);           break;
		case 2:	iscan = sscanf(word, "%Lf", &fmt->weight);     break;
		case 3:	
		  /* checks to see if 3rd number is an integer - if it is, assume its a pixel number */
		  iscan = sscanf(word, "%Lf", &temp_pixel);
		  fmt->pixel=(floorl(temp_pixel)-temp_pixel==0) ? (int)temp_pixel :0;
		  break;
		}
		if (iscan == 1) ird++;
		word += word_len;
	    } while (word && ird < 4);
	    if (ird < 2) {
		/* polygon format requires at least 2 integer arguments */
		if (strcmp(keyword, "polygon") == 0 || strcmp(keyword, "spolygon") == 0) {
		    WHERE;
		    fprintf(stderr, " expecting two integers <id> <number_of_caps> following keyword %s\n", keyword);
		    return(-1);
		/* default number of objects per line is variable (fmt->n = 0) */
		} else {
		    fmt->n = 0;
		}
	    }
	}
	/* format defines only one polygon */
	if (strcmp(keyword, "polygon") == 0 || strcmp(keyword, "spolygon") == 0) {
	    fmt->single = 1;
	/* if number of objects per line is explicitly set to 0,
	   interpret as allsky, and allow only one allsky polygon */
	} else if (ird >= ((strcmp(keyword, "edges") == 0)? 3 : 2) && fmt->n == 0) {
	    fmt->single = 1;
	/* format may apply to many polygons */
	} else {
	    fmt->single = 0;
	    /* for rectangles, require 1 rectangle per line */
	    if (strcmp(keyword, "rectangle") == 0) fmt->n = 1;
	}
	/* number of angles to read per object */
	if (strcmp(keyword, "circle") == 0) {
	    fmt->nn = 3;
	} else if (strcmp(keyword, "edges") == 0) {
	    fmt->nn = 2;
	} else if (strcmp(keyword, "rectangle") == 0) {
	    fmt->nn = 4;
	} else if (strcmp(keyword, "vertices") == 0) {
	    fmt->innve = 1;
	    fmt->nn = 2;
	}

    /* skip */
    } else if (strcmp(keyword, "skip") == 0) {
	ird = sscanf(*line_rest, skip_fmt, &fmt->skip);
	if (ird != 1) {
	    WHERE;
	    fprintf(stderr, " expecting integer following keyword skip\n");
	    return(-1);
	}

    /* end */
    } else if (strcmp(keyword, "end") == 0) {
	ird = sscanf(*line_rest, end_fmt, &fmt->end);
	if (ird != 1) {
	    WHERE;
	    fprintf(stderr, " expecting integer following keyword end\n");
	    return(-1);
	}
	file.end = fmt->end;

    /* angular units */
    } else if (strcmp(keyword, "unit") == 0) {
	ird = sscanf(*line_rest, unit_fmt, &fmt->inunitp);
	if (ird != 1) {
	    WHERE;
		fprintf(stderr, " expecting character (one of %s) following keyword unit\n", UNITS);
	    return(-1);
	} else if (!strchr(UNITS, fmt->inunitp)) {
	    WHERE;
	    fprintf(stderr, " unit %c must be one of %s\n", fmt->inunitp, UNITS);
	    return(-1);
	}

    } else if (strcmp(keyword, "graphics") == 0) {
	WHERE;
	fprintf(stderr, " %s format can only be written, not read\n", keyword);
	return(-1);

    /*
    } else if (strcmp(keyword, "area") == 0
	|| (strcmp(keyword, "id") == 0)
	|| (strcmp(keyword, "midpoint") == 0)
	|| (strcmp(keyword, "weight") == 0)) {
	WHERE;
	fprintf(stderr, " %s format can only be written, not read\n", keyword);
	return(-1);
    */

	/*parse pixelization information (from line like "pixelization 7s" for example) */
    } else if (strcmp(keyword, "pixelization") == 0) {

	ird = sscanf(*line_rest, pix_fmt, &res_max_temp, &scheme_temp);
	if (ird != 2) {
	    WHERE;
		fprintf(stderr, " expecting integer and character (one of %s) following keyword pixelization\n", SCHEMES);
	    return(-1);
	} else if (!strchr(SCHEMES, scheme_temp)) {
	    WHERE;
	    fprintf(stderr, " scheme %c must be one of %s\n", scheme_temp, SCHEMES);
	    return(-1);
	}
	
	/* if previously read file has already provided pixelization info, check for compatibility */
	if(pixelized>0){
	  if(scheme_temp!=scheme || res_max!=res_max_temp){
	    fprintf(stderr, " Pixelization %d%c incompatible with pixelization %d%c of previous file.  Pixelize files to be combined later using the same resolution and scheme.\n", res_max_temp, scheme_temp,res_max,scheme);
	    return(-1);
	  }
	  if(res_max==-1 || res_max_temp==-1){
	    fprintf(stderr, " Files pixelized adaptively cannot be combined later.  To use the adaptive pixelization feature, combine files before pixelization.\n");
	    return(-1);
	  }
	}
	/* if this is the first set of pixelization info encountered, set scheme and res_max to provided values */
	else{
	  if(res_max_temp==-1){
	    msg("Input polygons are pixelized adaptively using scheme %c\n",scheme_temp);
	  }
	  else{
	    msg("Input polygons are pixelized to resolution %d using scheme %c\n",res_max_temp,scheme_temp); 
	  }
	  scheme=scheme_temp;
	  res_max=res_max_temp;
	}
	pixelized++;
    } else if (strcmp(keyword, "snapped") == 0) {
      snapped++;
    } else if (strcmp(keyword, "balkanized") == 0) {
      balkanized++;
    }



    return(0);
}

/*------------------------------------------------------------------------------
   Extract word from string.

   Input: str = string from which word is to be extracted.
	  s = string of characters allowed in word.
	  in = 1 to treat s as allowed characters,
	       0 to treat s as NOT-allowed characters, i.e. as blanks.
  Output: *word_len = length of first allowed sequence in line.
  Return value: pointer to first allowed character of line,
		or null if there is no allowed sequence.
*/
char *get_word(char *str, const char *s, int in, size_t *word_len)
{
    char *first, *last;

    /* check for null string */
    if (!str) return(0x0);

    /* first allowed character */
    if (in) {
	for (first = str; *first && !strchr(s, *first); first++);
    } else {
	for (first = str; *first && strchr(s, *first); first++);
    }

    /* nothing in str */
    if (!*first) return(0x0);

    /* last allowed character */
    if (in) {
	for (last = first + 1; *last && strchr(s, *last); last++);
    } else {
	for (last = first + 1; *last && !strchr(s, *last); last++);
    }

    /* length of word */
    *word_len = last - first;

    /* pointer to word in str */
    return(first);
}

/*------------------------------------------------------------------------------
  Read integer number of thingys from 1st integer in line.

   Input: fmt = pointer to format structure.
	  lin = string containing the line.
  Return value: integer number of thingys,
		or 0 if error.
*/
int get_n(format *fmt, char *lin)
{
#define WARNMAX		3
    char *blank = " \t\n\r";
    char *plural = "s";
    static int warn = 0;
    char *ch;
    int n;

    /* find next non-blank character */
    for (ch = lin; *ch && strchr(blank, *ch); ch++);
    if (!*ch) return(0);
    /* check word is an integer, followed by at least one blank */
    while (*ch && strchr("0123456789", *ch)) ch++;
    if (!*ch) {
	WHERE;
	fprintf(stderr, " expecting more numbers after first integer on line\n");
	warn++;
	if (warn >= WARNMAX) {
	    fprintf(stderr, " We have a problem, Heuston.\n");
	    exit(1);
	}
	return(0);
    }
    if (!strchr(blank, *ch)) {
	WHERE;
	if (fmt->in[strlen(fmt->in) - 1] == 's') {
	    strcpy(plural, "");
	} else {
	    strcpy(plural, "s");
	}
	fprintf(stderr, " for variable length format, 1st number on line should be an integer;\n");
	fprintf(stderr, " to set fixed length format with <n> %s%s, use -%c<n>\n",
	    fmt->in, plural, fmt->in[0]);
	warn++;
	if (warn >= WARNMAX) {
	    fprintf(stderr, " We have a problem, Heuston.\n");
	    exit(1);
	}
	return(0);
    }

    /* read integer */
    sscanf(lin, "%d", &n);

    return(n);
}

/*------------------------------------------------------------------------------
  Read number of angles in line.

   Input: fmt = pointer to format structure.
	  lin = string containing the line.
	  nh = only first nh of fmt->nn angles may be in hms(RA) dms(Dec) format.
  Return value: integer number of angles.
*/
int get_nang(format *fmt, char *lin, int nh)
{
    char unit;
    char *next, *word;
    int i, iang, ird;
    long double angle;

    word = lin;
    iang = 0;
    ird = 0;
    /* read each angle till get to end of line */
    do {
	for (i = 0; i < fmt->nn; i++) {
	    unit = fmt->inunitp;
	    /* assume only first two angles per set can be in hms dms format */
	    if (i >= nh && fmt->inunitp == 'h') unit = 'd';
	    ird = rdangle(word, &next, unit, &angle);
	    if (ird != 1) {
		/* skip over `comment' line */
		if (iang == 0 || (iang == 1 && fmt->inunitp != 'h')) {
		    blankline++;
		    return (0);
		/* got to end of line with right number of angles */
		} else if (i == 0) {
		    break;
		/* got to end of line but the number of angle is wrong */
		} else {
		    WHERE;
		    fprintf(stderr, " expecting number of angles divisible by %d, found %d angles\n", fmt->nn, iang);
		    exit(1);
		}
	    }
	    /* increment angle count */
	    iang++;
	    /* point word to next word */
	    word = next;
	}
    } while (ird == 1);

    return(iang);
}

/*------------------------------------------------------------------------------
  Get polygon from data in specified format.

   Input: fmt = pointer to format structure.
  Return value: pointer to polyon.
*/
polygon *get_poly(format *fmt)
{
    int i;
    polygon *poly = 0x0;

    if (!fmt || !fmt->in) return(0x0);

    /* blank out first skip characters of line */
    for (i = 0; i < fmt->skip && file.line[i]; i++) file.line[i] = ' ';

    /* Region or polygon */
    if (strcmp(fmt->in, "Region") == 0
	|| strcmp(fmt->in, "polygon") == 0
	|| strcmp(fmt->in, "spolygon") == 0) {
	poly = rd_poly(fmt);

    /* healpix_weight */
    } else if (strcmp(fmt->in, "healpix_weight") == 0) {
        poly = rd_hpix(fmt);

    /* circle */
    } else if (strcmp(fmt->in, "circle") == 0) {
	poly = rd_circ(fmt);

    /* edges or vertices */
    } else if (strcmp(fmt->in, "edges") == 0
	|| strcmp(fmt->in, "vertices") == 0) {
	poly = rd_edge(fmt);

    /* rectangle */
    } else if (strcmp(fmt->in, "rectangle") == 0) {
	poly = rd_rect(fmt);
    }

    if (poly) {
	/* id number */
	poly->id = fmt->id;
        /* pixel number */
	poly->pixel = fmt->pixel;
	/* weight */
	poly->weight = fmt->weight;
    }

    return(poly);
}

/*------------------------------------------------------------------------------
  Read polygon defined by caps rp and cm as used by garea, gspher et al
  <rp_0> <rp_1> <rp_2> <cm>

   Input: fmt = pointer to format structure.
  Return value: pointer to polyon.
*/
polygon *rd_poly(format *fmt)
{
    int ip, ird;
    polygon *poly = 0x0;

    /* allocate memory for new polygon */
    poly = new_poly(fmt->n);
    if (!poly) {
	WHERE;
	fprintf(stderr, " failed to allocate memory for polygon of %d caps\n", fmt->n);
	return(0x0);
    }

    /* read caps */
    poly->np = fmt->n;
    for (ip = 0; ip < fmt->n; ip++) {
	/* read line of data */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) exit(1);
	/* EOF */
	if (ird == 0) {
	    WHERE;
	    fprintf(stderr, " unexpected EOF: expecting 4 reals\n");
	    free_poly(poly);
	    exit(1);
	}
	/* read rp and cm of cap from line */
	ird = sscanf(file.line, "%Lf %Lf %Lf %Lf",
	    &poly->rp[ip][0], &poly->rp[ip][1], &poly->rp[ip][2], &poly->cm[ip]);
	if (ird != 4) {
	    WHERE;
	    fprintf(stderr, " expecting 4 reals\n");
	    free_poly(poly);
	    exit(1);
	}
    }

    return(poly);
}

/*------------------------------------------------------------------------------
  Read polygon from HEALPix weights; in other words, simply read in the weights
  from the input file and use get_healpix_poly to assign each to the correct
  HEALPix polygon.

  Input: fmt = pointer to format structure.
  Return value: pointer to polygon.
*/

polygon *rd_hpix(format *fmt)
{
  int nweight, ird, nside;
  polygon *poly;
  long double weight;

  if(fmt->nweights == 0) return(0x0);

  if (fmt->n == 1) nweight = fmt->n;
  else {
    fprintf(stderr, "there should only be one weight per line\n");
    return(0x0);
  }

  /* too few weights in line*/
  if (nweight < 0) {
     WHERE;
     fprintf(stderr, " discarding polygon: supposedly line contains %d weights??\n", nweight);
     return(0x0);
  }

  /* allocate memory for new polygon */
  poly = new_poly(5);
  if (!poly) {
     WHERE;
     fprintf(stderr, " failed to allocate memory for polygon of 5 caps\n");
     return(0x0);
  }

  /* read weight from line */
  ird = sscanf(file.line, "%Lf", &weight);
  /* comment line at top of input file */
  if (ird != 1 && fmt->id == 0) {
     return(0x0);
  }
  /* blank/junk lines at EOF */
  if(ird != 1 && fmt->id > (fmt->nweights - 1)) {
    return(0x0);
  }
  /* blank/junk lines in middle of weights */
  if (ird != 1 && fmt->id < fmt->nweights) {
    return(0x0);
  }
  /* EOF */
  if (ird == -1) {
      WHERE;
      fprintf(stderr, " unexpected EOF: expecting 1 real\n");
      free_poly(poly);
      exit(1);
  }
  /* incorrect line */
  if (ird != 1) {
      WHERE;
      fprintf(stderr, " expecting 1 real\n");
      free_poly(poly);
      exit(1);
  }

  /* check whether fmt->nweights (i.e., total number of HEALPix polygons
     at this nside) has been exceeded; note that fmt->id starts at 0 */
  if(ird == 1 && fmt->id > (fmt->nweights - 1)) {
    fprintf(stderr, "error: the number of HEALPix weights in the input file is greater than %d\n", fmt->nweights);
    exit(1);
  }

  nside = get_nside(fmt->nweights);

  /* calculate correct HEALPix polygon */
  poly = get_healpix_poly(nside, fmt->id);
  if(!poly){
    fprintf(stderr, " error in calculating HEALPix polygon vertices\n");
    return(0x0);
  }
  poly->weight = weight;
  
  if(ird == 1 && (poly->cm[0] != 0 || poly->cm[1] != 0)) {
    fmt->id++;
    fmt->pixel = 0;
    fmt->weight = poly->weight;
  }

  return(poly);

}

/*------------------------------------------------------------------------------
  Read polygon from circles, each defined by
  <azimuth> <elevation> <radius>

   Input: fmt = pointer to format structure.
  Return value: pointer to polyon.
*/
polygon *rd_circ(format *fmt)
{
    char unit;
    char *next, *word;
    int i, iang, icirc, ird, ncirc;
    long double angle[3];
    polygon *poly = 0x0;

    /* point word to start of line */
    word = file.line;

    /* number of circles */
    if (fmt->n > 0 || fmt->single == 1) {
	ncirc = fmt->n;
    /* read number of circles from contents of line */
    } else {
	ncirc = get_nang(fmt, word, 2) / fmt->nn;
	if (ncirc == 0) return(0x0);
    }

    /* too few circles */
    if (ncirc < 0) {
	WHERE;
	fprintf(stderr, " discarding polygon: supposedly line contains %d circles??\n", ncirc);
	return(0x0);
    }

    /* allocate memory for new polygon */
    poly = new_poly(ncirc);
    if (!poly) {
	WHERE;
	fprintf(stderr, " failed to allocate memory for polygon of %d caps\n", ncirc);
	return(0x0);
    }

    /* read circles */
    iang = 0;
    poly->np = ncirc;
    for (icirc = 0; icirc < ncirc; icirc++) {
	/* read azimuth, elevation, radius of axis of circle from line */
	for (i = 0; i < 3; i++) {
	    unit = fmt->inunitp;
	    if (i == 2 && fmt->inunitp == 'h') unit = 'd';
	    ird = rdangle(word, &next, unit, &angle[i]);
	    if (ird != 1) {
		if (iang == 0 || (iang == 1 && fmt->inunitp != 'h')) {
		    blankline++;
		} else {
		    WHERE;
		    fprintf(stderr, " expecting %d, found %d angles\n",
			fmt->nn * ncirc, iang);
		    exit(1);
		}
		free_poly(poly);
		return(0x0);
	    }
	    /* scale angle to radians */
	    if (i == 1 && fmt->inunitp == 'h') unit = 'd';
	    scale(&angle[i], unit, 'r');
	    /* increment angle count */
	    iang++;
	    /* point word to next word */
	    word = next;
	}
	/* unnaturally small or large radius may indicate a problem */
	if (fabsl(angle[2]) > PI) {
	    WHERE;
	    if (icirc > 0) fprintf(stderr, " %d'th", icirc + 1);
	    fprintf(stderr, " circle has radius %.16Lg deg", places(angle[2] * 180./PI, 14));
	    if (angle[2] < -PI) {
		fprintf(stderr, " < 180 deg\n");
	    } else if (angle[2] > PI) {
		fprintf(stderr, " > 180 deg\n");
	    }
	}
	/* convert azimuth, elevation, radius to rp, cm */
	circ_to_rpcm(angle, poly->rp[icirc], &poly->cm[icirc]);
    }

    return(poly);
}

/*------------------------------------------------------------------------------
  Read polygon from vertices and edge-points
  <az_0> <el_0> <az_1> <el_1> ... <az_n> <el_n>

   Input: fmt = pointer to format structure.
  Return value: pointer to polyon.
*/
polygon *rd_edge(format *fmt)
{
/* whether to allow a polygon with a single cap to be specified with one edge */
#define ONEEDGESPECIAL	1
/* number of extra edges to allocate, to allow for expansion */
#define DNV		4
    static vertices *vert = 0x0;
    static int *ev = 0x0;

    const char *blank = " \t\n\r";
    char unit;
    char *next, *word;
    int anti, i, iang, iedge, iev, ird, iv, ive, ivert, nedge, nv, nve, nvert, reverse;
    long double *angle;
    long double az, el;
    polygon *poly = 0x0;

    /* point word to start of line */
    word = file.line;

    /* hack to deal with points winding left- not right-handedly */
    reverse = 0;
    /* check for 'r' which says to reverse order of points */
    while (*word && strchr(blank, *word)) word++;
    if (*word && *word == 'r') {
	reverse = 1;
	word++;
	/* require 'r' to be followed by a blank, to be safe */
	if (!*word || !strchr(blank, *word)) return(0x0);
    }

    /* no edges per line */
    if (fmt->innve == 0) return(0x0);

    /* number of edges */
    if (fmt->n > 0 || fmt->single == 1) {
	nedge = fmt->n;
    /* read number of edges from contents of line */
    } else {
	nedge = get_nang(fmt, word, fmt->nn) / (fmt->nn * fmt->innve);
	if (nedge == 0) return(0x0);
    }

    /* on input each edge contains fmt->innve points */
    nvert = fmt->innve * nedge;

    if (fmt->innve == 1) {
	/* impossible number of vertices */
	if (nvert < 3) {
	    WHERE;
	    fprintf(stderr, " discarding polygon");
	    if (nvert < 0) {
		fprintf(stderr, " supposedly with %d vertices??\n", nvert);
	    } else if (nvert == 0) {
		fprintf(stderr, ": a polygon with 0 vertices is ill-defined\n");
	    } else if (nvert == 1) {
		fprintf(stderr, ": a polygon with 1 vertex is ill-defined\n");
	    } else if (nvert == 2) {
		fprintf(stderr, ": a polygon with 2 vertices is null\n");
	    }
	    return(0x0);
	}
    } else {
	/* impossible number of edges */
	if (nedge <= 0) {
	    WHERE;
	    fprintf(stderr, " discarding polygon");
	    if (nedge < 0) {
		fprintf(stderr, " supposedly with %d edges??\n", nedge);
	    } else if (nedge == 0) {
		fprintf(stderr, ": a polygon with 0 edges is ill-defined\n");
	    }
	    return(0x0);
	}
    }

    /* only 1 point per edge */
    if (fmt->innve == 1) {
	nve = 1;
    /* keep only 2 points per edge */
    } else {
	nve = 2;
	/* keep extra point if only one edge but > 2 points per edge */
	if (ONEEDGESPECIAL && nedge == 1 && fmt->innve > 2) nve = 3;
    }
    nv = nve * nedge;

    /* ensure that vert contains enough space */
    if (!vert || vert->nvmax < nv) {
	free_vert(vert);
	vert = new_vert(nv + fmt->innve * DNV);
	if (!vert) {
	    WHERE;
	    fprintf(stderr, " failed to allocate memory for vertices structure of %d vertices\n", nv + fmt->innve * DNV);
	    return(0x0);
	}
	if (!ev) free(ev);
	ev = (int *)malloc(sizeof(int) * (nv + fmt->innve * DNV));
    }
    /* number of points in vertices structure */
    vert->nv = nv;

    /* read contents of line into vertices structure */
    iang = 0;
    ivert = 0;
    iv = -1;
    iev = 0;
    ev[iev] = 0;
    for (iedge = 0; iedge < nedge; iedge++) {
	for (ive = 0; ive < fmt->innve; ive++) {
	    for (i = 0; i < 2; i++) {
		unit = fmt->inunitp;
		angle = (i == 0)? &az : &el;
		ird = rdangle(word, &next, unit, angle);
		if (ird != 1) {
		    if (iang == 0 || (iang == 1 && fmt->inunitp != 'h')) {
			blankline++;
		    } else {
			WHERE;
			fprintf(stderr, "expecting %d, found %d angles on line %d of %s\n",
			    fmt->nn * fmt->innve * nedge, iang, file.line_number, file.name);
			exit(1);
		    }
		    return(0x0);
		}
		/* scale angle to radians */
		if (i == 1 && fmt->inunitp == 'h') unit = 'd';
		scale(angle, unit, 'r');
		/* increment angle count */
		iang++;
		/* point word at next word */
		word = next;
	    }
	    /* phase azimuth to previous */
	    if (iv >= 0) az -= rint((az - vert->v[iv].az) / TWOPI) * TWOPI;
	    /* pass angle to appropriate element of vertices structure */
	    if (ive == 0) {
		iv++;
		/* first point is vertex */
		vert->v[iv].az = az;
		vert->v[iv].el = el;
		/* intialize midpoint to equal vertex */
		if (nve >= 2) {
		    iv++;
		    vert->v[iv].az = az;
		    vert->v[iv].el = el;
		}
	    } else {
		/* update midpoint if previous midpoint was same as vertex */
		if (vert->v[iv - 1].az == vert->v[iv].az
		 && vert->v[iv - 1].el == vert->v[iv].el) {
		    vert->v[iv].az = az;
		    vert->v[iv].el = el;
		}
		/* extra point if want 3 points and this is last edge point */
		if (nve >= 3 && ive == fmt->innve - 1) {
		    iv++;
		    vert->v[iv].az = az;
		    vert->v[iv].el = el;
		}
	    }
	    /* number of points read so far */
	    ivert++;
	    /* number of points on this connected boundary */
	    ev[iev] = iv + 1;
	}

	/* point word at next non-blank character */
	while (*word && strchr(blank, *word)) word++;
	/* detect separate boundary on new line */
	if (iedge < nedge - 1 && *word == '\n') {
	    /* read line of data */
	    ird = rdline(&file);
	    /* serious error */
	    if (ird == -1) return(0x0);
	    /* EOF */
	    if (ird == 0) {
		WHERE;
		fprintf(stderr, "expecting %d, found %d angles on line %d of %s\n",
		    fmt->nn * fmt->innve * nedge, iang, file.line_number, file.name);
		exit(1);
	    }
	    /* disallow reverse hack */
	    if (reverse) {
		WHERE;
		fprintf(stderr, " reverse hack not supported for multi-boundary polygons\n");
		return(0x0);
	    }
	    /* point word to start of line */
	    word = file.line;
	    /* skip first skip characters of line */
	    for (i = 0; i < fmt->skip && *word; i++) word++;
	    /* increment boundary counter */
	    iev++;
	}
    }

    /* reverse order of points */
    if (reverse) {
	for (iv = 0; iv < nv/2; iv++) {
	    az = vert->v[iv].az;
	    vert->v[iv].az = vert->v[nv - 1 - iv].az;
	    vert->v[nv - 1 - iv].az = az;
	    el = vert->v[iv].el;
	    vert->v[iv].el = vert->v[nv - 1 - iv].el;
	    vert->v[nv - 1 - iv].el = el;
	}
    }

    /* make polygon with nedge boundaries */
    poly = new_poly(nedge);
    if (!poly) {
	WHERE;
	fprintf(stderr, " failed to allocate memory for polygon of %d caps\n", nedge);
	return(0x0);
    }

    /* convert vertices structure to polygon */
    edge_to_poly(vert, nve, ev, poly);

    /* check whether vertices are in antipodes */
    anti = antivert(vert, poly);

    /* vertices in antipodes */
    if (anti == 1) {
	WHERE;
	fprintf(stderr, " warning:");
	fprintf(stderr, " polygon %d may have its vertices ordered left- instead of right-handedly\n", fmt->id);
    }

    return(poly);
}

/*------------------------------------------------------------------------------
  Read polygon from rectangles, each defined by
  <azmin> <azmax> <elmin> <elmax>

   Input: fmt = pointer to format structure.
  Return value: pointer to polyon.
*/
polygon *rd_rect(format *fmt)
{
    char unit;
    char *next, *word;
    int i, iang, irect, ird, nrect;
    long double angle[4];
    polygon *poly = 0x0;
    polygon *poly_rect = 0x0;

    /* point word to start of line */
    word = file.line;

    /* number of rectangles */
    if (fmt->n > 0 || fmt->single == 1) {
	nrect = fmt->n;
    /* read number of rectangles from contents of line */
    } else {
	nrect = get_nang(fmt, word, fmt->nn) / fmt->nn;
	if (nrect == 0) return(0x0);
    }

    /* too few rectangles */
    if (nrect < 0) {
	WHERE;
	fprintf(stderr, " discarding polygon: supposedly line contains %d rectangles??\n", nrect);
	return(0x0);
    }

    /* allocate memory for new polygon */
    poly = new_poly(4 * nrect);
    if (!poly) {
	WHERE;
	fprintf(stderr, " failed to allocate memory for polygon of %d caps\n",
	    4 * nrect);
	return(0x0);
    }
    poly->np=0;

    /* read rectangles */
    iang = 0;
    for (irect = 0; irect < nrect; irect++) {
      /* read azmin, azmax, elmin, elmax of rectangle from line */
	for (i = 0; i < 4; i++) {
	    unit = fmt->inunitp;
	    ird = rdangle(word, &next, unit, &angle[i]);
	    if (ird != 1) {
		if (iang == 0 || (iang == 1 && fmt->inunitp != 'h')) {
		    blankline++;
		} else {
		    fprintf(stderr, " expecting %d, found %d angles on line %d of %s\n",
			fmt->nn * nrect, iang, file.line_number, file.name);
		}
		free_poly(poly);
		return(0x0);
	    }
	    /* scale angle to radians */
	    if (i >= 2 && fmt->inunitp == 'h') unit = 'd';
	    scale(&angle[i], unit, 'r');
	    /* increment angle count */
	    iang++;
	    /* point word to next word */
	    word = next;
	}
		    
	poly_rect=new_poly(4);
	if (!poly_rect) {
	  WHERE;
	  fprintf(stderr, " failed to allocate memory for polygon of 4 caps\n");
	  return(0x0);
	}

	rect_to_poly(angle, poly_rect);
	if(!poly_rect){
	  WHERE;
	  fprintf(stderr, " error in rect_to_poly: polygon is NULL.\n");
	  return(0x0);
	}
	poly_poly(poly,poly_rect,poly);
	free_poly(poly_rect);	
    }
    return(poly);
}
