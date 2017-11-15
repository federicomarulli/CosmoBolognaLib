/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "manglefn.h"

/*------------------------------------------------------------------------------
  Parse option <fopt> to -f<fopt> switch.

  Included inline to ensure uniform processing of arguments by all programs.
*/
int parse_fopt(void)
{
    int iscan, itr;
    char unit;
    char opt = 'f';
    char ins[16], outs[16];
    char *ch, *word;
    long double angle, az, el, azn, eln, azp;

    itr = 0;

    /* try -f<inframe>,<outframe> format */
    iscan = sscanf(fopt, "%[^ \t,]%*[ \t,]%[^ \t,]", ins, outs);
    if (iscan >= 1) {
	/* skip zeroth frame, which is "unknown" */
        fmt.inframe = strdictl(ins, &frames[1]) + 1;
	if (iscan == 1) {
	    fmt.outframe = fmt.inframe;
	} else if (iscan == 2) {
	    fmt.outframe = strdictl(outs, &frames[1]) + 1;
	}

	/* transformation angles */
	if (fmt.inframe > 0 && fmt.outframe > 0) {
	    /* pole */
	    az = 0.;	el = 90.;
	    /* az, el of old pole wrt new frame */
	    fframe_(&fmt.inframe, &az, &el, &fmt.outframe, &fmt.azp, &fmt.eln);
	    /* az, el of new pole wrt old frame */
	    fframe_(&fmt.outframe, &az, &el, &fmt.inframe, &fmt.azn, &fmt.eln);
	}
    }

    /* try -f<azn>,<eln>,<azp>[u] format */
    if (iscan < 1 || fmt.inframe <= 0 || fmt.outframe <= 0) {
	word = fopt;
	do {
	    /* look for angular unit in word */
	    for (ch = word; *ch && !strchr(UNITS, *ch); ch++);
	    fmt.trunit = (*ch)? *ch : 'd';
	    ch = word;
	    for (iscan = 0; iscan < 3; iscan++) {
		if (rdangle(ch, &ch, fmt.trunit, &angle) != 1) break;
		switch (iscan) {
		case 0:	azn = angle;	break;
		case 1:	eln = angle;	break;
		case 2:	azp = angle;	break;
		}
		while (*ch && strchr(" \t,", *ch)) ch++;
	    }
	    /* success */
	    if (iscan == 3) {
		/* check angular unit is 4th quantity in word */
		if (*ch && *ch != ':') {
		    iscan = sscanf(ch, "%c", &unit);
		    if (unit != fmt.trunit) {
			fprintf(stderr, "-%c%s: angular unit %c must be one of %s\n", opt, fopt, unit, UNITS);
			exit(1);
		    }
		    ch++;
		}
		/* convert to degrees */
		scale(&azn, fmt.trunit, 'd');
		scale(&eln, (fmt.trunit == 'h')? 'd' : fmt.trunit, 'd');
		scale(&azp, fmt.trunit, 'd');
		fmt.trunit = 'd';
		/* number of transformations */
		itr++;
		/* initialize transformation */
		if (itr == 1) {
		    fmt.azn = azn;
		    fmt.eln = eln;
		    fmt.azp = azp;
		/* accumulate transformation */
		} else {
		    azell_(&azn, &eln, &azp, &fmt.azp, &fmt.eln, &fmt.azn, &fmt.azn, &fmt.eln, &fmt.azp);
		}
		/* flag transformation is custom */
		fmt.outframe = -1;
	    /* failure */
	    } else {
		fmt.outframe = 0;
	    }
	    if (*ch && *ch == ':') word = ch + 1;
	} while (*ch && *ch == ':');

	if (fmt.outframe == 0) {
	    fprintf(stderr, "-%c%s: expecting -%c<inframe>[,<outframe>]\n", opt, fopt, opt);
	    fprintf(stderr, "Input and output angular frames <inframe> and <outframe> should be one of:\n");
	    for (fmt.inframe = 1; frames[fmt.inframe]; fmt.inframe++) {
		fprintf(stderr, " %s", frames[fmt.inframe]);
	    }
	    fprintf(stderr, "\n");
	    fprintf(stderr, "Alternatively, a general rotation may be specified by -%c<azn>,<eln>,<azp>[u]\n", opt);
	    fprintf(stderr, "<azn> = azimuth of new pole wrt old frame\n");
	    fprintf(stderr, "<eln> = elevation of new pole wrt old frame\n");
	    fprintf(stderr, "      = elevation of old pole wrt new frame\n");
	    fprintf(stderr, "<azp> = azimuth of old pole wrt new frame\n");
	    fprintf(stderr, "    u = r radians, d degrees, m arcmin, s arcsec, h hms(RA) & dms(Dec)\n");
	    fprintf(stderr, "A sequence of rotations -%c<azn_1>,<eln_1>,<azp_1>[u1]\n", opt);
	    fprintf(stderr, "            followed by -%c<azn_2>,<eln_2>,<azp_2>[u2] may be specified by\n", opt);
	    fprintf(stderr, " -%c<azn_1>,<eln_1>,<azp_1>[u1]:<azn_2>,<eln_2>,<azp_2>[u2]\n", opt);

	    exit(1);
	}
    }

    return(itr);
}
