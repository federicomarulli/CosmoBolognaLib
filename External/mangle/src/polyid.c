/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "inputfile.h"
#include "manglefn.h"
#include "defaults.h"

/* getopt options */
const char *optstr = "dqu:p:P:W";

/* allocate polygons as a global array */
polygon *poly_global[NPOLYSMAX];

/* declared in rdmask */
extern inputfile file;

/* local functions */
void	usage(void);
#ifdef	GCC
int	poly_ids(char *, char *, format *, int npoly, polygon *[npoly]);
#else
int	poly_ids(char *, char *, format *, int npoly, polygon *[/*npoly*/]);
#endif

/*------------------------------------------------------------------------------
  Main program.
*/
int main(int argc, char *argv[])
{
    int ifile, nfiles, npoly, npolys,i;
    polygon **poly;
    poly=poly_global;

    /* parse arguments */
    parse_args(argc, argv);

    /* at least two input and one output filename required as arguments */
    if (argc - optind < 3) {
	if (optind > 1 || argc - optind >= 1) {
	    fprintf(stderr, "%s requires at least 3 arguments: polygon_infile, azel_infile, and outfile\n", argv[0]);
	    usage();
	    exit(1);
	} else {
	    usage();
	    exit(0);
	}
    }

    msg("---------------- polyid ----------------\n");

    /* advise data format */
    advise_fmt(&fmt);

    /* read polygons */
    npoly = 0;
    nfiles = argc - 2 - optind;
    for (ifile = optind; ifile < optind + nfiles; ifile++) {
	npolys = rdmask(argv[ifile], &fmt, NPOLYSMAX - npoly, &poly[npoly]);
	if (npolys == -1) exit(1);
	npoly += npolys;
    }
    if (nfiles >= 2) {
	msg("total of %d polygons read\n", npoly);
    }
    if (npoly == 0) {
      msg("STOP\n");
      exit(0);
    }
    
    if (snapped==0) {
      msg("WARNING: 'snapped' keyword not found in all input files.\n");
      msg("Polygons should be snapped before performing other mangle operations.\n");
    }
    
    /* polygon id numbers */
    npolys = poly_ids(argv[argc - 2], argv[argc - 1], &fmt, npoly, poly);
    if (npolys == -1) exit(1);

  for(i=0;i<npoly;i++){
    free_poly(poly[i]);
  }
    return(0);
}

/*------------------------------------------------------------------------------
*/
void usage(void)
{
    printf("usage:\n");
    printf("polyid [-d] [-q] [-u<inunit>[,<outunit>]] [-p[+|-][<n>]] [-P[scheme][<p>][,<r>]] [-W] polygon_infile1 [polygon_infile2 ...] azel_infile outfile\n");
#include "usage.h"
}

/*------------------------------------------------------------------------------
*/
#include "parse_args.c"

/*------------------------------------------------------------------------------
  Id numbers of polygons containing az, el positions.
  The az, el positions are read from in_filename,
  and the results are written to out_filename.
  Implemented as interpretive read/write, to permit interactive behaviour.

   Input: in_filename = name of file to read from;
			"" or "-" means read from standard input.
	  out_filename = name of file to write to;
			"" or "-" means write to standard output.
	  fmt = pointer to format structure.
	  poly = array of pointers to polygons.
	  npoly = number of polygons in poly array.
  Return value: number of lines written,
		or -1 if error occurred.
*/
int poly_ids(char *in_filename, char *out_filename, format *fmt, int npoly, polygon *poly[/*npoly*/])
{
#define AZEL_STR_LEN	32
    char input[] = "input", output[] = "output";
    char *word, *next;
    char az_str[AZEL_STR_LEN], el_str[AZEL_STR_LEN];
    int i, idmax, idmin, idwidth, ird, len, nid, nids, nid0, nid2, np;
    int *id;
    long double *weight;
    azel v;
    char *out_fn;
    FILE *outfile;
    int *start;
    int *total;
    int *parent_pixels;
    int p, res, max_pixel, ier;
    max_pixel= poly[npoly-1]->pixel; 
    res_max=get_res(max_pixel, scheme);
    max_pixel=pixel_start(res_max+1,scheme);
    
    /* allocate memory for pixel info arrays start and total */ 
    printf("res_max=%d, max_pixel=%d\n",res_max,max_pixel);
    start = (int *) malloc(sizeof(int) * max_pixel);
    if (!start) {
      fprintf(stderr, "polyid: failed to allocate memory for %d integers\n", max_pixel);
      return(-1);
    }
    total = (int *) malloc(sizeof(int) * max_pixel);
    if (!total) {
      fprintf(stderr, "polyid: failed to allocate memory for %d integers\n", max_pixel);
      return(-1);
    }
    parent_pixels = (int *) malloc(sizeof(int) * (res_max+1));
    if (!parent_pixels) {
      fprintf(stderr, "polyid: failed to allocate memory for %d integers\n", res_max+1);
      return(-1);
    }

    /* build lists of starting indices of each pixel and total number of polygons in each pixel*/
    
    ier=pixel_list(npoly, poly, max_pixel, start, total);
    if (ier == -1) {
      fprintf(stderr, "poly_ids: error building pixel index lists\n");
      return(-1);
    }
    
    /* open in_filename for reading */
    if (!in_filename || strcmp(in_filename, "-") == 0) {
	file.file = stdin;
	file.name = input;
    } else {
	file.file = fopen(in_filename, "r");
	if (!file.file) {
	    fprintf(stderr, "cannot open %s for reading\n", in_filename);
	    return(-1);
	}
	file.name = in_filename;
    }
    file.line_number = 0;

    /* open out_filename for writing */
    if (!out_filename || strcmp(out_filename, "-") == 0) {
	outfile = stdout;
	out_fn = output;
    } else {
	outfile = fopen(out_filename, "w");
	if (!outfile) {
	    fprintf(stderr, "cannot open %s for writing\n", out_filename);
	    return(-1);
	}
	out_fn = out_filename;
    }

    /* advise angular units */
    msg("will take units of input az, el angles in %s to be ", file.name);
    switch (fmt->inunit) {
#include "angunit.h"
    }
    msg("\n");
    if (fmt->outunit != fmt->inunit) {
	msg("units of output az, el angles will be ");
	switch (fmt->outunit) {
#include "angunit.h"
	}
	msg("\n");
    }

    /* largest width of polygon id number */
    idmin = 0;
    idmax = 0;
    for (i = 0; i < npoly; i++) {
	if (!poly[i]) continue;
	if (poly[i]->id < idmin) idmin = poly[i]->id;
	if (poly[i]->id > idmax) idmax = poly[i]->id;
    }
    idmin = ((idmin < 0)? floorl(log10l((long double)-idmin)) + 2 : 1);
    idmax = ((idmax > 0)? floorl(log10l((long double)idmax)) + 1 : 1);
    idwidth = ((idmin > idmax)? idmin : idmax);

    /* write header */
    v.az = 0.;
    wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
    len = strlen(az_str);
    if (fmt->outunit == 'h') {
	sprintf(az_str, "az(hms)");
	sprintf(el_str, "el(dms)");
    } else {
	sprintf(az_str, "az(%c)", fmt->outunit);
	sprintf(el_str, "el(%c)", fmt->outunit);
    }
    fprintf(outfile, "%*s %*s", len, az_str, len, el_str);
    if (npoly > 0){
      if(polyid_weight==1){
	fprintf(outfile, " polygon_weights");
      }
      else{
	fprintf(outfile, " polygon_ids");	
      }
    }
    fprintf(outfile, "\n");

    /* interpretive read/write loop */
    np = 0;
    nid = 0;
    nids = 0;
    nid0 = 0;
    nid2 = 0;
    while (1) {
	/* read line */
	ird = rdline(&file);
	/* serious error */
	if (ird == -1) return(-1);
	/* EOF */
	if (ird == 0) break;

	/* read <az> */
	word = file.line;
	ird = rdangle(word, &next, fmt->inunit, &v.az);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* read <el> */
	word = next;
	ird = rdangle(word, &next, fmt->inunit, &v.el);
	/* skip header */
	if (ird != 1 && np == 0) continue;
	/* otherwise exit on unrecognized characters */
	if (ird != 1) break;

	/* convert az and el from input units to radians */
	scale_azel(&v, fmt->inunit, 'r');
	
	//find out what pixel the az el point is in at the maximum resolution
	p=which_pixel(v.az, v.el, res_max, scheme);
	//get the list of all the possible parent pixels
	get_parent_pixels(p, parent_pixels, scheme);

	nid=0;
	for(res=res_max;res>=0;res--){
	  p=parent_pixels[res];
	  //if this pixel isn't in the polygon list, go to next parent pixel
	  if(total[p]==0) continue;
	  // id numbers of the polygons containing position az, el 
	  nid = poly_id(total[p], &poly[start[p]], v.az, v.el, &id, &weight);
	}
	
	/* convert az and el from radians to output units */
	scale_azel(&v, 'r', fmt->outunit);

	/* write result */
	wrangle(v.az, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, az_str);
	wrangle(v.el, fmt->outunit, fmt->outprecision, AZEL_STR_LEN, el_str);
	fprintf(outfile, "%s %s", az_str, el_str);
	for (i = 0; i < nid; i++) {
	  if(polyid_weight==1){
	    fprintf(outfile, " %.18Lg", weight[i]);
	  } else{
	    fprintf(outfile, " %*d", idwidth, id[i]);
	  }
	}
	fprintf(outfile, "\n");
	fflush(outfile);

        /* increment counters of results */
	np++;
	nids += nid;
	if (nid == 0) {
	  nid0++;
	} else if (nid >= 2) {
	  nid2++;
	}
    }

    /* advise */
    if (nid0 > 0) msg("%d points were not inside any polygon\n", nid0);
    if (nid2 > 0) msg("%d points were inside >= 2 polygons\n", nid2);

    if (outfile != stdout) {
      if(polyid_weight==1){
	msg("polyid: %d weights at %d positions written to %s\n", nids, np, out_fn);
      } else {
	msg("polyid: %d id numbers at %d positions written to %s\n", nids, np, out_fn);
      }
    }
    
    free(start);
    free(total);
    free(parent_pixels);

    return(np);
}
