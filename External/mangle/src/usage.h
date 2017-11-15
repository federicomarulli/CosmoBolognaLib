/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------
  This is in-lined to ensure uniformity between all programs.
*/
    if (strchr(optstr, 'd')) printf("  -d\t\tadvise defaults and exit\n");

    if (strchr(optstr, 'q')) printf("  -q\t\texecute quietly\n");

    if (strchr(optstr, 'w')) printf("  -w<Wlmfile>\tname of file containing spherical harmonics\n");

    if (strchr(optstr, 'z')) printf("  -z<survey>\tname of survey, or of file containing list of weights\n");

    if (strchr(optstr, 'l')) printf("  -l<lmax>\tmaximum harmonic number\n");

    if (strchr(optstr, 'g')) printf("  -g<lsmooth>\tgaussian smoothing harmonic number (0 = default = no smooth)\n");

    if (strchr(optstr, 'c')) printf("  -c<seed>\tseed random number generator with integer <seed>\n");

    if (strchr(optstr, 'r')) printf("  -r<n>\t\tgenerate <n> random points\n");

    if (strchr(optstr, 'h')) printf("  -h\t\twrite only summary to output\n");

    if (strchr(optstr, 'S')) printf("  -S\t\tself-snap: snap edges only within each polygon\n");

    if (strchr(optstr, 'a')) printf("  -a<angle>[u]\tangle within which to snap cap axes together\n");

    if (strchr(optstr, 'b')) printf("  -b<angle>[u]\tangle within which to snap cap latitudes together\n");

    if (strchr(optstr, 't')) printf("  -t<angle>[u]\tangle within which to snap edge to circle\n");

    if (strchr(optstr, 'y')) printf("  -y<ratio>\tsnap edge to circle only if closer than <ratio> x edge length\n");

    if (strchr(optstr, 'm')) printf("  -m<angle>[u]\tangle within which close vertices are coincident\n");

    if (strchr(optstr, 'j')) printf("  -j[min][,max]\tkeep only polygons with weight in [min, max]\n");

    if (strchr(optstr, 'J')) printf("  -J[min][,max]\tkeep only polygons with id numbers in [min, max]\n");

    if (strchr(optstr, 'k')) printf("  -k[min][,max]\tkeep only polygons with area in [min, max] str\n");

    if (strchr(optstr, 'K')) printf("  -K[min][,max]\tkeep only polygons with pixel numbers in [min, max]\n");

    if (strchr(optstr, 'n')) printf("  -n\t\tintersect polygons of infile1 with those of same id in infile2\n");

    if (strchr(optstr, 's')) printf("  -s<n>\t\tskip first <n> characters of each line of polygon_infiles\n");

    if (strchr(optstr, 'e')) printf("  -e<n>\t\tread only to <n>'th character of each line (0 = no limit)\n");

    if (strchr(optstr, 'v')) {
      printf("  -vo|-vn|-vp\tassign old (o) or new (n) id numbers for output polygons\n");
      printf("             \tor use pixel (p) numbers as id numbers\n");
    }
    if (strchr(optstr, 'f')) {
	printf("  -f\t\tlist frames\n");
	printf("  -f<in>[,<ou>]\tinput, output angular frames\n");
	printf("  -f<azn>,<eln>,<azp>[u]\n\t\t<azn>,<eln> = azimuth, elevation of new pole wrt old frame\n\t\t<azp>,<eln> = azimuth, elevation of old pole wrt new frame\n");
    }

    if (strchr(optstr, 'u')) printf("  -u<in>[,<ou>]\tr radians, d degrees, m arcmin, s arcsec, h hms(RA) & dms(Dec)\n");

    if (strchr(optstr, 'p')) {
	printf("  -p[+|-][<n>]\t<n> digits after the decimal place in output angles\n");
	printf("            \toutput azimuths in: + [0, 2 pi) or - (-pi, pi]\n");
    
    }  
    if (strchr(optstr, 'P')) {
      printf("  -P[scheme][<p>][,<r>]\tpixelization scheme: s simple or d sdsspix\n");
      printf("                       \tpixelize to max resolution of <r>, with <p> polys per pixel\n");
    }
    if (strchr(optstr, 'B')) {
      printf("  -B[bmethod]\tmethod for combining weights in balkanize: l=last weight in polygon list,\n");
      printf("             \ta=add weights together, n=minimum weight, x=maximum weight\n");
    }
    if(strchr(optstr, 'U')) {
      printf("  -U\t\tremove pixelization info by unifying across whole mask\n"); 
    }
    if(strchr(optstr, 'W')) {
      printf("  -W\t\tprint weights in polyid output file rather than id numbers\n"); 
    }
    if (strchr(optstr, 'i')) printf("  -i<f>[<n>][u]\tread polygon_infile in format <f>, with <n> objects per line\n");

    if (strchr(optstr, 'o')) printf("  -o<f>[u]\twrite outfile in format <f>\n");

    if (strchr(optstr, 'o') || strchr(optstr, 'i')) {
	printf("  format <f>:\tc circle, e<i> edges, p polygon, r rectangle, R Region,\n");
	printf("             \ts spolygon, v vertices\n");
       	printf("             \tinput only: h healpix_weight\n");
    }
    if (strchr(optstr, 'o')) {
	printf("             \toutput only: a area, g<i> graphics, i id, m midpoint, w weight\n");
    }
    if (strchr(optstr, 'H')) {
      printf("  -H\t\twrite output file in healpix_weight format\n");
    }

    if (strchr(optstr, 'a') || strchr(optstr, 'b') || strchr(optstr, 't') || strchr(optstr, 'i') || strchr(optstr, 'o'))
	printf("  unit u:\tr radians, d degrees, m arcmin, s arcsec, h hms(RA) & dms(Dec)\n");

    printf("  - for infile\tmeans read from stdin\n");

    printf("  - for outfile\tmeans write to stdout\n");

    printf("mangle documentation is at http://space.mit.edu/home/tegmark/mangle/\n");
