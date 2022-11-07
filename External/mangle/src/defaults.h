/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
#include "defines.h"

/* LOCAL VARIABLES */

/* maximum harmonic */
static int lmax = LMAX;
/* smoothing parameters */
static long double lsmooth = LSMOOTH, esmooth = ESMOOTH;

/* name of file containing harmonics */
static char *Wlm_filename = 0x0;

/* name of survey */
static char *survey = 0x0;

/* option in -f<fopt> command line switch */
static char *fopt = 0x0;

/* seed for random number generator */
static unsigned int seed = SEED;
/* whether seed was read */
static int seed_read = 0;

/* number of random points to generate */
static int nrandom = NRANDOM;

/* write only summary to output */
static int summary = 0;

/* selfsnap = 0 snaps edges of polygons against edges of all other polygons;
   selfsnap = 1 snaps edges of polygons only against edges of the same polygon */
static int selfsnap = 0;

/* data format */
static format fmt = {
	0x0,		/* keyword defining the input data format */
	0x0,		/* default format of the output data */
	SKIP,		/* skip first skip characters of each line of data */
	END,		/* last character to read from line of data */
	0,		/* keyword does not define precisely one polygon */
	0,		/* the number of thingys defined by keyword */
	0,		/* the number of numbers per thingy */
	NVE,		/* the input number of points per edge */
	0,		/* controls interpetation of nve */
	NVE,		/* the output number of points per edge */
	0,		/* id number of current polygon */
	'o',		/* whether to use old or new id number */
        0,              /* default pixel number */
	1.,		/* weight of current polygon */
	INUNITP,	/* default unit of angles in input polygon data */
	OUTUNITP,	/* default unit of angles in output polygon data */
	0,		/* angular frame of input az, el data */
	0,		/* angular frame of output az, el data */
	INUNIT,		/* default unit of input az, el data */
	OUTUNIT,	/* default unit of output az, el data */
	-1,		/* digits after decimal point in output angles (-1 = automatic) */
	OUTPHASE,	/* '-' or '+' to make output azimuth in interval (-pi, pi] or [0, 2 pi) */
	AZN,		/* default			       */
	ELN,		/*	transformation		       */
	AZP,		/* 		between angular frames */
	TRUNIT,		/* unit of transformation angles */
	0,              /* default number of weights in healpix_weight input file */
};

/* GLOBAL VARIABLES */

/* default is to be verbose */
#ifdef DEBUG
int verbose = 2;
#else
int verbose = 1;
#endif

/* counter for input files read */
int infiles = 0;

/* tolerances */
long double axtol = AXTOL;		/* snap angle for axis */
char axunit = AXUNIT;		/* unit of snap angle for axis */
long double btol = BTOL;		/* snap angle for latitude */
char bunit = BUNIT;		/* unit of snap angle for latitude */
long double thtol = THTOL;		/* snap angle for edge */
char thunit = THUNIT;		/* unit of snap angle for edge */
long double ytol = YTOL;		/* edge to length tolerance */
long double mtol = MTOL;		/* tolerance angle for multiple intersections */
char munit = MUNIT;		/* unit of tolerance angle for multiple intersections */

/* whether min, max weight are turned on */
int is_weight_min = 0;
int is_weight_max = 0;
/* min, max weight to keep */
long double weight_min;
long double weight_max;

/* whether min, max area are turned on */
int is_area_min = 0;
int is_area_max = 0;
/* min, max area to keep */
long double area_min;
long double area_max;

/* whether min, max id are turned on */
int is_id_min = 0;
int is_id_max = 0;
/* min, max id to keep */
int id_min;
int id_max;

/* whether min, max pixel are turned on */
int is_pixel_min = 0;
int is_pixel_max = 0;
/* min, max pixel to keep */
int pixel_min;
int pixel_max;

/* whether to take intersection of polygons in input files */
int intersect = 0;

/* dictionary of keywords */
/* THE NAMES OF FORMATS MUST AGREE WITH THE INDICES IN defines.h */
char *keywords[] = {
    "area",
    "circle",
    "edges",
    "graphics",
    "healpix_weight",
    "id",
    "midpoint",
    "polygon",
    "rectangle",
    "Region",
    "spolygon",
    "vertices",
    "weight",
    "list",
    "pixelization",
    "skip",
    "end",
    "unit",
    "balkanized",
    "snapped",
    '\0'
};

/* dictionary of frames */
/* THE ORDER OF FRAMES MUST AGREE WITH THAT IN frames.par ! */
char *frames[] = {
    "unknown",
    "eqB1950",
    "eqJ2000",
    "galactic",
    "ecliptic",
    "ecliptic2k",
    "sdss",
    '\0'
};

/*pixelization info*/
int res_max=RES_MAX;                  /*maximum resolution allowed for pixelization*/
int polys_per_pixel=POLYS_PER_PIXEL;  /*level of pixelization: number of polygons allowed per pixel*/
                                      /*set polys_per_pixel=0 to pixelize everything to max resolution*/
char scheme=SCHEME;                   /*default pixelization scheme*/
int unpixelize=0;                     /*switch for whether unify should unpixelize or not*/
int pixelized=0;                      /*counter for pixelized input files */

int snapped=0;                      /*flag for whether files have been snapped */
int balkanized=0;                      /*flag for whether files have been balkanized */

/*balkanization method*/
char bmethod=BMETHOD;
//Default is for a balkanized polygon to inherit the weight of 
//the last overlapping polygon in the input polygons.  other options
//are to add them together (bmethod=a), or to take the minimum 
//(bmethod=n) or the maximum (bmethod=x) of the weights 

int polyid_weight=0;                     /*0= polyid prints id numbers, 1= polyid prints weights*/
