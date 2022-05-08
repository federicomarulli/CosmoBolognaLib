#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_integration.h>

#define FAILURE 0
#define SUCCESS 1

#define PI    3.14159265358979323846
#define TWOPI 6.283185307179586476925287

/* useful constants */
#define c     299792.458
#define G     6.67300e-11

/* WMAP5 Cosmology */
#define H0 72 /* all values definied with h = 0.72 */
#define Omega_M 0.258
#define Omega_L 0.742

//#define EPS 1.0e-10
#define INF   1.0e30
#define ODD   0
#define EVEN  1
#define LEAF  0
#define NODE  1
#define RADEC 0
#define CART  1

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN
#define SWAP(a,b) {swap = (a); (a) = (b); (b) = swap;}
#define SQUARE(a) ((a)*(a))

#define NFIELD    500
#define NCHAR     20
#define NVERTICES 100

#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)     atoi(array+NCHAR*(col-1))
#define getCharValue(array,col)    array+NCHAR*(col-1)
#define getLine(array,i)           array+NFIELD*NCHAR*i

/*----------------------------------------------------------------*
 *New types                                                       *
 *----------------------------------------------------------------*/

typedef struct Config
{
  char fileRegInName[1000];
  char fileCatInName[1000];
  char fileOutName[1000];
  char fileNofZName[1000];
  int nz, zrange;
  int nx,ny,format;
  int xcol,ycol;
  int coordType, constDen;
  size_t npart,seed;
  double min[2], max[2], zmin, zmax;
  int minDefinied[2];
  int maxDefinied[2];
  
  /* cosmology */
  double a[4];
} Config;


typedef struct Complex
{
    double re;
    double im;
} Complex;

typedef struct Polygon
{
  int N, id;
  // double x[NVERTICES];
  //double y[NVERTICES];
  //double xmin[2];
  //double xmax[2];
  double *x;
  double *y;
  double *xmin;
  double *xmax;
} Polygon;

typedef struct Node
{
  int type, *root, id, *polysAll, SplitDim;
  double SplitValue;
  size_t Nnodes, Npolys, NpolysAll;
  int *poly_id;
  void *Left, *Right;
} Node;

/*----------------------------------------------------------------*
 *Global variables                                                *
 *----------------------------------------------------------------*/

char   MYNAME[100];
size_t IDERR;
double EPS;

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

void testPython();

int mask2d(const Config *para);
int flagCat(const Config *para);
int randomCat(const Config *para);

/*----------------------------------------------------------------*
 *Initialization                                                   *
 *----------------------------------------------------------------*/

int readParameters(int argc, char **argv, Config *para);

/*----------------------------------------------------------------*
 *Utils - geometric                                               *
 *----------------------------------------------------------------*/

int insidePolygon(Polygon *p,int Npoly, double x0,double y0,double x,double y, int *poly_id);
int insidePolygonTree(Node *polyTree, double x0[2], double x[2], int *poly_id);
Polygon *readPolygonFile(FILE *fileIn, int *Npolys, Node *polyTree);
Node *readPolygonFileTree(FILE *fileIn, double xmin[2], double xmax[2]);
Node *createNode(Polygon *polys, size_t Npolys, double minArea, int SplitDim, double xmin[2], double xmax[2], int firstCall);
void free_Polygon(Polygon *polygon, size_t N);
void free_Node(Node *node);
void cpyPolygon(Polygon *a, Polygon *b);

/*----------------------------------------------------------------*
 *Utils - numeric                                                 *
 *----------------------------------------------------------------*/

double distComo(double z, const double a[4]);
double drdz(double x, void * params);
double dvdz(double z, const double a[4]);
double distAngSpher(const double RA1, double DEC1, double RA2, double DEC2);
void rotate(double x0, double y0, double x, double y, double *xrot, double *yrot, double angle, int spherical);
double determineMachineEpsilon();
size_t determineSize_tError();
gsl_rng *randomInitialize(size_t seed);
FILE *fopenAndCheck(const char *filename,char *mode);
int getStrings(char *line, char *strings, char *delimit, size_t *N);
void printCount(const size_t *count, const size_t *total,  const size_t step);
int checkFileExt(const char *s1, const char *s2);
int roundToNi(double a);
int compareDoubles(const void *a,const void *b);

/*----------------------------------------------------------------*
 *FITS                                                            *
 *----------------------------------------------------------------*/

void *readFits(const Config *para, int *bitpix, int *status, long naxes[2], double (**convert)(void *,long ));
double convertCHAR(void *table, long i);
double convertSHORT(void *table, long i);
double convertLONG(void *table, long i);
double convertFLOAT(void *table, long i);
double convertDOUBLE(void *table, long i);

#define BYTE_IMG      8
#define SHORT_IMG    16
#define LONG_IMG     32
#define LONGLONG_IMG 64 
#define FLOAT_IMG   -32 
#define DOUBLE_IMG  -64
