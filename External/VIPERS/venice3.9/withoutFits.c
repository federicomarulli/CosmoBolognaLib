#include "main.h"

void *readFits(const Config *para, int *bitpix, int *status, long naxes[2], double (**convert)(void *,long )){

  /* used for the compilation without the FITS supports. Return an error message*/
  
  fprintf(stderr,"\n\n                           ***ERROR***\n\n");
  fprintf(stderr,"FITS file detected. This binary version was not compiled with FITS support.\n");
  fprintf(stderr,"Please re-compile %s with FITS=yes and make sure to have\n",MYNAME); 
  fprintf(stderr,"the CFITSIO library installed: http://heasarc.gsfc.nasa.gov/fitsio/\n");
  fprintf(stderr,"Exiting...\n");
  exit(EXIT_FAILURE);
  
}

