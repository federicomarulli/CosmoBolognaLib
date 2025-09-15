// SWIG Interface to GlobalFunc

%module cblGlobalFunc

%include "Path.i"
%include "Kernel.i"
%import "Catalogue.i"

%{
#include "Data.h"
#include "Data1D.h"
#include "Data1D_collection.h"
#include "Data1D_extra.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
#include "Cosmology.h"
#include "LCDM.h"
#include "kLCDM.h"
#include "w0CDM.h"
#include "kw0CDM.h"
#include "w0waCDM.h"
#include "kw0waCDM.h"
#include "LnuCDM.h"
#include "kLnuCDM.h"
#include "w0nuCDM.h"
#include "kw0nuCDM.h"
#include "w0wanuCDM.h"
#include "kw0wanuCDM.h" 
#include "Chi2.h"
#include "CombinedPosterior.h"
#include "PosteriorParameters.h"
#include "LikelihoodParameters.h"
  
#include "GlobalFunc.h"
%}

%include "GlobalFunc.h"
