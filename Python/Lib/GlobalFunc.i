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
#include "CatalogueChainMesh.h"
#include "Chi2.h"
#include "CombinedPosterior.h"
#include "PosteriorParameters.h"
#include "LikelihoodParameters.h"
  
#include "GlobalFunc.h"
%}

%include "GlobalFunc.h"
