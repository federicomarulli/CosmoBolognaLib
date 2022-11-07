// SWIG Interface to LogNormal

%module cblLogNormal

%import "Path.i"
%include "Kernel.i"
%import "Catalogue.i"

%{
#include "Data.h"
#include "Data1D.h"
#include "Data1D_extra.h"
#include "Data1D_collection.h"
#include "TaperedCovarianceMatrix.h"
#include "Chi2.h"
#include "PosteriorParameters.h"
#include "LikelihoodParameters.h"
#include "CombinedPosterior.h"
#include "CatalogueChainMesh.h"
  
#include "LogNormal.h"
#include "LogNormalFull.h"
%}

%include "LogNormal.h"
%include "LogNormalFull.h"
