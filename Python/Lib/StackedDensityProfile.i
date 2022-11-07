// SWIG Interface to StackedDensityProfile

%module cblStackedDensityProfile

%import "Path.i"
%import "Kernel.i"
%import "Cosmology.i"
%import "Catalogue.i"
%import "Measure.i"

%{
#include "Data.h"
#include "Data1D_collection.h"
#include "Data1D_extra.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
#include "CatalogueChainMesh.h"
#include "Chi2.h"
#include "LikelihoodParameters.h"
#include "PosteriorParameters.h"
#include "CombinedPosterior.h"
  
#include "StackedDensityProfile.h"
%}

%include "StackedDensityProfile.h"
