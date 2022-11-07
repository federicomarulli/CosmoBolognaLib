// SWIG Interface to Modelling_DensityProfile

%module cblModelling_DensityProfile

%include "Path.i"
%include "Kernel.i"
%import "Cosmology.i"
%include "Modelling.i"
  
%{
#include "Data.h"
#include "Data1D.h"
#include "Data1D_collection.h"
#include "Data1D_extra.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
#include "Chi2.h"
#include "LikelihoodParameters.h"
#include "PosteriorParameters.h"
#include "CombinedPosterior.h"
  
#include "Modelling_DensityProfile.h"
%}

%include "Modelling_DensityProfile.h"
