// SWIG Interface to Cosmology

%module cblCosmology

%include "Path.i"
%import "Kernel.i"
%import "Data.i"
%import "Stat.i"

%shared_ptr(cbl::cosmology::Cosmology);

%{
#include "Data.h"
#include "Data1D.h"
#include "Data1D_extra.h"
#include "Data1D_collection.h"
#include "TaperedCovarianceMatrix.h"

#include "PosteriorParameters.h"
#include "LikelihoodParameters.h"
#include "Likelihood.h"
#include "Chi2.h"
#include "CombinedPosterior.h"
  
#include "EisensteinHu.h"
#include "CAMB.h"  
#include "Cosmology.h"
#include "HaloProfile.h"
#include "SuperSampleCovariance.h"
%}

%include "Likelihood.h"
%include "EisensteinHu.h"
%include "CAMB.h"
%include "Cosmology.h"
%include "HaloProfile.h"
%include "SuperSampleCovariance.h"

%template(CosmologicalParameterVector) std::vector<cbl::cosmology::CosmologicalParameter>;
