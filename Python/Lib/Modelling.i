// SWIG Interface to Modelling

%module cblModelling

%include "Path.i"
%include "Kernel.i"
%import "Cosmology.i"

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
  
#include "Modelling.h"
#include "CombinedModelling.h"
#include "Modelling_Distribution.h"
%}

%include "Modelling.h"
%include "CombinedModelling.h"
%include "Modelling_Distribution.h"

%template(ModellingPtrVector) std::vector<std::shared_ptr<cbl::modelling::Modelling>>;
