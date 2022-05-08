// SWIG Interface to Modelling_PowerSpectrumAngular

%module cblModelling_PowerSpectrumAngular

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
  
#include "ModelFunction_PowerSpectrum_Angular.h"
#include "Modelling_PowerSpectrum_Angular.h"
%}

%include "ModelFunction_PowerSpectrum_Angular.h"
%include "Modelling_PowerSpectrum_Angular.h"
