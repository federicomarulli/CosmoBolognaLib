// SWIG Interface to AngularPowerSpectrum

%module cblPowerSpectrumAngular

%import "Path.i"
%import "Kernel.i"
%import "TwoPointCorrelation.i"  

%shared_ptr(cbl::measure::PowerSpectrum_Angular);

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
#include "LikelihoodParameters.h"
#include "PosteriorParameters.h"
#include "CombinedPosterior.h"
#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D.h"
#include "TwoPointCorrelation1D_angular.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "TwoPointCorrelation2D.h"
#include "TwoPointCorrelation2D_cartesian.h"
#include "TwoPointCorrelation2D_polar.h"
#include "TwoPointCorrelation_projected.h"
#include "TwoPointCorrelation_deprojected.h"
#include "TwoPointCorrelation_multipoles_direct.h"
#include "TwoPointCorrelation_multipoles_integrated.h"
#include "TwoPointCorrelation_wedges.h"
#include "TwoPointCorrelation1D_filtered.h"
#include "TwoPointCorrelationCross.h"
#include "TwoPointCorrelationCross1D.h"
#include "TwoPointCorrelationCross1D_monopole.h"
  
#include "PowerSpectrum_Angular.h"
#include "TwoPointCorrelation1D_angular.h"
%}

%include "PowerSpectrum_Angular.h"
%include "TwoPointCorrelation1D_angular.h"
