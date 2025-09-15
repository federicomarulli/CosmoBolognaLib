// SWIG Interface to ThreePointCorrelation

%module cblThreePointCorrelation

%import "Path.i"
%import "Kernel.i"
%import "Catalogue.i"
  
 /*
%shared_ptr(cbl::measure::threept::ThreePointCorrelation);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_angular_connected);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_angular_reduced);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_comoving_connected);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_comoving_reduced);
 */
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
  
#include "Triplet.h"
#include "Triplet1D.h"
#include "Triplet2D.h"
#include "Measure.h"
#include "ThreePointCorrelation.h"
#include "ThreePointCorrelation_angular_connected.h"
#include "ThreePointCorrelation_angular_reduced.h"
#include "ThreePointCorrelation_comoving_connected.h"
#include "ThreePointCorrelation_comoving_reduced.h"
#include "ThreePointCorrelation_comoving_multipoles.h"
#include "ThreePointCorrelation_comoving_multipoles_single.h"
#include "ThreePointCorrelation_comoving_multipoles_all.h"
%}

%include "Triplet.h"
%include "Triplet1D.h"
%include "Triplet2D.h"
%include "Measure.h"
%include "ThreePointCorrelation.h"
%include "ThreePointCorrelation_angular_connected.h"
%include "ThreePointCorrelation_angular_reduced.h"
%include "ThreePointCorrelation_comoving_connected.h"
%include "ThreePointCorrelation_comoving_reduced.h"
%include "ThreePointCorrelation_comoving_multipoles.h"
%include "ThreePointCorrelation_comoving_multipoles_single.h"
%include "ThreePointCorrelation_comoving_multipoles_all.h"

