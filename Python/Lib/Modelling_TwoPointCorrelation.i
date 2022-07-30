// SWIG Interface to Modelling_TwoPointCorrelation

%module cblModelling_TwoPointCorrelation

%include "Path.i"
%include "Kernel.i"
%import "Cosmology.i"
%import "Measure.i"
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

#include "ModelFunction_TwoPointCorrelation1D_angular.h"
#include "ModelFunction_TwoPointCorrelation1D_filtered.h"
#include "ModelFunction_TwoPointCorrelation1D.h"
#include "ModelFunction_TwoPointCorrelation1D_monopole.h"
#include "ModelFunction_TwoPointCorrelation2D_cartesian.h"
#include "ModelFunction_TwoPointCorrelation2D.h"
#include "ModelFunction_TwoPointCorrelation2D_polar.h"
#include "ModelFunction_TwoPointCorrelation_deprojected.h"
#include "ModelFunction_TwoPointCorrelation.h"
#include "ModelFunction_TwoPointCorrelation_multipoles.h"
#include "ModelFunction_TwoPointCorrelation_projected.h"
#include "ModelFunction_TwoPointCorrelation_wedges.h"
#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation1D.h"
#include "Modelling_TwoPointCorrelation1D_monopole.h"
#include "Modelling_TwoPointCorrelation1D_angular.h"
#include "Modelling_TwoPointCorrelation1D_filtered.h"
#include "Modelling_TwoPointCorrelation2D.h"
#include "Modelling_TwoPointCorrelation2D_cartesian.h"
#include "Modelling_TwoPointCorrelation2D_polar.h"
#include "Modelling_TwoPointCorrelation_multipoles.h"
#include "Modelling_TwoPointCorrelation_projected.h"
#include "Modelling_TwoPointCorrelation_deprojected.h"
#include "Modelling_TwoPointCorrelation_wedges.h"
%}

%include "ModelFunction_TwoPointCorrelation1D_angular.h"
%include "ModelFunction_TwoPointCorrelation1D_filtered.h"
%include "ModelFunction_TwoPointCorrelation1D.h"
%include "ModelFunction_TwoPointCorrelation1D_monopole.h"
%include "ModelFunction_TwoPointCorrelation2D_cartesian.h"
%include "ModelFunction_TwoPointCorrelation2D.h"
%include "ModelFunction_TwoPointCorrelation2D_polar.h"
%include "ModelFunction_TwoPointCorrelation_deprojected.h"
%include "ModelFunction_TwoPointCorrelation.h"
%include "ModelFunction_TwoPointCorrelation_multipoles.h"
%include "ModelFunction_TwoPointCorrelation_projected.h"
%include "ModelFunction_TwoPointCorrelation_wedges.h"
%include "Modelling_TwoPointCorrelation.h"
%include "Modelling_TwoPointCorrelation1D.h"
%include "Modelling_TwoPointCorrelation1D_monopole.h"
%include "Modelling_TwoPointCorrelation1D_angular.h"
%include "Modelling_TwoPointCorrelation1D_filtered.h"
%include "Modelling_TwoPointCorrelation2D.h"
%include "Modelling_TwoPointCorrelation2D_cartesian.h"
%include "Modelling_TwoPointCorrelation2D_polar.h"
%include "Modelling_TwoPointCorrelation_multipoles.h"
%include "Modelling_TwoPointCorrelation_projected.h"
%include "Modelling_TwoPointCorrelation_deprojected.h"
%include "Modelling_TwoPointCorrelation_wedges.h"
