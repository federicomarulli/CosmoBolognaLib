// SWIG Interface to Modelling_ThreePointCorrelation

%module cblModelling_ThreePointCorrelation

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


#include "ModelFunction_ThreePointCorrelation_angular_connected.h"
#include "ModelFunction_ThreePointCorrelation_angular_reduced.h"
#include "ModelFunction_ThreePointCorrelation_comoving_connected.h"
#include "ModelFunction_ThreePointCorrelation_comoving_reduced.h"
#include "ModelFunction_ThreePointCorrelation.h"
#include "Modelling_ThreePointCorrelation.h"
#include "Modelling_ThreePointCorrelation_angular_connected.h"
#include "Modelling_ThreePointCorrelation_angular_reduced.h"
#include "Modelling_ThreePointCorrelation_comoving_connected.h"
#include "Modelling_ThreePointCorrelation_comoving_reduced.h"
%}

%include "ModelFunction_ThreePointCorrelation_angular_connected.h"
%include "ModelFunction_ThreePointCorrelation_angular_reduced.h"
%include "ModelFunction_ThreePointCorrelation_comoving_connected.h"
%include "ModelFunction_ThreePointCorrelation_comoving_reduced.h"
%include "ModelFunction_ThreePointCorrelation.h"
%include "Modelling_ThreePointCorrelation.h"
%include "Modelling_ThreePointCorrelation_angular_connected.h"
%include "Modelling_ThreePointCorrelation_angular_reduced.h"
%include "Modelling_ThreePointCorrelation_comoving_connected.h"
%include "Modelling_ThreePointCorrelation_comoving_reduced.h"
