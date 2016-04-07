// SWIG Interface to Modelling

%module cblModelling

%{
#include "ModelBias.h"
#include "Modelling.h"
#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation_monopole.h"
#include "Modelling_TwoPointCorrelation_deprojected.h"
#include "Modelling_TwoPointCorrelation_projected.h"
%}

%include "ModelBias.h"
%include "Modelling.h"
%include "Modelling_TwoPointCorrelation.h"
%include "Modelling_TwoPointCorrelation_monopole.h"
%include "Modelling_TwoPointCorrelation_deprojected.h"
%include "Modelling_TwoPointCorrelation_projected.h"

