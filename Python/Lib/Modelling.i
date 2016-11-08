// SWIG Interface to Modelling

%module cblModelling

%{
#include "Modelling.h"
#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation1D.h"
#include "Modelling_TwoPointCorrelation2D.h"
#include "Modelling_TwoPointCorrelation_monopole.h"
#include "Modelling_TwoPointCorrelation_deprojected.h"
#include "Modelling_TwoPointCorrelation_projected.h"
#include "Modelling_TwoPointCorrelation_cartesian.h"
%}

%include "Modelling.h"
%include "Modelling_TwoPointCorrelation.h"
%include "Modelling_TwoPointCorrelation1D.h"
%include "Modelling_TwoPointCorrelation2D.h"
%include "Modelling_TwoPointCorrelation_monopole.h"
%include "Modelling_TwoPointCorrelation_deprojected.h"
%include "Modelling_TwoPointCorrelation_projected.h"
%include "Modelling_TwoPointCorrelation_cartesian.h"

