// SWIG Interface to TwoPointCorrelation

%module cblTwoPointCorrelation

%{
#include "Pair.h"
#include "Pair1D.h"
#include "Pair2D.h"
#include "TwoPointCorrelation.h"
#include "TwoPointCorrelation1D.h"
#include "TwoPointCorrelation1D_angular.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "TwoPointCorrelation2D.h"
#include "TwoPointCorrelation2D_cartesian.h"
#include "TwoPointCorrelation2D_polar.h"
#include "TwoPointCorrelation_projected.h"
#include "TwoPointCorrelation_deprojected.h"
#include "TwoPointCorrelation_multipoles.h"
#include "TwoPointCorrelation_wedges.h"
#include "TwoPointCorrelation1D_filtered.h"
%}

%include "Pair.h"
%include "Pair1D.h"
%include "Pair2D.h"
%include "TwoPointCorrelation.h"
%include "TwoPointCorrelation1D.h"
%include "TwoPointCorrelation1D_angular.h"
%include "TwoPointCorrelation1D_monopole.h"
%include "TwoPointCorrelation2D.h"
%include "TwoPointCorrelation2D_cartesian.h"
%include "TwoPointCorrelation2D_polar.h"
%include "TwoPointCorrelation_projected.h"
%include "TwoPointCorrelation_deprojected.h"
%include "TwoPointCorrelation_multipoles.h"
%include "TwoPointCorrelation_wedges.h"
%include "TwoPointCorrelation1D_filtered.h"
