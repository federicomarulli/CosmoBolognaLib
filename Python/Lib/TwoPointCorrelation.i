// SWIG Interface to TwoPointCorrelation

%module cblTwoPointCorrelation

%shared_ptr(cosmobl::twopt::TwoPointCorrelation);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation1D);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation1D_angular);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation1D_monopole);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation2D);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation2D_cartesian);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation2D_polar);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation_projected);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation_deprojected);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation_multipoles);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation_wedges);
%shared_ptr(cosmobl::twopt::TwoPointCorrelation1D_filtered);

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
