// SWIG Interface to TwoPointCorrelation

%module cblTwoPointCorrelation

%import "Path.i"
%import "Kernel.i"
%import "Catalogue.i"

%shared_ptr(cbl::pairs::Pair);
%shared_ptr(cbl::pairs::Pair1D);
%shared_ptr(cbl::pairs::Pair1D_angular);
%shared_ptr(cbl::pairs::Pair1D_angular_lin);
%shared_ptr(cbl::pairs::Pair1D_angular_log);
%shared_ptr(cbl::pairs::Pair1D_comoving);
%shared_ptr(cbl::pairs::Pair1D_comoving_lin);
%shared_ptr(cbl::pairs::Pair1D_comoving_log);
%shared_ptr(cbl::pairs::Pair1D_comoving_multipoles);
%shared_ptr(cbl::pairs::Pair1D_comoving_multipoles_lin);
%shared_ptr(cbl::pairs::Pair1D_comoving_multipoles_log);
%shared_ptr(cbl::pairs::Pair2D);
%shared_ptr(cbl::pairs::Pair2D_comovingCartesian);
%shared_ptr(cbl::pairs::Pair2D_comovingCartesian_linlin);
%shared_ptr(cbl::pairs::Pair2D_comovingCartesian_linlog);
%shared_ptr(cbl::pairs::Pair2D_comovingCartesian_loglin);
%shared_ptr(cbl::pairs::Pair2D_comovingCartesian_loglog);
%shared_ptr(cbl::pairs::Pair2D_comovingPolar);
%shared_ptr(cbl::pairs::Pair2D_comovingPolar_linlin);
%shared_ptr(cbl::pairs::Pair2D_comovingPolar_linlog);
%shared_ptr(cbl::pairs::Pair2D_comovingPolar_loglin);
%shared_ptr(cbl::pairs::Pair2D_comovingPolar_loglog);

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
#include "CatalogueChainMesh.h"
  
#include "Pair.h"
#include "Pair1D.h"
#include "Pair2D.h"
#include "Measure.h"
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
%}

%include "Pair.h"
%include "Pair1D.h"
%include "Pair2D.h"
%include "Measure.h"
%include "TwoPointCorrelation.h"
%include "TwoPointCorrelation1D.h"
%include "TwoPointCorrelation1D_angular.h"
%include "TwoPointCorrelation1D_monopole.h"
%include "TwoPointCorrelation2D.h"
%include "TwoPointCorrelation2D_cartesian.h"
%include "TwoPointCorrelation2D_polar.h"
%include "TwoPointCorrelation_projected.h"
%include "TwoPointCorrelation_deprojected.h"
%include "TwoPointCorrelation_multipoles_direct.h"
%include "TwoPointCorrelation_multipoles_integrated.h"
%include "TwoPointCorrelation_wedges.h"
%include "TwoPointCorrelation1D_filtered.h"
%include "TwoPointCorrelationCross.h"
%include "TwoPointCorrelationCross1D.h"
%include "TwoPointCorrelationCross1D_monopole.h"
