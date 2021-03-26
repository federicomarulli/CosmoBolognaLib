// SWIG Interface to CovarianceMatrix

%module cblCovarianceMatrix

%shared_ptr(cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation1D_monopole);
%shared_ptr(cbl::measure::covmat::CovarianceMatrix_TwoPointCorrelation_multipoles);

%{
#include "CovarianceMatrix_TwoPointCorrelation1D_monopole.h"
#include "CovarianceMatrix_TwoPointCorrelation_multipoles.h"
%}

%include "CovarianceMatrix_TwoPointCorrelation1D_monopole.h"
%include "CovarianceMatrix_TwoPointCorrelation_multipoles.h"
