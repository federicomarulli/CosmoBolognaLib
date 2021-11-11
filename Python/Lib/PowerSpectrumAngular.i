// SWIG Interface to AngularPowerSpectrum

%module cblPowerSpectrumAngular

%shared_ptr(cbl::measure::PowerSpectrum_Angular);

%{
#include "PowerSpectrum_Angular.h"
#include "TwoPointCorrelation1D_angular.h"
%}

%include "PowerSpectrum_Angular.h"
%include "TwoPointCorrelation1D_angular.h"
