// SWIG Interface to Measure

%module cblMeasure

%import "Path.i"
%include "Kernel.i"
%import "Data.i"

%{
#include "Data1D_collection.h"
#include "Data1D_extra.h"
#include "Measure.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"
%}

%include "Measure.h"
