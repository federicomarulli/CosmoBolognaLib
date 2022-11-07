// SWIG Interface to Func

%module cblFunc

%import "Path.i"
%import "Kernel.i"
%import "FuncGrid.i"

%{
#include "Func.h"
#include "LegendrePolynomials.h"
#include "SphericalHarmonics_Coefficients.h"
%}

%include "Func.h"
%include "LegendrePolynomials.h"
%include "SphericalHarmonics_Coefficients.h"


