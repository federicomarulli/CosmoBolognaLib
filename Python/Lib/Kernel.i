// SWIG Interface to Kernel

%module cblKernel

%import "Path.i"

%{
#include "Kernel.h"
%}

%include "Kernel.h"

%template(nint) cbl::nint<int>;
%template(ndouble) cbl::nint<double>;
