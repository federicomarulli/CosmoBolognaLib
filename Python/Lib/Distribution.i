// SWIG Interface to Distribution

%module cblDistribution

%import "Path.i"
%import "Kernel.i"
%import "FuncGrid.i"
%import "Random.i"

%ignore *::operator[];
%shared_ptr(cbl::glob::Distribution);
%shared_ptr(cbl::glob::CombinedDistribution);

%{
#include "Distribution.h"
#include "CombinedDistribution.h"
%}

%include "Distribution.h"
%include "CombinedDistribution.h"
