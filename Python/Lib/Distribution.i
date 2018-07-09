// SWIG Interface to Distribution

%module cblDistribution

%shared_ptr(cbl::glob::Distribution);

%{
#include "Distribution.h"
%}

%include "Distribution.h"
