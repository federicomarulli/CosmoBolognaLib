// SWIG Interface to Distribution

%module cblDistribution

%shared_ptr(cbl::glob::Distribution);
%shared_ptr(cbl::glob::CombinedDistribution);

%{
#include "Distribution.h"
#include "CombinedDistribution.h"
%}

%include "Distribution.h"
%include "CombinedDistribution.h"
