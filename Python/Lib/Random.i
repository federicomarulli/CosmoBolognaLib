// SWIG Interface to RandomNumbers

%module cblRandomNumbers

%shared_ptr(cbl::glob::RandomNumbers);
%shared_ptr(cbl::glob::ConstantRandomNumbers);
%shared_ptr(cbl::glob::UniformRandomNumbers);
%shared_ptr(cbl::glob::PoissonRandomNumbers);
%shared_ptr(cbl::glob::NormalRandomNumbers);
%shared_ptr(cbl::glob::DiscreteRandomNumbers);
%shared_ptr(cbl::glob::DistributionRandomNumbers);
%shared_ptr(cbl::glob::CustomDistributionRandomNumbers);

%{
#include "RandomNumbers.h"
%}

%include "RandomNumbers.h"
