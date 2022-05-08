// SWIG Interface to Cosmology

%module cblCosmology

%shared_ptr(cbl::cosmology::Cosmology);

%{
#include "Cosmology.h"
%}

%include "Cosmology.h"

%template(CosmologicalParameterVector) std::vector<cbl::cosmology::CosmologicalParameter>;
