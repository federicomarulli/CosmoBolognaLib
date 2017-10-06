// SWIG Interface to Cosmology

%module cblCosmology

%shared_ptr(cosmobl::cosmology::Cosmology);

%{
#include "Cosmology.h"
%}

%include "Cosmology.h"

%template(CosmoParVector) std::vector<cosmobl::cosmology::CosmoPar>;
