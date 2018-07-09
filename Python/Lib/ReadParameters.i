// SWIG Interface to ReadParameters

%module cblReadParameters

%{
#include "ReadParameters.h"
%}

%include "ReadParameters.h"

%extend cbl::glob::ReadParameters
{
  %template(findString) find< std::string >;
  %template(findBool) find< bool >;
  %template(findInt) find< int >;
  %template(findLong) find< long >;
  %template(findFloat) find< float >;
  %template(findDouble) find< double >;

  %template(findVectorString) find_vector< std::string >;
  %template(findVectorBool) find_vector< bool >;
  %template(findVectorInt) find_vector< int >;
  %template(findVectorLong) find_vector< long >;
  %template(findVectorFloat) find_vector< float >;
  %template(findVectorDouble) find_vector< double >;
}
