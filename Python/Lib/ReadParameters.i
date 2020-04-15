// SWIG Interface to ReadParameters

%module cblReadParameters

%{
#include "ReadParameters.h"
#include "ParameterFile.h"
%}

%include "ReadParameters.h"
%include "ParameterFile.h"

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

%extend cbl::glob::ParameterFile {
  std::vector<std::string> & __getitem__(const std::string key) {
    return (*($self))[key];
  }
}
