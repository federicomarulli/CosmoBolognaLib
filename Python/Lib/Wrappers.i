// SWIG Interface to Wrapper

%module cblWrappers

%{
#include "CUBAwrapper.h"
#include "FITSwrapper.h"
  //#include "GSLfunction.h"
#include "GSLwrapper.h"
%}

%include "CUBAwrapper.h"
%include "FITSwrapper.h"
 //%include "GSLfunction.h"
%include "GSLwrapper.h"
