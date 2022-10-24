// SWIG Interface to Wrapper

%module cblWrappers

%import "Path.i"
%import "Kernel.i"

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
