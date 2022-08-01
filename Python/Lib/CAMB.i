// SWIG Interface to CAMB

%module cblCAMB

%import "Path.i"

%{
#include "CAMB.h"
%}

%include "CAMB.h"
