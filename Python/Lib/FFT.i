// SWIG Interface to FFTlog

%module cblFFT

%import "Path.i"
%import "Kernel.i"

%{
#include "FFTlog.h"
%}

%include "FFTlog.h"
