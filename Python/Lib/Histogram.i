// SWIG Interface to Histogram

%module cblHistogram

%shared_ptr(cbl::glob::Histogram);
%shared_ptr(cbl::glob::Histogram1D);
%shared_ptr(cbl::glob::Histogram2D);

%{
#include "Histogram.h"
%}

%include "Histogram.h"
