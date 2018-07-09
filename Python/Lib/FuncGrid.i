// SWIG Interface to FuncGrid

%module cblFuncGrid

%shared_ptr(cbl::glob::FuncGrid);
%shared_ptr(cbl::glob::FuncGrid2D);

%{
#include "FuncGrid.h"
%}

%include "FuncGrid.h"
