// SWIG Interface to Field3D

%module cblField

%import "Path.i"
%import "Kernel.i"

%{
#include "Field3D.h"
%}

%include "Field3D.h"
