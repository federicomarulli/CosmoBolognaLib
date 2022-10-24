// SWIG Interface to ChainMesh

%module cblChainMesh

%import "Path.i"
%import "Kernel.i"

%{
#include "ChainMesh.h"
%}

%include "ChainMesh.h"
