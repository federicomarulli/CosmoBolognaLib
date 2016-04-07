// SWIG Interface to ThreePointCorrelation

%module cblThreePointCorrelation

%{
#include "Triplet.h"
#include "Triplet1D.h"
#include "Triplet2D.h"
#include "ThreePointCorrelation.h"
#include "ThreePointCorrelation_angular_connected.h"
#include "ThreePointCorrelation_angular_reduced.h"
#include "ThreePointCorrelation_comoving_connected.h"
#include "ThreePointCorrelation_comoving_reduced.h"
%}

%include "Triplet.h"
%include "Triplet1D.h"
%include "Triplet2D.h"
%include "ThreePointCorrelation.h"
%include "ThreePointCorrelation_angular_connected.h"
%include "ThreePointCorrelation_angular_reduced.h"
%include "ThreePointCorrelation_comoving_connected.h"
%include "ThreePointCorrelation_comoving_reduced.h"
