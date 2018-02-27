// SWIG Interface to ThreePointCorrelation

%module cblThreePointCorrelation
 /*
%shared_ptr(cosmobl::measure::threept::ThreePointCorrelation);
%shared_ptr(cosmobl::measure::threept::ThreePointCorrelation_angular_connected);
%shared_ptr(cosmobl::measure::threept::ThreePointCorrelation_angular_reduced);
%shared_ptr(cosmobl::measure::threept::ThreePointCorrelation_comoving_connected);
%shared_ptr(cosmobl::measure::threept::ThreePointCorrelation_comoving_reduced);
 */
%{
#include "Triplet.h"
#include "Triplet1D.h"
#include "Triplet2D.h"
#include "Measure.h"
#include "ThreePointCorrelation.h"
#include "ThreePointCorrelation_angular_connected.h"
#include "ThreePointCorrelation_angular_reduced.h"
#include "ThreePointCorrelation_comoving_connected.h"
#include "ThreePointCorrelation_comoving_reduced.h"
%}

%include "Triplet.h"
%include "Triplet1D.h"
%include "Triplet2D.h"
%include "Measure.h"
%include "ThreePointCorrelation.h"
%include "ThreePointCorrelation_angular_connected.h"
%include "ThreePointCorrelation_angular_reduced.h"
%include "ThreePointCorrelation_comoving_connected.h"
%include "ThreePointCorrelation_comoving_reduced.h"
