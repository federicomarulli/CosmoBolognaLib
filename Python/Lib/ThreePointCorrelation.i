// SWIG Interface to ThreePointCorrelation

%module cblThreePointCorrelation
 /*
%shared_ptr(cbl::measure::threept::ThreePointCorrelation);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_angular_connected);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_angular_reduced);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_comoving_connected);
%shared_ptr(cbl::measure::threept::ThreePointCorrelation_comoving_reduced);
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
#include "ThreePointCorrelation_comoving_multipoles.h"
#include "ThreePointCorrelation_comoving_multipoles_single.h"
#include "ThreePointCorrelation_comoving_multipoles_all.h"
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
%include "ThreePointCorrelation_comoving_multipoles.h"
%include "ThreePointCorrelation_comoving_multipoles_single.h"
%include "ThreePointCorrelation_comoving_multipoles_all.h"

