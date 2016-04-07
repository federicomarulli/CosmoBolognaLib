// SWIG Interface to Stat

%module cblStat

%{
#include "Chain.h"
#include "Prior.h"
#include "Parameter.h"
#include "Model.h"
#include "Chi2.h"
#include "Likelihood.h"
%}

%include "Chain.h"
%include "Prior.h"
%include "Parameter.h"
%include "Model.h"
%include "Chi2.h"
%include "Likelihood.h"

