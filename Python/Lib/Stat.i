// SWIG Interface to Stat

%module cblStat

%shared_ptr(cosmobl::statistics::Chain);
%shared_ptr(cosmobl::statistics::Prior);
%shared_ptr(cosmobl::statistics::Posterior);
%shared_ptr(cosmobl::statistics::Parameter);
%shared_ptr(cosmobl::statistics::BaseParameter);
%shared_ptr(cosmobl::statistics::DerivedParameter);
%shared_ptr(cosmobl::statistics::Model);
%shared_ptr(cosmobl::statistics::Model1D);
%shared_ptr(cosmobl::statistics::Model2D);
%shared_ptr(cosmobl::statistics::LikelihoodParameters);
%shared_ptr(cosmobl::statistics::Likelihood);

%{
#include "Model.h"
#include "Model1D.h"
#include "Model2D.h"
#include "Chain.h"
#include "Prior.h"
#include "Posterior.h"
#include "Parameter.h"
#include "BaseParameter.h"
#include "DerivedParameter.h"
#include "LikelihoodParameters.h"
#include "Sampler.h"
#include "LikelihoodFunction.h"
#include "Likelihood.h"
%}

%include "Model.h"
%include "Model1D.h"
%include "Model2D.h"
%include "Chain.h"
%include "Prior.h"
%include "Posterior.h"
%include "Parameter.h"
%include "BaseParameter.h"
%include "DerivedParameter.h"
%include "LikelihoodParameters.h"
%include "Sampler.h"
%include "LikelihoodFunction.h"
%include "Likelihood.h"

%template(PriorVector) std::vector<cosmobl::statistics::Prior>;
%template(ParameterVector) std::vector<cosmobl::statistics::Parameter>;
%template(BaseParameterVector) std::vector<cosmobl::statistics::BaseParameter>;
%template(DerivedParameterVector) std::vector<cosmobl::statistics::DerivedParameter>;
