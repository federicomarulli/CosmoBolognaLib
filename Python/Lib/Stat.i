// SWIG Interface to Stat

%module cblStat

%include "stl.i"
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include <std_shared_ptr.i>

%shared_ptr(cbl::statistics::PriorDistribution);
%shared_ptr(cbl::statistics::PosteriorDistribution);
%shared_ptr(cbl::statistics::Prior);
%shared_ptr(cbl::statistics::ModelParameters);
%shared_ptr(cbl::statistics::LikelihoodParameters);
%shared_ptr(cbl::statistics::PosteriorParameters);
%shared_ptr(cbl::statistics::Model);
%shared_ptr(cbl::statistics::Model1D);
%shared_ptr(cbl::statistics::Model2D);
%shared_ptr(cbl::statistics::LikelihoodFunction);
%shared_ptr(cbl::statistics::Likelihood);
%shared_ptr(cbl::statistics::Chi2);
%shared_ptr(cbl::statistics::Sampler);
%shared_ptr(cbl::statistics::Posterior);
%shared_ptr(cbl::statistics::CombinedPosterior);

%{
#include "PriorDistribution.h"
#include "PosteriorDistribution.h"
#include "Prior.h"
#include "ModelParameters.h"
#include "LikelihoodParameters.h"
#include "PosteriorParameters.h"
#include "Model.h"
#include "Model1D.h"
#include "Model2D.h"
#include "LikelihoodFunction.h"
#include "Likelihood.h"
#include "Chi2.h"
#include "Sampler.h"
#include "Posterior.h"
#include "CombinedPosterior.h"

%}

%include "PriorDistribution.h"
%include "PosteriorDistribution.h"
%include "Prior.h"
%include "ModelParameters.h"
%include "LikelihoodParameters.h"
%include "PosteriorParameters.h"
%include "Model.h"
%include "Model1D.h"
%include "Model2D.h"
%include "LikelihoodFunction.h"
%include "Likelihood.h"
%include "Chi2.h"
%include "Sampler.h"
%include "Posterior.h"
%include "CombinedPosterior.h"

%template(ParameterTypeVector) std::vector<cbl::statistics::ParameterType>;
%template(PriorDistributionVector) std::vector<cbl::statistics::PriorDistribution>;
%template(PriorDistributionPtrVector) std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>>;
%template(PosteriorDistributionVector) std::vector<cbl::statistics::PosteriorDistribution>;
%template(PosteriorDistributionPtrVector) std::vector<std::shared_ptr<cbl::statistics::PosteriorDistribution>>;
%template(PosteriorPtrVector) std::vector<std::shared_ptr<cbl::statistics::Posterior>>;
