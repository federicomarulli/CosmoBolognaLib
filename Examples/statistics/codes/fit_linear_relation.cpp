// ==================================================
// Example code: how to model linear relation to data
// ==================================================

#include "Data1D.h"
#include "Likelihood.h"

using namespace cosmobl;
using namespace statistics;
using namespace random;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


double linear_relation (const double x, const shared_ptr<void> fixed_parameter, const vector<double> parameters)
{
  auto par = static_pointer_cast<vector<double>>(fixed_parameter);

  return parameters[0] * x + parameters[1] * (*par)[0];
}


int main () {

  // construct the dataset
  Data1D data("../input/data.dat");
  auto ptr_data = make_shared<Data1D>(data);

  // parameters of the model
  double A, B;

  // prior limits 
  double minA = -10., maxA = 10.;
  double minB = -10., maxB = 10.;

  // construct the free parameters A and B 
  Parameter parameterA(A, minA, maxA, false, "parameter A");
  Parameter parameterB(B, minB, maxB, false, "parameter B");

  // set one fixed parameter
  vector<double> fixed_parameter = {2.3};
  auto ptr_fixed_parameter = make_shared<vector<double>>(fixed_parameter);
  
  // construct the model
  Model1D model({parameterA, parameterB}, ptr_fixed_parameter, &linear_relation);
  auto ptr_model = make_shared<Model1D>(model);

  // define a gaussian likelihood
  Likelihood likelihood(ptr_data, ptr_model, LikelihoodType::_GaussianLikelihood_Error_);

  // minimize the logarithm of the Likelihood
  likelihood.minimize_LogLikelihood(true);

  // sample the likelihood
  int nchains = 20, chain_size = 10000;
  likelihood.sample_stretch_move(nchains, chain_size, 123123);

  // write the chain output
  likelihood.write_chain("../output/", "chains_linear_relation.dat");

  return 0;
}

 
