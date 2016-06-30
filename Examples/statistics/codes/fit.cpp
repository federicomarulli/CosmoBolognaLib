// ====================================================================
// Example code: how to model a set of data points with a generic model
// ====================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Likelihood.h"

using namespace cosmobl;
using namespace cosmology;
using namespace statistics;
using namespace random;
using namespace data;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


double model_function (const double x, const shared_ptr<void> modelInput, const vector<double> parameter)
{
  Cosmology &cosm = *static_pointer_cast<Cosmology>(modelInput);

  return parameter[0] * x + parameter[1] + parameter[2] * cosm.Omega_matter();
}


int main () {

  // construct the dataset (by reading an input file, in this example case)
  Data1D data(par::DirLoc+"../input/data.dat");
  auto ptr_data = make_shared<Data1D>(data);

  // set parameters of the model: A and B are free, C is fixed
  double A, B, C = 1.;

  // prior limits for A and B
  double minA = -10., maxA = 10.;
  double minB = -10., maxB = 10.;
  
  // construction the free parameters A and B
  Parameter parameterA(A, minA, maxA, _free_, "parameter A");
  Parameter parameterB(B, minB, maxB, _free_, "parameter B");

  // construction of the fixed parameter C
  Parameter parameterC(C, _fixed_, "parameter C");
  
  // set the stuff used to construct the model: here an object of class cosmology, just as an example 
  Cosmology cosmology;
  auto ptr_modelInput = make_shared<Cosmology>(cosmology);
  
  // construct the model
  Model1D model({parameterA, parameterB, parameterC}, ptr_modelInput, &model_function);
  auto ptr_model = make_shared<Model1D>(model);

  // define a gaussian likelihood
  Likelihood likelihood(ptr_data, ptr_model, LikelihoodType::_GaussianLikelihood_Error_);

  // minimize the logarithm of the Likelihood
  likelihood.minimize_LogLikelihood(true);

  // sample the likelihood
  int nchains = 20, chain_size = 10000;
  likelihood.sample_stretch_move(nchains, chain_size, 123123);

  // write the chain output
  likelihood.write_chain(par::DirLoc+"../output/", "chains_linear_relation.dat");

  return 0;
}

 
