// ====================================================================
// Example code: how to model a set of data points with a generic model
// ====================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Likelihood.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

double model_function (const double x, const shared_ptr<void> modelInput, const vector<double> parameter)
{
  cosmobl::cosmology::Cosmology &cosm = *static_pointer_cast<cosmobl::cosmology::Cosmology>(modelInput);

  return parameter[0] * x + parameter[1] + parameter[2] * cosm.Omega_matter();
}


int main () {

  try {
  
    // construct the dataset (by reading an input file, in this example case)
    cosmobl::data::Data1D data(cosmobl::par::DirLoc+"../input/data.dat");
    auto ptr_data = make_shared<cosmobl::data::Data1D>(data);

    // set parameters of the model: A and B are free, C is fixed
    double A, B, C = 1.;

    // prior limits for A and B
    double minA = -10., maxA = 10.;
    double minB = -10., maxB = 10.;
  
    // construction the free parameters A and B
    cosmobl::statistics::Parameter parameterA(A, minA, maxA, cosmobl::statistics::_free_, "parameter A");
    cosmobl::statistics::Parameter parameterB(B, minB, maxB, cosmobl::statistics::_free_, "parameter B");

    // construction of the fixed parameter C
    cosmobl::statistics::Parameter parameterC(C, cosmobl::statistics::_fixed_, "parameter C");
  
    // set the stuff used to construct the model: here an object of class cosmology, just as an example 
    cosmobl::cosmology::Cosmology cosmology;
    auto ptr_modelInput = make_shared<cosmobl::cosmology::Cosmology>(cosmology);
  
    // construct the model
    cosmobl::statistics::Model1D model({parameterA, parameterB, parameterC}, ptr_modelInput, &model_function);
    auto ptr_model = make_shared<cosmobl::statistics::Model1D>(model);

    // define a gaussian likelihood
    cosmobl::statistics::Likelihood likelihood(ptr_data, ptr_model, cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Error_);

    // minimize the logarithm of the Likelihood
    likelihood.minimize_LogLikelihood(true);

    // sample the likelihood
    int nchains = 20, chain_size = 10000;
    likelihood.sample_stretch_move(nchains, chain_size, 123123);

    // write the chain output
    likelihood.write_chain(cosmobl::par::DirLoc+"../output/", "chains_linear_relation.dat");

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }

  return 0;
}

 
