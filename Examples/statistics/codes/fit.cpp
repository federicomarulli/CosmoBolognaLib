// ====================================================================
// Example code: how to fit a set of data points with a generic model
// ====================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Likelihood.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

vector<double> model_function (const vector<double> x, const shared_ptr<void> modelInput, const vector<double> parameter)
{
  cosmobl::cosmology::Cosmology &cosm = *static_pointer_cast<cosmobl::cosmology::Cosmology>(modelInput);

  vector<double> model(x.size(), 0.);
  for (size_t i=0; i<x.size(); i++)
    model[i] = parameter[0] * x[i] + parameter[1] + parameter[2] * cosm.Omega_matter();

  return model;
}


int main () {

  try {
  
    // construct the dataset (by reading an input file, in this example case)
    cosmobl::data::Data1D data(cosmobl::par::DirLoc+"../input/data.dat");
    auto ptr_data = make_shared<cosmobl::data::Data1D>(data);

    // set parameters of the model: A and B are free, C is fixed
    double A = 1., B = 1., C = 1.;

    // prior limits for A and B
    double minA = -10., maxA = 10.;
    double minB = -10., maxB = 10.;
  
    // construction the free parameters A and B
    auto parameterA = make_shared<cosmobl::statistics::BaseParameter>(cosmobl::statistics::BaseParameter(minA, maxA, 24314, "parameter A"));
    auto parameterB = make_shared<cosmobl::statistics::BaseParameter>(cosmobl::statistics::BaseParameter(minB, maxB, 5635, "parameter B"));

    // construction of the fixed parameter C
    auto parameterC = make_shared<cosmobl::statistics::BaseParameter>(cosmobl::statistics::BaseParameter(C, "parameter C"));
  
    // set the stuff used to construct the model: here an object of class cosmology, just as an example 
    cosmobl::cosmology::Cosmology cosmology;
    auto ptr_modelInput = make_shared<cosmobl::cosmology::Cosmology>(cosmology);
  
    // construct the model
    cosmobl::statistics::Model1D model(&model_function, ptr_modelInput);
    auto ptr_model = make_shared<cosmobl::statistics::Model1D>(model);

    // define a gaussian likelihood
    cosmobl::statistics::Likelihood likelihood(ptr_data, ptr_model, {parameterA, parameterB, parameterC}, cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Error_);

    // minimize the logarithm of the Likelihood
    vector<double> guess = {A, B, C};

    likelihood.maximize(guess);

    // sample the likelihood
    int nwalkers=100, chain_size = 10000;
    double radius = 1.e-3;
    int seed = 43333213;
    likelihood.initialize_chains(chain_size, nwalkers, seed, guess, radius);
    likelihood.sample_stretch_move(seed);

    // write the chain output
    const int burn_in = 100;
    const int thin = 10;

    likelihood.show_results(burn_in, thin, seed);
    likelihood.write_results(cosmobl::par::DirLoc+"../output/", "chains_linear_relation.dat", burn_in, thin, seed);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
}

 
