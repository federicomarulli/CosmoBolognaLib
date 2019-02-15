// ========================================================================================
// Example code: how to perform a Bayesian fit to a set of data points with a generic model
// ========================================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Posterior.h"

using namespace std;

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


// =====================================================================


// this example model has 4 parameters: A, B, C, D; A and B are free parameters, C is fixed, D is a derived parameter

vector<double> model_function (const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter)
{
  // the object Cosmology, used in this example to compute Omega_matter
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);

  vector<double> model(x.size(), 0.);
  for (size_t i=0; i<x.size(); ++i)
    model[i] = parameter[0]*x[i]+parameter[1]+parameter[2]*cosm.Omega_matter(); // the model

  parameter[3] = parameter[0]+parameter[1]+parameter[2]; // parameter[3] is a derived parameter

  return model;
}


// =====================================================================


int main () {

  try {

    // --- set the input/output file/directories ---
    
    const string dir_input = cbl::par::DirLoc+"../input/";
    const string dir_output = cbl::par::DirLoc+"../output/";
    const string file_data = "data.dat";
    const string file_output_start = "model_starting_values.dat";
    const string file_output_bestfit = "model_bestfit.dat";
    

    // --- construct the dataset by reading an input file ---
    
    const cbl::data::Data1D data(dir_input+file_data);
    shared_ptr<cbl::data::Data> ptr_data = make_shared<cbl::data::Data1D>(data);

    
    // --- set the model to construct the likelihood ---
    
    // number of model parameters
    const int nparameters = 4;

    // names of the model parameters
    const vector<string> parNames = {"A", "B", "C", "D"};
    
    // starting values for A and B, and fixed values for C
    double valA = 1., valB = 1., valC = 1.; 
     
    // vector containing the 4 model parameters
    vector<cbl::statistics::ParameterType> parType(nparameters-1, cbl::statistics::ParameterType::_Base_);
    parType.emplace_back(cbl::statistics::ParameterType::_Derived_);   

    // set the stuff used to construct the model: here an object of class cosmology, just as an example 
    const cbl::cosmology::Cosmology cosmology;
    auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);

    // construct the model
    const cbl::statistics::Model1D model(&model_function, nparameters, parType, parNames, ptr_modelInput);
    auto ptr_model = make_shared<cbl::statistics::Model1D>(model);


    // --- construct and maximize the likelihood ---
    
    // define a Gaussian likelihood
    cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Error_);

    // write the model used to construct the likelihood, at the starting parameter values (A=1, B=1, C=1=fixed)
    const vector<double> start = {valA, valB};
    likelihood.parameters()->fix(2, valC); // fix the parameter C, i.e. parameter[2], to valC=1
    likelihood.write_model(dir_output, file_output_start, start);

    // limits for A and B
    const double minA = -10., maxA = 10.;
    const double minB = -10., maxB = 10.;
    const vector<vector<double>> limits = { {minA, maxA} , {minB, maxB} };
    
    // maximize the likelihood and write the output
    likelihood.maximize(start, limits);
    likelihood.write_model_at_bestfit(dir_output, file_output_bestfit);

    
    // --- construct the priors ---
    
    // (remember to give different seed to priors!)
    auto prior_A = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minA, maxA, 666));
    auto prior_B = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minB, maxB, 999));
    auto prior_C = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Constant_, valC));
    const vector<shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions = {prior_A, prior_B, prior_C};

    
    // --- construct, maximize and sample the posterior ---
    
    cbl::statistics::Posterior posterior(prior_distributions, likelihood, 696);

    // maximize the posterior
    //posterior.maximize({valA, valB});

    // sample the posterior (starting the MCMC chain from the maximum of the posterior to speed up the chain convergence)
    const int nwalkers = 10;
    const int chain_size = 5000;
    posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB});
    posterior.sample_stretch_move(2);

    // show the median MCMC values of the four parameters on screen
    cout << endl;
    for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
      cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;

    // show all the MCMC statistics on screen
    const int burn_in = 0;
    const int thin = 1;
    posterior.show_results(burn_in, thin);

    // store the chain ouputs
    posterior.write_results(cbl::par::DirLoc+"../output/", "chains_linear_relation", burn_in, thin);

    // store the best-fit model
    posterior.write_model_from_chain(cbl::par::DirLoc+"../output/", "model_from_chain.dat", {}, {}, burn_in, thin);

  }
  
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}


