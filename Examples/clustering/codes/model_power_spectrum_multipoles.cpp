// ======================================================================
// Example code: how to how to model the multipoles of the power spectrum
// ======================================================================

#include "w0wanuCDM.h"
#include "Data1D.h"
#include "Modelling_PowerSpectrum_multipoles.h"

using namespace std;


int main () {

  try {

    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------
  
    const double h = 0.67;
    const double OmegaB = 0.05;
    const double OmegaNu = 0.00143655;
    const double OmegaM = 0.32;
    const double OmegaRad = 0.0;
    const double As = 2.12605e-9;
    const double scalar_pivot = 0.05;
    const double n_spec = 0.96;
    const double tau = 0.0925;
    const double w0 = -1.;
    const double wa = 0.;
    const double massless_nu = 2.046;
    const int Nmassive_nu = 1;
    const bool unit = false;    
    auto cosmology = make_shared<cbl::cosmology::w0wanuCDM>(OmegaM, OmegaB, OmegaRad, h, As, scalar_pivot, n_spec, tau, w0, wa, OmegaNu, massless_nu, Nmassive_nu, unit);

  
    // ----------------------------------------------
    // ------------- reading the dataset ------------
    // ----------------------------------------------

    const string dir_input = "../input/";
    const string dir_output = "../output/";
    const string dir_chains = dir_output+"chains/";
    const string MK = "mkdir -p "+dir_output+" "+dir_chains; if (system(MK.c_str())) {}
    const string file_data = "Pl_Kaiser_z1.dat";
    
    const auto Powspec_dataset = make_shared<cbl::data::Data1D>(cbl::data::Data1D(dir_input+file_data)); 

    
    // ------------------------------------------------------
    // ------------- set the modelling structure ------------
    // ------------------------------------------------------
  
    cbl::modelling::powspec::Modelling_PowerSpectrum_multipoles model_multipoles(Powspec_dataset);

    const vector<double> zz{1.};
    const vector<int> ll{0, 2, 4};
    const vector<double> bias{1.46};
  
    model_multipoles.set_data_model(cosmology, zz, ll, bias, "EisensteinHu");
    
   
    // --------------------------------------------------------------
    // ------------- set the parameters and their priors ------------
    // --------------------------------------------------------------
    
    vector<std::string> CosmoPar = {"h", "n_spec"};
    cbl::statistics::PriorDistribution hh_prior, ns_prior;
    hh_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Gaussian_, {0.67, 0.01}, 0.4, 1.);
    ns_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Gaussian_, {0.96, 0.004}, 0.9, 1.);
    vector<cbl::statistics::PriorDistribution> CosmoPriors = {hh_prior, ns_prior}; 
    
    model_multipoles.set_model_multipoles_Kaiser(CosmoPar, CosmoPriors); 

    
    // ---------------------------------------------
    // ------------- set the likelihood ------------
    // ---------------------------------------------  

    model_multipoles.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);
    

    // -------------------------------------------
    // ------------- fit the data set ------------
    // -------------------------------------------

    const vector<double> start = {0.67, 0.97};
    model_multipoles.maximize_posterior(start);

    // retrieve and show the best-fit parameter values
    for (unsigned int i=0; i<model_multipoles.posterior()->parameters()->nparameters(); ++i)
      cout << "the best-fit value of " << model_multipoles.posterior()->parameters()->name(i) << " is " << model_multipoles.posterior()->parameters()->bestfit_value(i) << endl;

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}

