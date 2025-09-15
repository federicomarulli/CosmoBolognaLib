// ========================================================================================================================
// Example code: how to how to model the baryon acoustic oscillations in the monopole of the two-point correlation function
// ========================================================================================================================

#include "LCDM.h"
#include "Data1D.h"
#include "Modelling_TwoPointCorrelation1D_monopole.h"

using namespace std;


int main () {

  try {

    // ------------------------------------------------------------------------------------------
    // ---------------- set the cosmological parameters (as in Anderson et al. 2013) ------------
    // ------------------------------------------------------------------------------------------
  
    const double hh = 0.7;
    const double Omega_b = 0.0224/(hh*hh);
    const double Omega_CDM = 0.274; 
    const double n_s = 0.96;
    const double sigma8 = 0.8;

    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    cosmology->set_parameters({"Omega_CDM", "Omega_baryon", "n_spec", "sigma8"}, {Omega_CDM, Omega_b, n_s, sigma8});

  
    // ----------------------------------------------
    // ------------- reading the dataset ------------
    // ----------------------------------------------

    const string dir_input = "../input/";
    const string dir_output = "../output/";
    const string dir_chains = dir_output+"chains/";
    const string MK = "mkdir -p "+dir_output+" "+dir_chains; if (system(MK.c_str())) {}
    
    const string file_xi = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon.dat";
    const string file_cov = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon_cov.dat";

    const int skipped_lines = 1;
    
    const auto twop_dataset = make_shared<cbl::data::Data1D>(cbl::data::Data1D(file_xi, skipped_lines)); 
    twop_dataset->set_covariance(file_cov);
    
  
    // ------------------------------------------------------
    // ------------- set the modelling structure ------------
    // ------------------------------------------------------
  
    cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole model_twop(twop_dataset);
    
    
    // --------------------------------------------------------------------------
    // ------------- set the underlying dark matter clustering model ------------
    // --------------------------------------------------------------------------

    const double redshift = 0.57;
    const double sigmaNL = 8.;
  
    model_twop.set_data_model(cosmology, redshift, "CAMB", sigmaNL, sigmaNL, false, 2.);
    
  
    // -------------------------------------------------------
    // ------------- set the priors and the model ------------
    // -------------------------------------------------------

    cbl::statistics::PriorDistribution alpha_prior(cbl::glob::DistributionType::_Uniform_, 0.7, 1.3, 4133); 
    cbl::statistics::PriorDistribution BB_prior(cbl::glob::DistributionType::_Uniform_, 0., 5.);
    cbl::statistics::PriorDistribution A0_prior(cbl::glob::DistributionType::_Uniform_, -100., 100.);
    cbl::statistics::PriorDistribution A1_prior(cbl::glob::DistributionType::_Uniform_, -100., 100.);
    cbl::statistics::PriorDistribution A2_prior(cbl::glob::DistributionType::_Uniform_, -1000., 1000.);
  
    model_twop.set_model_BAO(alpha_prior, BB_prior, A0_prior, A1_prior, A2_prior);
    
  
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    const int chain_size = 1000;
    const int nwalkers = 100;

    const string chain_file = "Anderson_2013_CMASSDR11_monopole_prerecon_chains.dat";
  
    const double xmin = 40., xmax = 160.;
    model_twop.set_fit_range(xmin, xmax);
    
    model_twop.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Covariance_);

    model_twop.sample_posterior(chain_size, nwalkers);
    
    const int burn_in = 100; // discard the first 100 chain steps 
    const int thin = 10;     // take 1 step every 10
    model_twop.show_results(burn_in, thin);
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}

