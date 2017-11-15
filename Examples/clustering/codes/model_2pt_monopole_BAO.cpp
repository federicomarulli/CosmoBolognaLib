// ========================================================================================================================
// Example code: how to how to model the baryon acoustic oscillations in the monopole of the two-point correlation function
// ========================================================================================================================

#include "Data1D.h"
#include "Modelling_TwoPointCorrelation1D_monopole.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {

    // ------------------------------------------------------------------------------------------
    // ---------------- set the cosmological parameters (as in Anderson et al. 2013) ------------
    // ------------------------------------------------------------------------------------------
  
    const double hh = 0.7;
    const double Omega_b = 0.0224/(hh*hh);
    const double Omega_M = 0.274; 
    const double n_s = 0.96;
    const double sigma8 = 0.8;

    cosmobl::cosmology::Cosmology cosmology;
    cosmology.set_parameters({cosmobl::cosmology::CosmoPar::_Omega_matter_LCDM_, cosmobl::cosmology::CosmoPar::_Omega_baryon_, cosmobl::cosmology::CosmoPar::_hh_, cosmobl::cosmology::CosmoPar::_n_spec_, cosmobl::cosmology::CosmoPar::_sigma8_}, {Omega_M, Omega_b, hh, n_s, sigma8});

  
    // ----------------------------------------------
    // ------------- reading the dataset ------------
    // ----------------------------------------------

    const string dir_input = cosmobl::par::DirLoc+"../input/";
    const string dir_output = cosmobl::par::DirLoc+"../output/";
    const string dir_chains = dir_output+"chains/";
    const string MK = "mkdir -p "+dir_output+" "+dir_chains; if (system(MK.c_str())) {}
    
    const string file_xi = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon.dat";
    const string file_cov = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon_cov.dat";

    const int skipped_lines = 1;
    
    const auto twop_dataset = make_shared<cosmobl::data::Data1D>(cosmobl::data::Data1D(file_xi, skipped_lines)); 
    twop_dataset->set_covariance(file_cov);
    
  
    // ------------------------------------------------------
    // ------------- set the modelling structure ------------
    // ------------------------------------------------------
  
    cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole model_twop(twop_dataset);

    
    // --------------------------------------------------------------------------
    // ------------- set the underlying dark matter clustering model ------------
    // --------------------------------------------------------------------------

    const double redshift = 0.57;
    const double sigmaNL = 8.;
  
    model_twop.set_data_model(cosmology, redshift, "CAMB", sigmaNL, sigmaNL, false, 2.);

  
    // -------------------------------------------------------
    // ------------- set the priors and the model ------------
    // -------------------------------------------------------

    const cosmobl::statistics::Prior alpha_prior(cosmobl::glob::DistributionType::_UniformDistribution_, 0.7, 1.3); 
    const cosmobl::statistics::Prior BB_prior(cosmobl::glob::DistributionType::_UniformDistribution_, 0., 5.);
    const cosmobl::statistics::Prior A0_prior(cosmobl::glob::DistributionType::_UniformDistribution_, -100., 100.);
    const cosmobl::statistics::Prior A1_prior(cosmobl::glob::DistributionType::_UniformDistribution_, -100., 100.);
    const cosmobl::statistics::Prior A2_prior(cosmobl::glob::DistributionType::_UniformDistribution_, -1000., 1000.);
  
    model_twop.set_model_BAO_LinearPoint(alpha_prior, BB_prior, A0_prior, A1_prior, A2_prior);

  
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    const int chain_size = 2000;
    const int nwalkers = 200;
    const int seed = 4232;
    vector<double> starting_parameters = {0., 0., 0., 1., 0., 2., 0., 0., 0.};
    double radius = 1.e-3;

    const string chain_file = "Anderson_2013_CMASSDR11_monopole_prerecon_chains.dat";
  
    const double xmin = 40., xmax = 160.;
    model_twop.set_fit_range(xmin, xmax);
    
    model_twop.set_likelihood(cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Covariance_);
    
    model_twop.run_MCMC(chain_size, nwalkers, seed, starting_parameters, radius);
    
    model_twop.show_results(500, 10, 3213);
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

