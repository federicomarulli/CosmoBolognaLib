// ========================================================================================================================
// Example code: how to how to model the baryon acoustic oscillations in the monopole of the two-point correlation function
// ========================================================================================================================

#include "Data1D.h"
#include "Modelling_TwoPointCorrelation_monopole.h"

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
    cosmology.set_parameter({cosmobl::cosmology::CosmoPar::_Omega_matter_LCDM_, cosmobl::cosmology::CosmoPar::_Omega_baryon_, cosmobl::cosmology::CosmoPar::_hh_, cosmobl::cosmology::CosmoPar::_n_spec_, cosmobl::cosmology::CosmoPar::_sigma8_}, {Omega_M, Omega_b, hh, n_s, sigma8});

  
    // ----------------------------------------------
    // ------------- reading the dataset ------------
    // ----------------------------------------------

    const string dir_input = cosmobl::par::DirLoc+"../input/";
    const string dir_output = cosmobl::par::DirLoc+"../output/";
    const string dir_chains = dir_output+"chains/";

    const string file_xi = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon.dat";
    const string file_cov = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon_cov.dat";

    const string MK = "mkdir -p "+dir_output+" "+dir_chains; if (system(MK.c_str())) {}

    const auto twop_dataset = make_shared<cosmobl::data::Data1D>(cosmobl::data::Data1D(file_xi, true)); 
    twop_dataset->set_covariance(file_cov);

  
    // ------------------------------------------------------
    // ------------- set the modelling structure ------------
    // ------------------------------------------------------
  
    cosmobl::modelling::Modelling_TwoPointCorrelation_monopole model_twop(twop_dataset);
  

    // --------------------------------------------------------------------------
    // ------------- set the underlying dark matter clustering model ------------
    // --------------------------------------------------------------------------

    const double redshift = 0.57;
    const double sigmaNL = 8.;
  
    const vector<double> model_scales = cosmobl::linear_bin_vector(100, 1., 300.);
  
    model_twop.set_parameters_xiDM(model_scales, cosmology, redshift, "CAMB", sigmaNL);

  
    // -------------------------------------------------------
    // ------------- set the priors and the model ------------
    // -------------------------------------------------------

    const cosmobl::statistics::Prior alpha_prior(cosmobl::statistics::PriorType::_UniformPrior_, 0.7, 1.3); 
    const cosmobl::statistics::Prior bsigma8_prior(cosmobl::statistics::PriorType::_UniformPrior_, 1., 15.);
    const cosmobl::statistics::Prior A0_prior(cosmobl::statistics::PriorType::_UniformPrior_, -100., 100.);
    const cosmobl::statistics::Prior A1_prior(cosmobl::statistics::PriorType::_UniformPrior_, -100., 100.);
    const cosmobl::statistics::Prior A2_prior(cosmobl::statistics::PriorType::_UniformPrior_, -1000., 1000.);
  
    model_twop.set_model_BAO(alpha_prior, bsigma8_prior, A0_prior, A1_prior, A2_prior);

  
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    const int nChains = 100;
    const int chain_size = 2000;
    const double rmin = 40.;
    const double rmax = 200.;

    const string chain_file = "Anderson_2013_CMASSDR11_monopole_prerecon_chains.dat";
  
    model_twop.sample_likelihood(rmin, rmax, cosmobl::statistics::_GaussianLikelihood_Covariance_, nChains, chain_size, 43123, dir_chains, chain_file);
  
    const int npt = 200;
    const vector<double> model_scales_output = cosmobl::linear_bin_vector(npt, rmin, rmax);
    const string output_file = "Anderson_2013_CMASSDR11_monopole_prerecon_bestfit_model.dat";
    model_twop.write_model(model_scales_output, dir_chains, output_file);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

