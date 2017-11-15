// ====================================================================================================
// Example code: how to model the projected two-point correlation function to constrain the linear bias
// ====================================================================================================

#include "Modelling_TwoPointCorrelation_projected.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------------------------
    // ---------------- use default cosmological parameters and set sigma8 ------------
    // --------------------------------------------------------------------------------
  
    cosmobl::cosmology::Cosmology cosmology;
    cosmology.set_sigma8(0.8);

  
    // --------------------------------------------------------------------
    // ---------------- Input/Output files and directories ----------------
    // --------------------------------------------------------------------
  
    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";

    const string dir_output = cosmobl::par::DirLoc+"../output/";
    const string dir_pairs = dir_output+"pairs/";
    const string dir_covariance = dir_output+"covariance/";
  
    const string MK = "mkdir -p "+dir_output+" "+dir_pairs+" "+dir_covariance; if (system(MK.c_str())) {}

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------

    cout << "I'm reading the input catalogue..." << endl;

    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------

    const double N_R = 2.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters

    const double rpMin = 5.;     // minimum separation in the first dimension
    const double rpMax = 50.;    // maximum separation in the first dimension 
    const int nbins_D1 = 10;     // number of bins in the first dimension
    const double shift_D1 = 0.5; // spatial shift used to set the bin centre in the first dimension
    const double piMin = 0.;     // minimum separation in the second dimension
    const double piMax = 50.;    // maximum separation in the second dimension 
    const int nbins_D2 = 20;     // number of bins in the second dimension
    const double shift_D2 = 0.5; // spatial shift used to set the bin centre in the second dimension
    const double piMax_integral = 30.; // upper limit of the integral to compute the projected correlation function

  
    // measure the projected two-point correlation function and store the result

    const auto TwoP = cosmobl::measure::twopt::TwoPointCorrelation::Create(cosmobl::measure::twopt::TwoPType::_1D_projected_, catalogue, random_catalogue, cosmobl::_logarithmic_, rpMin, rpMax, nbins_D1, shift_D1, piMin, piMax, nbins_D2, shift_D2, piMax_integral);

    TwoP->measure(cosmobl::measure::ErrorType::_Poisson_, dir_output, {dir_output});
    TwoP->write(dir_output, "wp.dat");

  
    // -----------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------- model the projected two-point correlation function and estimate the linear bias ------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------------------

    // object used for modelling
    cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_projected model_twop(TwoP);

    // redshift of the snapshot of the simulation box, at which the model is computed
    const double redshift = 1.;
    
    // set the data used to construct the model (all the unspecified parameters will be set to default values) 
    model_twop.set_data_model(cosmology, redshift); 
  
    // set the model used to constrain the linear bias, imposing a flat prior on b*sigma8
    model_twop.set_model_linearBias({cosmobl::glob::DistributionType::_UniformDistribution_, 0.8, 3.});
    
    // set the scale limits use for fitting
    model_twop.set_fit_range(5., 20.);
    
    // set the likelihood
    model_twop.set_likelihood(cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Covariance_);

    // run the MCMC method to sample the posterior
    const int chain_size = 1000; // size of the chains
    const int nwalkers = 100; // number of chains
    model_twop.run_MCMC(chain_size, nwalkers, 32113); 

    // write the results on the standard output
    const int burn_in = 0.1*chain_size; // number of step to be excluded
    const int thin = 10;                // take 1 step every 10
    model_twop.show_results(burn_in, thin, 3213);
    
    // store the results on file
    const string filename = "model_projected";
    model_twop.write_results(dir_output, filename, burn_in, thin, 3213);
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

