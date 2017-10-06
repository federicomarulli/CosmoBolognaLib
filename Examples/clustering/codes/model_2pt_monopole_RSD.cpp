// ===============================================================================================
// Example code: how to model the monopole of the two-point correlation function in redshift space
// ===============================================================================================

#include "Modelling_TwoPointCorrelation1D_monopole.h"

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

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------

    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";
  
    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------------
    // ---------------- measure the monopole of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 0.1;  // minimum separation 
    const double rMax = 50.;  // maximum separation 
    const int nbins = 25;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre 

    const string dir = cosmobl::par::DirLoc+"../output/";
    const string file = "xi.dat";

  
    // measure the 2D Cartesian two-point correlation function and estimate Poissonian errors

    auto TwoP = cosmobl::measure::twopt::TwoPointCorrelation::Create(cosmobl::measure::twopt::TwoPType::_1D_monopole_, catalogue, random_catalogue, cosmobl::_linear_, rMin, rMax, nbins, shift);

    TwoP->measure(cosmobl::measure::ErrorType::_Poisson_, dir);
    TwoP->write(dir, file);


    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------- model the monopole two-point correlation function and estimate the f*sigma8 and b*sigma8 -----------------
    // ----------------------- (f*sigma8 and b*sigma8 will be degenerate, if the prior is uniform for both) -----------------------
    // ----------------------------------------------------------------------------------------------------------------------------

    // object used for modelling
    cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole model_twop(TwoP); 

    // mean redshift of the sample
    const double redshift = 1.; 

    // set the data used to construct the model
    model_twop.set_data_model(cosmology, redshift);
  
    // set the priors and the model
    const cosmobl::statistics::Prior fsigma8_prior {cosmobl::glob::DistributionType::_UniformDistribution_, 0., 2.}; // flat prior for the f*sigma8
    const cosmobl::statistics::Prior bsigma8_prior {cosmobl::glob::DistributionType::_UniformDistribution_, 0., 2.}; // flat prior for the b*sigma8
    model_twop.set_model_Kaiser(fsigma8_prior, bsigma8_prior);

  
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    // minimum and maxium scales used in the fit
    const double xmin = 10.;
    const double xmax = 40.;
    model_twop.set_fit_range(xmin, xmax);

    const int chain_size = 10000; // the size the chain lenght
    const int nwalkers = 200;     // the number of parallel walkers in the MCMC chains
    const int seed = 4232;        // the base seed for initialization

    // set the likelihood type
    model_twop.set_likelihood(cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Error_);

    // sample the likelihood
    model_twop.sample_likelihood(chain_size, nwalkers, seed);

    const int burn_in = 100; // discard the first 100 chain steps 
    const int thin = 10;     // take 1 step every 10

    // write the results on screen 
    model_twop.show_results(burn_in, thin, seed);

    // store the results in file
    model_twop.write_results(dir, "model_RSD", burn_in, thin, seed);
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

