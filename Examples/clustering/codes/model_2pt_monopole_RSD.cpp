// ===============================================================================================
// Example code: how to model the monopole of the two-point correlation function in redshift space
// ===============================================================================================

#include "Modelling_TwoPointCorrelation1D_monopole.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------------------------
    // ---------------- use default cosmological parameters and set sigma8 ------------
    // --------------------------------------------------------------------------------
  
    cbl::cosmology::Cosmology cosmology;
    cosmology.set_sigma8(0.8);

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------

    const std::string file_catalogue = cbl::par::DirLoc+"../input/cat.dat";
  
    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

  
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------------
    // ---------------- measure the monopole of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 10.;  // minimum separation 
    const double rMax = 30.;  // maximum separation 
    const int nbins = 10;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre 

    const std::string dir = cbl::par::DirLoc+"../output/";
    const std::string file = "xi.dat";

  
    // measure the monopole of the two-point correlation function and estimate Poissonian errors

    auto TwoP = cbl::measure::twopt::TwoPointCorrelation::Create(cbl::measure::twopt::TwoPType::_monopole_, catalogue, random_catalogue, cbl::BinType::_linear_, rMin, rMax, nbins, shift);

    TwoP->measure(cbl::measure::ErrorType::_Poisson_, dir);
    TwoP->write(dir, file);


    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------- model the monopole two-point correlation function and estimate the f*sigma8 and b*sigma8 -----------------
    // ----------------------- (f*sigma8 and b*sigma8 will be degenerate, if the prior is uniform for both) -----------------------
    // ----------------------------------------------------------------------------------------------------------------------------

    // object used for modelling
    cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole model_twop(TwoP); 

    // mean redshift of the sample
    const double redshift = 1.; 

    // set the data used to construct the model
    model_twop.set_data_model(cosmology, redshift);
  
    // set the priors and the model
    const cbl::statistics::PriorDistribution fsigma8_prior {cbl::glob::DistributionType::_Uniform_, 0., 2.}; // flat prior for the f*sigma8
    const cbl::statistics::PriorDistribution bsigma8_prior {cbl::glob::DistributionType::_Uniform_, 0., 2.}; // flat prior for the b*sigma8
    model_twop.set_model_Kaiser(fsigma8_prior, bsigma8_prior);

    
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    // minimum and maxium scales used in the fit
    const double xmin = 10.;
    const double xmax = 40.;
    model_twop.set_fit_range(xmin, xmax);

    // set the likelihood type
    model_twop.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);

    // maximise the posterior
    model_twop.maximize_posterior({1., 1.}, 10000, 1.e-5);
    
    // run the MCMC method to sample the posterior
    const int chain_size = 1000; // the size the chain lenght
    const int nwalkers = 10;     // the number of parallel walkers in the MCMC chains
    const int seed = 666;        // the base seed for initialization
    model_twop.sample_posterior(chain_size, nwalkers, seed);

    // write the results on screen
    const int burn_in = 100; // discard the first 100 chain steps 
    const int thin = 10;     // take 1 step every 10
    model_twop.show_results(burn_in, thin);

    // store the results in file
    model_twop.write_results(dir, "model_RSD", burn_in, thin);

    // store the best-fit model
    model_twop.write_model_from_chains(dir, "bestfit_model.dat", cbl::logarithmic_bin_vector(100, 0.1, 100.), burn_in, thin);
    
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

