// ===================================================================================================================
// Example code: how to model the Cartesian 2D two-point correlation function to constrain the linear bias and sigma12
// ===================================================================================================================

#include "Modelling_TwoPointCorrelation2D_cartesian.h"

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

    const double N_R = 3.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 5.;   // minimum separation 
    const double rMax = 50.;  // maximum separation 
    const int nbins = 10;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre 

    const std::string dir = cbl::par::DirLoc+"../output/";
    const std::string file = "xi2D.dat";

  
    // measure the 2D Cartesian two-point correlation function and estimate Poissonian errors

    const auto TwoP = cbl::measure::twopt::TwoPointCorrelation::Create(cbl::measure::twopt::TwoPType::_2D_Cartesian_, catalogue, random_catalogue, cbl::BinType::_linear_, rMin, rMax, nbins, shift, cbl::BinType::_linear_, rMin, rMax, nbins, shift);

    TwoP->measure(cbl::measure::ErrorType::_Poisson_, dir, {dir});
    TwoP->write(dir, file);
  

    // ----------------------------------------------------------------------------------------------------------------------------------
    // ----------------- model the Cartesian 2D two-point correlation function and estimate the linear bias and sigma12 ----------------- 
    // ----------------------------------------------------------------------------------------------------------------------------------

    cbl::modelling::twopt::Modelling_TwoPointCorrelation2D_cartesian model_twop(TwoP); // object used for modelling

  
    // flat prior for f*sigma8
    const std::vector<double> fsigma8_limits = {0., 1.}; 
    const cbl::statistics::PriorDistribution fsigma8_prior {cbl::glob::DistributionType::_Uniform_, fsigma8_limits[0], fsigma8_limits[1], 413414}; 
  
    // flat prior for b*sigma8
    const std::vector<double> bsigma8_limits = {0.8*cosmology.sigma8(), 3.*cosmology.sigma8()}; 
    const cbl::statistics::PriorDistribution bsigma8_prior {cbl::glob::DistributionType::_Uniform_, bsigma8_limits[0], bsigma8_limits[1],63656}; 
  
    // flat prior for sigma12
    //const std::vector<double> sigma12_limits = {1., 1000.}; 
    const cbl::statistics::PriorDistribution sigma12_prior {cbl::glob::DistributionType::_Constant_, 0.};// sigma12_limits[0], sigma12_limits[1], 5411}; 

    // mean redshift of the sample
    const double redshift = 1.;
    
    // set the data used to construct de model
    model_twop.set_data_model(cosmology, redshift);
    
    // set the model for the redshift-space 2D two-point correlation 
    model_twop.set_model_dispersionModel(fsigma8_prior, bsigma8_prior, sigma12_prior); 

    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    const double min = 1., max = 40.;
  
    model_twop.set_fit_range(min, max, min, max);
    
    const int chain_size = 100;
    const int nwalkers = 10;
    const int seed = 666;

    std::vector<double> starting_parameters = {0.5, 1.5, 100.};
    
    const std::string chain_file = "chain_cartesian_bias_sigma12.dat";
    
    model_twop.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);
    
    model_twop.sample_posterior(chain_size, nwalkers, seed);

    const int burn_in = 0;
    const int thin = 1; 
    model_twop.show_results(burn_in, thin);
    model_twop.write_chain(dir, chain_file, burn_in, thin);
    model_twop.write_model_from_chains(dir, "model", {}, {}, burn_in, thin);

  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

