// ================================================================================================================
// Example code: how to model the monopole of the two-point correlation function to constrain f*sigma8 and b*sigma8
// ================================================================================================================

#include "Modelling_TwoPointCorrelation_monopole.h"

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

    auto TwoP = cosmobl::twopt::TwoPointCorrelation::Create(cosmobl::twopt::TwoPType::_1D_monopole_, catalogue, random_catalogue, cosmobl::_linear_, rMin, rMax, nbins, shift);

    TwoP->measure(cosmobl::twopt::_Poisson_, dir);
    TwoP->write(dir, file);
  

    // ----------------------------------------------------------------------------------------------------------------------------
    // ----------------- model the monopole two-point correlation function and estimate the f*sigma8 and b*sigma8 -----------------
    // ----------------------------------------------------------------------------------------------------------------------------

    cosmobl::modelling::Modelling_TwoPointCorrelation_monopole model_twop(TwoP); // object used for modelling
  
  
    const vector<double> model_scales = cosmobl::linear_bin_vector(100, 0.1, 50.); // scales used in the fit
    const double redshift = 1.; // mean redshift of the sample
  
    model_twop.set_parameters_xiDM(model_scales, cosmology, redshift); // set the model parameters


    const cosmobl::statistics::Prior fsigma8_prior {cosmobl::statistics::PriorType::_UniformPrior_, 0., 2.}; // flat prior for the f*sigma8
    const cosmobl::statistics::Prior bsigma8_prior {cosmobl::statistics::PriorType::_UniformPrior_, 1.2, 1.3}; // flat prior for the b*sigma8
  
    model_twop.set_model_Kaiser(fsigma8_prior, bsigma8_prior); // set the priors 

  
    // sample the likelihood and store the chains
  
    const int nChains = 100; // number of chains
    const int chain_size = 1000; // size of the chains
    const string chain_file = "chain_monopole_fsigma8_bsigma8.dat";

    model_twop.sample_likelihood(10., 40., cosmobl::statistics::_GaussianLikelihood_Error_, nChains, chain_size, 32113, dir, chain_file); 

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

