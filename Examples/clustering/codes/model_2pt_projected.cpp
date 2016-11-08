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

    const double N_R = 1.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters

    const double rMin = 1.;            // minimum separation 
    const double rMax = 50.;           // maximum separation 
    const int nbins = 20;              // number of bins
    const double shift = 0.5;          // spatial shift used to set the bin centre 
    const double piMax_integral = 30.; // upper limit of the integral

  
    // measure the projected two-point correlation function and estimate Poissonian errors

    const auto TwoP = cosmobl::twopt::TwoPointCorrelation::Create(cosmobl::twopt::TwoPType::_1D_projected_, catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift, rMin, rMax, nbins, shift, piMax_integral);

    TwoP->measure(cosmobl::twopt::ErrorType::_Poisson_, dir_output, {dir_output});
    TwoP->write(dir_output, "wp.dat");

  
    // -----------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------- model the projected two-point correlation function and estimate the linear bias ------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------------------

    cosmobl::modelling::Modelling_TwoPointCorrelation_projected model_twop(TwoP); // object used for modelling

  
    const double redshift = 1.; // redshift of the sample
  
    const vector<double> rad_xiDM = cosmobl::logarithmic_bin_vector(200, 0.1, 50.); // scales used
  
    model_twop.set_parameters_xiDM(rad_xiDM, cosmology, redshift, "CAMB"); // set the model parameters
  

    const cosmobl::statistics::Prior bias_prior {cosmobl::statistics::PriorType::_UniformPrior_, 0.8, 3.}; // flat prior for the bias

    model_twop.set_model_linearBias(bias_prior); // set the prior on the bias
  

    // sample the likelihood and store the chains
  
    const int nChains = 100;     // number of chains
    const int chain_size = 1000; // size of the chains
    const double rpMin = 5.;     // minimum scale used for fitting
    const double rpMax = 20.;    // maximum scale used for fitting
    const string chain_file = "chain.dat";
  
    model_twop.sample_likelihood(rpMin, rpMax, cosmobl::statistics::_GaussianLikelihood_Error_, nChains, chain_size, 32113, dir_output, chain_file); 
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

