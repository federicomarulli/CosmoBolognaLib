// ========================================================================
// Example code: how to measure the multipoles of the two-point correlation
// function, using the "direct" method and the "integrated" method
// ========================================================================

#include "TwoPointCorrelation_multipoles_direct.h"
#include "TwoPointCorrelation_multipoles_integrated.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // -----------------------------------------------------------------
    // ---------------- use default cosmological parameters ------------
    // -----------------------------------------------------------------

    const cosmobl::cosmology::Cosmology cosmology {cosmobl::cosmology::_Planck15_};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";

    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 10.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the multipoles of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // output directory
    const string dir = cosmobl::par::DirLoc+"../output/";

    // binning parameters and output data

    const double rMin = 1.;   // minimum separation 
    const double rMax = 50.;  // maximum separation 
    const int nbins = 20;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre
  
  
    // measure the multipoles using direct estimator

    cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct TwoP_direct {catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift};
  
    TwoP_direct.measure(cosmobl::measure::ErrorType::_Poisson_, dir);
    

    // store the output data
  
    string file = "xil_direct.dat";

    TwoP_direct.write(dir, file);

    
    // measure the multipoles with integrated estimator 

    const double muMin = 0.;
    const double muMax = 1.;
    const int nbinsMu = 20;

    cosmobl::measure::twopt::TwoPointCorrelation_multipoles_integrated TwoP_integrated {catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift, muMin, muMax, nbinsMu, shift};
  
    TwoP_integrated.measure(cosmobl::measure::ErrorType::_Poisson_, dir);

    
    // store the output data
  
    file = "xil_integrated.dat";

    TwoP_integrated.write(dir, file);
  
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

