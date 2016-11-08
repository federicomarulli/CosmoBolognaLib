// =========================================================================
// Example code: how to measure the projected two-point correlation function
// =========================================================================

#include "TwoPointCorrelation1D.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // -----------------------------------------------------------------
    // ---------------- use default cosmological parameters ------------
    // -----------------------------------------------------------------
  
    const cosmobl::cosmology::Cosmology cosmology;

  
    // ---------------------------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) --------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------

    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";

    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 1.;           // minimum separation 
    const double rMax = 50.;          // maximum separation 
    const int nbins = 20;             // number of bins
    const double shift = 0.5;         // spatial shift used to set the bin centre 
    const double piMax_integral = 30; // upper limit of the integral

    const string dir = cosmobl::par::DirLoc+"../output/";
    const string file = "xi_projected.dat";
  

    // measure the projected two-point correlation function
  
    const auto TwoP = cosmobl::twopt::TwoPointCorrelation::Create(cosmobl::twopt::TwoPType::_1D_projected_, catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift, rMin, rMax, nbins, shift, piMax_integral);

    TwoP->measure(cosmobl::twopt::_Poisson_, dir);

    TwoP->write(dir, file);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

