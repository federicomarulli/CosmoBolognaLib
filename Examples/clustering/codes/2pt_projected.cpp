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

    const double N_R = 3.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rpMin = 5.;     // minimum separation in the first dimension
    const double rpMax = 50.;    // maximum separation in the first dimension 
    const int nbins_D1 = 10;     // number of bins in the first dimension
    const double shift_D1 = 0.5; // spatial shift used to set the bin centre in the first dimension
    const double piMin = 0.;     // minimum separation in the second dimension
    const double piMax = 50.;    // maximum separation in the second dimension 
    const int nbins_D2 = 20;     // number of bins in the second dimension
    const double shift_D2 = 0.5; // spatial shift used to set the bin centre in the second dimension
    
    const double piMax_integral = 30.; // upper limit of the integral
    
    const string dir = cosmobl::par::DirLoc+"../output/";
    const string file = "xi_projected.dat";
  

    // measure the projected two-point correlation function
  
    const auto TwoP = cosmobl::twopt::TwoPointCorrelation::Create(cosmobl::twopt::TwoPType::_1D_projected_, catalogue, random_catalogue, cosmobl::_logarithmic_, rpMin, rpMax, nbins_D1, shift_D1, piMin, piMax, nbins_D2, shift_D2, piMax_integral);

    TwoP->measure(cosmobl::twopt::_Poisson_, dir);

    TwoP->write(dir, file);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

