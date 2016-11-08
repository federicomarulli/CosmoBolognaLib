// =======================================================================
// Example code: how to measure the angular two-point correlation function 
// =======================================================================

#include "TwoPointCorrelation1D_angular.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // -----------------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec) --------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat2d.dat";
  
    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
    const int zbins = 10;  // number of bins used to compute the redshift distribution 
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_square_, catalogue, N_R, zbins};

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the angular of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters and output data

    const double angMin = 0.01;                // minimum angular separation 
    const double angMax = 1.;                  // maximum angular separation 
    const int nbins = 20;                      // number of bins
    const double shift = 0.5;                  // shift used to set the bin centre 
    const cosmobl::CoordUnits angularUnits = cosmobl::_degrees_; // angular units

    const string dir = cosmobl::par::DirLoc+"../output/";
    const string file = "xi_angular.dat";

  
    // measure the angular two-point correlation function
  
    cosmobl::twopt::TwoPointCorrelation1D_angular TwoP {catalogue, random_catalogue, cosmobl::_linear_, angMin, angMax, nbins, shift, angularUnits};
  
    TwoP.measure(cosmobl::twopt::_Poisson_, dir);

    TwoP.write(dir, file);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

