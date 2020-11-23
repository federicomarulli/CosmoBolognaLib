// =======================================================================
// Example code: how to measure the angular two-point correlation function 
// =======================================================================

#include "TwoPointCorrelation1D_angular.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


int main () {

  try {
  
    // -----------------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec) --------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    const std::string file_catalogue = cbl::par::DirLoc+"../input/cat2d.dat";
  
    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}};

  
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_square_, catalogue, N_R};
    
  
    // ------------------------------------------------------------------------------------
    // ---------------- measure the angular two-point correlation function ----------------
    // ------------------------------------------------------------------------------------

    // binning parameters and output data

    const double angMin = 0.01;                                              // minimum angular separation 
    const double angMax = 1.;                                                // maximum angular separation 
    const int nbins = 20;                                                    // number of bins
    const double shift = 0.5;                                                // shift used to set the bin centre 
    const cbl::CoordinateUnits angularUnits = cbl::CoordinateUnits::_degrees_; // angular units

    const std::string dir = cbl::par::DirLoc+"../output/";
    const std::string file = "xi_angular.dat";

  
    // measure the angular two-point correlation function and store the results
  
    cbl::measure::twopt::TwoPointCorrelation1D_angular TwoP {catalogue, random_catalogue, cbl::BinType::_linear_, angMin, angMax, nbins, shift, angularUnits};
  
    TwoP.measure(cbl::measure::ErrorType::_Poisson_, dir);

    TwoP.write(dir, file);

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

