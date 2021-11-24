// =========================================================================
// Example code: how to measure the projected two-point correlation function
// =========================================================================

#include "TwoPointCorrelation1D.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


int main () {

  try {
  
    // -----------------------------------------------------------------
    // ---------------- use default cosmological parameters ------------
    // -----------------------------------------------------------------
  
    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

  
    // ---------------------------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) --------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------

    const std::string file_catalogue = cbl::par::DirLoc+"../input/cat.dat";

    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 3.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};


    // --------------------------------------------------------------------------------------
    // ---------------- measure the projected two-point correlation function ----------------
    // --------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rpMin = 10.;     // minimum separation in the first dimension
    const double rpMax = 30.;    // maximum separation in the first dimension 
    const int nbins_D1 = 3;     // number of bins in the first dimension
    const double shift_D1 = 0.5; // spatial shift used to set the bin centre in the first dimension
    const double piMin = 0.;     // minimum separation in the second dimension
    const double piMax = 50.;    // maximum separation in the second dimension 
    const int nbins_D2 = 10;     // number of bins in the second dimension
    const double shift_D2 = 0.5; // spatial shift used to set the bin centre in the second dimension
    
    const double piMax_integral = 50.; // upper limit of the integral
    
    const std::string dir = cbl::par::DirLoc+"../output/";
    const std::string file = "xi_projected.dat";
  

    // measure the projected two-point correlation function
  
    const auto TwoP = cbl::measure::twopt::TwoPointCorrelation::Create(cbl::measure::twopt::TwoPType::_projected_, catalogue, random_catalogue, cbl::BinType::_logarithmic_, rpMin, rpMax, nbins_D1, shift_D1, piMin, piMax, nbins_D2, shift_D2, piMax_integral);

    TwoP->measure(cbl::measure::ErrorType::_Poisson_, dir);

    TwoP->write(dir, file);
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

