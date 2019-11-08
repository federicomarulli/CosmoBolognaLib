// =================================================================================================
// Example code: how to measure the angle-averaged two-point correlation function, i.e. the monopole 
// =================================================================================================

#include "TwoPointCorrelation1D_monopole.h"

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

    
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    const std::string file_catalogue = cbl::par::DirLoc+"../input/cat.dat";

    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};
    
  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};
    
  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the monopole of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 5.;   // minimum separation 
    const double rMax = 20.;  // maximum separation 
    const int nbins = 5;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre
  
    const std::string dir = cbl::par::DirLoc+"../output/";
    const std::string file = "xi.dat";

  
    // measure the monopole and compute Poisson errors 

    cbl::measure::twopt::TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, cbl::BinType::_logarithmic_, rMin, rMax, nbins, shift};

    TwoP.measure(cbl::measure::ErrorType::_Poisson_, dir);
    
  
    // store the output data
  
    TwoP.write(dir, file);
  
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

