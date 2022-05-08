// ========================================================================
// Example code: how to measure the multipoles of the two-point correlation
// function, using the "direct" method and the "integrated" method
// ========================================================================

#include "TwoPointCorrelation_multipoles_direct.h"
#include "TwoPointCorrelation_multipoles_integrated.h"

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

    const double N_R = 10.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the multipoles of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // output directory
    const std::string dir = cbl::par::DirLoc+"../output/";

    // binning parameters and output data

    const double rMin = 10.;   // minimum separation 
    const double rMax = 30.;  // maximum separation 
    const int nbins = 3;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre
  
  
    // measure the multipoles using direct estimator

    cbl::measure::twopt::TwoPointCorrelation_multipoles_direct TwoP_direct {catalogue, random_catalogue, cbl::BinType::_logarithmic_, rMin, rMax, nbins, shift};
  
    TwoP_direct.measure(cbl::measure::ErrorType::_Poisson_, dir);
    

    // store the output data
  
    std::string file = "xil_direct.dat";

    TwoP_direct.write(dir, file);
    
    
    // measure the multipoles with integrated estimator 

    const double muMin = 0.;
    const double muMax = 1.;
    const int nbinsMu = 20;

    cbl::measure::twopt::TwoPointCorrelation_multipoles_integrated TwoP_integrated {catalogue, random_catalogue, cbl::BinType::_logarithmic_, rMin, rMax, nbins, shift, muMin, muMax, nbinsMu, shift};
    
    TwoP_integrated.measure(cbl::measure::ErrorType::_Poisson_, dir);
    
    
    // store the output data
  
    file = "xil_integrated.dat";

    TwoP_integrated.write(dir, file);
  
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

