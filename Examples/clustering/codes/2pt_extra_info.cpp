// ==============================================================================================
// Example code: how to measure the angle-averaged two-point correlation function with extra info
// ==============================================================================================

#include "TwoPointCorrelation1D_monopole.h"
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

    const double N_R = 1.; // random/data ratio
  
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the monopole of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters and output data

    const double rMin = 1.;   // minimum separation 
    const double rMax = 50.;  // maximum separation 
    const int nbins = 20;     // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre

    string dir, file="xi.dat", cmd;
  
    // measure the monopole and compute Poisson errors 

    cosmobl::measure::twopt::TwoPointCorrelation1D_monopole TwoP0 {catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift, cosmobl::CoordUnits::_radians_, nullptr, true};
  
    dir = cosmobl::par::DirLoc+"../output/xi0/";
    cmd = "mkdir -p "+dir;
    if(system(cmd.c_str())) {}

    TwoP0.measure(cosmobl::measure::ErrorType::_Poisson_, dir);
  
    // store the output data
  
    TwoP0.write(dir, file);

    // measure the multipoles with direct method and compute Poisson errors 

    cosmobl::measure::twopt::TwoPointCorrelation_multipoles_direct TwoPld {catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift, cosmobl::CoordUnits::_radians_, nullptr, true};

    dir = cosmobl::par::DirLoc+"../output/xil_d/";
    cmd = "mkdir -p "+dir;
    if(system(cmd.c_str())) {}

    TwoPld.measure(cosmobl::measure::ErrorType::_Poisson_, dir);
    // store the output data

    TwoPld.write(dir, file);
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

