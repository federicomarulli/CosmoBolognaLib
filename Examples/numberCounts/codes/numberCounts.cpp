// ==============================================================================================
// Example code: how to measure the number counts of a catalogue, i.e. the redshift distribution 
// ==============================================================================================

#include "NumberCounts.h"

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

    
    // -------------------------------------------------------------------
    // ---------------- measure the redshift distribution ----------------
    // -------------------------------------------------------------------

    // binning parameters and output data

    const int nbin = 10;
    const string dir = cosmobl::par::DirLoc+"../output/";
    const string file = "redshift_distribution.dat";

    
    // measure the redshift distribution and compute Poisson errors

    cosmobl::measure::numbercounts::NumberCounts NC {catalogue, cosmobl::catalogue::Var::_Redshift_};

    NC.measure(cosmobl::measure::ErrorType::_Poisson_, nbin);

    
    // store the output data
  
    NC.write(dir, file);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

