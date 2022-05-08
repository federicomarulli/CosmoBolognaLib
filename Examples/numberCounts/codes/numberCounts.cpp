// ==============================================================================================
// Example code: how to measure the number counts of a catalogue, i.e. the redshift distribution 
// ==============================================================================================

#include "NumberCounts1D_Redshift.h"


int main () {

  try {

    // -----------------------------------------------------------------
    // ---------------- use default cosmological parameters ------------
    // -----------------------------------------------------------------

    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

    
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    const std::string file_catalogue = "../input/cat.dat";

    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

    
    // -------------------------------------------------------------------
    // ---------------- measure the redshift distribution ----------------
    // -------------------------------------------------------------------

    // binning parameters and output data

    const int nbin = 10;
    const std::string dir = "../output/";
    const std::string file = "redshift_distribution.dat";

    
    // measure the redshift distribution and compute Poisson errors

    cbl::measure::numbercounts::NumberCounts1D_Redshift NC {catalogue, nbin};

    NC.measure(cbl::measure::ErrorType::_Poisson_);

    
    // store the output data
  
    NC.write(dir, file);

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

