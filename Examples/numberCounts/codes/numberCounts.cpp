// ==============================================================================================
// Example code: how to measure the number counts of a catalogue, i.e. the redshift distribution 
// ==============================================================================================

#include "LCDM.h"
#include "NumberCounts1D_Redshift.h"

using namespace std;

int main () {

  try {

    // -----------------------------------------------------------------
    // ---------------- use default cosmological parameters ------------
    // -----------------------------------------------------------------

    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");

    
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    const string file_catalogue = "../input/cat.dat";

    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

    
    // -------------------------------------------------------------------
    // ---------------- measure the redshift distribution ----------------
    // -------------------------------------------------------------------

    // binning parameters and output data

    const int nbin = 10;
    const string dir = "../output/";
    const string file = "redshift_distribution.dat";

    
    // measure the redshift distribution and compute Poisson errors

    cbl::measure::numbercounts::NumberCounts1D_Redshift NC {catalogue, nbin};

    NC.measure(cbl::measure::ErrorType::_Poisson_);

    
    // store the output data
  
    NC.write(dir, file);

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}

