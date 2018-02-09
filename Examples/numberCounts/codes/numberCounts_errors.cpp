// =========================================================================================================================
// Example code: how to measure the number counts of a catalogue, i.e. the redshift distribution, computing Possonian errors
// =========================================================================================================================

#include "NumberCounts.h"
#include "GlobalFunc.h"

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

    cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

    
    // construct the sub-regions used for jackknife and bootstrap

    cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << endl;
    const double cellSize = 2;
    cosmobl::set_ObjectRegion_RaDec(catalogue, cellSize);

    
    // -------------------------------------------------------------------
    // ---------------- measure the redshift distribution ----------------
    // -------------------------------------------------------------------

    // binning parameters and output data

    const int nbin = 10;
    const string dir = cosmobl::par::DirLoc+"../output/";

    
    // measure the redshift distribution anc compute Poisson errors

    cosmobl::measure::numbercounts::NumberCounts NC {catalogue, cosmobl::catalogue::Var::_Redshift_};

    
    // measure the number counts and compute Poissonian errors 

    NC.measure(cosmobl::measure::ErrorType::_Poisson_, nbin);
    NC.write(dir, "redshift_distribution_Poisson.dat");

    
    // measure the number counts and compute Jackknife errors 

    NC.measure(cosmobl::measure::ErrorType::_Jackknife_, nbin);
    NC.write(dir, "redshift_distribution_Jackknife.dat");
    NC.write_covariance(dir, "redshift_distribution_Jackknife_covariance.dat");

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

