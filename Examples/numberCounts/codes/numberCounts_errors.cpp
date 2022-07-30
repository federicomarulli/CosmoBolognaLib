// =========================================================================================================================
// Example code: how to measure the number counts of a catalogue, i.e. the redshift distribution, computing Possonian errors
// =========================================================================================================================

#include "NumberCounts1D_Redshift.h"
#include "GlobalFunc.h"

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

    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

    
    // construct the sub-regions used for jackknife and bootstrap

    std::cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << std::endl;
    const int nCells_Ra = 3;
    const int nCells_Dec = 3;
    cbl::set_ObjectRegion_RaDec(catalogue, nCells_Ra, nCells_Dec);

    
    // -------------------------------------------------------------------
    // ---------------- measure the redshift distribution ----------------
    // -------------------------------------------------------------------

    // binning parameters and output data

    const int nbin = 10;
    const std::string dir = "../output/";

    
    // measure the redshift distribution anc compute Poisson errors

    cbl::measure::numbercounts::NumberCounts1D_Redshift NC {catalogue, nbin};

    
    // measure the number counts and compute Poissonian errors 

    NC.measure(cbl::measure::ErrorType::_Poisson_);
    NC.write(dir, "redshift_distribution_Poisson.dat");

    
    // measure the number counts and compute Jackknife errors 

    NC.measure(cbl::measure::ErrorType::_Jackknife_);
    NC.write(dir, "redshift_distribution_Jackknife.dat");
    NC.write_covariance(dir, "redshift_distribution_Jackknife_covariance.dat");

    // measure the number counts and compute Jackknife errors 

    NC.measure(cbl::measure::ErrorType::_Bootstrap_, dir+"Bootstrap/", 1000);
    NC.write(dir, "redshift_distribution_Bootstrap.dat");
    NC.write_covariance(dir, "redshift_distribution_Bootstrap_covariance.dat");
    
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

