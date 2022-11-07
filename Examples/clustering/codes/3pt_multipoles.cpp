// ============================================================================
// Example code: how to measure the three-point correlation function multipoles
// ============================================================================

#include "ThreePointCorrelation_comoving_multipoles.h"
#include "GlobalFunc.h"

int main () {

  try {
  
    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------

    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    std::string file_catalogue = "../input/cat.dat";

    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 2.; // random/data ratio
   
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};
  
    // -------------------------------------------------------------------------------
    // ---------------- measure the three-point correlation functions ----------------
    // -------------------------------------------------------------------------------

    // binning parameters

    const double rMin = 5.;  
    const double rMax = 10.;
    const double binSize = 2;
    const int nOrders = 11;
  
    // output data
  
    const std::string dir_output = "../output/";
    const std::string dir_triplets = "../output/triplets/";
    const std::string file_output = "zeta_multipoles.dat";
  
    // measure the connected three-point correlation functions legendre coefficients and write the output

    const auto ThreeP = cbl::measure::threept::ThreePointCorrelation::Create(catalogue, random_catalogue, rMin, rMax, binSize, nOrders);

    ThreeP->measure(cbl::measure::ErrorType::_None_, dir_triplets, {}, 0, true, false, 0.1);
  
    ThreeP->write(dir_output, file_output);

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

