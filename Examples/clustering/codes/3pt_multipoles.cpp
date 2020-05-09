// ========================================================================================
// Example code: how to measure the connected and reduced three-point correlation functions
// ========================================================================================

#include "ThreePointCorrelation_comoving_multipoles.h"
#include "GlobalFunc.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------

    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    std::string file_catalogue = cbl::par::DirLoc+"../input/cat.dat";

    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 1000.; // random/data ratio
   
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};

    // -------------------------------------------------------------------------------
    // ---------------- measure the three-point correlation functions ----------------
    // -------------------------------------------------------------------------------

    // binning parameters

    const double rMin = 5.;  
    const double rMax = 11.;
    const double binSize = 2;
    const int nOrders = 2;
  
    // output data
  
    const std::string dir_output = cbl::par::DirLoc+"../output/";
    std::string dir_triplets = cbl::par::DirLoc+"../output/triplets/Single/";
    std::string file_output = "zeta_multipoles_single.dat";
  
    // measure the connected three-point correlation functions legendre coefficients and write the output

    const auto ThreeP_Single = cbl::measure::threept::ThreePointCorrelation_comoving_multipoles::Create(catalogue, random_catalogue,
        9., 11., 9., 11., nOrders);

    ThreeP_Single->measure(cbl::measure::ErrorType::_None_, dir_triplets);
  
    ThreeP_Single->write(dir_output, file_output);

    // measure the connected three-point correlation functions legendre coefficients and write the output

    dir_triplets = cbl::par::DirLoc+"../output/triplets/All/";
    file_output = "zeta_multipoles_all.dat";

    const auto ThreeP_All = cbl::measure::threept::ThreePointCorrelation_comoving_multipoles::Create(catalogue, random_catalogue,
        rMin, rMax, binSize, nOrders);

    ThreeP_All->measure(cbl::measure::ErrorType::_None_, dir_triplets);
  
    ThreeP_All->write(dir_output, file_output);

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

