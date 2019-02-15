// ========================================================================================
// Example code: how to measure the connected and reduced three-point correlation functions
// ========================================================================================

#include "ThreePointCorrelation_comoving_reduced.h"
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

    const double N_R = 2.; // random/data ratio
   
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};
  
    // construct the sub-regions used for jackknife and bootstrap

    std::cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << std::endl;
    const int nx = 3, ny = 3, nz = 3;
    cbl::set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz);
  
    // -------------------------------------------------------------------------------
    // ---------------- measure the three-point correlation functions ----------------
    // -------------------------------------------------------------------------------

    // binning parameters

    const double side_s = 5.;  // 1st side of the triangle
    const double side_u = 2.;   // ratio between the 1st and 2nd sides of the triangle (u*s)
    const int nbins = 10;        // number of bins

  
    // output data
  
    const std::string dir_output = cbl::par::DirLoc+"../output/";
    const std::string dir_triplets = dir_output;
    const std::string dir_2pt = dir_output;
    const std::string file_output = "3ptJK.dat";

  
    // measure the connected and reduced three-point correlation functions and write the output

    const auto ThreeP = cbl::measure::threept::ThreePointCorrelation::Create(cbl::measure::threept::ThreePType::_comoving_connected_, catalogue, random_catalogue, cbl::triplets::TripletType::_comoving_side_, side_s, 2., side_s*side_u, 2., nbins);

    ThreeP->measure(cbl::measure::ErrorType::_None_, dir_triplets);
  
    ThreeP->write(dir_output, file_output);

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

