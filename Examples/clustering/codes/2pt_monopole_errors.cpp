// ============================================================================================================================
// Example code: how to measure the angle-averaged two-point correlation function, estimating the errors with different methods 
// ============================================================================================================================

#include "TwoPointCorrelation1D_monopole.h"
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

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------
  
    const double N_R = 3.; // random/data ratio
  
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};

  
    // construct the sub-regions used for jackknife and bootstrap

    std::cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << std::endl;
    const int nx = 3, ny = 3, nz = 3;
    cbl::set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz);

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the monopole of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters

    const double rMin = 10.;   // minimum separation 
    const double rMax = 30.;  // maximum separation 
    const int nbins = 3;      // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre
  
  
    // construct the object used to measure the two-point correlation function
  
    cbl::measure::twopt::TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, cbl::BinType::_logarithmic_, rMin, rMax, nbins, shift};

  
    // Input/Output directories
  
    const std::string dir_output = "../output/";
    const std::string dir_pairs = dir_output+"pairs/";

  
    // measure the monopole and compute Poissonian errors 
  
    TwoP.measure(cbl::measure::ErrorType::_Poisson_, dir_pairs);
    TwoP.write(dir_output, "xi_PoissonianErrors.dat");
  
  
    // measure the monopole and compute errors with jackknife (in cubic geometry)

    TwoP.measure(cbl::measure::ErrorType::_Jackknife_, dir_pairs);
    TwoP.write(dir_output, "xi_JackknifeErrors.dat");

  
    // measure the monopole and compute errors with bootstrap
  
    const int nM = 100; // number of mocks generated for bootstrap resampling
    TwoP.measure(cbl::measure::ErrorType::_Bootstrap_, dir_pairs, {dir_pairs}, "", nM);
    TwoP.write(dir_output, "xi_BootstrapErrors.dat");

    cbl::Print(TwoP.dataset()->covariance());

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

