// ============================================================================================================================
// Example code: how to measure the angle-averaged two-point correlation function, estimating the errors with different methods 
// ============================================================================================================================

#include "TwoPointCorrelation1D_monopole.h"
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

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------
  
    const double N_R = 1.; // random/data ratio
  
    cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the monopole of the two-point correlation function ----------------
    // --------------------------------------------------------------------------------------------

    // binning parameters

    const double rMin = 1.;   // minimum separation 
    const double rMax = 50.;  // maximum separation 
    const int nbins = 5;      // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre
  
  
    // construct the object used to measure the two-point correlation function
  
    cosmobl::measure::twopt::TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, cosmobl::_logarithmic_, rMin, rMax, nbins, shift};

  
    // Input/Output directories
  
    const string dir_output = cosmobl::par::DirLoc+"../output/";
    const string dir_pairs = dir_output+"pairs/";

  
    // construct the sub-regions used for jackknife and bootstrap

    cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << endl;
    const int nx = 3, ny = 3, nz = 3;
    cosmobl::set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz);

  
    // measure the monopole and compute Poissonian errors 
  
    TwoP.measure(cosmobl::measure::ErrorType::_Poisson_, dir_pairs);
    TwoP.write(dir_output, "xi_PoissonianErrors.dat");
  
  
    // measure the monopole and compute errors with jackknife (in cubic geometry)

    TwoP.measure(cosmobl::measure::ErrorType::_Jackknife_, dir_pairs);
    TwoP.write(dir_output, "xi_JackknifeErrors.dat");

  
    // measure the monopole and compute errors with bootstrap
  
    const int nM = 100; // number of mocks generated for bootstrap resampling
    TwoP.measure(cosmobl::measure::ErrorType::_Bootstrap_, dir_pairs, {dir_pairs}, "", nM);
    TwoP.write(dir_output, "xi_BootstrapErrors.dat");

    cosmobl::print(TwoP.dataset()->covariance());

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

