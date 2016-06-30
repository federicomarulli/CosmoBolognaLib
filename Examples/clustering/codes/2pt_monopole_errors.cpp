// ============================================================================================================================
// Example code: how to measure the angle-averaged two-point correlation function, estimating the errors with different methods 
// ============================================================================================================================

#include "TwoPointCorrelation1D_monopole.h"
#include "GlobalFunc.h"

using namespace cosmobl;
using namespace cosmology;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -----------------------------------------------------------------
  // ---------------- use default cosmological parameters ------------
  // -----------------------------------------------------------------

  Cosmology cosmology;
  
  
  // -----------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
  // -----------------------------------------------------------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";
  
  Catalogue catalogue {_Galaxy_, _observedCoordinates_, {file_catalogue}, cosmology};

  
  // --------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue (with cubic geometry) ----------------
  // --------------------------------------------------------------------------------------
  
  double N_R = 1.; // random/data ratio
  
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};

  
  // --------------------------------------------------------------------------------------------
  // ---------------- measure the monopole of the two-point correlation function ----------------
  // --------------------------------------------------------------------------------------------

  // binning parameters

  double rMin = 1.;   // minimum separation 
  double rMax = 50.;  // maximum separation 
  int nbins = 5;      // number of bins
  double shift = 0.5; // spatial shift used to set the bin centre
  
  
  // construct the object used to measure the two-point correlation function
  
  TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift};

  
  // Input/Output directories
  
  string dir_output = par::DirLoc+"../output/";
  string dir_pairs = dir_output+"pairs/";

  
  // construct the sub-regions used for jackknife and bootstrap

  cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << endl;
  int nx = 3, ny = 3, nz = 3;
  set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz);

  
  // measure the monopole and compute Poissonian errors 
  
  TwoP.measure(_Poisson_, dir_pairs);
  TwoP.write(dir_output, "xi_PoissonianErrors.dat");
  
  
  // measure the monopole and compute errors with jackknife (in cubic geometry)

  TwoP.measure(_Jackknife_, dir_pairs);
  TwoP.write(dir_output, "xi_JackknifeErrors.dat");

  
  // measure the monopole and compute errors with bootstrap
  
  int nM = 20; // number of mocks generated for bootstrap resampling
  TwoP.measure(_Bootstrap_, dir_pairs, {dir_pairs}, "", nM);
  TwoP.write(dir_output, "xi_BootstrapErrors.dat");

  print(TwoP.dataset()->covariance());
  
  return 0;

}

