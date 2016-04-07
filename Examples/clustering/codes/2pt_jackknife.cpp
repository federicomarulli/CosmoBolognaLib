// ====================================================================================
// How to measure the two-point correlation function and estimate errors with jackknife 
// ====================================================================================

#include "TwoPointCorrelation1D_monopole.h"
#include "GlobalFunc.h"

using namespace cosmobl;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -----------------------------------------------------------------
  // ---------------- use default cosmological parameters ------------
  // -----------------------------------------------------------------

  Cosmology cosmology;

  
  // --------------------------------------------------------------------
  // ---------------- Input/Output files and directories ----------------
  // --------------------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  string dir_output = par::DirLoc+"../output/";
  string dir_pairs = dir_output+"pairs/";
  string dir_covariance = dir_output+"covariance/";
  
  string MK = "mkdir -p "+dir_output+" "+dir_pairs+" "+dir_covariance; if (system(MK.c_str())) {};

  
  // ------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with polar coordinates: RA, Dec, redshift) ----------------
  // ------------------------------------------------------------------------------------------------------
  
  cout << "I'm reading the input catalogue..." << endl;

  Catalogue catalogue {_Galaxy_, {file_catalogue}, cosmology};

  
  // --------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue (with cubic geometry) ----------------
  // --------------------------------------------------------------------------------------
  
  double N_R = 1.; // random/object ratio
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};

  
  // --------------------------------------------------------------------------------------------
  // ---------------- measure the monopole of the two-point correlation function ----------------
  // --------------------------------------------------------------------------------------------

  // binning parameters

  double rMin = 1.;   // minimum separation 
  double rMax = 50.;  // maximum separation 
  int nbins = 20;     // number of bins
  double shift = 0.5; // spatial shift used to set the bin centre 

  
  // measure the monopole and compute errors with jackknife (in cubic geometry)

  int nx = 3, ny = 3, nz = 3;
  set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz);
  
  TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift};

  TwoP.measure(_Jackknife_, dir_output);
  
  TwoP.write(dir_output, "xi_Jackknife.dat");

  
  return 0;

}

