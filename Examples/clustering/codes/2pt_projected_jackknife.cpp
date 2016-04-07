// ==============================================================================================
// How to measure the projected two-point correlation function and estimate errors with jackknife
// ==============================================================================================

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

  
  // ----------------------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with polar coordinates: RA, Dec, redshift) --------------------------------
  // ----------------------------------------------------------------------------------------------------------------------

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

  double rMin = 1.;           // minimum separation 
  double rMax = 50.;          // maximum separation 
  int nbins = 20;             // number of bins
  double shift = 0.5;         // spatial shift used to set the bin centre 
  double piMax_integral = 30; // upper limit of the integral
  

  // measure the projected two-point correlation function and estimate errors with jackknife (with cubic geometry)

  int nx = 3, ny = 3, nz = 3;
  set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz);
  
  auto TwoP = TwoPointCorrelation::Create(TwoPType::_1D_projected_, catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift, rMin, rMax, nbins, shift, piMax_integral);

  TwoP->measure(_Jackknife_, dir_output);

  TwoP->write(dir_output, "xi_projected_jackknife.dat");

  return 0;
}

