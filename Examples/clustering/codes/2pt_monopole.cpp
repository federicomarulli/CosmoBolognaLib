// =================================================================================================
// Example code: how to measure the angle-averaged two-point correlation function, i.e. the monopole 
// =================================================================================================

#include "TwoPointCorrelation1D_monopole.h"

using namespace cosmobl;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -----------------------------------------------------------------
  // ---------------- use default cosmological parameters ------------
  // -----------------------------------------------------------------

  Cosmology cosmology;

  
  // ------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with polar coordinates: RA, Dec, redshift) ----------------
  // ------------------------------------------------------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  Catalogue catalogue {_Galaxy_, {file_catalogue}, cosmology};

  
  // --------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue (with cubic geometry) ----------------
  // --------------------------------------------------------------------------------------

  double N_R = 1.; // random/object ratio
  
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};

  
  // --------------------------------------------------------------------------------------------
  // ---------------- measure the monopole of the two-point correlation function ----------------
  // --------------------------------------------------------------------------------------------

  // binning parameters and output data

  double rMin = 1.;   // minimum separation 
  double rMax = 50.;  // maximum separation 
  int nbins = 20;     // number of bins
  double shift = 0.5; // spatial shift used to set the bin centre
  
  string dir = par::DirLoc+"../output/";
  string file = "xi.dat";

  
  // measure the monopole and compute Poisson errors 

  TwoPointCorrelation1D_monopole TwoP {catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift};
  
  TwoP.measure(_Poisson_, dir);

  
  // store the output data
  
  TwoP.write(dir, file);
  
  
  return 0;
}

