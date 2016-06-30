// =========================================================================
// Example code: how to measure the projected two-point correlation function
// =========================================================================

#include "TwoPointCorrelation1D.h"

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

  
  // ---------------------------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) --------------------------------
  // ---------------------------------------------------------------------------------------------------------------------------

  string file_catalogue = par::DirLoc+"../input/cat.dat";

  Catalogue catalogue {_Galaxy_, _observedCoordinates_, {file_catalogue}, cosmology};

  
  // --------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue (with cubic geometry) ----------------
  // --------------------------------------------------------------------------------------

  double N_R = 1.; // random/data ratio
  
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};


  // --------------------------------------------------------------------------------------
  // ---------------- measure the projected two-point correlation function ----------------
  // --------------------------------------------------------------------------------------

  // binning parameters and output data

  double rMin = 1.;           // minimum separation 
  double rMax = 50.;          // maximum separation 
  int nbins = 20;             // number of bins
  double shift = 0.5;         // spatial shift used to set the bin centre 
  double piMax_integral = 30; // upper limit of the integral

  string dir = par::DirLoc+"../output/";
  string file = "xi_projected.dat";
  

  // measure the projected two-point correlation function
  
  auto TwoP = TwoPointCorrelation::Create(TwoPType::_1D_projected_, catalogue, random_catalogue, _logarithmic_, rMin, rMax, nbins, shift, rMin, rMax, nbins, shift, piMax_integral);

  TwoP->measure(_Poisson_, dir);

  TwoP->write(dir, file);

  
  return 0;
}

