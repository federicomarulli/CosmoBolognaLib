// =======================================================================
// Example code: how to measure the angular two-point correlation function 
// =======================================================================

#include "TwoPointCorrelation1D_angular.h"

using namespace cosmobl;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // ---------------------------------------------------------------------------------------------------------------------------
  // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) --------------------------------
  // ---------------------------------------------------------------------------------------------------------------------------

  string file_catalogue = par::DirLoc+"../input/cat.dat";
  
  Catalogue catalogue {_Galaxy_, _observedCoordinates_, {file_catalogue}};

  
  // --------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue (with cubic geometry) ----------------
  // --------------------------------------------------------------------------------------

  double N_R = 1.; // random/data ratio
  
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};


  // --------------------------------------------------------------------------------------------
  // ---------------- measure the angular of the two-point correlation function ----------------
  // --------------------------------------------------------------------------------------------

  // binning parameters and output data

  double angMin = 1.;                  // minimum angular separation 
  double angMax = 50.;                 // maximum angular separation 
  int nbins = 20;                      // number of bins
  double shift = 0.5;                  // shift used to set the bin centre 
  CoordUnits angularUnits = _degrees_; // angular units

  string dir = par::DirLoc+"../output/";
  string file = "xi_projected.dat";
  

  // measure the angular two-point correlation function
  
  TwoPointCorrelation1D_angular TwoP {catalogue, random_catalogue, _linear_, angMin, angMax, nbins, shift, angularUnits};
  
  TwoP.measure(_Poisson_, dir);

  TwoP.write(dir, file);

  
  return 0;

}

