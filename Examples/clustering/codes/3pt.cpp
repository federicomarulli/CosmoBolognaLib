// ========================================================================================
// Example code: how to measure the connected and reduced three-point correlation functions
// ========================================================================================

#include "ThreePointCorrelation_comoving_reduced.h"

using namespace cosmobl;
using namespace cosmology;
using namespace catalogue;
using namespace threept;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // --------------------------------------------------------------
  // ---------------- set the cosmological parameters  ------------
  // --------------------------------------------------------------

  double OmegaM = 0.25;
  double Omega_b = 0.045;
  double Omega_nu = 0.;
  double massless_neutrinos = 3.04;
  int    massive_neutrinos = 0; 
  double OmegaL = 1.-OmegaM;
  double Omega_radiation = 0.;
  double hh = 0.73;
  double scalar_amp = 2.742e-9;
  double n_s = 1;
  double wa = 0.;
  double w0 = -1.;   

  Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s, w0, wa};

  
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

  
  // -------------------------------------------------------------------------------
  // ---------------- measure the three-point correlation functions ----------------
  // -------------------------------------------------------------------------------

  // binning parameters

  double side_s = 20.;  // 1st side of the triangle
  double side_u = 2.;   // ratio between the 1st and 2nd sides of the triangle (u*s)
  double perc = 0.0225; // tolerance
  int nbins = 15;       // number of bins

  
  // output data
  
  string dir_output = par::DirLoc+"../output/";
  string dir_triplets = dir_output;
  string dir_2pt = dir_output;
  string file_output = "3pt.dat";

  
  // measure the connected and reduced three-point correlation functions and write the output

  auto ThreeP = ThreePointCorrelation::Create(_comoving_reduced_, catalogue, random_catalogue, triplets::_comoving_theta_, side_s, side_u, perc, nbins);

  ThreeP->measure(dir_triplets, dir_2pt);
  
  ThreeP->write(dir_output, file_output, 1);

  
  return 0;
}

