// ==============================================================
// Example code: how to measure 2D two-point correlation function
// ==============================================================

#include "TwoPointCorrelation2D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace twopt;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {
  
  // -------------------------------------------------------------
  // ---------------- set the cosmological parameters ------------
  // -------------------------------------------------------------

  double OmegaM = 0.27;
  double Omega_b = 0.046;
  double Omega_nu = 0.;
  double massless_neutrinos = 3.04;
  int    massive_neutrinos = 0; 
  double OmegaL = 1.-OmegaM;
  double Omega_radiation = 0.;
  double hh = 0.7;
  double scalar_amp = 2.46e-9;
  double n_s = 0.96;
  double wa = 0.;
  double w0 = -1.;   
  
  Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s, w0, wa};
  
  
  // ----------------------------------------------------------
  // ---------------- read the input catalogue ----------------
  // ----------------------------------------------------------
  
  string file_catalogue = par::DirLoc+"../input/cat.dat";

  ifstream fin(file_catalogue.c_str());  checkIO(file_catalogue, 1);
 
  string line;
  double RA, DEC, RED;
  vector<double> ra, dec, redshift;
 
  while (getline(fin, line)) {
    stringstream ss(line);
    ss >> RA; ss >> DEC; ss >> RED;
    ra.push_back(RA);
    dec.push_back(DEC);
    redshift.push_back(RED);
  }
  
  fin.clear(); fin.close();
  
  Catalogue catalogue {_Galaxy_, _observedCoordinates_, ra, dec, redshift, cosmology};
  
  
  // --------------------------------------------------------------------------------------
  // ---------------- construct the random catalogue (with cubic geometry) ---------------- 
  // -------------------------------------------------------------------------------------- 

  double N_R = 1.; // random/data ratio
 
  Catalogue random_catalogue {_createRandom_box_, catalogue, N_R};
  
  
  // -----------------------------------------------------------------------------------------------
  // ---------------- measure the 2D two-point correlation in Cartesian coordinates ----------------
  // -----------------------------------------------------------------------------------------------

  // ----- output data ----- 

  string dir_pairs = par::DirLoc+"../output/";
  string dir_output = par::DirLoc+"../output/";

  
  // ----- measure the 2D correlation function in Cartesian coordinates, xi(rp,pi), and store the outputs ----- 
  
  double rpMin = 1.;     // minimum separation in the first dimension
  double rpMax = 50.;    // maximum separation in the first dimension 
  int nbins_D1 = 50;     // number of bins in the first dimension
  double shift_D1 = 0.5; // spatial shift used to set the bin centre in the first dimension
  double piMin = 1.;     // minimum separation in the second dimension
  double piMax = 50.;    // maximum separation in the second dimension 
  int nbins_D2 = 50;     // number of bins in the second dimension
  double shift_D2 = 0.5; // spatial shift used to set the bin centre in the second dimension
  
  // construct the object using a static factory
  auto xi2DCart = TwoPointCorrelation::Create(_2D_Cartesian_, catalogue, random_catalogue, _linear_, rpMin, rpMax, nbins_D1, shift_D1, _linear_, piMin, piMax, nbins_D2, shift_D2);

  // measure the 2D correlation function and compute Poisson errors
  xi2DCart->measure(_Poisson_, dir_pairs);
  xi2DCart->write(dir_output, "xi_rp_pi_linlin.dat");

  
  return 0;

}

