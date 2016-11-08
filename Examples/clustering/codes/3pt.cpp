// ========================================================================================
// Example code: how to measure the connected and reduced three-point correlation functions
// ========================================================================================

#include "ThreePointCorrelation_comoving_reduced.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------

    const double OmegaM = 0.25;
    const double Omega_b = 0.045;
    const double Omega_nu = 0.;
    const double massless_neutrinos = 3.04;
    const int    massive_neutrinos = 0; 
    const double OmegaL = 1.-OmegaM;
    const double Omega_radiation = 0.;
    const double hh = 0.73;
    const double scalar_amp = 2.742e-9;
    const double n_s = 1;
    const double wa = 0.;
    const double w0 = -1.;   

    cosmobl::cosmology::Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s, w0, wa};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";

    cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
   
    cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};

  
    // -------------------------------------------------------------------------------
    // ---------------- measure the three-point correlation functions ----------------
    // -------------------------------------------------------------------------------

    // binning parameters

    const double side_s = 20.;  // 1st side of the triangle
    const double side_u = 2.;   // ratio between the 1st and 2nd sides of the triangle (u*s)
    const double perc = 0.0225; // tolerance
    const int nbins = 15;       // number of bins

  
    // output data
  
    const string dir_output = cosmobl::par::DirLoc+"../output/";
    const string dir_triplets = dir_output;
    const string dir_2pt = dir_output;
    const string file_output = "3pt.dat";

  
    // measure the connected and reduced three-point correlation functions and write the output

    const auto ThreeP = cosmobl::threept::ThreePointCorrelation::Create(cosmobl::threept::_comoving_reduced_, catalogue, random_catalogue, cosmobl::triplets::_comoving_theta_, side_s, side_u, perc, nbins);

    ThreeP->measure(dir_triplets, dir_2pt);
  
    ThreeP->write(dir_output, file_output, 1);

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

