// ========================================================================================================================
// Example code: how to model the 2D two-point correlation function in Cartesian coordinates, xi(rp, pi), in redshift-space
// ========================================================================================================================

#include "TwoPointCorrelation2D.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // -------------------------------------------------------------
    // ---------------- set the cosmological parameters ------------
    // -------------------------------------------------------------

    const double OmegaM = 0.27;
    const double Omega_b = 0.046;
    const double Omega_nu = 0.;
    const double massless_neutrinos = 3.04;
    const int    massive_neutrinos = 0; 
    const double OmegaL = 1.-OmegaM;
    const double Omega_radiation = 0.;
    const double hh = 0.7;
    const double scalar_amp = 2.46e-9;
    const double n_s = 0.96;
    const double wa = 0.;
    const double w0 = -1.;   
  
    const cosmobl::cosmology::Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s, w0, wa};
  
  
    // ----------------------------------------------------------
    // ---------------- read the input catalogue ----------------
    // ----------------------------------------------------------
  
    const string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";

    const cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};
  
  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ---------------- 
    // -------------------------------------------------------------------------------------- 

    const double N_R = 1.; // random/data ratio
 
    const cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};
  
  
    // -----------------------------------------------------------------------------------------------
    // ---------------- measure the 2D two-point correlation in Cartesian coordinates ----------------
    // -----------------------------------------------------------------------------------------------

    // ----- output data ----- 

    const string dir_pairs = cosmobl::par::DirLoc+"../output/";
    const string dir_output = cosmobl::par::DirLoc+"../output/";

  
    // ----- measure the 2D correlation function in Cartesian coordinates, xi(rp,pi), and store the outputs ----- 
  
    const double rpMin = 1.;     // minimum separation in the first dimension
    const double rpMax = 50.;    // maximum separation in the first dimension 
    const int nbins_D1 = 50;     // number of bins in the first dimension
    const double shift_D1 = 0.5; // spatial shift used to set the bin centre in the first dimension
    const double piMin = 1.;     // minimum separation in the second dimension
    const double piMax = 50.;    // maximum separation in the second dimension 
    const int nbins_D2 = 50;     // number of bins in the second dimension
    const double shift_D2 = 0.5; // spatial shift used to set the bin centre in the second dimension
  
    // construct the object using a static factory
    const auto xi2DCart = cosmobl::twopt::TwoPointCorrelation::Create(cosmobl::twopt::_2D_Cartesian_, catalogue, random_catalogue, cosmobl::_linear_, rpMin, rpMax, nbins_D1, shift_D1, cosmobl::_linear_, piMin, piMax, nbins_D2, shift_D2);

    // measure the 2D correlation function and compute Poisson errors
    xi2DCart->measure(cosmobl::twopt::_Poisson_, dir_pairs);
    xi2DCart->write(dir_output, "xi_rp_pi_linlin.dat");

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}

