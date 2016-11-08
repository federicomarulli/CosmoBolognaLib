// ============================================
// Example code: how to estimate f(z)*sigma8(z) 
// ============================================

#include "Cosmology.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {
  
    // -------------------------------------------------------------
    // ---------------- set the cosmological parameters ------------
    // -------------------------------------------------------------

    double OmegaM = 0.25;
    double Omega_b = 0.045;
    double Omega_nu = 0.;
    double massless_neutrinos = 3.04;
    int massive_neutrinos = 0; 
    double OmegaL = 1.-OmegaM;
    double Omega_radiation = 0.;
    double hh = 0.73;
    double scalar_amp = 2.742e-9;
    double n_s = 1;
    double wa = 0.;
    double w0 = -1.;   

    cosmobl::cosmology::Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s, w0, wa};


    // ----------------------------------------------------------------
    // ---------------- estimate f*sigma8 at z=1 with CAMB ------------
    // ----------------------------------------------------------------

    double redshift = 1.;
    string method = "CAMB";

    double fs8 = cosmology.fsigma8(redshift, method);

    cout << "f*sigma8(z=" << redshift << ") = " << fs8 << endl;
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
} 

