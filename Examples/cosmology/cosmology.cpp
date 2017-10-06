// =============================================
// Example code: how to set a cosmological model
// =============================================

#include "Cosmology.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {

    // ------------------------------------------------------
    // ---------------- using default parameters ------------
    // ------------------------------------------------------
    
    cosmobl::cosmology::Cosmology cosmo1;

    
    // ---------------------------------------------------------------------------
    // ---------------- using one of the built-in cosmological models ------------
    // ---------------------------------------------------------------------------
    
    cosmobl::cosmology::Cosmology cosmo2 {cosmobl::cosmology::_Planck15_};

    
    // -----------------------------------------------------------------
    // ---------------- setting the cosmological parameters ------------
    // -----------------------------------------------------------------

    const double OmegaM = 0.25;
    const double Omega_b = 0.045;
    const double Omega_nu = 0.;
    const double massless_neutrinos = 3.04;
    const int massive_neutrinos = 0; 
    const double OmegaL = 1.-OmegaM;
    const double Omega_radiation = 0.;
    const double hh = 0.73;
    const double scalar_amp = 2.742e-9;
    const double scalar_pivot = 0.05;
    const double n_s = 1;
    const double wa = 0.;
    const double w0 = -1.;   

    cosmobl::cosmology::Cosmology cosmo3 {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, scalar_pivot, n_s, w0, wa};
    
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

