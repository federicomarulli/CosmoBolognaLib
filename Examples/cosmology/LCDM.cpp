// ==================================================
// Example code: how to set a LCDM cosmological model
// ==================================================

#include "LCDM.h"

using namespace std;


int main () {

  try {
    
    // --------------------------------------------------------------------------
    // ----- using one of the built-in set of cosmological model parameters -----
    // --------------------------------------------------------------------------
    
    cbl::cosmology::LCDM cosmology("Planck18");

    cout << endl << "Planck18:" << endl;
    cosmology.print_parameters();
    
    
    // ------------------------------------------------------------------
    // ----- setting directly all the cosmological model parameters -----
    // ------------------------------------------------------------------

    const double OmegaM = 0.25;
    const double Omega_b = 0.045;
    const double Omega_radiation = 0.;
    const double little_h = 0.73;
    const double scalar_amp = 2.742e-9;
    const double scalar_pivot = 0.05;
    const double n_s = 1.;
    const double tau = 0.089;
    const bool unit = true;

    cbl::cosmology::LCDM cosmology_test(OmegaM, Omega_b, Omega_radiation, little_h, scalar_amp, scalar_pivot, n_s, tau, unit);

    cout << endl << "Test parameters:" << endl;
    cosmology_test.print_parameters();
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
} 

