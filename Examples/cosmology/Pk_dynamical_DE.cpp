// ================================================================================
// Example code: how to compute the matter power spectrum and growth rate with CAMB
// ================================================================================

#include "Cosmology.h"

using namespace std;

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;

int main () {

  try {
    
    // -----------------------------------------------------------------
    // ------------------ set a cosmological model  -----------------
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
    // use non-standard values of wa and w0
    const double wa = 0.3;
    const double w0 = -1.1;

    cbl::cosmology::Cosmology cosm {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, scalar_pivot, n_s, w0, wa};

    
    // -----------------------------------------------------------------
    // --------- compute power spectrum and growth factor  -------------
    // -----------------------------------------------------------------    

    // select 2 different redshifts
    const double redshift1 = 0.5;
    const double redshift2 = 2.;

    // select a scale at which the power spectrum is computed
    const double kk = 10.;

    // compute the ratio of the power spectrum computed a the two redshifts
    const double ratio_Pk = cosm.Pk_matter(kk, "CAMB", false, redshift2)/cosm.Pk_matter(kk, "CAMB", false, redshift1);
    cout << endl << "P(z=" << redshift2 << ")/P(z=" << redshift1 << ") = " << ratio_Pk << endl;

    // compare the latter value with the squared value of the growth factor normalised at z=0
    const double GF2 = pow(cosm.DN(redshift2, redshift1), 2.);
    cout << "[D(z=" << redshift2 << ")/D(z=" << redshift1 << ")]^2 = " << GF2 << endl << endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

