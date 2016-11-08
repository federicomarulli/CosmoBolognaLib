// ===================================================================================================================
// Example code: how to convert redshifts into comoving distances -- testing the performances of two different methods 
// ===================================================================================================================

#include "Cosmology.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    if (system ("clear")) {}

  
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

    cosmobl::cosmology::Cosmology cosmology {OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, n_s};


    // ---------------------------------------------------------------------------------------
    // ---------------- set the redshifts to be converted into comoving distances ------------
    // ---------------------------------------------------------------------------------------

    int step = 1e6;
    double z_min = 0.1;
    double z_max = 2.;
    vector<double> redshift = cosmobl::linear_bin_vector<double>(step, z_min, z_max);

  
    // --------------------------------------------------------------------------------------------------------
    // ---------------- convert the redshifts into comoving distances using a standard integration ------------
    // --------------------------------------------------------------------------------------------------------
  
    cout << "Computing comoving distances with a standard integration..." << endl;

    time_t start, end;
    double diffT, D_C;
  
    time (&start);

    for (auto &&zz : redshift) 
      D_C = cosmology.D_C(zz);
  
    time (&end);
    double T1 = difftime(end, start);
    cout << "The computetional time is: " << T1 << " sec" << endl << endl;
 
  
    // -----------------------------------------------------------------------------------------------------
    // ---------------- convert the redshifts into comoving distances using a smart integration ------------
    // -----------------------------------------------------------------------------------------------------

    cout << "Computing comoving distances with a faster method..." << endl;

    time (&start);

    for (auto &&zz : redshift) 
      D_C = cosmology.D_C_LCDM(zz);
  
    time (&end);
    double T2 = difftime(end, start);
    cout << "The computational time is: " << T2 << " sec" << endl << endl;

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; }
  
  return 0;
}
