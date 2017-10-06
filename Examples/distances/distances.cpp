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

    cosmobl::cosmology::Cosmology cosmology {cosmobl::cosmology::_Planck15_};


    // ---------------------------------------------------------------------------------------
    // ---------------- set the redshifts to be converted into comoving distances ------------
    // ---------------------------------------------------------------------------------------

    int step = 5e4;
    double z_min = 0.1;
    double z_max = 2.;
    vector<double> redshift = cosmobl::linear_bin_vector<double>(step, z_min, z_max);

  
    // --------------------------------------------------------------------------------------------------------
    // ---------------- convert the redshifts into comoving distances using a standard integration ------------
    // --------------------------------------------------------------------------------------------------------
  
    cout << "Computing comoving distances with a standard integration..." << endl;

    time_t start, end;
    double D_C;
    
    time (&start);

    for (auto &&zz : redshift) 
      D_C = cosmology.D_C(zz);

    cout << "D_C = " << D_C << endl;
    
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

    cout << "D_C = " << D_C << endl;
  
    time (&end);
    double T2 = difftime(end, start);
    cout << "The computational time is: " << T2 << " sec" << endl << endl;

  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
