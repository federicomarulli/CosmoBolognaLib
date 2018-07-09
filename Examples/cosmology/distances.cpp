// ===================================================================================================================
// Example code: how to convert redshifts into comoving distances -- testing the performances of two different methods 
// ===================================================================================================================

#include "Cosmology.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


int main () {

  try {
  
    if (system ("clear")) {}

  
    // -------------------------------------------------------------
    // ---------------- set the cosmological parameters ------------
    // -------------------------------------------------------------

    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};


    // ---------------------------------------------------------------------------------------
    // ---------------- set the redshifts to be converted into comoving distances ------------
    // ---------------------------------------------------------------------------------------

    int step = 10000;
    double z_min = 0.1;
    double z_max = 2.;
    std::vector<double> redshift = cbl::linear_bin_vector<double>(step, z_min, z_max);

  
    // --------------------------------------------------------------------------------------------------------
    // ---------------- convert the redshifts into comoving distances using a standard integration ------------
    // --------------------------------------------------------------------------------------------------------
  
    std::cout << "Computing comoving distances with a standard integration..." << std::endl;

    time_t start, end;
    double D_C;
    
    time (&start);

    for (auto &&zz : redshift) 
      D_C = cosmology.D_C(zz);

    std::cout << "D_C = " << D_C << std::endl;
    
    time (&end);
    double T1 = difftime(end, start);
    std::cout << "The computetional time is: " << T1 << " sec" << std::endl << std::endl;
 
  
    // -----------------------------------------------------------------------------------------------------
    // ---------------- convert the redshifts into comoving distances using a smart integration ------------
    // -----------------------------------------------------------------------------------------------------

    std::cout << "Computing comoving distances with a faster method..." << std::endl;

    time (&start);

    for (auto &&zz : redshift) 
      D_C = cosmology.D_C_LCDM(zz);

    std::cout << "D_C = " << D_C << std::endl;
  
    time (&end);
    double T2 = difftime(end, start);
    std::cout << "The computational time is: " << T2 << " sec" << std::endl << std::endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
