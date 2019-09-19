// ============================================
// Example code: how to estimate f(z)*sigma8(z) 
// ============================================

#include "Cosmology.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;

int main () {

  try {
  
    // -------------------------------------------------------------
    // ---------------- set the cosmological parameters ------------
    // -------------------------------------------------------------
    
    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};


    // ----------------------------------------------------------------
    // ---------------- estimate f*sigma8 at z=1 with CAMB ------------
    // ----------------------------------------------------------------

    double redshift = 1.;
    std::string method = "CAMB";

    double fs8 = cosmology.fsigma8(redshift, method);

    std::cout << "f*sigma8(z=" << redshift << ") = " << fs8 << std::endl;
    
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

