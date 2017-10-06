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

    cosmobl::cosmology::Cosmology cosmology {cosmobl::cosmology::_Planck15_};


    // ----------------------------------------------------------------
    // ---------------- estimate f*sigma8 at z=1 with CAMB ------------
    // ----------------------------------------------------------------

    double redshift = 1.;
    string method = "CAMB";

    double fs8 = cosmology.fsigma8(redshift, method);

    cout << "f*sigma8(z=" << redshift << ") = " << fs8 << endl;
  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

