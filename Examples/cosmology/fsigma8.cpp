// ============================================
// Example code: how to estimate f(z)*sigma8(z) 
// ============================================

#include "Cosmology.h"

int main () {

  try {
  
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // ---------------------------------------------------------
    
    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck18_};


    // ----------------------------------------------------------------
    // ---------------- estimate f*sigma8 at z=1 with CAMB ------------
    // ----------------------------------------------------------------

    double redshift = 1.;
    std::string method = "CAMB";

    double fs8 = cosmology.fsigma8(redshift, method);

    std::cout << "f*sigma8(z=" << redshift << ") = " << fs8 << std::endl;

    cbl::Beep("the linear growth rate times sigma8, at redshift "+cbl::conv(redshift, cbl::par::fDP0)+", is equal to "+cbl::conv(fs8, cbl::par::fDP1));
  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

