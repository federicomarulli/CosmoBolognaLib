// ====================================================
// Example code: how to estimate the void size function
// ====================================================

#include "LCDM.h"
#include "SizeFunction.h"

using namespace std;

int main () {
    
  try {

    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");


    // -------------------------------------------------------------
    // ---------------- estimate the void size function ------------
    // -------------------------------------------------------------
    
    cbl::cosmology::SizeFunction sf(cosmology);
    
    const double SF = sf.size_function(1., 1., "SvdW", 1.);

    cout << "The void size function at R=1 Mpc/h and z=1, predicted by the Sheth and van de Weygaert (2004) model is: " << SF << endl;
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
