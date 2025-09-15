// ====================================================
// Example code: how to estimate the halo mass function
// ====================================================

#include "LCDM.h"
#include "MassFunction.h"

using namespace std;

int main () {
    
  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    
    // -------------------------------------------------------------
    // ---------------- estimate the halo mass function ------------
    // -------------------------------------------------------------
    
    cbl::cosmology::MassFunction mf(cosmology);
    
    const double MF = mf.mass_function(1.e13, 1., "PS", "EisensteinHu");
				       
    cout << "The halo mass function at M=1e13 Msun/h and z=1, predicted by the Press & Schechter model is:" << MF << endl;

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
