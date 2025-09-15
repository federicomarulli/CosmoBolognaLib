// ===========================================
// Example code: how to estimate the halo bias
// ===========================================

#include "LCDM.h"
#include "Bias.h"

using namespace std;

int main () {
    
  try {

    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");

    
    // ----------------------------------------------------
    // ---------------- estimate the halo bias ------------
    // ----------------------------------------------------

    cbl::cosmology::Bias bias(cosmology);
    
    const double bias_halo = bias.bias_halo(1.e13, 1., "ST99", "EisensteinHu");
    
    cout << "The halo bias at M=1e13 Msun/h and z=1, predicted by the Sheth & Tormen (1999) model, is:" << bias_halo << endl;
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
