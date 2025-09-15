// =============================================================
// Example code: how to estimate the sound horizon at drag epoch 
// =============================================================

#include "LCDM.h"
#include "BAO.h"

using namespace std;

int main () {
    
  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------

    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    
    // ----------------------------------------------------------------------
    // ---------------- estimate the sound horizon at drag epoch ------------
    // ----------------------------------------------------------------------

    cbl::cosmology::BAO bao(cosmology);

    const double rs = bao.rs("EisensteinHu");

    cout << "The sound horizon at drag epoch is: " << rs << endl;
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
