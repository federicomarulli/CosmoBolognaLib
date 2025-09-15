// ==========================================
// Example code: how to estimate the skewness
// ==========================================

#include "LCDM.h"
#include "NG.h"

using namespace std;

int main () {
    
  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
     
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");

    
    // ---------------------------------------------------
    // ---------------- estimate the skewness ------------
    // ---------------------------------------------------

    cbl::cosmology::NG ng(cosmology);
    
    const double SK = ng.skewness(1.e13, "EisensteinHu");
 
    cout << "The skewness at M=1e13 Msun/h is:" << SK << endl << endl;
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
