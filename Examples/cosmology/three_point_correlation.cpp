// ======================================================
// Example code: how to estimate the function Gamma_{123} 
// ======================================================

#include "LCDM.h"
#include "3PCF.h"

using namespace std;

int main () {

  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // ---------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");

    
    // ----------------------------------------------------------------
    // ---------------- estimate the function Gamma_{123}  ------------
    // ----------------------------------------------------------------
  
    cbl::cosmology::ThreePointCorrelation threept(cosmology);
    
    const double Gamma = threept.Gamma_3PCF(1., 1., 1., {1., 1,}, {1., 1.});
    
    cout << "The non-local contribution to the three-point correlation function at r1=r2=1 Mpc/h is: " << Gamma << endl;
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
