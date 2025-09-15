// ============================================
// Example code: how to estimate f(z)*sigma8(z) 
// ============================================

#include "LCDM.h"
#include "RSD.h"

using namespace std;

int main () {

  try {
  
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    

    // ------------------------------------------------------
    // ---------------- estimate f*sigma8 at z=1 ------------
    // ------------------------------------------------------

    cbl::cosmology::RSD rsd(cosmology);
    
    const double redshift = 1.;
    const string method = "EisensteinHu";

    const double fs8 = rsd.fsigma8(redshift, method);

    cout << "f*sigma8(z=" << redshift << ") = " << fs8 << endl;

    cbl::Beep("the linear growth rate times sigma8, at redshift "+cbl::conv(redshift, cbl::par::fDP0)+", is equal to "+cbl::conv(fs8, cbl::par::fDP1));
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
} 

