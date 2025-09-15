// ==============================================================================
// Example code: how to compute the z-space matter monopole from the Kaiser model
// ==============================================================================

#include "LCDM.h"
#include "PkXizSpace.h"

using namespace std;

int main () {
    
  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    cosmology->set_parameter("sigma8",0.8111);
    
    
    // ----------------------------------------------------------------
    // ---------------- estimate the matter power spectrum ------------
    // ----------------------------------------------------------------
    
    cbl::cosmology::PkXizSpace PXS(cosmology);

    cout << "xi0(r=1 Mpc/h, fsigma8=1, bsigma8=1, z=1) = " << PXS.xi0_Kaiser(1., 1., 1., "EisensteinHu", 1.) << endl;
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
