// ===========================================================================
// Example code: how to compute the real-space matter nonlinear power spectrum
// ===========================================================================

#include "LCDM.h"
#include "PkXiNonLinear.h"

using namespace std;

int main () {
    
  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    
    // ----------------------------------------------------------------
    // ---------------- estimate the matter power spectrum ------------
    // ----------------------------------------------------------------
    
    cbl::cosmology::PkXiNonLinear PXN(cosmology);

    vector<double> kk = cbl::logarithmic_bin_vector(100, 1.e-4, 500.);

    cout << "P_delta-delta(k=1e-4, z=1) = " << PXN.Pk_DeltaDelta(kk, 1., "EisensteinHu")[0] << endl;
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
