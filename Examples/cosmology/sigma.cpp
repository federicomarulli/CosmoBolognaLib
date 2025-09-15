// ===============================================
// Example code: how to estimate the mass variance
// ===============================================

#include "LCDM.h"
#include "Sigma.h"

using namespace std;

int main () {
    
  try {

    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
     
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    

    // --------------------------------------------------------
    // ---------------- estimate the mass variance ------------
    // --------------------------------------------------------

    cbl::cosmology::Sigma SG(cosmology);
    
    const vector<double> mass = {1.e15};
    const double redshift = 1.;
    
    const double sigma2M = SG.sigma2M(mass, "EisensteinHu", redshift)[0];
    
    cout << "The mass variance at M=1e15 Msun and z=1 is: " << sigma2M << endl;

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
