// =======================================================================
// Example code: how to estimate the matter two-point correlation function
// =======================================================================

#include "LCDM.h"
#include "PkXi.h"

using namespace std;

int main () {
    
  try {

    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
     
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    cosmology->set_parameter("sigma8",0.8);

    
    // --------------------------------------------------------------------------------
    // ---------------- estimate the matter two-point correlation function ------------
    // --------------------------------------------------------------------------------
   
    cbl::cosmology::PkXi PX(cosmology);
    
    const double rr = 50.;
    const double redshift = 1.;
    
    const double xi_matter = PX.xi_matter(rr, "EisensteinHu", false, redshift);
    
    cout << "The matter two-point correlation function at r=50 Mpc/h and z=1 is: " << xi_matter << endl;

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
