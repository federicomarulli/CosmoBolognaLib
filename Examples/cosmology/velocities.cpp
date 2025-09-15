// ====================================================
// Example code: how to estimate the Cosmic Mach Number
// ====================================================

#include "LCDM.h"
#include "Velocities.h"

using namespace std;

int main () {
    
  try {

    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
     
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");

    
    // -------------------------------------------------------------
    // ---------------- estimate the Cosmic Mach Number ------------
    // -------------------------------------------------------------
   
    cbl::cosmology::Velocities vel(cosmology);
    
    const double CMN = vel.CMN(2., 2., "CAMB", 1.);
    
    cout << "The Cosmic Mach Number at r=2 Mpc/h and z=1 is: " << CMN << endl;

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
