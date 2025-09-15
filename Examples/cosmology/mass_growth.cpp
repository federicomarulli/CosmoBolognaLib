// ===========================================================================
// Example code: how to compute the cumulative formation redshift distribution
// ===========================================================================

#include "LCDM.h"
#include "MassGrowth.h"

using namespace std;

int main () {
    
  try {
    
    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");

    
    // -------------------------------------------------------------------------------------
    // ---------------- estimate the cumulative formation redshift distribution ------------
    // -------------------------------------------------------------------------------------
    
    cbl::cosmology::MassGrowth mg(cosmology);
    
    const double ww = 1.; // rescaled variable w as in Lacey and Coles (1993)
    const double ff = 0.5; // assembled fraction
    const string author = "GTS"; // Giocoli et al. (2012)
				       
    cout << "The cumulative formation redshift distribution at w=1 and f=0.5, predicted by Giocoli et al. (2012) model is: " << mg.cumPw(ww, ff, author) << endl;
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }

  return 0;
}
