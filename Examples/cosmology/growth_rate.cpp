// ======================================================================
// Example code: how to compute the matter power spectrum and growth rate
// ======================================================================

#include "kw0wanuCDM.h"
#include "PkXi.h"

using namespace std;

int main () {

  try {
    
    // --------------------------------------------
    // --------- set a cosmological model ---------  
    // --------------------------------------------

    auto cosmology = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

    cosmology->set_parameter("w0",-1.1);
    cosmology->set_parameter("wa",0.3);
    
    cbl::cosmology::PkXi PX(cosmology);
    
    
    // -----------------------------------------------------------------------
    // --------- compute the growth factor in two different ways -------------
    // -----------------------------------------------------------------------  

    // select two different redshifts
    const double redshift1 = 0.5;
    const double redshift2 = 2.;

    // select a scale at which the power spectra are computed
    const double kk = 10.;

    // compute the ratio of the power spectra at the two redshifts
    const double ratio_Pk = PX.Pk_matter({kk}, "CAMB", false, redshift2)[0]/PX.Pk_matter({kk}, "CAMB", false, redshift1)[0];
    cout << endl << "P(z=" << redshift2 << ")/P(z=" << redshift1 << ") = " << ratio_Pk << endl;

    // compare the latter value with the squared value of the growth factor normalised at z=0
    const double GF2 = pow(cosmology->DN(redshift2, redshift1), 2.);
    cout << "[D(z=" << redshift2 << ")/D(z=" << redshift1 << ")]^2 = " << GF2 << endl << endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

