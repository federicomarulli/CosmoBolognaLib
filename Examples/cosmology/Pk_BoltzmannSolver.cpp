// =================================================================================================
// Example code: how to compute the matter power spectrum with CAMB and CLASS at different redshifts
// =================================================================================================

#include "LCDM.h"
#include "PkXi.h"

using namespace std;

int main () {

  try {
    
    // -----------------------------------------------------------------
    // ------------------ set a cosmological model  -----------------
    // -----------------------------------------------------------------

    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    
    // -----------------------------------------------------------------
    // --------- compute power spectrum with CAMB and CLASS  -----------
    // -----------------------------------------------------------------    

    cbl::cosmology::PkXi PX(cosmology);
    
    // choose a vector of redshifts
    const vector<double> redshifts = {0., 0.5, 1., 1.5};

    // choose a vector of scales
    const vector<double> kk = cbl::logarithmic_bin_vector(100, 0.001, 1.);

    // choose the linear power spectrum
    const bool do_NL = false;

    // compute the power spectra with different Boltzmann solver and compare them
    vector<vector<double>> Pk_CAMB = PX.Pk_matter(kk, "CAMB", do_NL, redshifts);
    vector<vector<double>> Pk_CLASS = PX.Pk_matter(kk, "CLASS", do_NL, redshifts);

    for (size_t ii=0; ii<redshifts.size(); ii++)
      cout << "At redshift z = "+cbl::conv(redshifts[ii], cbl::par::fDP2)+", at scale k = "+cbl::conv(kk[50], cbl::par::fDP2)+" Mpc/h, the relative percentage difference is: "+cbl::conv((Pk_CAMB[ii][50]-Pk_CLASS[ii][50])/Pk_CAMB[ii][50]*100., cbl::par::fDP2) << "%"<< endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
} 

