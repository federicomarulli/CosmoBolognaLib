// =============================================================
// Example code: how to set and use different cosmologica models
// =============================================================

#include "Cosmologies.h"

using namespace std;


int main () {

  try {
    
    // ---------------------------------------------------------------------
    // ----- set various cosmological models, with Planck18 parameters -----
    // ---------------------------------------------------------------------

    // flat Lambda-cold dark matter model
    auto LCDM = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    // flat Lambda-cold dark matter model with massive neutrinos
    auto LnuCDM = make_shared<cbl::cosmology::LnuCDM>("Planck18");

    // flat w-cold dark matter model
    auto w0CDM = make_shared<cbl::cosmology::w0CDM>("Planck18");
    
    // non-flat Lambda-cold dark matter model
    auto kLCDM = make_shared<cbl::cosmology::kLCDM>("Planck18");

    // flat w-cold dark matter model with massive neutrinos
    auto w0nuCDM = make_shared<cbl::cosmology::w0nuCDM>("Planck18");
    
    // flat w0wa-cold dark matter model
    auto w0waCDM = make_shared<cbl::cosmology::w0waCDM>("Planck18");
    
    // non-flat Lambda-cold dark matter model with massive neutrinos
    auto kLnuCDM = make_shared<cbl::cosmology::kLnuCDM>("Planck18");

    // non-flat w-cold dark matter model
    auto kw0CDM = make_shared<cbl::cosmology::kw0CDM>("Planck18");

    // flat w0wa-cold dark matter model with massive neutrinos
    auto w0wanuCDM = make_shared<cbl::cosmology::w0wanuCDM>("Planck18");

    // non-flat w-cold dark matter model with massive neutrinos
    auto kw0nuCDM = make_shared<cbl::cosmology::kw0nuCDM>("Planck18");
    
    // non-flat w0wa-cold dark matter 
    auto kw0waCDM = make_shared<cbl::cosmology::kw0waCDM>("Planck18");
    
    // non-flat w0wa-cold dark matter model with massive neutrinos
    auto kw0wanuCDM = make_shared<cbl::cosmology::kw0wanuCDM>("Planck18");

    // vector containing all the different cosmological models
    vector<shared_ptr<cbl::cosmology::Cosmology>> cosmology = {LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM};

    
    // -------------------------------------------------------------------------------------
    // ----- print the model names and the values of the cosmological model parameters -----
    // -------------------------------------------------------------------------------------

    for (size_t i=0; i<cosmology.size(); ++i) {
      cout << endl << cosmology[i]->model() << endl;
      cosmology[i]->print_parameters();
    }
    
  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
} 

