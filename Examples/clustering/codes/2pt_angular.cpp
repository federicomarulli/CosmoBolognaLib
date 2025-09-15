// =======================================================================
// Example code: how to measure the angular two-point correlation function 
// =======================================================================

#include "LCDM.h"
#include "TwoPointCorrelation1D_angular.h"

using namespace std;


int main () {

  try {
    
    // --------------------------------------------------------------------
    // ---------------- set the cosmological model to Planck18 ------------
    // --------------------------------------------------------------------

    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");


    // -----------------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec) --------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    const string file_catalogue = "../input/cat2d.dat";
    
    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};
    

  
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_square_, catalogue, N_R, 10, cosmology};
    
  
    // ------------------------------------------------------------------------------------
    // ---------------- measure the angular two-point correlation function ----------------
    // ------------------------------------------------------------------------------------

    // binning parameters and output data

    const double angMin = 0.01;                                              // minimum angular separation 
    const double angMax = 1.;                                                // maximum angular separation 
    const int nbins = 20;                                                    // number of bins
    const double shift = 0.5;                                                // shift used to set the bin centre 
    const cbl::CoordinateUnits angularUnits = cbl::CoordinateUnits::_degrees_; // angular units

    const string dir = "../output/";
    const string file = "xi_angular.dat";

  
    // measure the angular two-point correlation function and store the results
  
    cbl::measure::twopt::TwoPointCorrelation1D_angular TwoP {catalogue, random_catalogue, cbl::BinType::_linear_, angMin, angMax, nbins, shift, angularUnits};
  
    TwoP.measure(cbl::measure::ErrorType::_Poisson_, dir);

    TwoP.write(dir, file);

  }

  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}

