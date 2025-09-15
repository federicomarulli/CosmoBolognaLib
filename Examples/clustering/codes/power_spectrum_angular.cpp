// =======================================================
// Example code: how to measure the angular power spectrum 
// =======================================================

#include "LCDM.h"
#include "PowerSpectrum_Angular.h"
#include "GlobalFunc.h"

using namespace std;


int main () {
  
  try {

    // --------------------------------------------------------------------
    // ---------------- set the cosmological model to Planck18 ------------
    // --------------------------------------------------------------------

    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    
    
    // -----------------------------------------------------------------------------------------------------------------
    // -------------------- read the input catalogue (in observed coordinates, RA, Dec, Redshift) ----------------------
    // -----------------------------------------------------------------------------------------------------------------

    string file_catalogue = "../input/catalogue.csv";   
    vector<cbl::catalogue::Var> attribute = {cbl::catalogue::Var::_RA_, cbl::catalogue::Var::_Dec_, cbl::catalogue::Var::_Redshift_};
    vector<int> column = {1, 2, 3};
    const cbl::CoordinateUnits angularUnits=cbl::CoordinateUnits::_radians_;
    
    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, attribute, column, {file_catalogue}, 0, 1, 1, cosmology, angularUnits};
    catalogue = catalogue.sub_catalogue(cbl::catalogue::Var::_Redshift_, 0.9, 1.1);
												 
    
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------
       
    const double N_R = 3.; // random/data ratio
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_square_, catalogue, N_R, 10, cosmology};
    
    cout << "I'm constructing the sub-regions used for jackknife and bootstrap..." << endl;
    const int nRa = 3, nDec = 3;
    cbl::set_ObjectRegion_RaDec(catalogue, random_catalogue, nRa, nDec);

    
    // --------------------------------------------------------------------------------------------------
    // ------------------ angular power spectrum: binning parameters and output data ---------------------
    // --------------------------------------------------------------------------------------------------
    
    double ell_min = 18.; // minimum multipole, starts from 1 in fast estimator
    double ell_max = 26.; // maximum multipole
    const int Nell = 8;   // number of multipoles
    const string dir_output = "../output/";

    
    // --------------------------------------------------------------------------------------------------
    // ------ METHOD 1: fast estimator, using angular two-point correlation function --------------------
    // --------------------------------------------------------------------------------------------------

    cbl::measure::angularpk::AngularEstimator estimator = cbl::measure::angularpk::AngularEstimator::_Fast_;
    const string dir_correlation_output = "../../clustering/output/"; 
    const string file_correlation_output = "xi_angular.dat"; 
    const cbl::CoordinateUnits correlation_angular_Units = cbl::CoordinateUnits::_arcminutes_;
    cbl::BinType power_spectrum_bintype = cbl::BinType::_linear_;                                  

    const string dir_correlation_input = "";  // input directory (leave empty if you want to measure w(theta))
    const string file_correlation_input = ""; // input file (leave empty if you want to measure w(theta))
    const int n_lines_header = 1;
    const cbl::CoordinateUnits inputUnits = cbl::CoordinateUnits::_arcminutes_;  // angular units of the input file
       
    const double angMin = 1;     // minimum angular separation 
    const double angMax = 1000;  // maximum angular separation 
    const int nbins = 100;       // number of bins
    const double shift = 0.5;    // shift used to set the bin centre 
    cbl::BinType correlation_bintype = cbl::BinType::_logarithmic_;

    
    // measure the angular power spectrum and store the results

    cbl::measure::angularpk::PowerSpectrum_angular Pow{catalogue, random_catalogue, ell_min, ell_max, Nell, power_spectrum_bintype, angMin, angMax, nbins, shift, correlation_angular_Units, cbl::BinType::_linear_};
    
    Pow.measure(estimator, dir_correlation_input, file_correlation_input, n_lines_header, inputUnits, correlation_bintype, cbl::measure::ErrorType::_Jackknife_, dir_correlation_output, file_correlation_output);
    
    Pow.write(dir_output, "Cl_fast.dat");

    
    // --------------------------------------------------------------------------------------------------
    // ------------------------- METHOD 2: spherical armonic estimator ----------------------------------
    // --------------------------------------------------------------------------------------------------

    estimator = cbl::measure::angularpk::AngularEstimator::_SphericalArmonic_; 
    cbl::measure::angularpk::PowerSpectrum_angular Pow2{catalogue, ell_min, ell_max};
    
    Pow2.measure(estimator);
    Pow2.write(dir_output, "Cl_spherical_harmonic.dat");
  }
  
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}
