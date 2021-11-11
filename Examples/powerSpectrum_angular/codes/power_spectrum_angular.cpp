// =======================================================================
// Example code: how to measure the angular power spectrum 
// =======================================================================

#include "PowerSpectrum_Angular.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;

int main () {
  
  try {

    // --------------------------------------------------------------------
    // ---------------- set the cosmological model to Planck15 ------------
    // --------------------------------------------------------------------
    
    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

    
    // -----------------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec) --------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    
    std::string file_catalogue = cbl::par::DirLoc+"../input/data_0.4_0.6";
    std::vector<cbl::catalogue::Var> attribute = {cbl::catalogue::Var::_RA_, cbl::catalogue::Var::_Dec_};
    std::vector<int> column = {1,2};
    const cbl::CoordinateUnits angularUnits=cbl::CoordinateUnits::_degrees_;
    
    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, attribute, column, {file_catalogue}, 0, 1, 1, cosmology, angularUnits};

    
    // ----------------------------------------------------------------
    // ---------------- construct the random catalogue ----------------
    // ----------------------------------------------------------------
    
    const double N_R = 2.; // random/data ratio
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_square_, catalogue, N_R};

    
    // --------------------------------------------------------------------------------------------------
    // ------ angular two-point correlation function: binning parameters and output data ----------------
    // --------------------------------------------------------------------------------------------------
    
    const double angMin = 0.02;                                                  // minimum angular separation 
    const double angMax = 30;                                                    // maximum angular separation 
    const int nbins = 20;                                                        // number of bins
    const double shift = 0.5;                                                    // shift used to set the bin centre 
    const std::string dir_correlation_output=cbl::par::DirLoc+"../../clustering/output/";       //output directory
    const std::string file_correlation_output="xi_angular.dat";                  //output file
    const cbl::CoordinateUnits correlation_angular_Units = cbl::CoordinateUnits::_arcminutes_;  // angular units

    
    // --------------------------------------------------------------------------------------------------
    // ------------------ angular power spectrum: binning parameters and input data ---------------------
    // --------------------------------------------------------------------------------------------------
    
    cbl::measure::angularpk::Estimator estimator =cbl::measure::angularpk::Estimator::_Fast_;          //or _SphericalArmonic_
    const double ell_min=100;                                                    //minimum multipole
    const double ell_max=1000;                                                   //maximum multipole
    const int Nell=9;                                                            //number of multipoles
    const std::string dir_correlation_input="";                                  //input directory (leave empty if you want to measure w(theta))
    const std::string file_correlation_input="";                                 //input file (leave empty if you want to measure w(theta))
    const int n_lines_header=1;
    const cbl::CoordinateUnits inputUnits = cbl::CoordinateUnits::_arcminutes_;  // angular units of the input file

    const std::string dir_output = cbl::par::DirLoc+"../output/";
    const std::string file_output = "Cell.dat";

    
    // measure the angular power spectrum and store the results

    cbl::measure::angularpk::PowerSpectrum_angular Pow{catalogue, random_catalogue, ell_min, ell_max, Nell, cbl::BinType::_logarithmic_, angMin, angMax, nbins, shift, correlation_angular_Units};
    
    Pow.measure(estimator, dir_correlation_input, file_correlation_input, n_lines_header, inputUnits, cbl::BinType::_logarithmic_, cbl::measure::ErrorType::_Poisson_,  dir_correlation_output, file_correlation_output);
    
    Pow.write(dir_output, file_output);
 
  }
  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}
