// =====================================================
// Example code: how to model the angular power spectrum 
// =====================================================

#include "LCDM.h"
#include "Data1D.h"
#include "ModelFunction_PowerSpectrum_Angular.h"
#include "Modelling_PowerSpectrum_Angular.h"

using namespace std;


int main () {

  try {
    
    // ---------------------------------------------------
    // ---------------- set the Cosmology ----------------
    // ---------------------------------------------------
    
    auto cosmology = make_shared<cbl::cosmology::LCDM>("Planck18");
    cosmology->set_parameter("sigma8", 0.83);

    
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
    
    string file_catalogue = "../input/catalogue.csv";
    vector<cbl::catalogue::Var> attribute = {cbl::catalogue::Var::_RA_, cbl::catalogue::Var::_Dec_, cbl::catalogue::Var::_Redshift_};
    vector<int> column = {1, 2, 3};
    const cbl::CoordinateUnits angularUnits = cbl::CoordinateUnits::_radians_;
    
    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, attribute, column, {file_catalogue}, 0, 1, 1, cosmology, angularUnits};  

    
    // -----------------------------------------------------------------------------------------------------------
    // ------------------- construct the dataset by reading an input file ----------------------------------------
    // -----------------------------------------------------------------------------------------------------------

    const string file_measure = "../input/Cl_spherical_harmonic.dat";
    const int skipped_lines = 1; 
    const auto pow_dataset = make_shared<cbl::data::Data1D>(cbl::data::Data1D(file_measure, skipped_lines));
    
    
    // ------------------------------------------------------
    // ------------- set the modelling structure ------------
    // ------------------------------------------------------

    cbl::modelling::angularpk::Modelling_PowerSpectrum_angular model_pow(pow_dataset);
    const vector<double> dN_par = {0.5, 1., -0.2, 3.}; // parameters of the dN/dz distribution, not normalized
    model_pow.set_data_model(cosmology, catalogue.Min(cbl::catalogue::Var::_Redshift_),catalogue.Max(cbl::catalogue::Var::_Redshift_), "EisensteinHu", true, 0, 0.001, 100., dN_par);
    
    
    // -------------------------------------------------------
    // ------------- set the priors and the model ------------
    // -------------------------------------------------------

    vector<std::string> CosmoPar = {"Omega_CDM", "sigma8"};
    cbl::statistics::PriorDistribution Omega_CDM_prior, sigma8_prior;
    
    Omega_CDM_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Constant_, cosmology->Omega_CDM());
    sigma8_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Constant_, cosmology->sigma8());
    
    vector<cbl::statistics::PriorDistribution> CosmoPriors = {Omega_CDM_prior, sigma8_prior}; 
    cbl::statistics::PriorDistribution bias_prior(cbl::glob::DistributionType::_Uniform_, 0., 100.); 
    cbl::statistics::PriorDistribution shot_noise_prior(cbl::glob::DistributionType::_Constant_, 0.); 

    model_pow.set_model(CosmoPar, CosmoPriors, bias_prior, shot_noise_prior); 
    
    const double xmin = 10., xmax = 30.; 
    model_pow.set_fit_range(xmin, xmax);
    vector<double> ell = cbl::linear_bin_vector(xmax-xmin, xmin, xmax);
    
    model_pow.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);
    

    // -------------------------------------------
    // ------------- fit the data set ------------
    // -------------------------------------------

    model_pow.maximize_posterior({1.}, 10000, 1.e-4);
    
    const string file_output = "model_bestfit.dat";
    model_pow.write_model_at_bestfit("../output/", file_output, ell);
  
  }
  
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}
  
