// ========================================================================================
// Example code: how to model the angular power spectrum 
// ========================================================================================

#include "Data1D.h"
#include "ModelFunction_PowerSpectrum_Angular.h"
#include "Modelling_PowerSpectrum_Angular.h"

using namespace std;
using namespace cbl;


// =====================================================================


int main () {

  try {
    
    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};
    cosmology.set_sigma8(0.83);
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
    
    std::string file_catalogue = "../input/catalogue.csv";
    std::vector<cbl::catalogue::Var> attribute = {cbl::catalogue::Var::_RA_, cbl::catalogue::Var::_Dec_, cbl::catalogue::Var::_Redshift_};
    std::vector<int> column = {1,2,3};
    const cbl::CoordinateUnits angularUnits=cbl::CoordinateUnits::_radians_;
    
    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, attribute, column, {file_catalogue}, 0, 1, 1, cosmology, angularUnits};

    // ------------------------------------------------------------------------------------------------------------------------------------------
    // ---- insert mean of the offset and slope of the normalized redshift distribution (or measure them from catalogue number count!) ----------
    // ------------------------------------------------------------------------------------------------------------------------------------------
    
    const double mean_boffset = 208.84611; 
    const double mean_bslope = -175.71603;    

    // -----------------------------------------------------------------------------------------------------------
    // ------------------- construct the dataset by reading an input file ----------------------------------------
    // -----------------------------------------------------------------------------------------------------------

    const std::string file_measure = "../input/Cl_spherical_harmonic.dat";
    const int skipped_lines = 1;   //line header
    const auto pow_dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(file_measure, skipped_lines));
   
    // ------------------------------------------------------
    // ------------- set the modelling structure ------------
    // ------------------------------------------------------
    
    cbl::modelling::angularpk::Modelling_PowerSpectrum_angular model_pow(pow_dataset);
    model_pow.set_data_model(cosmology, catalogue.Min(cbl::catalogue::Var::_Redshift_),catalogue.Max(cbl::catalogue::Var::_Redshift_), "EisensteinHu", true, 0, 0.001, 100., {mean_boffset, mean_bslope});
    // -------------------------------------------------------
    // ------------- set the priors and the model ------------
    // -------------------------------------------------------

    std::vector<cbl::cosmology::CosmologicalParameter> CosmoPar={cbl::cosmology::CosmologicalParameter::_Omega_matter_LCDM_, cbl::cosmology::CosmologicalParameter::_sigma8_};
    cbl::statistics::PriorDistribution Omega_matter_LCDM_prior, sigma8_prior;

    Omega_matter_LCDM_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Constant_, cosmology.Omega_matter());
    sigma8_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Constant_, cosmology.sigma8());
    
    std::vector<cbl::statistics::PriorDistribution> CosmoPriors = {Omega_matter_LCDM_prior, sigma8_prior}; //,OmegaB_prior
    cbl::statistics::PriorDistribution bias_prior(cbl::glob::DistributionType::_Uniform_, 0., 100.); 

    model_pow.set_model_limber(CosmoPar, CosmoPriors, bias_prior);
    
    const double xmin = 10., xmax = 30.; 
    model_pow.set_fit_range(xmin, xmax);
    std::vector<double> ell= cbl::linear_bin_vector(xmax-xmin, xmin, xmax);
    
    model_pow.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Covariance_);
    // maximise the posterior
    std::vector<double> start= {1.};
    model_pow.maximize_posterior(start, 10000, 1.e-5);
    
    const string file_output = "model_bestfit.dat";
    model_pow.write_model_at_bestfit("../output/",file_output, ell);
  
  }
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  
  return 0;
}
  
