// ========================================================================================
// Example code: how to model the multipoles of the two-point correlation in redshift space
// ========================================================================================

#include "TwoPointCorrelation_multipoles_direct.h"
#include "TwoPointCorrelation_multipoles_integrated.h"
#include "Modelling_TwoPointCorrelation_multipoles.h"

int main () {

  try {
  
    // --------------------------------------------------------------------
    // ---------------- set the cosmological model to Planck15 ------------
    // --------------------------------------------------------------------

    const cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    const std::string file_catalogue = "../input/cat.dat";

    const cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 10.; // random/data ratio
  
    const cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};

  
    // --------------------------------------------------------------------------------------------
    // ---------------- measure the multipoles of the two-point correlation function --------------
    // --------------------------------------------------------------------------------------------

    // output directory
    const std::string dir = "../output/";

    // binning parameters and output data
    const double rMin = 10.;  // minimum separation 
    const double rMax = 30.;  // maximum separation 
    const int nbins = 5;      // number of bins
    const double shift = 0.5; // spatial shift used to set the bin centre
  
    // measure the multipoles (using the direct estimator)
    cbl::measure::twopt::TwoPointCorrelation_multipoles_direct TwoP_direct {catalogue, random_catalogue, cbl::BinType::_logarithmic_, rMin, rMax, nbins, shift};
    TwoP_direct.measure(cbl::measure::ErrorType::_Poisson_, dir);
    
    // store the output data
    TwoP_direct.write(dir, "xil_direct.dat");

    
    // ------------------------------------------------------------------------------------------
    // ---------------- model the multipoles of the two-point correlation function --------------
    // ------------------------------------------------------------------------------------------

    // object used for modelling and set the data used to construct the model
    auto ptr_TwoP = std::make_shared<cbl::measure::twopt::TwoPointCorrelation_multipoles_direct>(TwoP_direct);
    cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles model_multipoles(ptr_TwoP);

    // set the priors
    const cbl::statistics::PriorDistribution alpha_perpendicular_prior {cbl::glob::DistributionType::_Constant_, 1.}; // constant prior 
    const cbl::statistics::PriorDistribution alpha_parallel_prior {cbl::glob::DistributionType::_Constant_, 1.};      // constant prior 
    const cbl::statistics::PriorDistribution SigmaNL_perpendicular_prior {cbl::glob::DistributionType::_Constant_, 0.}; // constant prior 
    const cbl::statistics::PriorDistribution SigmaNL_parallel_prior {cbl::glob::DistributionType::_Constant_, 0.};  // constant prior 
    const cbl::statistics::PriorDistribution fsigma8_prior {cbl::glob::DistributionType::_Uniform_, 0., 2.}; // Uniform prior for the f*sigma8
    const cbl::statistics::PriorDistribution bsigma8_prior {cbl::glob::DistributionType::_Uniform_, 0., 2.}; // Uniform prior for the b*sigma8
    const cbl::statistics::PriorDistribution SigmaS_prior {cbl::glob::DistributionType::_Uniform_, 0., 10.}; // Uniform prior for the SigmaS

    // provide the cosmological model and the redshift    
    const double redshift = 1.;
    model_multipoles.set_data_model(cosmology, redshift);

    // set the scale range for the fit    
    const double xmin = 10., xmax = 50.;
    model_multipoles.set_fit_range(xmin, xmax, 3);

    // set the model for the full-shape analyses fo the clustering multipoles    
    model_multipoles.set_model_fullShape_DeWiggled(alpha_perpendicular_prior, alpha_parallel_prior, SigmaNL_perpendicular_prior, SigmaNL_parallel_prior, fsigma8_prior, bsigma8_prior, SigmaS_prior);
    
    // set the likelihood type
    model_multipoles.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);

    // maximise the posterior
    const std::vector<double> start = {0.2, 0.2, 1.};
    model_multipoles.maximize_posterior(start);

    // retrieve and show the best-fit parameter values
    for (unsigned int i=0; i<model_multipoles.posterior()->parameters()->nparameters(); ++i)
      std::cout << "the best-fit value of " << model_multipoles.posterior()->parameters()->name(i) << " is " << model_multipoles.posterior()->parameters()->bestfit_value(i) << std::endl;

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

