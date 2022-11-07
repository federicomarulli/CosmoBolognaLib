// ==================================================================
// Example code: how to measure and model the void size function,
// extracting constraints on the cosmological parameters of the model
// ==================================================================

#include "Cosmology.h"
#include "Catalogue.h"
#include "Posterior.h"
#include "NumberCounts1D_Size.h"
#include "Modelling_NumberCounts1D_Size.h"

using namespace std;


int main () {

  try {

    // --- set the input/output file/directories ---
    const std::string dir = "../output/";
    const string file_output_start = "model_starting_values.dat";
    const string file_output_bestfit = "model_bestfit.dat";


    // ------------------------------------------
    // ----- load the input void catalogue ------
    // ------------------------------------------

    // ASCII void catalogue 
    std::string file_voids_in = "../input/cleaned_void_catalogue.out";

    // std::vector containing the variable name list to read from file
    std::vector<cbl::catalogue::Var> var_names_voids = {cbl::catalogue::Var::_X_, cbl::catalogue::Var::_Y_, cbl::catalogue::Var::_Z_, cbl::catalogue::Var::_Radius_};
    
    // std::vector containing the columns corresponding to each attribute
    std::vector<int> columns_voids = {1, 2, 3, 4};
    
    // catalogue constructor
    cbl::catalogue::Catalogue void_catalogue {cbl::catalogue::ObjectType::_Void_, cbl::CoordinateType::_comoving_, var_names_voids, columns_voids, {file_voids_in}, 0};

    // ------------------------------------------------------------------
    // ---------------- measure the void size function ------------------
    // ------------------------------------------------------------------

    // binning parameters and o
    const int nbin = 8;
    const double rmin = 23.5;
    const double rmax = 35.5;
    const double shift = 0.5;
    const double vol = 480.*480.*480.; // cut borders
    
    // measure the void number counts and compute Poisson errors
    cbl::measure::numbercounts::NumberCounts1D_Size NC {void_catalogue, nbin, rmin, rmax, shift, cbl::glob::HistogramType::_dn_dlnV_, vol};
    shared_ptr<cbl::measure::numbercounts::NumberCounts1D> ptr_NC = make_shared<cbl::measure::numbercounts::NumberCounts1D>(NC);

    
    // measure the number counts and compute Poissonian errors 
    ptr_NC->measure(cbl::measure::ErrorType::_Poisson_, "../output/NumberCounts.out", 0, 3213);
    ptr_NC->write(dir, "Size_distribution_Poisson.dat");


    // ------------------------------------------------------------------
    // ------------------ model the void size function ------------------
    // ------------------------------------------------------------------

    
    // --- set the model to construct the likelihood ---

    // names of the cosmological model parameters
    const vector<string> cosmo_parNames = {"sigma8", "Omega_matter"};
    
    // starting values
    double val_s8 = 0.809;
    double val_Om = 0.2711;

    // set the cosmological parameters
    const double OmegaM = 0.2711;
    const double Omega_b = 0.0451;
    const double Omega_nu = 0.;
    const double massless_neutrinos = 3.04;
    const int massive_neutrinos = 0;
    const double OmegaL = 0.7289;
    const double Omega_radiation = 0.;
    const double hh = 0.703;
    const double scalar_amp = 2.194e-9;
    const double scalar_pivot = 0.05;
    const double n_s =  0.96;
    const double w0 = -1.;
    const double wa = 0.;
    cbl::cosmology::Cosmology cosmology(OmegaM, Omega_b, Omega_nu, massless_neutrinos, massive_neutrinos, OmegaL, Omega_radiation, hh, scalar_amp, scalar_pivot, n_s, w0, wa);
    cosmology.set_sigma8(0.809);   
    const double redshift = 0.00;
    
    vector<cbl::cosmology::CosmologicalParameter> ModelCosmoPar = cbl::cosmology::CosmologicalParameterCast(cosmo_parNames);
    
    // set the limits for the cosmological parameters
    const double min_s8 = 0.6, max_s8 = 0.9;
    const double min_Om = 0.1, max_Om = 0.4;
    const vector<vector<double>> limits = {{min_s8, max_s8}, {min_Om, max_Om}};

    // --- construct the priors ---
    auto prior_s8 = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, min_s8, max_s8, 666);
    auto prior_Om = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, min_Om, max_Om, 777);
    const vector<cbl::statistics::PriorDistribution> prior_cosmoPar = {prior_s8, prior_Om};  

    vector<double> size; // vector containing the void radii
    for (size_t i=0; i<ptr_NC->dataset()->xx().size(); i++){
      size.push_back(ptr_NC->dataset()->xx()[i]);}

    // create the void size function model using the measured void number counts
    cbl::modelling::numbercounts::Modelling_NumberCounts1D_Size model_SF(ptr_NC);

    // set the theoretical model and the parameter priors
    double val_beff = 1.1218389;               // the tracer bias
    double val_slp = 0.85676;                  // the slope of the linear relation
    double val_offs = 0.42153;                 // the offset of the linear relation
    double del_vNL = -0.7;                     // the linear underdensity threshold
    double del_c = cosmology.deltac(redshift); // the linear overdensity threshold 
    model_SF.set_data_model_SF(cosmology, size, redshift, "Vdn", val_beff, val_slp, val_offs, del_vNL, del_c);
    model_SF.set_model_NumberCounts_cosmology(ModelCosmoPar, prior_cosmoPar);

    // set the likelihood
    model_SF.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);
       
    // maximize the posterior
    model_SF.maximize_posterior({val_s8, val_Om});
    
  }

  catch(cbl::glob::Exception &exc) {cerr << exc.what() << endl; exit(1); }

  return 0;
}
