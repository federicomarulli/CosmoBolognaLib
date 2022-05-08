// ========================================================================================
// Example code: how to measure the connected and reduced three-point correlation functions
// ========================================================================================

#include "Modelling_ThreePointCorrelation_comoving_reduced.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
std::string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------

    cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    std::string file_catalogue = cbl::par::DirLoc+"../input/cat.dat";

    cbl::catalogue::Catalogue catalogue {cbl::catalogue::ObjectType::_Galaxy_, cbl::CoordinateType::_observed_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 2.; // random/data ratio
   
    cbl::catalogue::Catalogue random_catalogue {cbl::catalogue::RandomType::_createRandom_box_, catalogue, N_R};

  
    // -------------------------------------------------------------------------------
    // ---------------- measure the three-point correlation functions ----------------
    // -------------------------------------------------------------------------------

    // binning parameters

    const double side_s = 20.;  // 1st side of the triangle
    const double side_u = 2.;   // ratio between the 1st and 2nd sides of the triangle (u*s)
    const double perc = 0.0225; // tolerance
    const int nbins = 5;       // number of bins

  
    // output data
  
    const std::string dir_output = cbl::par::DirLoc+"../output/";
    const std::string dir_triplets = dir_output;
    const std::string dir_2pt = dir_output;
    const std::string file_output = "3pt.dat";

  
    // measure the connected and reduced three-point correlation functions and write the output

    const auto ThreeP = cbl::measure::threept::ThreePointCorrelation::Create(cbl::measure::threept::ThreePType::_comoving_reduced_, catalogue, random_catalogue, cbl::triplets::TripletType::_comoving_theta_, side_s, side_u, perc, nbins);

    ThreeP->measure(dir_triplets, dir_2pt);
  
    ThreeP->write(dir_output, file_output, 1);


    // --------------------------------------------------------------
    // ---------------- read Q dark matter (DEMNUNI) ----------------
    // --------------------------------------------------------------

    const std::string file_Q = cbl::par::DirLoc+"../input/zeta_lin_DM_z1.1_u2s5.00.dat";
    std::ifstream fin(file_Q); cbl::checkIO(fin, file_Q);

    double theta, Q, err;
    std::vector<double> Q_DM;
    while (fin >> theta >> Q >> err) 
        Q_DM.emplace_back(Q);

    fin.clear(); fin.close();
    

    // ------------------------------------------------------------------------------------------------------------
    // ----------------- model the reduced three-point correlation function to estimate b1 and b2 -----------------
    // ------------------------------------------------------------------------------------------------------------

    // object used for modelling
    cbl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced model_threep(ThreeP); 

    // set the data used to construct the model
    model_threep.set_data_model(Q_DM);
  
    // set the priors and the model
    const cbl::statistics::PriorDistribution b1_prior {cbl::glob::DistributionType::_Uniform_, 0., 5., 43142}; // flat prior for b1
    const cbl::statistics::PriorDistribution b2_prior {cbl::glob::DistributionType::_Uniform_, -3., 3., 4342}; // flat prior for b2
    model_threep.set_model_nonlinear_localbias(b1_prior, b2_prior);

  
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    // minimum and maxium scales used in the fit
    const double theta_min = 0.;
    const double theta_max = 1.;
    model_threep.set_fit_range(theta_min, theta_max);

    const int chain_size = 1000; // the size the chain lenght
    const int nwalkers = 10;     // the number of parallel walkers in the MCMC chains
    const int seed = 666;        // the base seed for initialization

    // set the likelihood type
    model_threep.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_);

    // run the MCMC method to sample the posterior
    model_threep.sample_posterior(chain_size, nwalkers, seed);

    const int burn_in = 100; // discard the first 100 chain steps 
    const int thin = 10;     // take 1 step every 10
    
    // write the results on screen 
    model_threep.show_results(burn_in, thin, seed);

    // store the results in file
    model_threep.write_results(dir_output, "model_Q", burn_in, thin);

  }

  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

