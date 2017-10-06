// ========================================================================================
// Example code: how to measure the connected and reduced three-point correlation functions
// ========================================================================================

#include "Modelling_ThreePointCorrelation_comoving_reduced.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;


int main () {

  try {
  
    // --------------------------------------------------------------
    // ---------------- set the cosmological parameters  ------------
    // --------------------------------------------------------------

    cosmobl::cosmology::Cosmology cosmology {cosmobl::cosmology::_Planck15_};

  
    // -----------------------------------------------------------------------------------------------------------
    // ---------------- read the input catalogue (with observed coordinates: R.A., Dec, redshift) ----------------
    // -----------------------------------------------------------------------------------------------------------
  
    string file_catalogue = cosmobl::par::DirLoc+"../input/cat.dat";

    cosmobl::catalogue::Catalogue catalogue {cosmobl::catalogue::_Galaxy_, cosmobl::_observedCoordinates_, {file_catalogue}, cosmology};

  
    // --------------------------------------------------------------------------------------
    // ---------------- construct the random catalogue (with cubic geometry) ----------------
    // --------------------------------------------------------------------------------------

    const double N_R = 1.; // random/data ratio
   
    cosmobl::catalogue::Catalogue random_catalogue {cosmobl::catalogue::_createRandom_box_, catalogue, N_R};

  
    // -------------------------------------------------------------------------------
    // ---------------- measure the three-point correlation functions ----------------
    // -------------------------------------------------------------------------------

    // binning parameters

    const double side_s = 5.;  // 1st side of the triangle
    const double side_u = 2.;   // ratio between the 1st and 2nd sides of the triangle (u*s)
    const double perc = 0.0225; // tolerance
    const int nbins = 20;       // number of bins

  
    // output data
  
    const string dir_output = cosmobl::par::DirLoc+"../output/";
    const string dir_triplets = dir_output;
    const string dir_2pt = dir_output;
    const string file_output = "3pt.dat";

  
    // measure the connected and reduced three-point correlation functions and write the output

    const auto ThreeP = cosmobl::measure::threept::ThreePointCorrelation::Create(cosmobl::measure::threept::_comoving_reduced_, catalogue, random_catalogue, cosmobl::triplets::_comoving_theta_, side_s, side_u, perc, nbins);

    ThreeP->measure(dir_triplets, dir_2pt);
  
    ThreeP->write(dir_output, file_output, 1);


    // --------------------------------------------------------------
    // ---------------- read Q dark matter (DEMNUNI) ----------------
    // --------------------------------------------------------------

    const string file_Q = cosmobl::par::DirLoc+"../input/zeta_lin_DM_z1.1_u2s5.00.dat";
    ifstream fin(file_Q); cosmobl::checkIO(fin, file_Q);

    double theta, Q, err;
    vector<double> Q_DM;
    while (fin >> theta >> Q >> err) 
        Q_DM.emplace_back(Q);

    fin.clear(); fin.close();

    // ------------------------------------------------------------------------------------------------------------
    // ----------------- model the reduced three-point correlation function to estimate b1 and b2 -----------------
    // ------------------------------------------------------------------------------------------------------------

    // object used for modelling
    cosmobl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced model_threep(ThreeP); 

    // set the data used to construct the model
    model_threep.set_data_model(Q_DM);
  
    // set the priors and the model
    const cosmobl::statistics::Prior b1_prior {cosmobl::glob::DistributionType::_UniformDistribution_, 0., 5.}; // flat prior for b1
    const cosmobl::statistics::Prior b2_prior {cosmobl::glob::DistributionType::_UniformDistribution_, -3., 3.}; // flat prior for b2
    model_threep.set_model_nonlinear_localbias(b1_prior, b2_prior);

  
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    // minimum and maxium scales used in the fit
    //const double theta_min = 0.;
    //const double theta_max = 1.;
    //model_threep.set_fit_range(theta_min, theta_max);

    const int chain_size = 10000; // the size the chain lenght
    const int nwalkers = 200;     // the number of parallel walkers in the MCMC chains
    const int seed = 4232;        // the base seed for initialization

    // set the likelihood type
    model_threep.set_likelihood(cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Error_);

    // sample the likelihood
    model_threep.sample_likelihood(chain_size, nwalkers, seed);

    const int burn_in = 100; // discard the first 100 chain steps 
    const int thin = 10;     // take 1 step every 10

    // write the results on screen 
    model_threep.show_results(burn_in, thin, seed);

    // store the results in file
    model_threep.write_results(dir_output, "model_Q", burn_in, thin, seed);





  }

  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }
  
  return 0;
}

