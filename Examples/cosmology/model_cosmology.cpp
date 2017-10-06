// ==================================================
// Example code: how to model cosmological parameters
// Data from Addison et al. 2013
// ==================================================

#include "Data1D.h"
#include "Modelling_Cosmology.h"

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cosmobl::par::DirCosmo = DIRCOSMO, cosmobl::par::DirLoc = DIRL;

int main () {

  try {

    // ---------------------------------------------------------------------------
    // ---------------- using one of the built-in cosmological models ------------
    // ---------------------------------------------------------------------------

    cosmobl::cosmology::Cosmology cosmology {cosmobl::cosmology::_Planck15_, "LCDM", false};


    // ----------------------------------------------
    // ---------------- Read the dataset ------------
    // ---------------------------------------------- 

    string data_file = cosmobl::par::DirCosmo+"/External/Data/BAO/BAO_Addison2013.dat";
    string covariance_file = cosmobl::par::DirCosmo+"/External/Data/BAO/BAO_Addison2013_covariance.dat";

    auto data = make_shared<cosmobl::data::Data1D>(cosmobl::data::Data1D(data_file, 1));
    data->set_covariance(covariance_file, 2, 1);


    // ------------------------------------------------
    // ---------------- Read the data type ------------
    // ------------------------------------------------

    vector<string> data_type;
    ifstream fin(data_file.c_str());
    string line;
    
    // skip header
    getline(fin,line);

    while(getline(fin, line)) {
      double A;
      string dt;
      stringstream ss(line);
      ss >> A >> A >> A >> dt;
      data_type.push_back(dt);
    }
    fin.clear(); fin.close();


    // ------------------------------------------------------------------------------------
    // ---------------- modelling cosmological parameters from BAO information ------------
    // ------------------------------------------------------------------------------------

    cosmobl::modelling::cosmology::Modelling_Cosmology modelCosmo(data, data_type);

    
    // ---------------------------------------------------------------------
    // ---------------- define cosmological parameters to model ------------
    // ---------------------------------------------------------------------

    modelCosmo.set_fiducial_cosmology(cosmology);

    vector<cosmobl::cosmology::CosmoPar> Cpar = {cosmobl::cosmology::CosmoPar::_Omega_matter_LCDM_, cosmobl::cosmology::CosmoPar::_H0_, cosmobl::cosmology::CosmoPar::_rs_};

    cosmobl::statistics::Prior OmegaM_prior(cosmobl::glob::DistributionType::_UniformDistribution_, 0.1, 0.5);
    cosmobl::statistics::Prior H0_prior(cosmobl::glob::DistributionType::_UniformDistribution_, 50, 100);
    cosmobl::statistics::Prior rs_prior(cosmobl::glob::DistributionType::_GaussianDistribution_, {cosmology.rs_CAMB(), 10}, 130, 180);

    modelCosmo.set_cosmological_parameters(Cpar, {OmegaM_prior, H0_prior, rs_prior});

    
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    const int chain_size = 2000;
    const int nwalkers = 200;
    const int seed = 4232;
    vector<double> starting_parameters = {cosmology.Omega_matter(), cosmology.H0(), 150.};
    double radius = 1.e-3;

    modelCosmo.set_likelihood(cosmobl::statistics::LikelihoodType::_GaussianLikelihood_Covariance_);
    modelCosmo.sample_likelihood(chain_size, nwalkers, seed, starting_parameters, radius);
    modelCosmo.show_results(500, 10, 3213);
  }

  
  catch(cosmobl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
} 

