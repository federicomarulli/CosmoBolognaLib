// ==================================================
// Example code: how to model cosmological parameters
// Data from Addison et al. 2013
// ==================================================

#include "LCDM.h"
#include "Data1D.h"
#include "Modelling_Cosmology.h"

int main () {

  try {

    // --------------------------------------------------------
    // ---------------- set the cosmological model ------------
    // --------------------------------------------------------

    auto cosmology = std::make_shared<cbl::cosmology::LCDM>("Planck18");


    // ----------------------------------------------
    // ---------------- read the dataset ------------
    // ---------------------------------------------- 

    cbl::Path path;
    const std::string data_file = path.DirCosmo()+"/External/Data/BAO/BAO_Addison2013.dat";
    const std::string covariance_file = path.DirCosmo()+"/External/Data/BAO/BAO_Addison2013_covariance.dat";

    auto data = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(data_file, 1));
    data->set_covariance(covariance_file, 2, 1);


    // ------------------------------------------------
    // ---------------- read the data type ------------
    // ------------------------------------------------

    std::vector<std::string> data_type;
    std::ifstream fin(data_file.c_str());
    std::string line;
    
    // skip header
    getline(fin,line);

    while(getline(fin, line)) {
      double A;
      std::string dt;
      std::stringstream ss(line);
      ss >> A >> A >> A >> dt;
      data_type.push_back(dt);
    }
    fin.clear(); fin.close();


    // ------------------------------------------------------------------------------------
    // ---------------- modelling cosmological parameters from BAO information ------------
    // ------------------------------------------------------------------------------------

    cbl::modelling::cosmo::Modelling_Cosmology modelCosmo(data, data_type);

    
    // ---------------------------------------------------------------------
    // ---------------- define cosmological parameters to model ------------
    // ---------------------------------------------------------------------

    modelCosmo.set_fiducial_cosmology(cosmology);


    std::vector<std::string> Cpar = {"Omega_matter_LCDM", "H0", "rs"};

    cbl::statistics::PriorDistribution OmegaM_prior(cbl::glob::DistributionType::_Uniform_, 0.1, 0.5, 5452);
    cbl::statistics::PriorDistribution H0_prior(cbl::glob::DistributionType::_Uniform_, 50, 100, 6764);
    cbl::statistics::PriorDistribution rs_prior(cbl::glob::DistributionType::_Gaussian_, {cosmology->sound_horizon_at_drag_epoch(), 10}, 130, 180, 5645);

    modelCosmo.set_cosmological_parameters(Cpar, {OmegaM_prior, H0_prior, rs_prior});

    
    // ----------------------------------------------------------------------
    // ------------- run chains and write output chain and model ------------
    // ----------------------------------------------------------------------

    const int chain_size = 100;
    const int nwalkers = 10;
    const int seed = 4232;
    std::vector<double> starting_parameters = {cosmology->Omega_matter(), cosmology->Hubble(0.), 150.};

    modelCosmo.set_fit_range(cbl::Min(data->xx()), cbl::Max(data->xx()));
    
    modelCosmo.set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Covariance_);
    
    modelCosmo.sample_posterior(chain_size, nwalkers, seed);

    modelCosmo.show_results();
  }

  
  catch(cbl::glob::Exception &exc) { std::cerr << exc.what() << std::endl; exit(1); }

  return 0;
} 

