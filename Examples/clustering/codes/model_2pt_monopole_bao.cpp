// ========================================================================================================================
// Example code: how to how to model the baryon acoustic oscillations in the monopole of the two-point correlation function
// ========================================================================================================================

#include "Data1D.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation_monopole.h"
#include "GlobalFunc.h"

using namespace cosmobl;
using namespace cosmology;
using namespace catalogue;
using namespace statistics;
using namespace data;
using namespace twopt;
using namespace modelling;

string par::DirCosmo = DIRCOSMO, par::DirLoc = DIRL;


int main () {

  // ------------------------------------------------------------------------------------------
  // ---------------- set the cosmological parameters (as in Anderson et al. 2013) ------------
  // ------------------------------------------------------------------------------------------
  
  double hh = 0.7;
  double Omega_b = 0.0224/(hh*hh);
  double Omega_M = 0.274; 
  double n_s = 0.96;
  double sigma8 = 0.8;

  Cosmology cosmology;
  cosmology.set_parameter({CosmoPar::_Omega_matter_LCDM_, CosmoPar::_Omega_baryon_, CosmoPar::_hh_, CosmoPar::_n_spec_, CosmoPar::_sigma8_}, {Omega_M, Omega_b, hh, n_s, sigma8});

  double redshift = 0.57;
  cout << "Dv at redshift = " << redshift << " is " << cosmology.D_V(redshift)/hh << " Mpc" << endl;

  
  // ----------------------------------------------
  // ------------- reading the dataset ------------
  // ----------------------------------------------

  string dir_input = par::DirLoc+"../input/";
  string dir_output = par::DirLoc+"../output/";
  string dir_chains = dir_output+"chains/";

  string file_xi = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon.dat";
  string file_cov = dir_input+"Anderson_2013_CMASSDR11_monopole_prerecon_cov.dat";

  string MK = "mkdir -p "+dir_output+" "+dir_chains; if (system(MK.c_str())) {}

  auto twop_dataset = make_shared<Data1D> (Data1D(file_xi, 1)); 
  twop_dataset->set_covariance(file_cov);

  
  // ------------------------------------------------------
  // ------------- set the modelling structure ------------
  // ------------------------------------------------------
  
  auto model_twop = Modelling_TwoPointCorrelation::Create(twop_dataset, twopt::TwoPType::_1D_monopole_);
  model_twop->set_data(twop_dataset);
  

  // --------------------------------------------------------------------------
  // ------------- set the underlying dark matter clustering model ------------
  // --------------------------------------------------------------------------

  vector<double> model_scales = linear_bin_vector(100, 1., 300.);
  double sigmaNL = 8;
  model_twop->set_parameters_twop_DM(model_scales, cosmology, redshift, "CAMB", sigmaNL);

  
  // -------------------------------------------------------
  // ------------- set the priors and the model ------------
  // -------------------------------------------------------

  statistics::Prior alpha_prior(PriorType::_UniformPrior_, 0.7,1.3);
  statistics::Prior B_prior (PriorType::_UniformPrior_, 1,15);
  statistics::Prior A0_prior(PriorType::_UniformPrior_, -100, 100);
  statistics::Prior A1_prior(PriorType::_UniformPrior_, -100, 100);
  statistics::Prior A2_prior(PriorType::_UniformPrior_, -1000, 1000);

  model_twop->set_model_AP_isotropic(alpha_prior, B_prior, A0_prior, A1_prior, A2_prior);

  
  // ----------------------------------------------------------------------
  // ------------- run chains and write output chain and model ------------
  // ----------------------------------------------------------------------

  int nChains = 100;
  int chain_size = 2000;
  double rmin = 40;
  double rmax = 200;

  string chain_file = "Anderson_2013_CMASSDR11_monopole_prerecon_chains.dat";
  model_twop->sample_likelihood(rmin, rmax, statistics::_GaussianLikelihood_Covariance_, nChains, chain_size, 43123, dir_chains, chain_file);

  int npt = 200;
  vector<double> model_scales_output = linear_bin_vector(npt, rmin, rmax);
  string output_file = "Anderson_2013_CMASSDR11_monopole_prerecon_bestfit_model.dat";
  model_twop->write_model(model_scales_output, dir_chains, output_file);

  return 0;
}

