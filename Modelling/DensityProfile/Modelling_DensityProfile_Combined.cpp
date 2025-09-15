/********************************************************************
 *  Copyright (C) 2024 by Federico Marulli and Giorgio Lesci        *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   * 
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/**
 *  @file
 *  Modelling/DensityProfile/Modelling_DensityProfile_Combined.cpp
 *
 *  @brief Methods of the class Modelling_DensityProfile_Combined
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_DensityProfile_Combined, i.e. the common functions to model
 *  the galaxy cluster surface density profiles
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "Modelling_DensityProfile_Combined.h"
#include "Data1D.h"

using namespace std;

using namespace cbl;


// ===========================================================================================

cbl::modelling::densityprofilecomb::Modelling_DensityProfile_Combined::Modelling_DensityProfile_Combined (const std::vector<std::shared_ptr<cbl::data::Data>> dataset, const std::string probe, const std::vector<std::vector<double>> zl_bin_edges_for_N_zs, const std::vector<std::vector<double>> zs, const std::vector<std::vector<std::vector<double>>> N_zs, const std::vector<std::vector<double>> purity_background)
{
  if ( (dataset.size() != zl_bin_edges_for_N_zs.size()) || (dataset.size() != zs.size()) || (dataset.size() != N_zs.size()) || (dataset.size() != purity_background.size()) )
    ErrorCBL("Wrong size of the input vectors!", "Modelling_DensityProfile_Combined", "Modelling_DensityProfile_Combined.cpp");
  
  // Build the Modelling_DensityProfile objects
  std::vector<double> test_xx;
  
  for (size_t i=0; i<dataset.size(); i++) {
    m_data_model.single_object.emplace_back(cbl::modelling::densityprofile::Modelling_DensityProfile(dataset[i], probe, zl_bin_edges_for_N_zs[i], zs[i], N_zs[i], purity_background[i]));

    if (i==0)
      test_xx = dataset[i]->xx();
    else
      if (test_xx != dataset[i]->xx())
	ErrorCBL("For each object, the measurements must be derived in the same radial points!", "Modelling_DensityProfile_Combined", "Modelling_DensityProfile_Combined.cpp");
  }
  

  // ---------------- Define the combined data set ----------------

  // x and data points
  std::vector<double> xx, data;
  for (size_t i=0; i<dataset.size(); i++) {
    for (size_t j=0; j<dataset[i]->xx().size(); j++) {
      xx.emplace_back(dataset[i]->xx(j));
      data.emplace_back(dataset[i]->data(j));
    }
  }
  
  // Covariance matrix
  std::vector<std::vector<double>> covariance (xx.size(), std::vector<double>(xx.size(),0.));
  int add = 0;
  
  for (size_t i=0; i<dataset.size(); i++) {
    for (size_t j=0; j<dataset[i]->xx().size(); j++) {
      for (size_t k=0; k<dataset[i]->xx().size(); k++) {
	covariance[j+add][k+add] = dataset[i]->covariance(j,k);
      }
    }
    add += dataset[i]->xx().size();
  }

  // Data set
  m_data = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(xx, data, covariance));
}

// ===========================================================================================
void cbl::modelling::densityprofilecomb::Modelling_DensityProfile_Combined::write_model_from_combined_chain (const std::string output_dir, const std::string output_file, const int start, const int thin, const std::vector<double> xx)
{
  // Build a vector containing the chain values of the first common parameters
  const int chain_size = this->posterior_parameters()->chain_size();
  const int nwalkers = this->posterior_parameters()->chain_nwalkers();
  
  std::vector<std::vector<double>> chain_values (m_data_model.Cpar.size()+20, std::vector<double>(chain_size * nwalkers));
  for (size_t p=0; p<chain_values.size(); p++)
    for (int i=0; i<chain_size; i++)
      for (int j=0; j<nwalkers; j++)
        chain_values[p][i * nwalkers + j] = this->posterior_parameters()->chain_value(p, i, j);

  // Check
  if ( chain_values.size() + m_data_model.single_object.size() + 2 != this->posterior_parameters()->nparameters() )
    ErrorCBL("Mistake in the new implementation. This function works only for a well defined number of repeated parameters.", "write_model_from_combined_chain", "Modelling_DensityProfile_Combined.cpp");
  
  // Write the models
  std::vector<double> xvec = (xx.size() > 0) ? xx : cbl::different_elements(m_data_fit->xx());
  
  for (size_t p=0; p<m_data_model.single_object.size(); p++) {

    std::vector<double> purity_values(chain_size * nwalkers);
    for (int i=0; i<chain_size; i++)
      for (int j=0; j<nwalkers; j++)
        purity_values[i * nwalkers + j] = this->posterior_parameters()->chain_value(chain_values.size() + p, i, j);
    
    std::vector<std::vector<double>> chain_values_single = chain_values;
    chain_values_single.emplace_back(purity_values);
    chain_values_single.emplace_back(purity_values); // emplace back whatever values for the splashback radius, as it is a derived parameter
    chain_values_single.emplace_back(purity_values); // emplace back whatever values for the splashback radius, as it is a derived parameter

    m_data_model.single_object[p].m_set_posterior(666);
    m_data_model.single_object[p].posterior_parameters()->set_chain_values(chain_values_single, nwalkers);
    m_data_model.single_object[p].likelihood()->get_m_model()->write_from_chains(output_dir, std::to_string(p)+"_"+output_file, xvec, start, thin);
    
  }
}

// ===========================================================================================

void cbl::modelling::densityprofilecomb::Modelling_DensityProfile_Combined::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const std::vector<std::vector<double>> redshift_obs, const std::vector<std::vector<double>> mass_proxy_obs, const std::vector<std::vector<double>> Pz_params, const std::vector<std::vector<double>> Pproxy_params, const double redshift_pivot, const double proxy_pivot, const double logM_base, const double mass_pivot, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  if ( (m_data_model.single_object.size() != redshift_obs.size()) || (m_data_model.single_object.size() != mass_proxy_obs.size()) || (m_data_model.single_object.size() != Pz_params.size()) || (m_data_model.single_object.size() != Pproxy_params.size()) )
    ErrorCBL("Wrong size of the input vectors!", "set_data_model", "Modelling_DensityProfile_Combined.cpp");

  m_data_model.cosmology = move(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units

  // Set the data model for each Modelling_DensityProfile object
  for (size_t i=0; i<m_data_model.single_object.size(); i++)    
    m_data_model.single_object[i].set_data_model(cosmology, profile_author, two_halo, halo_def, Delta, redshift_obs[i], mass_proxy_obs[i], Pz_params[i], Pproxy_params[i], redshift_pivot, proxy_pivot, logM_base, mass_pivot, bias_author, method_Pk, interp_type);

  // Derive the redshift and mass proxy limits
  std::vector<double> local_z_min, local_z_max, local_proxy_min, local_proxy_max;

  for (size_t i=0; i<m_data_model.single_object.size(); i++) {
    local_z_min.emplace_back(cbl::Min(m_data_model.single_object[i].data_model().redshifts));
    local_z_max.emplace_back(cbl::Max(m_data_model.single_object[i].data_model().redshifts));
    local_proxy_min.emplace_back(cbl::Min(m_data_model.single_object[i].data_model().mass_proxies));
    local_proxy_max.emplace_back(cbl::Max(m_data_model.single_object[i].data_model().mass_proxies));
  }

  m_data_model.z_min = cbl::Min(local_z_min);
  m_data_model.z_max = cbl::Max(local_z_max);
  m_data_model.mass_proxy_min = cbl::Min(local_proxy_min);
  m_data_model.mass_proxy_max = cbl::Max(local_proxy_max);
}

// ===========================================================================================

void cbl::modelling::densityprofilecomb::Modelling_DensityProfile_Combined::set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution c0_prior, const statistics::PriorDistribution cM_prior, const statistics::PriorDistribution cz_prior, const statistics::PriorDistribution f_off0_prior, const statistics::PriorDistribution f_offM_prior, const statistics::PriorDistribution f_offz_prior, const statistics::PriorDistribution sigma_off0_prior, const statistics::PriorDistribution sigma_offM_prior, const statistics::PriorDistribution sigma_offz_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<statistics::PriorDistribution> purity_prior)
{
  if (purity_prior.size() != m_data_model.single_object.size())
    ErrorCBL("Wrong number of input purity priors!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile_Combined.cpp");
  
  for (size_t i=0; i<m_data_model.single_object.size(); i++) {
    
    m_data_model.single_object[i].set_model_DensityProfile_cosmology(cosmo_param, cosmo_prior, z_evo, Rt_prior, c0_prior, cM_prior, cz_prior, f_off0_prior, f_offM_prior, f_offz_prior, sigma_off0_prior, sigma_offM_prior, sigma_offz_prior, anisotropic_boost_prior, orientation_boost_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior, purity_prior[i]);

    // set the likelihood for writing the models' 68% confidence levels
    m_data_model.single_object[i].set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Covariance_, {});

  }
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size() + 22 + (int)(purity_prior.size());
  const int n_derivedPars = 2;
  
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[Par_type.size()-1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-2] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  
  if (m_data_model.single_object[0].data_model().profile_author == "NFW_trunc")
    param_prior[cosmo_param.size()] = Rt_prior;
  else {
    param_prior[cosmo_param.size()] = cbl::statistics::PriorDistribution {cbl::glob::DistributionType::_Constant_, 0.};
    cbl::WarningMsgCBL("I set the truncation factor equal to zero and constant, since it is not used in the "+m_data_model.single_object[0].data_model().profile_author+" model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  }
  
  Par_string[cosmo_param.size()+1] = "c0";
  param_prior[cosmo_param.size()+1] = c0_prior;
  Par_string[cosmo_param.size()+2] = "cM";
  param_prior[cosmo_param.size()+2] = cM_prior;
  Par_string[cosmo_param.size()+3] = "cz";
  param_prior[cosmo_param.size()+3] = cz_prior;
  
  Par_string[cosmo_param.size()+4] = "f_off0";
  param_prior[cosmo_param.size()+4] = f_off0_prior;
  Par_string[cosmo_param.size()+5] = "f_offM";
  param_prior[cosmo_param.size()+5] = f_offM_prior;
  Par_string[cosmo_param.size()+6] = "f_offz";
  param_prior[cosmo_param.size()+6] = f_offz_prior;
  Par_string[cosmo_param.size()+7] = "sigma_off0";
  param_prior[cosmo_param.size()+7] = sigma_off0_prior;
  Par_string[cosmo_param.size()+8] = "sigma_offM";
  param_prior[cosmo_param.size()+8] = sigma_offM_prior;
  Par_string[cosmo_param.size()+9] = "sigma_offz";
  param_prior[cosmo_param.size()+9] = sigma_offz_prior;

  Par_string[cosmo_param.size()+10] = "AB_fact";
  param_prior[cosmo_param.size()+10] = anisotropic_boost_prior;
  Par_string[cosmo_param.size()+11] = "OB_fact";
  param_prior[cosmo_param.size()+11] = orientation_boost_prior;
  
  Par_string[cosmo_param.size()+12] = "alpha";
  param_prior[cosmo_param.size()+12] = alpha_prior;
  Par_string[cosmo_param.size()+13] = "beta";
  param_prior[cosmo_param.size()+13] = beta_prior;
  Par_string[cosmo_param.size()+14] = "gamma";
  param_prior[cosmo_param.size()+14] = gamma_prior;
  Par_string[cosmo_param.size()+15] = "scatter0";
  param_prior[cosmo_param.size()+15] = scatter0_prior;
  Par_string[cosmo_param.size()+16] = "scatterM";
  param_prior[cosmo_param.size()+16] = scatterM_prior;
  Par_string[cosmo_param.size()+17] = "scatterM_exponent";
  param_prior[cosmo_param.size()+17] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+18] = "scatterz";
  param_prior[cosmo_param.size()+18] = scatterz_prior;
  Par_string[cosmo_param.size()+19] = "scatterz_exponent";
  param_prior[cosmo_param.size()+19] = scatterz_exponent_prior;

  for (size_t i=0; i<purity_prior.size(); i++) {
    Par_string[cosmo_param.size()+20+i] = "purity_"+cbl::conv(i, cbl::par::fINT);
    param_prior[cosmo_param.size()+20+i] = purity_prior[i];
  }

  Par_string[nParams - 2] = "Rsp";
  Par_string[nParams - 1] = "Rsp_rDelta";
  
  // input data used to construct the model
  auto inputs = make_shared<STR_ProfileComb_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  m_data_model.Par_type = Par_type;
  m_data_model.Par_string = Par_string;

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&combined_model_density_scaling_relation_evolving_concentration_offcentering, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofilecomb::combined_model_density_scaling_relation_evolving_concentration_offcentering (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  const std::vector<double> radius_unique = cbl::different_elements(radius);
  
  // structure contaning the required input data
  shared_ptr<STR_ProfileComb_data_model> pp = static_pointer_cast<STR_ProfileComb_data_model>(inputs);
  
  // redefine the cosmology
  std::shared_ptr<cosmology::Cosmology> cosmology = pp->cosmology->clone();
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmology->set_parameter(pp->Cpar[i], parameter[i]);

  // halo parameters
  const double Rt = parameter[pp->Cpar.size()];
  const double c0 = parameter[pp->Cpar.size()+1];
  const double cM = parameter[pp->Cpar.size()+2];
  const double cz = parameter[pp->Cpar.size()+3];
  const double f_off0 = parameter[pp->Cpar.size()+4];
  const double f_offM = parameter[pp->Cpar.size()+5];
  const double f_offz = parameter[pp->Cpar.size()+6];
  const double sigma_off0 = parameter[pp->Cpar.size()+7];
  const double sigma_offM = parameter[pp->Cpar.size()+8];
  const double sigma_offz = parameter[pp->Cpar.size()+9];
  const double AB_fact = parameter[pp->Cpar.size()+10];
  const double OB_fact = parameter[pp->Cpar.size()+11];

  // scaling relation parameters  
  const double alpha = parameter[pp->Cpar.size()+12];
  const double beta = parameter[pp->Cpar.size()+13];
  const double gamma = parameter[pp->Cpar.size()+14];
  const double scatter0 = parameter[pp->Cpar.size()+15];
  const double scatterM = parameter[pp->Cpar.size()+16];
  const double scatterM_exp = parameter[pp->Cpar.size()+17];
  const double scatterz = parameter[pp->Cpar.size()+18];
  const double scatterz_exp = parameter[pp->Cpar.size()+19];

  // purity
  std::vector<double> purity (pp->single_object.size(), 0.);
  for (size_t i=0; i<pp->single_object.size(); i++) {
    
    purity[i] = parameter[pp->Cpar.size() + 20 + i];    
    if ( (purity[i] < 0) || (purity[i] > 1) )
      ErrorCBL("The purity must range between 0 and 1!", "combined_model_density_scaling_relation_evolving_concentration_offcentering", "Modelling_DensityProfile_Combined.cpp");
    
  }
  

  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Interpolate DeltaSigma and Sigma as a function of redshift and mass
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  auto dm = pp->single_object[0].data_model();
  
  // Define the (optimal) redshift vector
  const std::vector<double> redshift_vector = cbl::linear_bin_vector(std::max(4, (int)((pp->z_max-pp->z_min)/0.1)), pp->z_min, pp->z_max);

  // Define the (optimal) mass vector
  const std::vector<double> min_max = modelling::densityprofileaux::min_max_logMass(pp->z_min, pp->z_max, pp->mass_proxy_min, pp->mass_proxy_max, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmology, dm.scaling_relation, dm.redshift_pivot, dm.proxy_pivot, dm.logM_base);

  const double M_min = pow(dm.logM_base, min_max[0]) * dm.mass_pivot;
  const double M_max = pow(dm.logM_base, min_max[1]) * dm.mass_pivot;
  const std::vector<double> Mass_vector = cbl::logarithmic_bin_vector(std::max(5, (int)((log10(M_max)-log10(M_min))*5)), M_min, M_max);

  // Interpolate the profiles
  std::vector<std::vector<cbl::glob::FuncGrid2D>> interp_vec = cbl::modelling::densityprofileaux::interpolated_profiles (radius_unique, redshift_vector, Mass_vector, cosmology, dm.profile_author, dm.two_halo, dm.halo_def, dm.Delta, AB_fact, dm.bias_author, dm.method_Pk, dm.conc_scaling_relation_func, c0, cM, cz, sigma_off0, sigma_offM, sigma_offz, Rt, dm.redshift_pivot, dm.mass_pivot);
  std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_cen = interp_vec[0];
  std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_off = interp_vec[1];
  std::vector<cbl::glob::FuncGrid2D> Sigma_interp_cen = interp_vec[2];
  std::vector<cbl::glob::FuncGrid2D> Sigma_interp_off = interp_vec[3];


  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the model for each Modelling_DensityProfile instance
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  std::vector<double> profile;

  for (size_t i=0; i<pp->single_object.size(); i++) {

    auto dm_single = pp->single_object[i].data_model();

    // If the probe is g, interpolate Sigma_crit
    cbl::glob::FuncGrid inv_Scrit_interp;
    if (dm_single.probe == "g")
      inv_Scrit_interp = modelling::densityprofileaux::Sigma_crit_factor_interp(redshift_vector, dm_single.zl_bin_edges_for_N_zs, dm_single.zs, dm_single.N_zs, cosmology); 

    // Compute the profile  
    std::vector<double> profile_single = cbl::modelling::densityprofileaux::compute_model_density_scaling_relation(dm_single.redshifts, dm_single.mass_proxies, radius_unique, dm_single.probe, DeltaSigma_interp_cen, DeltaSigma_interp_off, Sigma_interp_cen, Sigma_interp_off, inv_Scrit_interp, cosmology, f_off0, f_offM, f_offz, OB_fact, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, purity[i], dm_single.scaling_relation, dm_single.redshift_pivot, dm_single.proxy_pivot, dm_single.mass_pivot, dm_single.logM_base, dm_single.purity_background, dm_single.N_zs_index);

    // Emplace
    for (size_t j=0; j<profile_single.size(); j++)
      profile.emplace_back(profile_single[j]);
    
  }
  
  // At the moment, the splashback radius derivation is not implemented for this model.
  // We extract a random number between -1 and 0 to avoid problems in the posterior derivation.
  srand(time(0));
  cbl::random::UniformRandomNumbers extract(-1, 0, rand());
  parameter[parameter.size()-2] = extract();
  parameter[parameter.size()-1] = extract();

  return profile;
}
