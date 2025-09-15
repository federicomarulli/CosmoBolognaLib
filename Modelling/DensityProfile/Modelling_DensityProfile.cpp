/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  Modelling/DensityProfile/Modelling_DensityProfile.cpp
 *
 *  @brief Methods of the class Modelling_DensityProfile
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_DensityProfile, i.e. the common functions to model
 *  the galaxy cluster surface density profiles
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "Modelling_DensityProfile.h"
#include "Modelling_NumberCounts.h"
#include "Data1D.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


cbl::modelling::densityprofile::Modelling_DensityProfile::Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> dataset, const std::string probe, const std::vector<double> zl_bin_edges_for_N_zs, const std::vector<double> zs, const std::vector<std::vector<double>> N_zs, const std::vector<double> purity_background)
{  
  m_data = dataset->dataset();
  m_mass_is_derived = false;
  m_set_probe_model(probe, zl_bin_edges_for_N_zs, zs, N_zs, purity_background);  
}


// ===========================================================================================

cbl::modelling::densityprofile::Modelling_DensityProfile::Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const std::string probe, const std::vector<double> zl_bin_edges_for_N_zs, const std::vector<double> zs, const std::vector<std::vector<double>> N_zs, const std::vector<double> purity_background)
{
  m_data = dataset;
  m_mass_is_derived = false;
  m_set_probe_model(probe, zl_bin_edges_for_N_zs, zs, N_zs, purity_background);
}


// ===========================================================================================


void cbl::modelling::densityprofile::Modelling_DensityProfile::m_set_probe_model (const std::string probe, const std::vector<double> zl_bin_edges_for_N_zs, const std::vector<double> zs, const std::vector<std::vector<double>> N_zs, const std::vector<double> purity_background)
{
  // Set the function returning the probe's model
  if ( (probe != "DeltaSigma") && (probe != "g") )
    ErrorCBL("Wrong declaration of the probe!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  m_data_model.probe = probe;

  // checks regarding the source redshift distributions
  if ( (zs.size() == 0) && (N_zs.size() == 0) && (probe == "g") )
    ErrorCBL("If the reduced shear is the probe, the source redshift distribution must be provided!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");
  
  if ( (zl_bin_edges_for_N_zs.size() - 1 != N_zs.size()) && (probe == "g") )
    ErrorCBL("The size of zl_bin_edges_for_N_zs must be equal to the one of N_zs plus 1!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  if ( (purity_background.size() != N_zs.size()) && (probe == "g") )
    ErrorCBL("The sizes of purity_background and N_zs must be equal!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  for (size_t i=0; i<N_zs.size(); i++)
    if ( (zs.size() != N_zs[i].size()) && (probe == "g") )
      ErrorCBL("The sizes of zs and the ones of the N_zs elements must be equal!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  for (size_t i=0; i<purity_background.size(); i++)
    if ( ( (purity_background[i] < 0.) || (purity_background[i] > 1.) ) && (probe == "g") )
      ErrorCBL("The elements of purity_background must be defined between 0 and 1!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  // Set the parameters
  m_data_model.zs = zs;
  m_data_model.N_zs = N_zs;
  m_data_model.zl_bin_edges_for_N_zs = zl_bin_edges_for_N_zs;
  m_data_model.purity_background = purity_background;

  // interpolate P_background as a function of lens redshift
  if (zl_bin_edges_for_N_zs.size()>2) {
    
    std::vector<double> zl_bin_centres(zl_bin_edges_for_N_zs.size() - 1);
    for (size_t i=0; i<zl_bin_edges_for_N_zs.size() - 1; i++)
      zl_bin_centres[i] = 0.5 * (zl_bin_edges_for_N_zs[i] + zl_bin_edges_for_N_zs[i+1]);
  
    m_data_model.purity_background_interp = cbl::glob::FuncGrid(zl_bin_centres, m_data_model.purity_background, "Spline");
    
  }
}


// ===========================================================================================


void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const double redshift, const double logM_base, const double mass_pivot, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{ 
  m_data_model.cosmology = move(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units

  m_data_model.profile_author = profile_author;
  m_data_model.two_halo = two_halo;
  m_data_model.halo_def = halo_def;
  m_data_model.Delta = Delta;
  
  m_data_model.redshift = redshift;
  m_data_model.logM_base = logM_base;
  m_data_model.mass_pivot = mass_pivot;

  m_data_model.bias_author = bias_author;
  m_data_model.method_Pk = method_Pk;
  m_data_model.interp_type = interp_type;

  // source redshift distribution
  if ( (m_data_model.N_zs.size() != 1) && (m_data_model.probe == "g") )
    ErrorCBL("Exactly one source redshift distribution must be provided in this case! "+cbl::conv(m_data_model.N_zs.size(), cbl::par::fINT)+" were provided.", "set_data_model", "Modelling_DensityProfile.cpp");
}


// ===========================================================================================


void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const std::vector<double> redshift_obs, const std::vector<double> mass_proxy_obs, const std::vector<double> Pz_params, const std::vector<double> Pproxy_params, const double redshift_pivot, const double proxy_pivot, const double logM_base, const double mass_pivot, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  if ( redshift_obs.size() != mass_proxy_obs.size() )
    ErrorCBL("Redshift and mass proxy vectors must have the same size!", "set_data_model", "Modelling_DensityProfile.cpp");

  m_data_model.min_extracted = 1.e-3; // minimum value an extracted quantity can assume

  // define the skew-normal PDF
  std::function<double(double, std::shared_ptr<void>, std::vector<double>)> skewnorm = 
    [](double x, std::shared_ptr<void> pp, std::vector<double> params) -> double {

      (void) pp;
      
      double mu = params[0];
      double sigma = params[1];
      double alpha = params[2];
		   
      double normpdf = 1. / ( sigma * sqrt(2. * cbl::par::pi) ) * exp( - pow(x - mu, 2.) / (2. * pow(sigma, 2)) ); 
      double normcdf = 0.5 * ( 1. + erf( alpha * (x - mu) / sigma / sqrt(2.) ) );
    
      return 2 * normpdf * normcdf;
	
    };

  // extract redshift values from P(z_tr|z_ob,proxy_ob)
  if (Pz_params.size() != 3)
    ErrorCBL("Pz_params must be a vector with size equal to 3.", "set_data_model", "Modelling_DensityProfile.cpp");

  if ( (Pz_params[0] == 0.) && (Pz_params[1] == 0.) ) {
    
    m_data_model.redshifts = redshift_obs;

  }
  else {

    cbl::random::CustomDistributionRandomNumbers rand (skewnorm, NULL, Pz_params, 666, -1., 1.);

    m_data_model.redshifts.resize(redshift_obs.size());

    for (size_t i=0; i<m_data_model.redshifts.size(); i++)
      m_data_model.redshifts[i] = std::max(m_data_model.min_extracted, rand() * (1. + redshift_obs[i]) + redshift_obs[i]);

  }

  // extract mass proxy values from P(proxy_tr|proxy_ob,z_ob)
  if (Pproxy_params.size() != 3)
    ErrorCBL("Pproxy_params must be a vector with size equal to 3.", "set_data_model", "Modelling_DensityProfile.cpp");

  if ( (Pproxy_params[0] == 0.) && (Pproxy_params[1] == 0.) ) {
    
    m_data_model.mass_proxies = mass_proxy_obs;

  }
  else {

    cbl::random::CustomDistributionRandomNumbers rand (skewnorm, NULL, Pproxy_params, 666, -1., 1.);

    m_data_model.mass_proxies.resize(mass_proxy_obs.size());

    for (size_t i=0; i<m_data_model.mass_proxies.size(); i++)
      m_data_model.mass_proxies[i] = std::max(m_data_model.min_extracted, rand() * mass_proxy_obs[i] + mass_proxy_obs[i]);

  }
  
  // set the rest of the model parameters
  m_data_model.cosmology = move(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units

  m_data_model.profile_author = profile_author;
  m_data_model.two_halo = two_halo;
  m_data_model.halo_def = halo_def;
  m_data_model.Delta = Delta;
  
  m_data_model.redshift_pivot = redshift_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  
  m_data_model.logM_base = logM_base;
  m_data_model.mass_pivot = mass_pivot;

  m_data_model.bias_author = bias_author;
  m_data_model.method_Pk = method_Pk;
  m_data_model.interp_type = interp_type;

  m_mass_is_derived = true;

  // Build a dummy dataset for the scaling relation Modelling object, useful only to avoid internal errors
  std::vector<double> dummy_vec = {1.};
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(dummy_vec, dummy_vec, dummy_vec));
  
  // Build the scaling relation object
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation (dataset);
  m_data_model.scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);

  (m_data_model.scaling_relation)->set_data_model(cosmology, {m_data_model.redshifts}, redshift_pivot, proxy_pivot, logM_base);

  // assign a background sample purity to each cluster, based on its observed redshift
  if ( (m_data_model.probe == "g") && (m_data_model.zl_bin_edges_for_N_zs.size() != 2) && (redshift_obs.size() == 1) )
    ErrorCBL("If only one cluster is present in the sample, only one source redshift distribution must be provided!", "set_data_model", "Modelling_DensityProfile.cpp");
  
  m_data_model.N_zs_index.resize(redshift_obs.size(),-1);
  
  if (m_data_model.zl_bin_edges_for_N_zs.size() != 0) {
  
    for (size_t j=0; j<redshift_obs.size(); j++)
      for (size_t i=0; i<m_data_model.zl_bin_edges_for_N_zs.size() - 1; i++)
	if ( (redshift_obs[j] >= m_data_model.zl_bin_edges_for_N_zs[i]) && (redshift_obs[j] < m_data_model.zl_bin_edges_for_N_zs[i+1]) )
	  m_data_model.N_zs_index[j] = i;
	else if (redshift_obs[j] == m_data_model.zl_bin_edges_for_N_zs[m_data_model.zl_bin_edges_for_N_zs.size()])
	  m_data_model.N_zs_index[j] = m_data_model.zl_bin_edges_for_N_zs.size() - 2;

    for (size_t j=0; j<redshift_obs.size(); j++)
      if (m_data_model.N_zs_index[j] == -1)
	ErrorCBL("I could not find a redshift bin for at least one cluster!", "set_data_model", "Modelling_DensityProfile.cpp");

  } else {

    for (size_t j=0; j<redshift_obs.size(); j++)
      m_data_model.N_zs_index[j] = 0;
    
  }
}


// ===========================================================================================


void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const double redshift_min, const double redshift_max, const double mass_proxy_min, const double mass_proxy_max, const double z_error, const double proxy_rel_error, const double redshift_pivot, const double proxy_pivot, const double logM_base, const double mass_pivot, const std::vector<double> redshift_points_completeness, const std::vector<double> proxy_points_completeness, const std::vector<std::vector<double>> completeness, const std::string bias_author, const std::string method_Pk, std::string interp_type, const std::string MF_author)
{  
  m_data_model.cosmology = move(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  
  if ( (redshift_points_completeness.size() > 0) && (proxy_points_completeness.size() > 0) )
    m_data_model.completeness_interp = cbl::glob::FuncGrid2D(redshift_points_completeness, proxy_points_completeness, completeness, "Cubic");  
  else {    
    std::vector<double> dummy_values = cbl::linear_bin_vector(10, 0., 10000.);
    std::vector<std::vector<double>> dummy_completeness(dummy_values.size(), std::vector<double>(dummy_values.size(), 1.));
    m_data_model.completeness_interp = cbl::glob::FuncGrid2D(dummy_values, dummy_values, dummy_completeness, "Cubic");
  }
  
  m_data_model.profile_author = profile_author;
  m_data_model.two_halo = two_halo;
  m_data_model.halo_def = halo_def;
  m_data_model.Delta = Delta;
  
  m_data_model.z_min = redshift_min;
  m_data_model.z_max = redshift_max;
  m_data_model.mass_proxy_min = mass_proxy_min;
  m_data_model.mass_proxy_max = mass_proxy_max;
  
  m_data_model.z_error = z_error;
  this->set_P_proxy(0., 0., 0., proxy_rel_error, 0., 0.);
  
  m_data_model.redshift_pivot = redshift_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  
  m_data_model.logM_base = logM_base;
  m_data_model.mass_pivot = mass_pivot;

  m_data_model.bias_author = bias_author;
  m_data_model.method_Pk = method_Pk;
  m_data_model.interp_type = interp_type;

  m_data_model.MF_author = MF_author;
  m_data_model.isTheoretical_MF = true;
  m_data_model.k_min = 1.e-4;
  m_data_model.k_max = 100;
  m_data_model.kk = logarithmic_bin_vector(500, m_data_model.k_min, m_data_model.k_max);
  m_data_model.norm = -1;
  m_data_model.store_output = false;
  m_data_model.output_root = "test";
  m_data_model.prec = 1.e-4;
  m_data_model.Mass_vector = logarithmic_bin_vector(200, 1.e10, 1.e16);

  m_mass_is_derived = true;

  // Build a dummy dataset for the scaling relation Modelling object, useful only to avoid internal errors
  std::vector<double> dummy_vec = {1.};
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(dummy_vec, dummy_vec, dummy_vec));
  
  // Build the scaling relation object
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation (dataset);
  m_data_model.scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);

  (m_data_model.scaling_relation)->set_data_model(cosmology, {0.}, redshift_pivot, proxy_pivot, logM_base);
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy (const double A_mu, const double B_mu, const double C_mu, const double A_sigma, const double B_sigma, const double C_sigma)
{
  if (m_model != NULL)
    ErrorCBL("This function must be called before the definition of the model!", "set_P_proxy", "Modelling_DensityProfile.cpp");
  
  m_data_model.Plambda_A_mu = A_mu;
  m_data_model.Plambda_B_mu = B_mu;
  m_data_model.Plambda_C_mu = C_mu;
  m_data_model.Plambda_A_sigma = A_sigma;
  m_data_model.Plambda_B_sigma = B_sigma;
  m_data_model.Plambda_C_sigma = C_sigma;
  
  m_data_model.Plambda_mean_fc = [] (const double proxy_tr, const double z_tr, const double A_mu, const double B_mu, const double C_mu)
				 {
				   return proxy_tr + A_mu * proxy_tr * exp(- proxy_tr * (B_mu + C_mu * z_tr));
				 };

  m_data_model.Plambda_std_fc = [] (const double proxy_tr, const double z_tr, const double A_sigma, const double B_sigma, const double C_sigma)
				 {
				   return A_sigma * proxy_tr * exp(- proxy_tr * (B_sigma + C_sigma * z_tr));
				 };
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology_DK14 (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution be_prior, const statistics::PriorDistribution se_prior, const statistics::PriorDistribution alpha_0_prior, const statistics::PriorDistribution alpha_nu_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_0_prior)
{
  if ( m_data_model.profile_author != "DK14" )
    ErrorCBL("This function can be used only if the DK14 profile is assumed!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");

  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+15; // The total number of parameters is given by the cosmological ones + 15
  const int n_derivedPars = 2;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[Par_type.size()-1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-2] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  if (cosmo_param.size() != cosmo_prior.size())
    ErrorCBL("cosmo_param and cosmo_prior must have the same size!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  for (size_t i=0; i<cosmo_param.size(); i++) {
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  param_prior[cosmo_param.size()] = Rt_prior;
  Par_string[cosmo_param.size()+1] = "concentration";
  param_prior[cosmo_param.size()+1] = concentration_prior;
  Par_string[cosmo_param.size()+2] = "logM";
  param_prior[cosmo_param.size()+2] = logM_prior;
  Par_string[cosmo_param.size()+3] = "f_off";
  param_prior[cosmo_param.size()+3] = f_off_prior;
  Par_string[cosmo_param.size()+4] = "sigma_off";
  param_prior[cosmo_param.size()+4] = sigma_off_prior;

  Par_string[cosmo_param.size()+5] = "AB_fact";
  param_prior[cosmo_param.size()+5] = anisotropic_boost_prior;
  Par_string[cosmo_param.size()+6] = "OB_fact";
  param_prior[cosmo_param.size()+6] = orientation_boost_prior;

  Par_string[cosmo_param.size()+7] = "b_e";
  param_prior[cosmo_param.size()+7] = be_prior;
  Par_string[cosmo_param.size()+8] = "s_e";
  param_prior[cosmo_param.size()+8] = se_prior;

  Par_string[cosmo_param.size()+9] = "alpha_0";
  param_prior[cosmo_param.size()+9] = alpha_0_prior;
  Par_string[cosmo_param.size()+10] = "alpha_nu";
  param_prior[cosmo_param.size()+10] = alpha_nu_prior;

  Par_string[cosmo_param.size()+11] = "beta";
  param_prior[cosmo_param.size()+11] = beta_prior;
  Par_string[cosmo_param.size()+12] = "gamma_0";
  param_prior[cosmo_param.size()+12] = gamma_0_prior;

  Par_string[cosmo_param.size()+13] = "Rsp";
  Par_string[cosmo_param.size()+14] = "Rsp_rDelta";

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (m_data_model.cosmology, m_data_model.two_halo, m_data_model.redshift, 2., 1.e14, m_data_model.Delta, m_data_model.profile_author, m_data_model.halo_def, 0., true, false, 0., 0., 0., 0., 0., m_data_model.bias_author, m_data_model.method_Pk);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) {(void)halo_profile; return conc;};

  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  m_data_model.Par_type = Par_type;
  m_data_model.Par_string = Par_string;

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density_DK14, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology_DK14 (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution be_prior, const statistics::PriorDistribution se_prior, const statistics::PriorDistribution alpha_0_prior, const statistics::PriorDistribution alpha_nu_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_0_prior)
{
  if ( m_data_model.profile_author != "DK14" )
    ErrorCBL("This function can be used only if the DK14 profile is assumed!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");

  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+15; // The total number of parameters is given by the cosmological ones + 15
  const int n_derivedPars = 3;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[cosmo_param.size()+1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-2] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  if (cosmo_param.size() != cosmo_prior.size())
    ErrorCBL("cosmo_param and cosmo_prior must have the same size!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  param_prior[cosmo_param.size()] = Rt_prior;
  
  Par_string[cosmo_param.size()+1] = "concentration";
  
  Par_string[cosmo_param.size()+2] = "logM";
  param_prior[cosmo_param.size()+1] = logM_prior;
  Par_string[cosmo_param.size()+3] = "f_off";
  param_prior[cosmo_param.size()+2] = f_off_prior;
  Par_string[cosmo_param.size()+4] = "sigma_off";
  param_prior[cosmo_param.size()+3] = sigma_off_prior;

  Par_string[cosmo_param.size()+5] = "AB_fact";
  param_prior[cosmo_param.size()+4] = anisotropic_boost_prior;
  Par_string[cosmo_param.size()+6] = "OB_fact";
  param_prior[cosmo_param.size()+5] = orientation_boost_prior;

  Par_string[cosmo_param.size()+7] = "b_e";
  param_prior[cosmo_param.size()+6] = be_prior;
  Par_string[cosmo_param.size()+8] = "s_e";
  param_prior[cosmo_param.size()+7] = se_prior;

  Par_string[cosmo_param.size()+9] = "alpha_0";
  param_prior[cosmo_param.size()+8] = alpha_0_prior;
  Par_string[cosmo_param.size()+10] = "alpha_nu";
  param_prior[cosmo_param.size()+9] = alpha_nu_prior;

  Par_string[cosmo_param.size()+11] = "beta";
  param_prior[cosmo_param.size()+10] = beta_prior;
  Par_string[cosmo_param.size()+12] = "gamma_0";
  param_prior[cosmo_param.size()+11] = gamma_0_prior;

  Par_string[cosmo_param.size()+13] = "Rsp";
  Par_string[cosmo_param.size()+14] = "Rsp_rDelta";

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (m_data_model.cosmology, m_data_model.two_halo, m_data_model.redshift, cM_author, 1.e14, m_data_model.Delta, m_data_model.profile_author, m_data_model.halo_def, 0., true, false, 0., 0., 0., 0., 0., m_data_model.bias_author, m_data_model.method_Pk);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) { (void)conc; return halo_profile.concentration(); };

  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  m_data_model.Par_type = Par_type;
  m_data_model.Par_string = Par_string;

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density_DK14, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior)
{
  if ( (m_data_model.profile_author != "NFW") && (m_data_model.profile_author != "NFW_trunc") && (m_data_model.profile_author != "Einasto") )
    ErrorCBL("You cannot use this set_model if the "+m_data_model.profile_author+" profile is assumed! Check the other set_model functions in this class. This error message was implemented to avoid the use of DK14 profile along with this set_model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+9; // The total number of parameters is given by the cosmological ones + 9, since the density profile has 9 parameters (Rt, conc, logM, f_off, sigma_off, AB_fact, OB_fact, Rsp, Rsp_rDelta)
  const int n_derivedPars = 2;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[Par_type.size()-1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-2] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  if (cosmo_param.size() != cosmo_prior.size())
    ErrorCBL("cosmo_param and cosmo_prior must have the same size!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  if (m_data_model.profile_author == "NFW_trunc")
    param_prior[cosmo_param.size()] = Rt_prior;
  else {
    param_prior[cosmo_param.size()] = cbl::statistics::PriorDistribution {cbl::glob::DistributionType::_Constant_, 0.};
    cbl::WarningMsgCBL("I set the truncation factor equal to zero and constant, since it is not used in the "+m_data_model.profile_author+" model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  }
  Par_string[cosmo_param.size()+1] = "concentration";
  param_prior[cosmo_param.size()+1] = concentration_prior;
  Par_string[cosmo_param.size()+2] = "logM";
  param_prior[cosmo_param.size()+2] = logM_prior;
  Par_string[cosmo_param.size()+3] = "f_off";
  param_prior[cosmo_param.size()+3] = f_off_prior;
  Par_string[cosmo_param.size()+4] = "sigma_off";
  param_prior[cosmo_param.size()+4] = sigma_off_prior;

  Par_string[cosmo_param.size()+5] = "AB_fact";
  param_prior[cosmo_param.size()+5] = anisotropic_boost_prior;
  Par_string[cosmo_param.size()+6] = "OB_fact";
  param_prior[cosmo_param.size()+6] = orientation_boost_prior;

  Par_string[cosmo_param.size()+7] = "Rsp";
  Par_string[cosmo_param.size()+8] = "Rsp_rDelta";

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (m_data_model.cosmology, m_data_model.two_halo, m_data_model.redshift, 2., 1.e14, m_data_model.Delta, m_data_model.profile_author, m_data_model.halo_def, 0., true, false, 0., 0., 0., 0., 0., m_data_model.bias_author, m_data_model.method_Pk);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) {(void)halo_profile; return conc;};

  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  m_data_model.Par_type = Par_type;
  m_data_model.Par_string = Par_string;

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior)
{
  if ( (m_data_model.profile_author != "NFW") && (m_data_model.profile_author != "NFW_trunc") && (m_data_model.profile_author != "Einasto") )
    ErrorCBL("You cannot use this set_model if the "+m_data_model.profile_author+" profile is assumed! Check the other set_model functions in this class. This error message was implemented to avoid the use of DK14 profile along with this set_model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+9;
  const int n_derivedPars = 3;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[cosmo_param.size()+1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-2] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  if (cosmo_param.size() != cosmo_prior.size())
    ErrorCBL("cosmo_param and cosmo_prior must have the same size!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  if (m_data_model.profile_author == "NFW_trunc")
    param_prior[cosmo_param.size()] = Rt_prior;
  else {
    param_prior[cosmo_param.size()] = cbl::statistics::PriorDistribution {cbl::glob::DistributionType::_Constant_, 0.};
    cbl::WarningMsgCBL("I set the truncation factor equal to zero and constant, since it is not used in the "+m_data_model.profile_author+" model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  }
  
  Par_string[cosmo_param.size()+1] = "concentration";

  Par_string[cosmo_param.size()+2] = "logM";
  param_prior[cosmo_param.size()+1] = logM_prior;
  Par_string[cosmo_param.size()+3] = "f_off";
  param_prior[cosmo_param.size()+2] = f_off_prior;
  Par_string[cosmo_param.size()+4] = "sigma_off";
  param_prior[cosmo_param.size()+3] = sigma_off_prior;

  Par_string[cosmo_param.size()+5] = "AB_fact";
  param_prior[cosmo_param.size()+4] = anisotropic_boost_prior;
  Par_string[cosmo_param.size()+6] = "OB_fact";
  param_prior[cosmo_param.size()+5] = orientation_boost_prior;

  Par_string[cosmo_param.size()+7] = "Rsp";
  Par_string[cosmo_param.size()+8] = "Rsp_rDelta";

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (m_data_model.cosmology, m_data_model.two_halo, m_data_model.redshift, cM_author, 1.e14, m_data_model.Delta, m_data_model.profile_author, m_data_model.halo_def, 0., true, false, 0., 0., 0., 0., 0., m_data_model.bias_author, m_data_model.method_Pk);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) { (void)conc; return halo_profile.concentration(); };
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  m_data_model.Par_type = Par_type;
  m_data_model.Par_string = Par_string;

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution c0_prior, const statistics::PriorDistribution cM_prior, const statistics::PriorDistribution cz_prior, const statistics::PriorDistribution f_off0_prior, const statistics::PriorDistribution f_offM_prior, const statistics::PriorDistribution f_offz_prior, const statistics::PriorDistribution sigma_off0_prior, const statistics::PriorDistribution sigma_offM_prior, const statistics::PriorDistribution sigma_offz_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution purity_prior)
{
  if ( (m_data_model.profile_author != "NFW") && (m_data_model.profile_author != "NFW_trunc") && (m_data_model.profile_author != "Einasto") )
    ErrorCBL("You cannot use this set_model if the "+m_data_model.profile_author+" profile is assumed! Check the other set_model functions in this class. This error message was implemented to avoid the use of DK14 profile along with this set_model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  if (m_mass_is_derived == false)
    ErrorCBL("If the mass is derived from the scaling relation, you must use the correct set_data_model!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+23; // The total number of parameters is given by the cosmological ones, + 21 base parameters and 2 derived parameters (Rsp, Rsp_rDelta)
  const int n_derivedPars = 2;
  
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[Par_type.size()-1] = statistics::ParameterType::_Derived_;
  Par_type[Par_type.size()-2] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  if (cosmo_param.size() != cosmo_prior.size())
    ErrorCBL("cosmo_param and cosmo_prior must have the same size!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  if (m_data_model.profile_author == "NFW_trunc")
    param_prior[cosmo_param.size()] = Rt_prior;
  else {
    param_prior[cosmo_param.size()] = cbl::statistics::PriorDistribution {cbl::glob::DistributionType::_Constant_, 0.};
    cbl::WarningMsgCBL("I set the truncation factor equal to zero and constant, since it is not used in the "+m_data_model.profile_author+" model.", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
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

  Par_string[cosmo_param.size()+20] = "purity";
  param_prior[cosmo_param.size()+20] = purity_prior;

  Par_string[cosmo_param.size()+21] = "Rsp";
  Par_string[cosmo_param.size()+22] = "Rsp_rDelta";

  // Set the scaling relation object
  (m_data_model.scaling_relation)->set_model_MassObservableRelation_cosmology (z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
  (m_data_model.scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {});

  // Set the function returning the concentration
  m_data_model.conc_scaling_relation_func = [] (const double c0, const double cM, const double cz, const double logM, const double logz) {
					      const double logc = c0 + cM*logM + cz*logz;
					      return pow(10, logc);
					    };
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  m_data_model.Par_type = Par_type;
  m_data_model.Par_string = Par_string;

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density_scaling_relation_evolving_concentration_offcentering, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::Modelling_DensityProfile::mass_expected_value (const std::vector<double> z_ob, const std::vector<double> mass_proxy_ob, const std::shared_ptr<cosmology::Cosmology> cosmology, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp)
{

  if (m_data_model.isTheoretical_MF == false)
    ErrorCBL("This function is designed only for the case of expected values weighted over the theoretical mass function.", "mass_expected_value", "Modelling_DensityProfile.cpp");;

  // Define the redshift and mass vectors used to interpolate the following cosmological quantities
  std::vector<std::vector<double>> int_limits = cbl::modelling::numbercounts::get_integration_limits
    (
     cosmology, cbl::Min(z_ob), cbl::Max(z_ob), cbl::Min(mass_proxy_ob), cbl::Max(mass_proxy_ob), m_data_model.z_error,
     m_data_model.Plambda_mean_fc, m_data_model.Plambda_A_mu, m_data_model.Plambda_B_mu, m_data_model.Plambda_C_mu,
     m_data_model.Plambda_std_fc, m_data_model.Plambda_A_sigma, m_data_model.Plambda_B_sigma, m_data_model.Plambda_C_sigma,
     m_data_model.scaling_relation, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp,
     m_data_model.redshift_pivot, m_data_model.proxy_pivot, m_data_model.mass_pivot, m_data_model.logM_base
     );

  double z_min = int_limits[2][0] * 0.99;
  double z_max = int_limits[2][1] * 1.01;

  double M_min = pow(m_data_model.logM_base, int_limits[0][0]) * m_data_model.mass_pivot * 0.99;
  double M_max = pow(m_data_model.logM_base, int_limits[0][1]) * m_data_model.mass_pivot * 1.01;
  
  // interpolate sigmaM and its derivative
  cosmology::PkXi PX(cosmology);
  std::vector<double> Pk = PX.Pk_matter(m_data_model.kk, m_data_model.method_Pk, false, 0., m_data_model.store_output, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, cbl::par::defaultString, true);
  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (cbl::logarithmic_bin_vector(500, M_min, M_max), cosmology, m_data_model.kk, Pk, "Spline", m_data_model.k_max);

  // interpolate the growth factor
  const std::vector<double> z_for_DN = cbl::linear_bin_vector(100, z_min, z_max);
  std::vector<double> DN (z_for_DN.size(), 0.);
  for (size_t i=0; i<z_for_DN.size(); i++)
    DN[i] = cosmology->DN(z_for_DN[i]);
  cbl::glob::FuncGrid interp_DN (z_for_DN, DN, "Spline");

  // mass expected values
  return cbl::modelling::densityprofileaux::compute_mass_expected_value (z_ob, mass_proxy_ob, m_data_model.MF_author, m_data_model.halo_def, m_data_model.Delta, m_data_model.z_error, m_data_model.Plambda_mean_fc, m_data_model.Plambda_A_mu, m_data_model.Plambda_B_mu, m_data_model.Plambda_C_mu, m_data_model.Plambda_std_fc, m_data_model.Plambda_A_sigma, m_data_model.Plambda_B_sigma, m_data_model.Plambda_C_sigma, interp[0], interp[1], interp_DN, cosmology, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, m_data_model.completeness_interp, m_data_model.scaling_relation, m_data_model.redshift_pivot, m_data_model.proxy_pivot, m_data_model.mass_pivot, m_data_model.logM_base);
  
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_density (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{  
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  std::shared_ptr<cosmology::Cosmology> cosmology = pp->cosmology->clone();

  // redefine the HaloProfile object
  cbl::cosmology::HaloProfile halo_profile = *pp->halo_profile;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmology->set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  halo_profile.set_cosmology(cosmology, false);
  
  const double Rt = parameter[pp->Cpar.size()];
  const double f_off = parameter[pp->Cpar.size()+3];
  const double sigma_off = parameter[pp->Cpar.size()+4];
  const double AB_fact = parameter[pp->Cpar.size()+5];
  const double OB_fact = parameter[pp->Cpar.size()+6];

  const double mass = pow(pp->logM_base, parameter[pp->Cpar.size()+2])*pp->mass_pivot;
  
  halo_profile.set_trunc_fact(Rt);
  halo_profile.set_mass(mass * (1+OB_fact));
  halo_profile.set_f_off(f_off);
  halo_profile.set_sigma_off(sigma_off);
  halo_profile.set_AB_fact(AB_fact);
  
  halo_profile.set_concentration(pp->conc_func(parameter[pp->Cpar.size()+1], halo_profile));  
  parameter[pp->Cpar.size()+1] = halo_profile.concentration();

  // Compute the splashback radius, if possible, and the profile
  if (pp->two_halo) {
    halo_profile.update_2halo();
    const double Rsp = cbl::modelling::densityprofileaux::splashback_radius(halo_profile.rho(radius), radius);
    const double overdensity =  cosmology->rho_crit(pp->redshift) * halo_profile.Delta();
    const double rDelta = pow( 3. * (mass * (1+OB_fact)) / (4. * cbl::par::pi * overdensity), 1./3. );
    parameter[parameter.size()-2] = Rsp;
    parameter[parameter.size()-1] = Rsp / rDelta;
  }
  else {
    // We extract a random number between -1 and 0 to avoid problems in the posterior derivation.
    srand(time(0));
    cbl::random::UniformRandomNumbers extract(-1, 0, rand());
    parameter[parameter.size()-2] = extract();
    parameter[parameter.size()-1] = extract();
  }
  
  if (pp->probe == "DeltaSigma")
    return halo_profile.DeltaSigma(radius);
  else {
    std::vector<double> profile = halo_profile.g(radius, pp->zs, pp->N_zs[0]);    
    for (size_t i=0; i<profile.size(); i++)
      profile[i] = profile[i] * pp->purity_background[0];
    return profile;
  }
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_density_DK14 (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  std::shared_ptr<cosmology::Cosmology> cosmology = pp->cosmology->clone();

  // redefine the HaloProfile object
  cbl::cosmology::HaloProfile halo_profile = *pp->halo_profile;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmology->set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  halo_profile.set_cosmology(cosmology, false);

  const double Rt = parameter[pp->Cpar.size()];
  const double f_off = parameter[pp->Cpar.size()+3];
  const double sigma_off = parameter[pp->Cpar.size()+4];
  const double AB_fact = parameter[pp->Cpar.size()+5];
  const double OB_fact = parameter[pp->Cpar.size()+6];
  const double b_e = parameter[pp->Cpar.size()+7];
  const double s_e = parameter[pp->Cpar.size()+8];
  const double alpha_0 = parameter[pp->Cpar.size()+9];
  const double alpha_nu = parameter[pp->Cpar.size()+10];
  const double beta = parameter[pp->Cpar.size()+11];
  const double gamma_0 = parameter[pp->Cpar.size()+12];

  const double mass = pow(pp->logM_base, parameter[pp->Cpar.size()+2])*pp->mass_pivot;
  
  halo_profile.set_trunc_fact(Rt);
  halo_profile.set_mass(mass * (1+OB_fact));
  halo_profile.set_f_off(f_off);
  halo_profile.set_sigma_off(sigma_off);
  halo_profile.set_b_e(b_e);
  halo_profile.set_s_e(s_e);
  halo_profile.set_AB_fact(AB_fact);
  halo_profile.set_alpha(alpha_0, alpha_nu);
  halo_profile.set_beta(beta);
  halo_profile.set_gamma_0(gamma_0);
  
  halo_profile.set_concentration(pp->conc_func(parameter[pp->Cpar.size()+1], halo_profile));  
  parameter[pp->Cpar.size()+1] = halo_profile.concentration();

  if (pp->two_halo)
    ErrorCBL("If the DK14 profile is assumed, the 2-halo term cannot be included.", "model_density_DK14", "Modelling_DensityProfile.cpp");

  // Compute the profile and the splashback radius
  const double Rsp = cbl::modelling::densityprofileaux::splashback_radius(halo_profile.rho(radius), radius);
  const double overdensity =  cosmology->rho_crit(pp->redshift) * halo_profile.Delta();
  const double rDelta = pow( 3. * (mass * (1+OB_fact)) / (4. * cbl::par::pi * overdensity), 1./3. );
  parameter[parameter.size()-2] = Rsp;
  parameter[parameter.size()-1] = Rsp / rDelta;
  
  if (pp->probe == "DeltaSigma")
    return halo_profile.DeltaSigma(radius);
  else {
    std::vector<double> profile = halo_profile.g(radius, pp->zs, pp->N_zs[0]);
    for (size_t i=0; i<profile.size(); i++)
      profile[i] = profile[i] * pp->purity_background[0];
    return profile;
  }
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_density_scaling_relation_evolving_concentration_offcentering (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);
  
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
  const double purity = parameter[parameter.size()-3];
  if ( (purity < 0) || (purity > 1) )
    ErrorCBL("The purity must range between 0 and 1!", "model_density_scaling_relation_evolving_concentration_offcentering", "Modelling_DensityProfile.cpp");


  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Interpolate DeltaSigma and Sigma as a function of redshift and mass
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  std::vector<double> redshift_vector, Mass_vector;
  
  if (pp->isTheoretical_MF) {

    // Define the (optimal) redshift and mass vectors
    std::vector<std::vector<double>> int_limits = cbl::modelling::numbercounts::get_integration_limits
    (
     cosmology, pp->z_min, pp->z_max, pp->mass_proxy_min, pp->mass_proxy_max, pp->z_error,
     pp->Plambda_mean_fc, pp->Plambda_A_mu, pp->Plambda_B_mu, pp->Plambda_C_mu, pp->Plambda_std_fc, pp->Plambda_A_sigma, pp->Plambda_B_sigma, pp->Plambda_C_sigma,
     pp->scaling_relation, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, pp->redshift_pivot, pp->proxy_pivot, pp->mass_pivot, pp->logM_base
     );

    double z_min = int_limits[2][0] * 0.99;
    double z_max = int_limits[2][1] * 1.01;
    redshift_vector = cbl::linear_bin_vector(std::max(4, (int)((z_max-z_min)/0.1)), z_min, z_max);

    double M_min = pow(pp->logM_base, int_limits[0][0]) * pp->mass_pivot * 0.99;
    double M_max = pow(pp->logM_base, int_limits[0][1]) * pp->mass_pivot * 1.01;
    Mass_vector = cbl::logarithmic_bin_vector(std::min( 20, std::max(5, (int)((log10(M_max)-log10(M_min))*5)) ), M_min, M_max);

  } else {

    // Define the (optimal) redshift vector
    double z_min = cbl::Min(pp->redshifts);
    double z_max = cbl::Max(pp->redshifts);
    redshift_vector = cbl::linear_bin_vector(std::max(4, (int)((z_max-z_min)/0.1)), z_min, z_max);
    
    // Define the (optimal) mass vector
    double proxy_min = cbl::Min(pp->mass_proxies);
    double proxy_max = cbl::Max(pp->mass_proxies);
    std::vector<double> min_max = modelling::densityprofileaux::min_max_logMass(z_min, z_max, proxy_min, proxy_max, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmology, pp->scaling_relation, pp->redshift_pivot, pp->proxy_pivot, pp->logM_base);

    double M_min = pow(pp->logM_base, min_max[0]) * pp->mass_pivot;
    double M_max = pow(pp->logM_base, min_max[1]) * pp->mass_pivot;
    Mass_vector = cbl::logarithmic_bin_vector(std::min( 20, std::max(5, (int)((log10(M_max)-log10(M_min))*5)) ), M_min, M_max);
    
  }

  // Interpolate the profiles
  std::vector<std::vector<cbl::glob::FuncGrid2D>> interp_vec = cbl::modelling::densityprofileaux::interpolated_profiles (radius, redshift_vector, Mass_vector, cosmology, pp->profile_author, pp->two_halo, pp->halo_def, pp->Delta, AB_fact, pp->bias_author, pp->method_Pk, pp->conc_scaling_relation_func, c0, cM, cz, sigma_off0, sigma_offM, sigma_offz, Rt, pp->redshift_pivot, pp->mass_pivot);
  std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_cen = interp_vec[0];
  std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_off = interp_vec[1];
  std::vector<cbl::glob::FuncGrid2D> Sigma_interp_cen = interp_vec[2];
  std::vector<cbl::glob::FuncGrid2D> Sigma_interp_off = interp_vec[3];


  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // If the probe is g, interpolate Sigma_crit
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  cbl::glob::FuncGrid inv_Scrit_interp;
  if (pp->probe == "g")
    inv_Scrit_interp = modelling::densityprofileaux::Sigma_crit_factor_interp(redshift_vector, pp->zl_bin_edges_for_N_zs, pp->zs, pp->N_zs, cosmology);

  
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Compute the profile and the splashback radius
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  std::vector<double> profile (radius.size(), 0.);
  
  if (pp->isTheoretical_MF) {

    // compute the power spectrum
    cosmology::PkXi PX(cosmology);
    std::vector<double> Pk = PX.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, cbl::par::defaultString, true);

    // interpolate sigmaM and its derivative
    const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (cbl::logarithmic_bin_vector(200, cbl::Min(Mass_vector), cbl::Max(Mass_vector)), cosmology, pp->kk, Pk, "Spline", pp->k_max);

    // interpolate the growth factor
    const std::vector<double> z_for_DN = cbl::linear_bin_vector(10, cbl::Min(redshift_vector), cbl::Max(redshift_vector));
    std::vector<double> DN (z_for_DN.size(), 0.);
    for (size_t i=0; i<z_for_DN.size(); i++)
      DN[i] = cosmology->DN(z_for_DN[i]);
    cbl::glob::FuncGrid DN_interp (z_for_DN, DN, "Spline");

    // interpolation range for the background purity
    const int vec_size = (int)(pp->zl_bin_edges_for_N_zs.size());
    std::vector<double> P_bkg_interp_range;

    if (vec_size > 0)
      P_bkg_interp_range = {0.5 * (pp->zl_bin_edges_for_N_zs[0] + pp->zl_bin_edges_for_N_zs[1]), 0.5 * (pp->zl_bin_edges_for_N_zs[vec_size-2] + pp->zl_bin_edges_for_N_zs[vec_size-1])};
    
    // compute the profile
    profile = cbl::modelling::densityprofileaux::compute_model_density_scaling_relation(pp->MF_author, pp->halo_def, pp->Delta, pp->z_min, pp->z_max, pp->mass_proxy_min, pp->mass_proxy_max, pp->z_error, pp->Plambda_mean_fc, pp->Plambda_A_mu, pp->Plambda_B_mu, pp->Plambda_C_mu, pp->Plambda_std_fc, pp->Plambda_A_sigma, pp->Plambda_B_sigma, pp->Plambda_C_sigma, interp[0], interp[1], DN_interp, radius, pp->probe, DeltaSigma_interp_cen, DeltaSigma_interp_off, Sigma_interp_cen, Sigma_interp_off, inv_Scrit_interp, cosmology, f_off0, f_offM, f_offz, OB_fact, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, purity, pp->completeness_interp, pp->scaling_relation, pp->redshift_pivot, pp->proxy_pivot, pp->mass_pivot, pp->logM_base, pp->purity_background_interp, P_bkg_interp_range);
    
  } else
    profile = cbl::modelling::densityprofileaux::compute_model_density_scaling_relation(pp->redshifts, pp->mass_proxies, radius, pp->probe, DeltaSigma_interp_cen, DeltaSigma_interp_off, Sigma_interp_cen, Sigma_interp_off, inv_Scrit_interp, cosmology, f_off0, f_offM, f_offz, OB_fact, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, purity, pp->scaling_relation, pp->redshift_pivot, pp->proxy_pivot, pp->mass_pivot, pp->logM_base, pp->purity_background, pp->N_zs_index);
  
  // At the moment, the splashback radius derivation is not implemented for this model.
  // We extract a random number between -1 and 0 to avoid problems in the posterior derivation.
  srand(time(0));
  cbl::random::UniformRandomNumbers extract(-1, 0, rand());
  parameter[parameter.size()-2] = extract();
  parameter[parameter.size()-1] = extract();

  return profile;
}

