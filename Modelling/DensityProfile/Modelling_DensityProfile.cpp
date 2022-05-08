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
#include "Data1D.h"

using namespace std;

using namespace cbl;


// ===========================================================================================

cbl::modelling::densityprofile::Modelling_DensityProfile::Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> profile, const std::string profile_author, const bool _2halo, const std::string halo_def, const double Delta)
{
  m_data = profile->dataset();
  m_profile_author = profile_author;
  m_2halo = _2halo;
  m_mass_is_derived = false;
  m_halo_def = halo_def;
  m_Delta = Delta;
}


// ===========================================================================================

cbl::modelling::densityprofile::Modelling_DensityProfile::Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const std::string profile_author, const bool _2halo, const std::string halo_def, const double Delta)
{
  m_data = dataset;
  m_profile_author = profile_author;
  m_2halo = _2halo;
  m_mass_is_derived = false;
  m_halo_def = halo_def;
  m_Delta = Delta;
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const cosmology::Cosmology cosmology, const double redshift, const double contrast, const double logM_base, const double mass_pivot, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  
  m_data_model.redshift = redshift;
  m_data_model.contrast = contrast;
  m_data_model.logM_base = logM_base;
  m_data_model.mass_pivot = mass_pivot;

  m_data_model.bias_author = bias_author;
  m_data_model.method_Pk = method_Pk;
  m_data_model.interp_type = interp_type;
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const cosmology::Cosmology cosmology, const double redshift, const double mass_proxy, const double redshift_pivot, const double proxy_pivot, const double contrast, const double logM_base, const double mass_pivot, const double Nclusters, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  
  m_data_model.redshift = redshift;
  m_data_model.mass_proxy = mass_proxy;
  m_data_model.redshift_pivot = redshift_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  
  m_data_model.contrast = contrast;
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

  (m_data_model.scaling_relation)->set_data_model(cosmology, {redshift}, redshift_pivot, proxy_pivot, logM_base, {Nclusters});
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior)
{  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+7; // The total number of parameters is given by the cosmological ones + 7, since the density profile has 7 parameters (Rt, conc, logM, f_off, sigma_off, AB_fact, OB_fact)

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
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

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (*(m_data_model.cosmology), m_data_model.redshift, 0., 0., m_Delta, m_profile_author, m_halo_def, 0., true, false, 0., 0.);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) {(void)halo_profile; return conc;};

  // Set the function returning the 2-halo term
  if (m_2halo)
    m_data_model.two_halo_func = [] (const std::vector<double> radius, cbl::cosmology::HaloProfile halo_profile, const std::string bias_author, const std::string method_Pk, const std::string interp_type) {
				   return halo_profile.DeltaSigma_2h(radius, bias_author, method_Pk, interp_type);};
  else
    m_data_model.two_halo_func = [] (const std::vector<double> radius, cbl::cosmology::HaloProfile halo_profile, const std::string bias_author, const std::string method_Pk, const std::string interp_type) {
				   (void)halo_profile; (void)bias_author; (void)method_Pk; (void)interp_type;
				   std::vector<double> res (radius.size(), 0.);
				   return res;};

  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior)
{  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+7; // The total number of parameters is given by the cosmological ones + 7, since the density profile has 6 base parameters (Rt, logM, f_off, sigma_off, AB_fact, OB_fact) and 1 derived parameter (conc)
  const int n_derivedPars = 1;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[cosmo_param.size()+1] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
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

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (*(m_data_model.cosmology), m_data_model.redshift, cM_author, 0., m_Delta, m_profile_author, m_halo_def, 0., true, false, 0., 0.);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) { (void)conc; return halo_profile.concentration(); };

  // Set the function returning the 2-halo term
  if (m_2halo)
    m_data_model.two_halo_func = [] (const std::vector<double> radius, cbl::cosmology::HaloProfile halo_profile, const std::string bias_author, const std::string method_Pk, const std::string interp_type) {
				   return halo_profile.DeltaSigma_2h(radius, bias_author, method_Pk, interp_type);};
  else
    m_data_model.two_halo_func = [] (const std::vector<double> radius, cbl::cosmology::HaloProfile halo_profile, const std::string bias_author, const std::string method_Pk, const std::string interp_type) {
				   (void)halo_profile; (void)bias_author; (void)method_Pk; (void)interp_type;
				   std::vector<double> res (radius.size(), 0.);
				   return res;};
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior)
{
  if (m_mass_is_derived == false)
    ErrorCBL("If the mass is derived from the scaling relation, you must use the correct set_data_model!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+14; // The total number of parameters is given by the cosmological ones, + 14 base parameters.
  
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "Rt";
  param_prior[cosmo_param.size()] = Rt_prior;
  Par_string[cosmo_param.size()+1] = "concentration";
  param_prior[cosmo_param.size()+1] = concentration_prior;
  
  Par_string[cosmo_param.size()+2] = "f_off";
  param_prior[cosmo_param.size()+2] = f_off_prior;
  Par_string[cosmo_param.size()+3] = "sigma_off";
  param_prior[cosmo_param.size()+3] = sigma_off_prior;

  Par_string[cosmo_param.size()+4] = "AB_fact";
  param_prior[cosmo_param.size()+4] = anisotropic_boost_prior;
  Par_string[cosmo_param.size()+5] = "OB_fact";
  param_prior[cosmo_param.size()+5] = orientation_boost_prior;
  
  Par_string[cosmo_param.size()+6] = "alpha";
  param_prior[cosmo_param.size()+6] = alpha_prior;
  Par_string[cosmo_param.size()+7] = "beta";
  param_prior[cosmo_param.size()+7] = beta_prior;
  Par_string[cosmo_param.size()+8] = "gamma";
  param_prior[cosmo_param.size()+8] = gamma_prior;
  Par_string[cosmo_param.size()+9] = "scatter0";
  param_prior[cosmo_param.size()+9] = scatter0_prior;
  Par_string[cosmo_param.size()+10] = "scatterM";
  param_prior[cosmo_param.size()+10] = scatterM_prior;
  Par_string[cosmo_param.size()+11] = "scatterM_exponent";
  param_prior[cosmo_param.size()+11] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+12] = "scatterz";
  param_prior[cosmo_param.size()+12] = scatterz_prior;
  Par_string[cosmo_param.size()+13] = "scatterz_exponent";
  param_prior[cosmo_param.size()+13] = scatterz_exponent_prior;

  // Build the HaloProfile object
  cosmology::HaloProfile halo_profile (*(m_data_model.cosmology), m_data_model.redshift, 0., 0., m_Delta, m_profile_author, m_halo_def, 0., true, false, 0., 0.);
  m_data_model.halo_profile = make_shared<cosmology::HaloProfile>(halo_profile);

  // Set the scaling relation object
  (m_data_model.scaling_relation)->set_model_MassObservableRelation_cosmology (z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
  (m_data_model.scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {});

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::cosmology::HaloProfile halo_profile) {(void)halo_profile; return conc;};

  // Set the function returning the 2-halo term
  if (m_2halo)
    m_data_model.two_halo_func_fast = [] (const std::vector<double> radius, cbl::cosmology::HaloProfile halo_profile, const double bias, const std::string method_Pk, const std::string interp_type) {
				   return halo_profile.DeltaSigma_2h(radius, bias, method_Pk, interp_type);};
  else
    m_data_model.two_halo_func_fast = [] (const std::vector<double> radius, cbl::cosmology::HaloProfile halo_profile, const double bias, const std::string method_Pk, const std::string interp_type) {
				   (void)halo_profile; (void)bias; (void)method_Pk; (void)interp_type;
				   std::vector<double> res (radius.size(), 0.);
				   return res;};
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_density_scaling_relation, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_density (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the HaloProfile object
  cbl::cosmology::HaloProfile halo_profile = *pp->halo_profile;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  halo_profile.set_cosmology(cosmo);

  const double mass = pow(pp->logM_base, parameter[pp->Cpar.size()+2])*pp->mass_pivot;

  halo_profile.set_trunc_fact(parameter[pp->Cpar.size()]);
  halo_profile.set_mass(mass * (1+parameter[pp->Cpar.size()+6]));
  halo_profile.set_f_off(parameter[pp->Cpar.size()+3]);
  halo_profile.set_sigma_off(parameter[pp->Cpar.size()+4]);
  
  halo_profile.set_concentration(pp->conc_func(parameter[pp->Cpar.size()+1], halo_profile));  
  parameter[pp->Cpar.size()+1] = halo_profile.concentration();

  // Compute DeltaSigma
  std::vector<double> one_halo = halo_profile.DeltaSigma(radius);  
  std::vector<double> two_halo = pp->two_halo_func(radius, halo_profile, pp->bias_author, pp->method_Pk, pp->interp_type);

  std::vector<double> total_profile (radius.size());  
  for (size_t j=0; j<radius.size(); j++)    
    total_profile[j] = one_halo[j] + two_halo[j] * (1.+parameter[pp->Cpar.size()+5]);

  return total_profile;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_density_scaling_relation (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);
  
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the HaloProfile object
  cbl::cosmology::HaloProfile halo_profile = *pp->halo_profile;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  auto cosmo_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmo);
  
  // set the cluster parameters
  halo_profile.set_cosmology(cosmo);

  halo_profile.set_trunc_fact(parameter[pp->Cpar.size()]);
  halo_profile.set_f_off(parameter[pp->Cpar.size()+2]);
  halo_profile.set_sigma_off(parameter[pp->Cpar.size()+3]);

  // set the concentration
  halo_profile.set_concentration(pp->conc_func(parameter[pp->Cpar.size()+1], halo_profile));
  parameter[pp->Cpar.size()+1] = halo_profile.concentration();

  // 
  std::vector<double> scalRel_pars;
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    scalRel_pars.push_back(parameter[i]);
  for (int i=6; i<14; i++)
    scalRel_pars.emplace_back(parameter[pp->Cpar.size()+i]);
  
  const double alpha = scalRel_pars[scalRel_pars.size()-8];
  const double beta = scalRel_pars[scalRel_pars.size()-7];
  const double gamma = scalRel_pars[scalRel_pars.size()-6];
  const double scatter0 = scalRel_pars[scalRel_pars.size()-5];
  const double scatterM = scalRel_pars[scalRel_pars.size()-4];
  const double scatterM_exp = scalRel_pars[scalRel_pars.size()-3];
  const double scatterz = scalRel_pars[scalRel_pars.size()-2];
  const double scatterz_exp = scalRel_pars[scalRel_pars.size()-1];

  const double AB_fact = parameter[pp->Cpar.size()+4];
  const double OB_fact = parameter[pp->Cpar.size()+5];

  // Define the integrand
  double sqrt_Nclusters = sqrt( (pp->scaling_relation)->data_model().Nclusters[0] );
  
  cbl::glob::FuncGrid DeltaSigma_Rj_interp;
  std::shared_ptr<void> ptr;
  
  auto integrand = [&] (const double x)
		   {
		     double mass = pow(pp->logM_base,x)*pp->mass_pivot;
		       
		     // Compute P(M|lambda,z)
		     double log_lambda = log(pp->mass_proxy/pp->proxy_pivot)/log(pp->logM_base);
		     double log_f_z = log( (pp->scaling_relation)->data_model().fz(pp->redshift, pp->redshift_pivot, cosmo_ptr) )/log(pp->logM_base);
		       
		     double mean = alpha + beta*log_lambda + gamma*log_f_z;
		     double scatter_intr = std::abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp)) / sqrt_Nclusters;
		     double P_M__lambda_z = cbl::gaussian(x, ptr, {mean,scatter_intr});

		     // Compute the halo excess surface profile
		     double DeltaSigma = DeltaSigma_Rj_interp(mass);
		       
		     return DeltaSigma * P_M__lambda_z;
		   };

  // Find the minimum and maximum masses, given the parameters of the scaling relation and the intrinsic scatter
  double log_lambda = log(pp->mass_proxy/pp->proxy_pivot)/log(pp->logM_base);
  double log_f_z = log( (pp->scaling_relation)->data_model().fz(pp->redshift, pp->redshift_pivot, cosmo_ptr) )/log(pp->logM_base);

  double logM = alpha + beta*log_lambda + gamma*log_f_z;

  double scatter = std::abs( scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp) ) / sqrt_Nclusters;

  double min_logM = logM-3.5*scatter;
  double max_logM = logM+3.5*scatter;

  // »»»»»»»»»»»»»»»»»»
  // Compute DeltaSigma
  // »»»»»»»»»»»»»»»»»»

  const std::vector<double> Mass_vector = cbl::logarithmic_bin_vector(5, pow(pp->logM_base,min_logM)*pp->mass_pivot * (1-std::abs(OB_fact)), pow(pp->logM_base,max_logM)*pp->mass_pivot * (1+std::abs(OB_fact)));

  // For the halo bias in the two-halo term, compute:
  // - Delta
  // - growth factor
  // - the interpolated sigmaM
  double Delta = halo_profile.Delta() / cosmo.OmegaM(pp->redshift);
  double DN = cosmo.DN(pp->redshift);

  std::vector<double> sigmaM (Mass_vector.size(), 0.);
  for (size_t i=0; i<sigmaM.size(); i++)
    sigmaM[i] = sqrt( cosmo.sigma2M(Mass_vector[i], pp->method_Pk, 0., false, "test", "Linear", 100.) );
  cbl::glob::FuncGrid sigmaM_interp (Mass_vector, sigmaM, "Spline");

  // Compute DeltaSigma between maximum and minimum masses, and interpolate it for a given radius
  std::vector<std::vector<double>> total_profile_4interp(Mass_vector.size(), std::vector<double>(radius.size()));

  // 0) Since DeltaSigma_2h / bias_halo is constant, compute it only once
  const std::vector<double> normalised_2h = pp->two_halo_func_fast(radius, halo_profile, 1., pp->method_Pk, pp->interp_type);

  // 1) For a given mass, compute DeltaSigma at all radii
  for (size_t i=0; i<total_profile_4interp.size(); i++) {

    double mass = Mass_vector[i] * (1+OB_fact);
    halo_profile.set_mass(mass);

    double bias = cosmo.bias_halo(mass, sigmaM_interp(mass), pp->redshift, DN, pp->bias_author, false, par::defaultString, "Linear", Delta, -1, -1, 0.001, 100., 1.e-2, pp->method_Pk);
    
    std::vector<double> one_halo = halo_profile.DeltaSigma(radius);
    std::vector<double> two_halo = normalised_2h;

    for (size_t j=0; j<total_profile_4interp[i].size(); j++)
      total_profile_4interp[i][j] = one_halo[j] + two_halo[j] * bias * (1.+AB_fact);
    
  }

  std::vector<double> total_profile(radius.size());
  
  for (size_t j=0; j<radius.size(); j++) {

    // 2) Given a radius, interpolate DeltaSigma as a function of mass
    std::vector<double> DeltaSigma_Rj_4interp(Mass_vector.size());
    
    for (size_t i=0; i<Mass_vector.size(); i++)
      DeltaSigma_Rj_4interp[i] = total_profile_4interp[i][j];

    DeltaSigma_Rj_interp = cbl::glob::FuncGrid (Mass_vector, DeltaSigma_Rj_4interp, "Spline");

    // Integrate
    total_profile[j] = wrapper::gsl::GSL_integrate_qag(integrand, min_logM, max_logM);
    
  }
  
  return total_profile;
}

// ===========================================================================================

