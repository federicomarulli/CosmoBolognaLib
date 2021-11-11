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

cbl::modelling::densityprofile::Modelling_DensityProfile::Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> profile, const std::string profile_author, const bool _2halo, const std::string halo_def)
{
  m_data = profile->dataset();
  m_profile_author = profile_author;
  m_2halo = _2halo;
  m_mass_is_derived = false;

  if (m_profile_author != "NFW")
    ErrorCBL("Wrong profile author declaration!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  m_halo_def = halo_def;

  if (m_halo_def != "200" && m_halo_def != "500")
    ErrorCBL("The code can only deal with critical densities!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");
}


// ===========================================================================================

cbl::modelling::densityprofile::Modelling_DensityProfile::Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const std::string profile_author, const bool _2halo, const std::string halo_def)
{
  m_data = dataset;
  m_profile_author = profile_author;
  m_2halo = _2halo;
  m_mass_is_derived = false;

  if (m_profile_author != "NFW")
    ErrorCBL("Wrong profile author declaration!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");

  m_halo_def = halo_def;

  if (m_halo_def != "200" && m_halo_def != "500")
    ErrorCBL("The code can only deal with critical densities!", "Modelling_DensityProfile", "Modelling_DensityProfile.cpp");
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const double redshift, const double contrast, const double trunc_fact, const double logM_base, const double mass_pivot, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  m_data_model.cluster = make_shared<catalogue::Cluster>(cluster);
  
  m_data_model.redshift = redshift;
  m_data_model.contrast = contrast;
  m_data_model.trunc_fact = trunc_fact;
  m_data_model.logM_base = logM_base;
  m_data_model.mass_pivot = mass_pivot;

  m_data_model.bias_author = bias_author;
  m_data_model.method_Pk = method_Pk;
  m_data_model.interp_type = interp_type;
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const double redshift, const double mass_proxy, const double redshift_pivot, const double proxy_pivot, const double contrast, const double trunc_fact, const double logM_base, const double mass_pivot, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  m_data_model.cluster = make_shared<catalogue::Cluster>(cluster);
  
  m_data_model.redshift = redshift;
  m_data_model.mass_proxy = mass_proxy;
  m_data_model.redshift_pivot = redshift_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  
  m_data_model.contrast = contrast;
  m_data_model.trunc_fact = trunc_fact;
  m_data_model.logM_base = logM_base;
  m_data_model.mass_pivot = mass_pivot;

  m_data_model.bias_author = bias_author;
  m_data_model.method_Pk = method_Pk;
  m_data_model.interp_type = interp_type;

  m_mass_is_derived = true;

  // Build a dummy dataset for the scaling relation Modelling object, useful only to avoid internal errors
  std::vector<double> dummy_vec = {1., 1.};
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(dummy_vec, dummy_vec, dummy_vec));
  
  // Build the scaling relation object
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation (dataset);
  m_data_model.scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);

  (m_data_model.scaling_relation)->set_data_model(cosmology, cluster, {redshift}, redshift_pivot, proxy_pivot, logM_base);
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior)
{  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+4; // The total number of parameters is given by the cosmological ones + 4, since the density profile has 4 parameters (conc, logM, f_off, sigma_off)

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "concentration";
  param_prior[cosmo_param.size()] = concentration_prior;
  Par_string[cosmo_param.size()+1] = "logM";
  param_prior[cosmo_param.size()+1] = logM_prior;
  Par_string[cosmo_param.size()+2] = "f_off";
  param_prior[cosmo_param.size()+2] = f_off_prior;
  Par_string[cosmo_param.size()+3] = "sigma_off";
  param_prior[cosmo_param.size()+3] = sigma_off_prior;

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::catalogue::Cluster cluster) {(void)cluster; return conc;};

  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  if (m_2halo == false)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated, nParams, Par_type, Par_string, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated_2halo, nParams, Par_type, Par_string, inputs));
}

// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior)
{  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+4; // The total number of parameters is given by the cosmological ones + 4, since the density profile has 3 base parameters (logM, f_off, sigma_off) and 1 derived parameter (conc)
  const int n_derivedPars = 1;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[0] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "concentration";
  
  Par_string[cosmo_param.size()+1] = "logM";
  param_prior[cosmo_param.size()] = logM_prior;
  Par_string[cosmo_param.size()+2] = "f_off";
  param_prior[cosmo_param.size()+1] = f_off_prior;
  Par_string[cosmo_param.size()+3] = "sigma_off";
  param_prior[cosmo_param.size()+2] = sigma_off_prior;

  // Set the profile (and the concentration-mass relation) for the Cluster object
  m_data_model.cluster->set_profile(*m_data_model.cosmology, 0., 0., cM_author, m_profile_author, m_halo_def);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::catalogue::Cluster cluster) { (void)conc; return cluster.concentration_from_mass(); };
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  if (m_2halo == false)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated, nParams, Par_type, Par_string, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated_2halo, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior)
{
  if (m_mass_is_derived == false)
    ErrorCBL("If the mass is derived from the scaling relation, you must use the correct set_data_model!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+12; // The total number of parameters is given by the cosmological ones, + 11 base parameters, and 1 derived (the mass).
  const int n_derivedPars = 1;
  
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[nParams-1] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "concentration";
  param_prior[cosmo_param.size()] = concentration_prior;
  Par_string[cosmo_param.size()+1] = "f_off";
  param_prior[cosmo_param.size()+1] = f_off_prior;
  Par_string[cosmo_param.size()+2] = "sigma_off";
  param_prior[cosmo_param.size()+2] = sigma_off_prior;
  Par_string[cosmo_param.size()+3] = "alpha";
  param_prior[cosmo_param.size()+3] = alpha_prior;
  Par_string[cosmo_param.size()+4] = "beta";
  param_prior[cosmo_param.size()+4] = beta_prior;
  Par_string[cosmo_param.size()+5] = "gamma";
  param_prior[cosmo_param.size()+5] = gamma_prior;
  Par_string[cosmo_param.size()+6] = "scatter0";
  param_prior[cosmo_param.size()+6] = scatter0_prior;
  Par_string[cosmo_param.size()+7] = "scatterM";
  param_prior[cosmo_param.size()+7] = scatterM_prior;
  Par_string[cosmo_param.size()+8] = "scatterM_exponent";
  param_prior[cosmo_param.size()+8] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+9] = "scatterz";
  param_prior[cosmo_param.size()+9] = scatterz_prior;
  Par_string[cosmo_param.size()+10] = "scatterz_exponent";
  param_prior[cosmo_param.size()+10] = scatterz_exponent_prior;

  Par_string[cosmo_param.size()+11] = "logM";

  // Set the scaling relation object
  (m_data_model.scaling_relation)->set_model_MassObservableRelation_cosmology (z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
  (m_data_model.scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {});

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::catalogue::Cluster cluster) {(void)cluster; return conc;};
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  if (m_2halo == false)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated_scaling_relation, nParams, Par_type, Par_string, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated_2halo_scaling_relation, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const std::string cM_author, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior)
{
  if (m_mass_is_derived == false)
    ErrorCBL("If the mass is derived from the scaling relation, you must use the correct set_data_model!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+12; // The total number of parameters is given by the cosmological ones, + 10 base parameters, and 2 derived (concentration and mass).
  const int n_derivedPars = 2;
  
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[0] = statistics::ParameterType::_Derived_;
  Par_type[nParams-1] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams-n_derivedPars);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "concentration";

  Par_string[cosmo_param.size()+1] = "f_off";
  param_prior[cosmo_param.size()] = f_off_prior;
  Par_string[cosmo_param.size()+2] = "sigma_off";
  param_prior[cosmo_param.size()+1] = sigma_off_prior;
  Par_string[cosmo_param.size()+3] = "alpha";
  param_prior[cosmo_param.size()+2] = alpha_prior;
  Par_string[cosmo_param.size()+4] = "beta";
  param_prior[cosmo_param.size()+3] = beta_prior;
  Par_string[cosmo_param.size()+5] = "gamma";
  param_prior[cosmo_param.size()+4] = gamma_prior;
  Par_string[cosmo_param.size()+6] = "scatter0";
  param_prior[cosmo_param.size()+5] = scatter0_prior;
  Par_string[cosmo_param.size()+7] = "scatterM";
  param_prior[cosmo_param.size()+6] = scatterM_prior;
  Par_string[cosmo_param.size()+8] = "scatterM_exponent";
  param_prior[cosmo_param.size()+7] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+9] = "scatterz";
  param_prior[cosmo_param.size()+8] = scatterz_prior;
  Par_string[cosmo_param.size()+10] = "scatterz_exponent";
  param_prior[cosmo_param.size()+9] = scatterz_exponent_prior;

  Par_string[cosmo_param.size()+11] = "logM";

  // Set the scaling relation object
  (m_data_model.scaling_relation)->set_model_MassObservableRelation_cosmology (z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
  (m_data_model.scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {});

  // Set the profile (and the concentration-mass relation) for the Cluster object
  m_data_model.cluster->set_profile(*m_data_model.cosmology, 0., 0., cM_author, m_profile_author, m_halo_def);

  // Set the function returning the concentration
  m_data_model.conc_func = [] (const double conc, cbl::catalogue::Cluster cluster) { (void)conc; return cluster.concentration_from_mass(); };
  
  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  if (m_2halo == false)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated_scaling_relation, nParams, Par_type, Par_string, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_NFW_truncated_2halo_scaling_relation, nParams, Par_type, Par_string, inputs));
}
				  

// ===========================================================================================

double cbl::modelling::densityprofile::NFW_truncated (cbl::cosmology::Cosmology cosm, cbl::catalogue::Cluster cluster, const double r, const double redshift, const double contrast, const double trunc_fact)
{
  // **** Equations ****
  //
  // Let us take Aguena et al. 2021 as the reference. The NFW profile is expressed as:
  //
  // rho(r) = rho_s / (r/r_s) / (1+r/r_s),    (1)
  //
  // where both the scale radius r_s and the characteristic density rho_s are functions of the cluster’s mass, M, and the concentration, c. Under the spherical overdensity mass definition
  // used in this code, r_s is expressed as:
  //
  // r_s = 1/c * (3*M/(4*pi*Delta*rho_bg))^(1/3),    (2)
  //
  // where in our case 
  //
  // Delta*rho_bg = Delta*rho_crit.    (3)
  //
  // In Eq. (1), rho_s is the characteristic density and it is expressed as follows:
  //
  // rho_s = c^3 * Delta*rho_bg / (3*I(c)),    (4)
  //
  // where
  //
  // I(c) = int_0^c dx x^2 rho_hat(x),    (5)
  //
  // where x=r/r_s and rho_hat(x) = rho(x*r_s)/rho_s is the dimensionless density. From Eq. (2), we can write Eq. (4) as
  //
  // rho_s = M / (4*pi*I(c)).    (6)
  //
  // **** End of equations ****

  
  // Calculate the radius enclosing a given mass, in the cosmology considered
  double H = cosm.HH(redshift)/cosm.HH(0.)*100./3.0857e+19; // in sec^-1
  double rho_crit = 3*H*H/8/cbl::par::pi/6.6732e-8/1.98892e+33*3.0857e+24*3.0857e+24*3.0857e+24; // in h^2*M_sun/Mpc^3
  double r_s_times_c = pow( 3*cluster.mass()/(4*cbl::par::pi*contrast*rho_crit), 1./3. ); // Eq. (2) multiplied by the concentration

  double r_s = r_s_times_c/cluster.concentration(); // Eq. (2)
  double tau = trunc_fact*cluster.concentration();

  std::function<double(double)> F = [] (double x){
				      // both operations must be generalised to complex numbers when x<1
				      // in that case, both num and den are imaginary, and their ratio is real again
				      std::complex<double> num = std::acos(std::complex<double>(1./x,0.));
				      std::complex<double> den = std::sqrt(std::complex<double>(x*x-1.,0.));  
				      std::complex<double> res = num/den;						 
				      return std::abs(real(res));
				    };
  std::function<double(double)> L = [&tau] (double x){
				      return log(x/(sqrt(x*x+tau*tau)+tau));
				    };

  std::function<double(double)> integrand_I = [&tau,&r_s] (double rad){ // Eq. (5)
						   double x = rad/r_s;
						   double rho_hat = 1./x/(1+x)/(1+x) * (tau*tau/(tau*tau+x*x))*(tau*tau/(tau*tau+x*x)); // in (M_sun/h)/((Mpc/h)^3)
						   return rad*rad*rho_hat;
						 };
  
  double a_meantau = 0, b_meantau = 0, inv_mis_scale = 1./cluster.sigma_off()/cluster.sigma_off(); 
  std::function<double(double)> tau_factor1 = [&a_meantau,&b_meantau,&inv_mis_scale] (double s){
						return (exp(-0.5*(s*s+a_meantau)*inv_mis_scale))/sqrt(b_meantau-a_meantau-s*s);
					      };
  std::function<double(double)> tau_factor2 = [&a_meantau,&b_meantau,&inv_mis_scale] (double s){
						return (exp(-0.5*(b_meantau-s*s)*inv_mis_scale))/sqrt(b_meantau-a_meantau-s*s);
					      };


  // Define the variable tau_for_misc
  double tau_for_misc = 0;
  
  std::function<double(double)> meantau_factor = [&tau_for_misc,&cluster,&a_meantau,&b_meantau,&inv_mis_scale,&tau_factor1,&tau_factor2] (double rad){
						   a_meantau = (tau_for_misc-rad)*(tau_for_misc-rad);
						   b_meantau = (tau_for_misc+rad)*(tau_for_misc+rad);
						   double c_meantau = (a_meantau+b_meantau)*0.5;
						   double upper_limit = std::min(sqrt(c_meantau-a_meantau),sqrt(25*cluster.sigma_off()*cluster.sigma_off()-a_meantau));
						   double meantau_factor = rad*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor1, 0, upper_limit);
						   if (c_meantau < 25*cluster.sigma_off()*cluster.sigma_off())
						     {
						       double int_left = std::max(0.,sqrt(b_meantau-25*cluster.sigma_off()*cluster.sigma_off()));
						       double int_right = sqrt(b_meantau-c_meantau);						     
						       meantau_factor += rad*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor2, int_left, int_right);
						     }
						   return meantau_factor;
						 };
  
  // Calculate rho_s
  double I_c = cbl::wrapper::gsl::GSL_integrate_cquad(integrand_I,0.,r_s_times_c);
  double rho_s = cluster.mass()/(4.*cbl::par::pi*I_c); // Eq. (6)

  std::function<double(double)> sigma_cen_fc = [&r_s,&tau,&rho_s,&F,&L] (double r){
						 double x = r/r_s;
						 double tausq = tau*tau; double tausq_xsq = tausq+x*x; double tau4th = tausq*tausq;
						 double prefact = rho_s*r_s*tau4th/(tausq+1.)/(tausq+1.)/(tausq+1.);
						 double a = 2*(tausq+1.)/(x*x-1.)*(1-F(x));
						 double b = 8*F(x);
						 double c = (tau4th-1.)/tausq/tausq_xsq;
						 double d = -cbl::par::pi*(4*tausq_xsq+tausq+1.)/pow(tausq_xsq,1.5);
						 double e = (tausq*(tau4th-1.)+(tausq_xsq)*(3*tau4th-6*tausq-1))*L(x)/(tau*tau*tau)/pow(tausq_xsq,1.5);    
						 return prefact*(a+b+c+d+e)/1.e+12; // in (M_sun/h)/((pc/h)^2) = h*M_sun/pc^2
					       };
  std::function<double(double)> sigma_cen_int = [&sigma_cen_fc] (double rad){
						  return sigma_cen_fc(rad)*rad;
						};
  std::function<double(double)> meansigma_off_int = [&tau_for_misc,&r,&cluster,&sigma_cen_fc,&meantau_factor] (double tau){
						      tau_for_misc = tau;
						      double int_left = std::max(0.,tau-5*cluster.sigma_off());
						      double int_right = std::min(r,tau+5*cluster.sigma_off());
						      return tau * sigma_cen_fc(tau) * cbl::wrapper::gsl::GSL_integrate_cquad(meantau_factor, int_left, int_right);
						    };
  std::function<double(double)> sigma_off_int = [&r,&cluster,&a_meantau,&b_meantau,&inv_mis_scale,&tau_factor1,&tau_factor2,&sigma_cen_fc] (double tau){
						  a_meantau = (tau-r)*(tau-r);
						  b_meantau = (tau+r)*(tau+r);
						  double c_meantau = (a_meantau+b_meantau)*0.5;
						  double upper_limit = std::min(sqrt(c_meantau-a_meantau),sqrt(25*cluster.sigma_off()*cluster.sigma_off()-a_meantau));
						  double sigma_off_int = r*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor1, 0., upper_limit);
						  if (c_meantau < 25*cluster.sigma_off()*cluster.sigma_off())
						    {
						      double int_left = std::max(0.,sqrt(b_meantau-25*cluster.sigma_off()*cluster.sigma_off()));
						      double int_right = sqrt(b_meantau-c_meantau);						     
						      sigma_off_int += r*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor2, int_left, int_right);
						    }
						  return sigma_off_int * tau * sigma_cen_fc(tau);
						};
  
  double sigma_cen = sigma_cen_fc(r);
  double meansigma_cen = cbl::wrapper::gsl::GSL_integrate_cquad(sigma_cen_int,0.,r) * 2./r/r;

  double deltasigma_cen = meansigma_cen - sigma_cen;
  
  double deltasigma_1h = 0, meansigma_off = 0, integrated_sigma_off = 0, deltasigma_off = 0;
  if (cluster.f_off() <= 1.e-5 || cluster.sigma_off() <= 1.e-5)
    deltasigma_1h = deltasigma_cen;
  else {
    meansigma_off = cbl::wrapper::gsl::GSL_integrate_cquad(meansigma_off_int,0.,r+5.*cluster.sigma_off()) *2./cbl::par::pi/cluster.sigma_off()/cluster.sigma_off()/r/r;
    integrated_sigma_off = cbl::wrapper::gsl::GSL_integrate_cquad(sigma_off_int,std::max(0., r-5*cluster.sigma_off()),r+5.*cluster.sigma_off()) /r/cbl::par::pi/cluster.sigma_off()/cluster.sigma_off();
    deltasigma_off = meansigma_off - integrated_sigma_off;
    deltasigma_1h = (1-cluster.f_off())*deltasigma_cen+cluster.f_off()*deltasigma_off;
  }
    
  return deltasigma_1h;
}

// ===========================================================================================

double cbl::modelling::densityprofile::two_halo (cbl::cosmology::Cosmology cosm, const double radius, const double mass, const double redshift, const double contrast, const std::string bias_author, const std::string method_Pk, std::string interp_type)
{
  double delta_bkg = cosm.Delta_vir(contrast, redshift);
  double Bias      = cosm.bias_halo(mass, redshift, bias_author, method_Pk, false, "test", interp_type, delta_bkg, -1., -1, 0.001, 100); 
  double Rho_m     = cosm.rho_m(redshift);
  double Dl        = cosm.D_A(redshift);
  
  double theta = radius/cosm.D_A(redshift); // in radians

  // Integrand function of the integral defining the 2-halo density profile
  std::function<double(double)> Func = [&redshift,&Dl,&method_Pk,&theta,&cosm] (double l)    
				       {
					 l = pow(10,l);
					 
					 double X  = theta*l;
					 double kl = l/((1+redshift)*Dl);         
					 double Pk = cosm.Pk_matter(kl, method_Pk, false, redshift, false);
					 double J2 = cbl::j2(X);    
					 return l*J2*Pk * l;
				       };

  const double integrand = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(1.e-6), log10(1.e6)) * log(10);
  
  return 1.e-12 * Bias*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+redshift,3)) * integrand; // h*M_sun/pc^2
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_NFW_truncated (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_cosmology(cosmo);
  
  cluster.set_mass(pow(pp->logM_base, parameter[pp->Cpar.size()+1])*pp->mass_pivot);
  cluster.set_f_off(parameter[pp->Cpar.size()+2]);
  cluster.set_sigma_off(parameter[pp->Cpar.size()+3]);

  cluster.set_concentration(pp->conc_func(parameter[pp->Cpar.size()], cluster));
  parameter[pp->Cpar.size()] = cluster.concentration();

  std::vector<double> density_profile(radius.size());
  
  for (size_t j=0; j<radius.size(); j++)
    density_profile[j] = NFW_truncated(cosmo, cluster, radius[j], pp->redshift, pp->contrast, pp->trunc_fact);

  return density_profile;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_NFW_truncated_2halo (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_cosmology(cosmo);
  
  cluster.set_mass(pow(pp->logM_base, parameter[pp->Cpar.size()+1])*pp->mass_pivot);
  cluster.set_f_off(parameter[pp->Cpar.size()+2]);
  cluster.set_sigma_off(parameter[pp->Cpar.size()+3]);

  cluster.set_concentration(pp->conc_func(parameter[pp->Cpar.size()], cluster));
  parameter[pp->Cpar.size()] = cluster.concentration();

  std::vector<double> density_profile(radius.size());
  
  for (size_t j=0; j<radius.size(); j++)
    density_profile[j] = NFW_truncated(cosmo, cluster, radius[j], pp->redshift, pp->contrast, pp->trunc_fact) + two_halo(cosmo, radius[j], cluster.mass(), pp->redshift, pp->contrast, pp->bias_author, pp->method_Pk, pp->interp_type);

  return density_profile;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_NFW_truncated_scaling_relation (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);
  
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_cosmology(cosmo);
  
  cluster.set_f_off(parameter[pp->Cpar.size()+1]);
  cluster.set_sigma_off(parameter[pp->Cpar.size()+2]);

  // compute the mass from the scaling relation, and set its value in the Cluster object
  std::vector<double> scalRel_pars;
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    scalRel_pars.push_back(parameter[i]);
  for (int i=3; i<11; i++)
    scalRel_pars.emplace_back(parameter[pp->Cpar.size()+i]);
  
  double logMass = (pp->scaling_relation)->likelihood()->get_m_model()->operator()(pp->mass_proxy, scalRel_pars); 
  parameter[pp->Cpar.size()+11] = logMass;

  double Mpiv = pp->mass_pivot;
  cluster.set_mass(pow(pp->logM_base, logMass)*Mpiv);

  // set the concentration
  cluster.set_concentration(pp->conc_func(parameter[pp->Cpar.size()], cluster));
  parameter[pp->Cpar.size()] = cluster.concentration();

  // Compute the density profile
  std::vector<double> density_profile(radius.size());
  
  for (size_t j=0; j<radius.size(); j++)
    density_profile[j] = NFW_truncated(cosmo, cluster, radius[j], pp->redshift, pp->contrast, pp->trunc_fact);
  
  return density_profile;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::model_NFW_truncated_2halo_scaling_relation (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_cosmology(cosmo);
  
  cluster.set_f_off(parameter[pp->Cpar.size()+1]);
  cluster.set_sigma_off(parameter[pp->Cpar.size()+2]);

  // compute the mass from the scaling relation, and set its value in the Cluster object
  std::vector<double> scalRel_pars;
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    scalRel_pars.push_back(parameter[i]);
  for (int i=3; i<11; i++)
    scalRel_pars.emplace_back(parameter[pp->Cpar.size()+i]);
  
  double logMass = (pp->scaling_relation)->likelihood()->get_m_model()->operator()(pp->mass_proxy, scalRel_pars);
  parameter[pp->Cpar.size()+11] = logMass;
  
  double Mpiv = pp->mass_pivot;
  cluster.set_mass(pow(pp->logM_base, logMass)*Mpiv);

  // set the concentration
  cluster.set_concentration(pp->conc_func(parameter[pp->Cpar.size()], cluster));
  parameter[pp->Cpar.size()] = cluster.concentration();

  // Compute the density profile
  std::vector<double> density_profile(radius.size());
  
  for (size_t j=0; j<radius.size(); j++)
    density_profile[j] = NFW_truncated(cosmo, cluster, radius[j], pp->redshift, pp->contrast, pp->trunc_fact) + two_halo(cosmo, radius[j], cluster.mass(), pp->redshift, pp->contrast, pp->bias_author, pp->method_Pk, pp->interp_type);

  return density_profile;
}

// ===========================================================================================
