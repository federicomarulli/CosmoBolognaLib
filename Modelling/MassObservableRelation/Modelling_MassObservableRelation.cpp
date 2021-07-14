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
 *  Modelling/MassObservableRelation/Modelling_MassObservableRelation.cpp
 *
 *  @brief Methods of the class Modelling_MassObservableRelation
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_MassObservableRelation, i.e. the common functions to model
 *  the galaxy cluster mass - mass proxy relation
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "Modelling_MassObservableRelation.h"
#include "ModelFunction_NumberCounts.h"

using namespace std;

using namespace cbl;


// ===========================================================================================

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const std::vector<double> redshift, const double redshift_pivot, const double proxy_or_mass_pivot, const double log_base)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  m_data_model.cluster = make_shared<catalogue::Cluster>(cluster);
  
  m_data_model.redshift = redshift;
  m_data_model.redshift_pivot = redshift_pivot;
  
  m_data_model.proxy_or_mass_pivot = proxy_or_mass_pivot;
  m_data_model.log_base = log_base;
}



// ===========================================================================================


void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const std::vector<double> z_edges, const std::vector<std::vector<double>> proxy_edges, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double log_base, const std::string method_Pk, const bool store_output, const int norm, const double Delta, const bool isDelta_vir, const std::string model_MF, const double area_degrees, const double prec)
{
  if (proxy_edges.size() != z_edges.size()-1)
    ErrorCBL("The number of proxy edges sets must match the number of redshift bins!","set_data_model","Modelling_MassObservableRelation.cpp");
  for (size_t i=0; i<proxy_edges.size(); i++)
    for (size_t j=0; j<proxy_edges[i].size(); j++)
      if (proxy_edges[i][j] <= 0.)
	ErrorCBL("The values of the proxy edges cannot be <= 0.","set_data_model","Modelling_MassObservableRelation.cpp");
  
  m_data_model.edges_proxy = proxy_edges;
  m_data_model.edges_z = z_edges;

  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cluster = make_shared<catalogue::Cluster>(cluster);
  
  m_data_model.redshift_pivot = z_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  m_data_model.mass_pivot = mass_pivot;
  m_data_model.log_base = log_base;
  
  m_data_model.method_Pk = method_Pk;
  m_data_model.k_min = 1.e-4;
  m_data_model.k_max = 100.;
  m_data_model.step = 500;
  m_data_model.kk = logarithmic_bin_vector(500, 1.e-4, 100.);
  m_data_model.norm = norm;
  
  m_data_model.store_output = store_output;
  m_data_model.output_root = "test";
  m_data_model.file_par = par::defaultString;

  m_data_model.isDelta_Vir = isDelta_vir;
  m_data_model.Delta = Delta;
  m_data_model.model_MF = model_MF;

  m_data_model.Mass_vector = logarithmic_bin_vector(200, 1.e10, 1.e16);

  m_data_model.prec = prec;

  m_data_model.area_rad = area_degrees*pow(par::pi/180., 2);

  m_data_model.is_sigma8_free = false;
}


// ===========================================================================================


void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_model_MassObservableRelation_cosmology (const std::string scalrel_z_evo, const std::string z_error_type, const std::string proxy_error_type, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const statistics::PriorDistribution proxy_bias_prior, const statistics::PriorDistribution z_error_prior, const statistics::PriorDistribution proxy_error_prior)
{  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+12;

  vector<statistics::ParameterType> Par_type(nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string(nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++) {
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_param_prior[i];
  }

  // Set the names and priors for the mass-observable relation parameters, and for P(lambda|z)
  Par_string[cosmo_param.size()] = "alpha";
  param_prior[cosmo_param.size()] = alpha_prior;
  Par_string[cosmo_param.size()+1] = "beta";
  param_prior[cosmo_param.size()+1] = beta_prior;
  Par_string[cosmo_param.size()+2] = "gamma";
  param_prior[cosmo_param.size()+2] = gamma_prior;
  Par_string[cosmo_param.size()+3] = "scatter0";
  param_prior[cosmo_param.size()+3] = scatter0_prior;
  Par_string[cosmo_param.size()+4] = "scatterM";
  param_prior[cosmo_param.size()+4] = scatterM_prior;
  Par_string[cosmo_param.size()+5] = "scatterM_exponent";
  param_prior[cosmo_param.size()+5] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+6] = "scatterz";
  param_prior[cosmo_param.size()+6] = scatterz_prior;
  Par_string[cosmo_param.size()+7] = "scatterz_exponent";
  param_prior[cosmo_param.size()+7] = scatterz_exponent_prior;
  Par_string[cosmo_param.size()+8] = "z_bias";
  param_prior[cosmo_param.size()+8] = z_bias_prior;
  Par_string[cosmo_param.size()+9] = "proxy_bias";
  param_prior[cosmo_param.size()+9] = proxy_bias_prior;
  Par_string[cosmo_param.size()+10] = "z_error";
  param_prior[cosmo_param.size()+10] = z_error_prior;
  Par_string[cosmo_param.size()+11] = "proxy_error";
  param_prior[cosmo_param.size()+11] = proxy_error_prior;

  // Set the functional form for the redshift evolution in the scaling relation 
  if (scalrel_z_evo == "E_z")
    m_data_model.fz = &Fz_Ez;
  else if (scalrel_z_evo == "direct")
    m_data_model.fz = &Fz_direct;
  else
    cbl::ErrorCBL("Error in the input parameter scalrel_z_evo: no such a possibility for f(z)!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp");

  // Set the error types for redshift and mass proxy
  if (z_error_type == "relative")
    m_data_model.z_error = &AbsoluteFromRelativeError;
  else if (z_error_type == "absolute")
    m_data_model.z_error = &ReturnAbsoluteError;
  else
    cbl::ErrorCBL("Error in the input parameter z_error_type: choose between \"relative\" and \"absolute\"!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp");

  if (proxy_error_type == "relative")
    m_data_model.proxy_error = &AbsoluteFromRelativeError;
  else if (proxy_error_type == "absolute")
    m_data_model.proxy_error = &ReturnAbsoluteError;
  else
    cbl::ErrorCBL("Error in the input parameter proxy_error_type: choose between \"relative\" and \"absolute\"!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp");

  // input data used to construct the model
  auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_mass_from_counts, nParams, Par_type, Par_string, inputs));
  
}



// ===========================================================================================

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_model_MassObservableRelation_cosmology (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<double> z_eff_err, const std::vector<double> proxy_eff_err)
{
  m_data_model.Cpar = {};
  
  const size_t nParams = 8;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors for the mass-observable relation parameters
  Par_string[0] = "alpha";
  param_prior[0] = alpha_prior;
  Par_string[1] = "beta";
  param_prior[1] = beta_prior;
  Par_string[2] = "gamma";
  param_prior[2] = gamma_prior;
  Par_string[3] = "scatter0";
  param_prior[3] = scatter0_prior;
  Par_string[4] = "scatterM";
  param_prior[4] = scatterM_prior;
  Par_string[5] = "scatterM_exponent";
  param_prior[5] = scatterM_exponent_prior;
  Par_string[6] = "scatterz";
  param_prior[6] = scatterz_prior;
  Par_string[7] = "scatterz_exponent";
  param_prior[7] = scatterz_exponent_prior;

  // set the redshift evolution function
  m_data_model.fz = &Fz_direct;  

  // set prior
  m_set_prior(param_prior);

  // set the errors on z and proxy/mass, if provided
  // and construct the model
  if (z_eff_err.size() == 0 && proxy_eff_err.size() == 0) {
    auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_scaling_relation, nParams, Par_type, Par_string, inputs));
  }
  else if (z_eff_err.size() > 0 && proxy_eff_err.size() > 0) {
    if (z_eff_err.size() != m_data_model.redshift.size() || proxy_eff_err.size() != m_data_model.redshift.size())
      ErrorCBL("The vectors of errors on redshift and proxy must have the same size of the redshift and proxy vectors!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp"); 
    m_data_model.z_eff_err = z_eff_err;
    m_data_model.proxy_eff_err = proxy_eff_err;
    auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_scaling_relation_with_errors, nParams, Par_type, Par_string, inputs));
  } else
    ErrorCBL("You must set (or not) both the vectors of errors on redshift and proxy!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp");
}



// ===========================================================================================

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_model_MassObservableRelation_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<double> z_eff_err, const std::vector<double> proxy_eff_err)
{
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+8; // The total number of parameters is given by the cosmological ones + 8

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the mass-observable relation parameters
  Par_string[cosmo_param.size()] = "alpha";
  param_prior[cosmo_param.size()] = alpha_prior;
  Par_string[cosmo_param.size()+1] = "beta";
  param_prior[cosmo_param.size()+1] = beta_prior;
  Par_string[cosmo_param.size()+2] = "gamma";
  param_prior[cosmo_param.size()+2] = gamma_prior;
  Par_string[cosmo_param.size()+3] = "scatter0";
  param_prior[cosmo_param.size()+3] = scatter0_prior;
  Par_string[cosmo_param.size()+4] = "scatterM";
  param_prior[cosmo_param.size()+4] = scatterM_prior;
  Par_string[cosmo_param.size()+5] = "scatterM_exponent";
  param_prior[cosmo_param.size()+5] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+6] = "scatterz";
  param_prior[cosmo_param.size()+6] = scatterz_prior;
  Par_string[cosmo_param.size()+7] = "scatterz_exponent";
  param_prior[cosmo_param.size()+7] = scatterz_exponent_prior;

  // set the redshift evolution function
  m_data_model.fz = &Fz_Ez;

  // set prior
  m_set_prior(param_prior);
  
  // set the errors on z and proxy/mass, if provided
  // and construct the model
  if (z_eff_err.size() == 0 && proxy_eff_err.size() == 0) {
    auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_scaling_relation, nParams, Par_type, Par_string, inputs));
  }
  else if (z_eff_err.size() > 0 && proxy_eff_err.size() > 0) {
    if (z_eff_err.size() != m_data_model.redshift.size() || proxy_eff_err.size() != m_data_model.redshift.size())
      ErrorCBL("The vectors of errors on redshift and proxy must have the same size of the redshift and proxy vectors!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp"); 
    m_data_model.z_eff_err = z_eff_err;
    m_data_model.proxy_eff_err = proxy_eff_err;
    auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_scaling_relation_with_errors, nParams, Par_type, Par_string, inputs));
  } else
    ErrorCBL("You must set (or not) both the vectors of errors on redshift and proxy!","set_data_model","Modelling_MassObservableRelation.cpp");
}


// ===========================================================================================


double cbl::modelling::massobsrel::scaling_relation (cbl::catalogue::Cluster cluster, const double div_proxy_or_mass, const double f_z)
{  
  return cluster.alpha_scaling_rel() + cluster.beta_scaling_rel() * div_proxy_or_mass + cluster.gamma_scaling_rel() * f_z;
}


// ===========================================================================================


double cbl::modelling::massobsrel::mass_from_counts (double (*fz)(std::vector<double>, std::shared_ptr<void>), double (*z_error)(std::vector<double>), double (*proxy_error)(std::vector<double>), const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, cbl::cosmology::Cosmology cosmology, cbl::catalogue::Cluster cluster, const double Area, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const cbl::glob::FuncGrid interp_sigmaM, const  cbl::glob::FuncGrid interp_DlnsigmaM, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base)
{  
  double fact = (cosmology.unit()) ? 1 : cosmology.hh();
  std::shared_ptr<void> pp;
  auto cosmology_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmology);

  double the_minM = 0;
  
  // M_tot integrand
  auto integrand_Mtot = [&the_minM,&z_error,&proxy_error,&fz,&cosmology_ptr,&cosmology,&cluster,&log_base,&proxy_pivot,&z_pivot,&mass_pivot,&proxy_min,&proxy_max,&redshift_min,&redshift_max,&interp_sigmaM,&interp_DlnsigmaM,&store_output,&fact,&model_MF,&Delta,&Area,&isDelta_vir,&pp] (const std::vector<double> x)
    {
      double Delta_ = (isDelta_vir) ? cosmology.Delta_vir(Delta, x[1]) : Delta;
      double Mass = x[0]*mass_pivot*pow(log_base,the_minM);
      double normM = x[0]*pow(log_base,the_minM);

      // Compute P(lambda|M,z)
      double log_M = log(normM)/log(log_base);
      double log_f_z = log( fz({x[1], z_pivot}, cosmology_ptr) )/log(log_base);

      double mean = cluster.alpha_scaling_rel() + cluster.beta_scaling_rel()*log_M + cluster.gamma_scaling_rel()*log_f_z + log(proxy_pivot)/log(log_base);
      double sigma = cluster.scatter0_scaling_rel() + cluster.scatterM_scaling_rel()*pow(log_M, cluster.scatterM_exponent_scaling_rel()) + cluster.scatterz_scaling_rel()*pow(log_f_z, cluster.scatterz_exponent_scaling_rel());
      double P_lambda__M_z = (cbl::gaussian(log(x[2])/log(log_base), pp, {mean,sigma})) / (x[2]*log(log_base));
      
      // Compute the integrals of P(z|z) and P(lambda|lambda)
      double mean_Pz = x[1] + cluster.zbias() * (1+x[1]);
      double int_P_z = 0.5 * ( erf( (redshift_max - mean_Pz) / (sqrt(2)*z_error({cluster.zerror(), redshift_max})) ) - erf( (redshift_min - mean_Pz) / (sqrt(2)*z_error({cluster.zerror(), redshift_min})) ) );
      double mean_Plambda = x[2] + cluster.proxybias() * (x[2]);
      double int_P_lambda = 0.5 * ( erf( (proxy_max - mean_Plambda) / (sqrt(2)*proxy_error({cluster.proxyerror(), proxy_max})) ) - erf( (proxy_min - mean_Plambda) / (sqrt(2)*proxy_error({cluster.proxyerror(), proxy_min})) ) );
      
      return normM * cosmology.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), x[1], model_MF, store_output, cbl::par::defaultString, Delta_)*Area*cosmology.dV_dZdOmega(x[1], true) * P_lambda__M_z * int_P_z * int_P_lambda;
    };  
  
  // Pure counts integrand
  auto integrand_counts = [&the_minM,&z_error,&proxy_error,&fz,&cosmology_ptr,&cosmology,&cluster,&log_base,&proxy_pivot,&z_pivot,&mass_pivot,&proxy_min,&proxy_max,&redshift_min,&redshift_max,&interp_sigmaM,&interp_DlnsigmaM,&store_output,&fact,&model_MF,&Delta,&Area,&isDelta_vir,&pp] (const std::vector<double> x)
    {
      double Delta_ = (isDelta_vir) ? cosmology.Delta_vir(Delta, x[1]) : Delta;
      double Mass = x[0]*mass_pivot*pow(log_base,the_minM);
      double normM = x[0]*pow(log_base,the_minM);

      // Compute P(lambda|M,z)
      double log_M = log(normM)/log(log_base);
      double log_f_z = log( fz({x[1], z_pivot}, cosmology_ptr) )/log(log_base);

      double mean = cluster.alpha_scaling_rel() + cluster.beta_scaling_rel()*log_M + cluster.gamma_scaling_rel()*log_f_z + log(proxy_pivot)/log(log_base);
      double sigma = cluster.scatter0_scaling_rel() + cluster.scatterM_scaling_rel()*pow(log_M, cluster.scatterM_exponent_scaling_rel()) + cluster.scatterz_scaling_rel()*pow(log_f_z, cluster.scatterz_exponent_scaling_rel());
      double P_lambda__M_z = (cbl::gaussian(log(x[2])/log(log_base), pp, {mean,sigma})) / (x[2]*log(log_base));
      
      // Compute the integrals of P(z|z) and P(lambda|lambda)
      double mean_Pz = x[1] + cluster.zbias() * (1+x[1]);
      double int_P_z = 0.5 * ( erf( (redshift_max - mean_Pz) / (sqrt(2)*z_error({cluster.zerror(), redshift_max})) ) - erf( (redshift_min - mean_Pz) / (sqrt(2)*z_error({cluster.zerror(), redshift_min})) ) );
      double mean_Plambda = x[2] + cluster.proxybias() * (x[2]);
      double int_P_lambda = 0.5 * ( erf( (proxy_max - mean_Plambda) / (sqrt(2)*proxy_error({cluster.proxyerror(), proxy_max})) ) - erf( (proxy_min - mean_Plambda) / (sqrt(2)*proxy_error({cluster.proxyerror(), proxy_min})) ) );
      
      return cosmology.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), x[1], model_MF, store_output, cbl::par::defaultString, Delta_)*Area*cosmology.dV_dZdOmega(x[1], true) * P_lambda__M_z * int_P_z * int_P_lambda;
    };
  
  // -------------------------------------------------------------

  // Find the minimum and maximum masses, given the parameters of the scaling relation
  double scatter = cluster.scatter0_scaling_rel();
  
  double log_lambda_min = std::max(log(proxy_min/proxy_pivot)/log(log_base) - 4*scatter, log(0.0001/proxy_pivot)/log(log_base));
  double log_lambda_max = log(proxy_max/proxy_pivot)/log(log_base) + 4*scatter;
  double log_f_z_min = log( fz({redshift_min, z_pivot}, cosmology_ptr) )/log(log_base);
  double log_f_z_max = log( fz({redshift_max, z_pivot}, cosmology_ptr) )/log(log_base);

  double M1 = (- cluster.alpha_scaling_rel() + log_lambda_min - cluster.gamma_scaling_rel()*log_f_z_min) / cluster.beta_scaling_rel();
  double M2 = (- cluster.alpha_scaling_rel() + log_lambda_max - cluster.gamma_scaling_rel()*log_f_z_min) / cluster.beta_scaling_rel();
  double M3 = (- cluster.alpha_scaling_rel() + log_lambda_min - cluster.gamma_scaling_rel()*log_f_z_max) / cluster.beta_scaling_rel();
  double M4 = (- cluster.alpha_scaling_rel() + log_lambda_max - cluster.gamma_scaling_rel()*log_f_z_max) / cluster.beta_scaling_rel();

  double min1 = std::min(M1, M2);
  double min2 = std::min(min1, M3);
  double minM_ = std::min(min2, M4);
  double max1 = std::max(M1, M2);
  double max2 = std::max(max1, M3);
  double maxM_ = std::max(max2, M4);

  double minM = std::max(minM_, log(1.e10/mass_pivot)/log(log_base));
  double maxM = std::min(maxM_, log(1.e16/mass_pivot)/log(log_base));

  the_minM = minM;

  
  // Define the integral limits
  int integral_dimension=3;
  std::vector<std::vector<double>> integration_limits(integral_dimension);
  integration_limits[0] = {pow(log_base,minM)/pow(log_base,minM), pow(log_base,maxM)/pow(log_base,minM)};
  integration_limits[1] = {std::max(redshift_min - 3.5*z_error({cluster.zerror(), redshift_min}), 0.), redshift_max + 3.5*z_error({cluster.zerror(), redshift_max})};
  integration_limits[2] = {std::max(proxy_min - 3.5*proxy_error({cluster.proxyerror(), proxy_min}), 0.00001), proxy_max + 3.5*proxy_error({cluster.proxyerror(), proxy_max})};

  // Compute the M_tot integral
  cbl::wrapper::cuba::CUBAwrapper CW (integrand_Mtot, integral_dimension);
  double nc = CW.IntegrateVegas(integration_limits,false);

  // Divide it by the counts
  cbl::wrapper::cuba::CUBAwrapper CW2 (integrand_counts, integral_dimension);
  double counts = CW2.IntegrateVegas(integration_limits,false);

  if (integration_limits[0][0] < integration_limits[0][1] && counts > 0) {
    nc /= counts;
    return log(nc)/log(log_base);
  } else
    return cbl::par::defaultDouble;
}


// ===========================================================================================

std::vector<double> cbl::modelling::massobsrel::model_scaling_relation (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_alpha_scaling_rel(parameter[pp->Cpar.size()]);
  cluster.set_beta_scaling_rel(parameter[pp->Cpar.size()+1]);
  cluster.set_gamma_scaling_rel(parameter[pp->Cpar.size()+2]);
  cluster.set_scatter0_scaling_rel(parameter[pp->Cpar.size()+3]);
  cluster.set_scatterM_scaling_rel(parameter[pp->Cpar.size()+4]);
  cluster.set_scatterM_exponent_scaling_rel(parameter[pp->Cpar.size()+5]);
  cluster.set_scatterz_scaling_rel(parameter[pp->Cpar.size()+6]);
  cluster.set_scatterz_exponent_scaling_rel(parameter[pp->Cpar.size()+7]);

  auto cosmo_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmo);
  std::vector<double> res(proxy_or_mass.size());
  for (size_t j=0; j<proxy_or_mass.size(); j++) {
    double log1 = log(proxy_or_mass[j]/pp->proxy_or_mass_pivot)/log(pp->log_base);
    double log2 = log(pp->fz({pp->redshift[j], pp->redshift_pivot}, cosmo_ptr))/log(pp->log_base);
    res[j] = scaling_relation(cluster, log1, log2);
  }

  return res;
}


// ===========================================================================================

std::vector<double> cbl::modelling::massobsrel::model_scaling_relation_with_errors (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_alpha_scaling_rel(parameter[pp->Cpar.size()]);
  cluster.set_beta_scaling_rel(parameter[pp->Cpar.size()+1]);
  cluster.set_gamma_scaling_rel(parameter[pp->Cpar.size()+2]);
  cluster.set_scatter0_scaling_rel(parameter[pp->Cpar.size()+3]);
  cluster.set_scatterM_scaling_rel(parameter[pp->Cpar.size()+4]);
  cluster.set_scatterM_exponent_scaling_rel(parameter[pp->Cpar.size()+5]);
  cluster.set_scatterz_scaling_rel(parameter[pp->Cpar.size()+6]);
  cluster.set_scatterz_exponent_scaling_rel(parameter[pp->Cpar.size()+7]);

  auto cosmo_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmo);

  std::shared_ptr<void> gauss_ptr;
  
  std::vector<double> res(proxy_or_mass.size());
  for (size_t j=0; j<proxy_or_mass.size(); j++) {
    
    std::function<double(double)> P_proxy = [&cosmo_ptr,&pp,&gauss_ptr,&j,&proxy_or_mass] (const double xx) { return log(xx/pp->proxy_or_mass_pivot)/log(pp->log_base)*cbl::gaussian(xx, gauss_ptr, {proxy_or_mass[j],pp->proxy_eff_err[j]}); };
    std::function<double(double)> P_z = [&cosmo_ptr,&pp,&gauss_ptr,&j] (const double xx) { return log(pp->fz({xx,pp->redshift_pivot}, cosmo_ptr))/log(pp->log_base)*cbl::gaussian(xx, gauss_ptr, {pp->redshift[j],pp->z_eff_err[j]}); };
    
    double int1 = cbl::wrapper::gsl::GSL_integrate_qag(P_proxy, std::max(proxy_or_mass[j]-3.5*pp->proxy_eff_err[j],0.), proxy_or_mass[j]+3.5*pp->proxy_eff_err[j], 1.e-5, 1000, 6);
    double int2 = cbl::wrapper::gsl::GSL_integrate_qag(P_z, std::max(pp->redshift[j]-3.5*pp->z_eff_err[j],0.), pp->redshift[j]+3.5*pp->z_eff_err[j], 1.e-5, 1000, 6);
    res[j] = scaling_relation(cluster, int1, int2);
    
  }

  return res;
}


// ===========================================================================================


std::vector<double> cbl::modelling::massobsrel::model_mass_from_counts (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)proxy;
  
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  // set the cluster parameters
  cluster.set_alpha_scaling_rel(parameter[pp->Cpar.size()]);
  cluster.set_beta_scaling_rel(parameter[pp->Cpar.size()+1]);
  cluster.set_gamma_scaling_rel(parameter[pp->Cpar.size()+2]);
  cluster.set_scatter0_scaling_rel(parameter[pp->Cpar.size()+3]);
  cluster.set_scatterM_scaling_rel(parameter[pp->Cpar.size()+4]);
  cluster.set_scatterM_exponent_scaling_rel(parameter[pp->Cpar.size()+5]);
  cluster.set_scatterz_scaling_rel(parameter[pp->Cpar.size()+6]);
  cluster.set_scatterz_exponent_scaling_rel(parameter[pp->Cpar.size()+7]);
  cluster.set_zbias(parameter[pp->Cpar.size()+8]);
  cluster.set_proxybias(parameter[pp->Cpar.size()+9]);
  cluster.set_zerror(parameter[pp->Cpar.size()+10]);
  cluster.set_proxyerror(parameter[pp->Cpar.size()+11]);

  // compute the power spectrum
  std::vector<double> Pk = cosmo.Pk_DM(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (pp->Mass_vector, cosmo, pp->kk, Pk, "Spline", pp->k_max);

  // Compute the counts
  std::vector<double> number_counts;

  for (size_t j=0; j<pp->edges_proxy.size(); j++)
    for (size_t k=0; k<pp->edges_proxy[j].size()-1; k++)
      number_counts.emplace_back( cbl::modelling::massobsrel::mass_from_counts(pp->fz, pp->z_error, pp->proxy_error, pp->edges_z[j], pp->edges_z[j+1], pp->edges_proxy[j][k], pp->edges_proxy[j][k+1], cosmo, cluster, pp->area_rad, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_Vir, interp[0], interp[1], pp->proxy_pivot, pp->redshift_pivot, pp->mass_pivot, pp->log_base) );

  return number_counts;
}


// ===========================================================================================


double cbl::modelling::massobsrel::Fz_Ez (const std::vector<double> x, const std::shared_ptr<void> cosmo)
{
  cbl::cosmology::Cosmology cosmology = *std::static_pointer_cast<cbl::cosmology::Cosmology>(cosmo);
  return cosmology.HH(x[0])/cosmology.HH(x[1]);
}


// ===========================================================================================


double cbl::modelling::massobsrel::Fz_direct (const std::vector<double> x, const std::shared_ptr<void> cosmo)
{
  (void)cosmo;
  return (1+x[0])/(1+x[1]);
}


// ===========================================================================================


double cbl::modelling::massobsrel::ReturnAbsoluteError (const std::vector<double> x)
{
  return x[0];
}


// ===========================================================================================


double cbl::modelling::massobsrel::AbsoluteFromRelativeError (const std::vector<double> x)
{
  return x[0]*x[1];
}
