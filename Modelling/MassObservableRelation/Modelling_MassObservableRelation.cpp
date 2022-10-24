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

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_data_model (const cosmology::Cosmology cosmology, const std::vector<double> redshift, const double redshift_pivot, const double proxy_pivot, const double log_base)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  
  m_data_model.redshift = redshift;
  m_data_model.redshift_pivot = redshift_pivot;
  
  m_data_model.proxy_pivot = proxy_pivot;
  m_data_model.log_base = log_base;
}


// ===========================================================================================


void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_data_model (const cosmology::Cosmology cosmology, const std::vector<double> redshift, const double redshift_pivot, const double proxy_pivot, const double log_base, const std::vector<double> Nclusters)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  
  m_data_model.redshift = redshift;
  m_data_model.redshift_pivot = redshift_pivot;
  
  m_data_model.proxy_pivot = proxy_pivot;
  m_data_model.log_base = log_base;

  m_data_model.Nclusters = Nclusters;
  if (Nclusters.size() != redshift.size())
    ErrorCBL("The vectors Nclusters and redshift must have the same size!","set_data_model","Modelling_MassObservableRelation.cpp");

  m_isSet_Nclusters = true;
}


// ===========================================================================================


void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_model_MassObservableRelation_cosmology (const std::string z_evo, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior)
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
  if (z_evo == "E_z")
    m_data_model.fz = [] (const double z, const double z_piv, const std::shared_ptr<void> cosmo) {cbl::cosmology::Cosmology cosmology = *std::static_pointer_cast<cbl::cosmology::Cosmology>(cosmo); return cosmology.HH(z)/cosmology.HH(z_piv);};
  else if (z_evo == "direct")
    m_data_model.fz = [] (const double z, const double z_piv, const std::shared_ptr<void> cosmo) {(void)cosmo; return (1+z)/(1+z_piv);};
  else
    ErrorCBL("Wrong redshift evolution!","set_model_MassObservableRelation_cosmology","Modelling_MassObservableRelation.cpp"); 

  // set prior
  m_set_prior(param_prior);
  
  // set the errors on z and proxy, if provided,
  // and construct the model
  auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);
  if (m_isSet_Nclusters)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_scaling_relation_int, nParams, Par_type, Par_string, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_scaling_relation, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================


double cbl::modelling::massobsrel::scaling_relation (const double alpha, const double beta, const double gamma, const double div_proxy, const double f_z)
{  
  return alpha + beta * div_proxy + gamma * f_z;
}


// ===========================================================================================


std::vector<double> cbl::modelling::massobsrel::model_scaling_relation (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  auto cosmo_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmo);
  std::vector<double> res(proxy.size());
  for (size_t j=0; j<proxy.size(); j++) {
    double log1 = log(proxy[j]/pp->proxy_pivot)/log(pp->log_base);
    double log2 = log(pp->fz(pp->redshift[j], pp->redshift_pivot, cosmo_ptr))/log(pp->log_base);
    res[j] = scaling_relation(parameter[pp->Cpar.size()], parameter[pp->Cpar.size()+1], parameter[pp->Cpar.size()+2], log1, log2);
  }

  return res;
}


// ===========================================================================================

std::vector<double> cbl::modelling::massobsrel::model_scaling_relation_int (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  auto cosmo_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmo);
  
  // define the integrand

  double dummy_Nclusters, dummy_proxy_eff, dummy_z_eff;
  std::shared_ptr<void> ptr;
  
  auto integrand = [&] (const std::vector<double> x)
		   {		     
		     // Compute P(M|lambda,z)
		     double log_lambda = log(dummy_proxy_eff/pp->proxy_pivot)/log(pp->log_base);
		     double log_f_z = log( pp->fz(dummy_z_eff, pp->redshift_pivot, cosmo_ptr) )/log(pp->log_base);
      
		     double mean = parameter[pp->Cpar.size()] + parameter[pp->Cpar.size()+1]*log_lambda + parameter[pp->Cpar.size()+2]*log_f_z;
		     double sigma = parameter[pp->Cpar.size()+3] + parameter[pp->Cpar.size()+4]*pow(log_lambda, parameter[pp->Cpar.size()+5]) + parameter[pp->Cpar.size()+6]*pow(log_f_z, parameter[pp->Cpar.size()+7]);
		     double P_M__lambda_z = (cbl::gaussian(x[0], ptr, {mean,sigma/dummy_Nclusters}));
      
		     return  x[0] * P_M__lambda_z;
		   };

  // compute the model  
  std::vector<double> model(proxy.size());
  
  for (size_t j=0; j<proxy.size(); j++) {
    dummy_proxy_eff = proxy[j];
    dummy_z_eff = pp->redshift[j];
    dummy_Nclusters = pp->Nclusters[j];

    double logLambda = log(proxy[j]/pp->proxy_pivot)/log(pp->log_base);
    double logFz = log(pp->fz(pp->redshift[j], pp->redshift_pivot, cosmo_ptr))/log(pp->log_base);

    double logM = parameter[pp->Cpar.size()] + parameter[pp->Cpar.size()+1]*logLambda + parameter[pp->Cpar.size()+2]*logFz;
    double sigma_intr = parameter[pp->Cpar.size()+3] + parameter[pp->Cpar.size()+4]*pow(logLambda, parameter[pp->Cpar.size()+5]) + parameter[pp->Cpar.size()+6]*pow(logFz, parameter[pp->Cpar.size()+7]);
    sigma_intr = sigma_intr/pp->Nclusters[j];

    std::vector<std::vector<double>> integration_limits(1);
    integration_limits[0] = {logM-3.5*sigma_intr, logM+3.5*sigma_intr};
    
    cbl::wrapper::cuba::CUBAwrapper CW (integrand, (int)(integration_limits.size()));
    model[j] = CW.IntegrateVegas(integration_limits,false);    
  }

  return model;
}
