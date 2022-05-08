/********************************************************************
 *  Copyright (C) 2021 by Giorgio Lesci and Federico Marulli        *
 *  giorgio.lesci2@unibo.it, federico.marulli3@unibo.it             *
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
 *  Modelling/NumberCounts/Modelling_NumberCounts1D_MassProxy.cpp
 *
 *  @brief Methods of the class Modelling_NumberCounts1D_MassProxy
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_NumberCounts1D_MassProxy, i.e. the functions to model
 *  the number counts as a function of a mass proxy
 *
 *  @authors Giorgio Lesci, Federico Marulli
 *
 *  @authors giorgio.lesci2@unibo.it, federico.marulli3@unibo.it
 */


#include "ModelFunction_NumberCounts1D_MassProxy.h"
#include "Modelling_NumberCounts1D_MassProxy.h"

using namespace std;

using namespace cbl;

using namespace modelling::numbercounts;


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_MassProxy::set_model_NumberCounts_cosmology (const std::string scalrel_z_evo, const std::string z_error_type, const std::string proxy_error_type, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const statistics::PriorDistribution proxy_bias_prior, const statistics::PriorDistribution z_error_prior, const statistics::PriorDistribution proxy_error_prior, const std::vector<statistics::PriorDistribution> Plambda_prior)
{
  if (Plambda_prior.size() != 3)
    cbl::ErrorCBL("Plambda_prior must have size = 3.","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");
  
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+15;

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
  Par_string[cosmo_param.size()+12] = "Plambda_a";
  param_prior[cosmo_param.size()+12] = Plambda_prior[0];
  Par_string[cosmo_param.size()+13] = "Plambda_b";
  param_prior[cosmo_param.size()+13] = Plambda_prior[1];
  Par_string[cosmo_param.size()+14] = "Plambda_c";
  param_prior[cosmo_param.size()+14] = Plambda_prior[2];

  // Set the functional form for the redshift evolution in the scaling relation 
  if (scalrel_z_evo == "E_z")
    m_data_model.fz = &fz_Ez;
  else if (scalrel_z_evo == "direct")
    m_data_model.fz = &fz_direct;
  else
    cbl::ErrorCBL("Error in the input parameter scalrel_z_evo: no such a possibility for f(z)!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");

  // Set the error types for redshift and mass proxy
  if (z_error_type == "relative")
    m_data_model.z_error = &absolute_from_relative_error;
  else if (z_error_type == "absolute")
    m_data_model.z_error = &return_absolute_error;
  else
    cbl::ErrorCBL("Error in the input parameter z_error_type: choose between \"relative\" and \"absolute\"!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");

  if (proxy_error_type == "relative")
    m_data_model.proxy_error = &absolute_from_relative_error;
  else if (proxy_error_type == "absolute")
    m_data_model.proxy_error = &return_absolute_error;
  else
    cbl::ErrorCBL("Error in the input parameter proxy_error_type: choose between \"relative\" and \"absolute\"!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");

  // input data used to construct the model
  m_data_model.transfer_func = &no_transfer;
  auto inputs = make_shared<STR_NC_data_model>(m_data_model);

  // input data used to construct the transfer function
  m_data_model.transfer_func = &bias_transfer;
  auto inputs2 = make_shared<STR_NC_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  switch (m_HistogramType) {
    
  case (glob::HistogramType::_N_V_):
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_proxy, nParams, Par_type, Par_string, inputs));
    m_transfer_func = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_proxy, nParams, Par_type, Par_string, inputs2));
    break;
    
  default:
    ErrorCBL("Only counts can be modelled! Set _N_V_ as the histogram type.", "set_model_NumberCounts_cosmology", "Modelling_NumberCounts1D_MassProxy.cpp");
    
  }
}



// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_MassProxy::set_model_NumberCounts_cosmology (const std::string scalrel_z_evo, const std::string z_error_type, const std::string proxy_error_type, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const statistics::PriorDistribution proxy_bias_prior, const statistics::PriorDistribution z_error_prior, const statistics::PriorDistribution proxy_error_prior)
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
    m_data_model.fz = &fz_Ez;
  else if (scalrel_z_evo == "direct")
    m_data_model.fz = &fz_direct;
  else
    cbl::ErrorCBL("Error in the input parameter scalrel_z_evo: no such a possibility for f(z)!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");

  // Set the error types for redshift and mass proxy
  if (z_error_type == "relative")
    m_data_model.z_error = &absolute_from_relative_error;
  else if (z_error_type == "absolute")
    m_data_model.z_error = &return_absolute_error;
  else
    cbl::ErrorCBL("Error in the input parameter z_error_type: choose between \"relative\" and \"absolute\"!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");

  if (proxy_error_type == "relative")
    m_data_model.proxy_error = &absolute_from_relative_error;
  else if (proxy_error_type == "absolute")
    m_data_model.proxy_error = &return_absolute_error;
  else
    cbl::ErrorCBL("Error in the input parameter proxy_error_type: choose between \"relative\" and \"absolute\"!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");

  // input data used to construct the model
  m_data_model.transfer_func = &no_transfer;
  auto inputs = make_shared<STR_NC_data_model>(m_data_model);

  // input data used to construct the transfer function
  m_data_model.transfer_func = &bias_transfer;
  auto inputs2 = make_shared<STR_NC_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  switch (m_HistogramType) {
    
  case (glob::HistogramType::_N_V_):
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_proxy_classic, nParams, Par_type, Par_string, inputs));
    m_transfer_func = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_proxy_classic, nParams, Par_type, Par_string, inputs2));
    break;
    
  default:
    ErrorCBL("Only counts can be modelled! Set _N_V_ as the histogram type.", "set_model_NumberCounts_cosmology", "Modelling_NumberCounts1D_MassProxy.cpp");
    
  }
}
