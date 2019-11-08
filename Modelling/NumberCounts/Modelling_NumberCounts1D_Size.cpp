/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  Modelling/NumberCounts/Modelling_NumberCounts1D_Size.cpp
 *
 *  @brief Methods of the class Modelling_NumberCounts1D_Size
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_NumberCounts1D_Size, i.e. the functions to model
 *  the size number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts1D_Size.h"
#include "Modelling_NumberCounts1D_Size.h"

using namespace std;

using namespace cbl;

using namespace modelling::numbercounts;


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_Size::set_model_NumberCounts_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_params, const std::vector<statistics::PriorDistribution> cosmo_param_priors)
{
  m_data_model_SF.Cpar = cosmo_params;

  const size_t nParams = cosmo_params.size();
  checkDim(cosmo_params, nParams, "cosmoPar_prior");

  vector<statistics::ParameterType> cosmoPar_type(nParams, statistics::ParameterType::_Base_);
  vector<string> cosmoPar_string(nParams);

  for (size_t i=0; i<nParams; i++)
    cosmoPar_string[i] = CosmologicalParameter_name(cosmo_params[i]);

  // input data used to construct the model
  auto inputs = make_shared<STR_NCSF_data_model>(m_data_model_SF);

  // set prior
  m_set_prior(cosmo_param_priors);

  // construct the model
  if (m_HistogramType==glob::HistogramType::_dn_dlnV_)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&size_function_model, nParams, cosmoPar_type, cosmoPar_string, inputs));
  else  ErrorCBL("no such a variable in the list!", "set_model_NumberCounts_cosmology", "Modelling_NumberCounts1D_Size.cpp");
  
}

// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_Size::set_model_NumberCounts_cosmology_and_bias (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_params, const std::vector<statistics::PriorDistribution> cosmo_param_priors, const std::vector<statistics::PriorDistribution> bias_priors)
{
  // set the model parameters
  const size_t nParams = cosmo_params.size()+3;
  checkDim(cosmo_params, cosmo_params.size(), "cosmoPar&bias_prior");

  vector<statistics::ParameterType> Par_type(nParams, statistics::ParameterType::_Base_);

  vector<string> Par_name(nParams);
  vector<statistics::PriorDistribution> priors(nParams);
  
  for (size_t i=0; i<cosmo_params.size(); i++){
    Par_name[i] = CosmologicalParameter_name(cosmo_params[i]);
    priors[i] = cosmo_param_priors[i];}

  Par_name[nParams-3] = "b_eff";
  priors[nParams-3] = bias_priors[0];
  Par_name[nParams-2] = "b_slope";
  priors[nParams-2] = bias_priors[1];
  Par_name[nParams-1] = "b_offset";
  priors[nParams-1] = bias_priors[2];

  // add the cosmological parameters
  m_data_model_SF.Cpar = cosmo_params;

  // input data used to construct the model
  auto inputs = make_shared<STR_NCSF_data_model>(m_data_model_SF);

  // set prior
  m_set_prior(priors);

  // construct the model
  if (m_HistogramType==glob::HistogramType::_dn_dlnV_)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&size_function_model, nParams, Par_type, Par_name, inputs));
  else  ErrorCBL("no such a variable in the list!", "set_model_NumberCounts_cosmology_and_bias", "Modelling_NumberCounts1D_Size.cpp");
  
}
