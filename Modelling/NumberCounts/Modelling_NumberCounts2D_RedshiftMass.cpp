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
 *  Modelling/NumberCounts/Modelling_NumberCounts2D_RedshiftMass.cpp
 *
 *  @brief Methods of the class Modelling_NumberCounts2D_RedshiftMass
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_NumberCounts2D_RedshiftMass, i.e. the functions to model
 *  the mass-redshift number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts2D_RedshiftMass.h"
#include "Modelling_NumberCounts2D_RedshiftMass.h"

using namespace std;

using namespace cbl;

using namespace modelling::numbercounts;


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts2D_RedshiftMass::set_model_NumberCounts_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior)
{
  std::vector<cbl::cosmology::CosmologicalParameter> param = cosmo_param;
  
  // check if sigma8 is a free parameter
  if (find(param.begin(), param.end(), cosmology::CosmologicalParameter::_sigma8_)!=param.end())
    m_data_model.is_sigma8_free = true;
  
  if (!m_data_model.is_sigma8_free)
    param.push_back(cosmology::CosmologicalParameter::_sigma8_);
  
  m_data_model.Cpar = param;

  const size_t nParams = param.size();
  checkDim(param, nParams, "cosmoPar_prior");

  vector<statistics::ParameterType> cosmoPar_type((m_data_model.is_sigma8_free) ? nParams : nParams-1, statistics::ParameterType::_Base_);
  if (!m_data_model.is_sigma8_free)
    cosmoPar_type.push_back(statistics::ParameterType::_Derived_);
  
  vector<string> cosmoPar_string(nParams);

  for (size_t i=0; i<nParams; i++)
    cosmoPar_string[i] = CosmologicalParameter_name(param[i]);
    
  // input data used to construct the model
  auto inputs = make_shared<STR_NC_data_model>(m_data_model);
  
  // set prior
  m_set_prior(cosmo_param_prior);
  
  // construct the model
  switch (m_HistogramType) {
  case (glob::HistogramType::_N_V_):
    m_model = make_shared<statistics::Model2D>(statistics::Model2D(&number_counts_redshift_mass, nParams, cosmoPar_type, cosmoPar_string, inputs));
    break;
  case (glob::HistogramType::_n_V_):
    m_model = make_shared<statistics::Model2D>(statistics::Model2D(&number_density_redshift_mass, nParams, cosmoPar_type, cosmoPar_string, inputs));
    break;
  case (glob::HistogramType::_dn_dV_):
    m_model = make_shared<statistics::Model2D>(statistics::Model2D(&mass_function_redshift_mass, nParams, cosmoPar_type, cosmoPar_string, inputs));
    break;
  default:
    ErrorCBL("no such a variable in the list!", "set_model_NumberCounts_cosmology", "Modelling_NumberCounts2D_RedshiftMass.cpp");
  }
}
