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


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_MassProxy::set_model_NumberCounts_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const std::vector<statistics::PriorDistribution> Plambda_prior)
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
  Par_string[cosmo_param.size()+9] = "Plambda_a";
  param_prior[cosmo_param.size()+9] = Plambda_prior[0];
  Par_string[cosmo_param.size()+10] = "Plambda_b";
  param_prior[cosmo_param.size()+10] = Plambda_prior[1];
  Par_string[cosmo_param.size()+11] = "Plambda_c";
  param_prior[cosmo_param.size()+11] = Plambda_prior[2];

  // input data used to construct the model
  auto inputs = make_shared<STR_NC_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  switch (m_HistogramType) {
  case (glob::HistogramType::_N_V_):
    if (m_data_model.scalrel_z_evo == "E_z")
      m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_proxy_Ez, nParams, Par_type, Par_string, inputs));
    else if (m_data_model.scalrel_z_evo == "direct")
      m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_proxy_zDirect, nParams, Par_type, Par_string, inputs));
    else
      cbl::ErrorCBL("No such a possibility for f(z)!","set_model_NumberCounts_cosmology","Modelling_NumberCounts1D_MassProxy.cpp");
    break;
  default:
    ErrorCBL("Only counts can be modelled! Set _N_V_ as the histogram type.", "set_model_NumberCounts_cosmology", "Modelling_NumberCounts1D_MassProxy.cpp");
  }
}
