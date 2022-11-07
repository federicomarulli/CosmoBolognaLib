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
 *  Modelling/NumberCounts/Modelling_NumberCounts1D_Mass.cpp
 *
 *  @brief Methods of the class Modelling_NumberCounts1D_Mass
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_NumberCounts1D_Mass, i.e. the functions to model
 *  the mass number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts1D_Mass.h"
#include "Modelling_NumberCounts1D_Mass.h"

using namespace std;

using namespace cbl;

using namespace modelling::numbercounts;


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_Mass::set_data_model_snapshot (const cosmology::Cosmology cosmology, const double redshift, const std::string method_Pk, const double k_min, const double k_max, const int step, const int norm, const double Delta, const bool isDelta_critical, const std::string model_MF, const double Volume, const double Mass_min, const double Mass_max, const int Mass_step, const double prec)
{
  m_data_model.isSnapshot = true;

  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.redshift = redshift;
  m_data_model.method_Pk = method_Pk;
  m_data_model.k_min = k_min;
  m_data_model.k_max = k_max;
  m_data_model.step = step;
  m_data_model.kk = logarithmic_bin_vector(step, k_min, k_max);
  m_data_model.norm = norm;
  
  m_data_model.output_root = "test";
  m_data_model.file_par = par::defaultString;

  m_data_model.isDelta_critical = isDelta_critical;
  m_data_model.Delta = Delta;
  m_data_model.model_MF = model_MF;

  m_data_model.Volume = Volume;

  m_data_model.Mass_min = Mass_min;
  m_data_model.Mass_max = Mass_max;
  m_data_model.Mass_step = Mass_step;
  m_data_model.Mass_vector = logarithmic_bin_vector(200, 1.e10, 1.e16);

  m_data_model.prec = prec;
}


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts1D_Mass::set_model_NumberCounts_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior)
{
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size();

  vector<statistics::ParameterType> cosmoPar_type(nParams, statistics::ParameterType::_Base_);
  vector<string> cosmoPar_string(nParams);

  for (size_t i=0; i<nParams; i++)
    cosmoPar_string[i] = CosmologicalParameter_name(cosmo_param[i]);

  // input data used to construct the model
  auto inputs = make_shared<STR_NC_data_model>(m_data_model);

  // set prior
  m_set_prior(cosmo_param_prior);

  // construct the model
  switch (m_HistogramType) {

  case (glob::HistogramType::_N_V_):
    if (m_data_model.isSnapshot == true)
      m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_mass_snapshot, nParams, cosmoPar_type, cosmoPar_string, inputs));
    else
      m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_counts_mass, nParams, cosmoPar_type, cosmoPar_string, inputs));
    break;
  case (glob::HistogramType::_n_V_):
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&number_density_mass, nParams, cosmoPar_type, cosmoPar_string, inputs));
    break;
  case (glob::HistogramType::_dn_dV_):
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&mass_function_mass, nParams, cosmoPar_type, cosmoPar_string, inputs));
    break;
  default:
    ErrorCBL("no such a variable in the list!", "set_model_NumberCounts_cosmology", "Modelling_NumberCounts1D_Mass.cpp");
  }
}
