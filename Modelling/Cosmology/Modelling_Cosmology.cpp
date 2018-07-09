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
 *  Modelling/Cosmology/Modelling_Cosmology.cpp
 *
 *  @brief Methods of the class Modelling_Cosmology
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_Cosmology, i.e. the common functions to
 *  model cosmological measurements
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_Cosmology.h"
#include "Modelling_Cosmology_DistancePrior.h"

using namespace std;

using namespace cbl;


// ============================================================================================


shared_ptr<cbl::modelling::cosmology::CMB_DistancePrior> cbl::modelling::cosmology::CMB_DistancePrior::Create (const string distance_prior_name)
{
  if (distance_prior_name == "Aubourg15_Planck15")
    return move(unique_ptr<cbl::modelling::cosmology::Aubourg15_Planck15>(new cbl::modelling::cosmology::Aubourg15_Planck15()));
  else if (distance_prior_name == "Aubourg15_WMAP09")
    return move(unique_ptr<cbl::modelling::cosmology::Aubourg15_WMAP09>(new cbl::modelling::cosmology::Aubourg15_WMAP09()));
  else ErrorCBL("Error in cbl::modelling::cosmology::CMB_DistancePrior::Create of Modelling_Cosmology_DistancePrior.h: no such type of CMB_DistancePrior!");

  return NULL;
}


// ============================================================================================


cbl::modelling::cosmology::Modelling_Cosmology::Modelling_Cosmology(const shared_ptr<cbl::data::Data> dataset, const vector<string> data_type)
{
  m_data = dataset;
  m_data_type = data_type;
}


// ============================================================================================


void cbl::modelling::cosmology::Modelling_Cosmology::set_fiducial_cosmology(const cbl::cosmology::Cosmology cosmology)
{
  m_cosmology = move(make_shared<cbl::cosmology::Cosmology>(cosmology));
}


// ============================================================================================


void cbl::modelling::cosmology::Modelling_Cosmology::set_cosmological_parameters(const vector<cbl::cosmology::CosmologicalParameter> cosmoPar_name, const vector<cbl::statistics::PriorDistribution> cosmoPar_prior, const string distance_prior, const vector<string> external_dataset)
{
  (void)distance_prior;
  (void)external_dataset;

  const size_t nParams = cosmoPar_name.size();
  checkDim(cosmoPar_prior, nParams, "cosmoPar_prior");

  vector<statistics::ParameterType> cosmoPar_type(nParams, statistics::ParameterType::_Base_);
  vector<string> cosmoPar_string(nParams);

  for (size_t i=0; i<nParams; i++)
    cosmoPar_string[i] = CosmologicalParameter_name(cosmoPar_name[i]);

  m_data_model.cosmology = m_cosmology;
  m_data_model.Cpar = cosmoPar_name;
  m_data_model.data_type = m_data_type;

  if(distance_prior != par::defaultString){
    m_data_model.distance_prior = cbl::modelling::cosmology::CMB_DistancePrior::Create(distance_prior);
    m_data_fit = cbl::data::join_dataset({m_data, m_data_model.distance_prior->dataset()});
    m_fit_range = true;
  }

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_cosmology>(m_data_model);

  // set the priors
  m_set_prior(cosmoPar_prior);

  // construct the model
  if(distance_prior != par::defaultString)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&cbl::modelling::cosmology::cosmological_measurements_model_CMB_DistancePrior, nParams, cosmoPar_type, cosmoPar_string, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&cbl::modelling::cosmology::cosmological_measurements_model, nParams, cosmoPar_type, cosmoPar_string, inputs));
}
