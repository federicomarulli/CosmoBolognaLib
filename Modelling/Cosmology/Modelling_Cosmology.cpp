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

using namespace cosmobl;

// ============================================================================================


shared_ptr<cosmobl::modelling::cosmology::CMB_DistancePrior> cosmobl::modelling::cosmology::CMB_DistancePrior::Create (const string distance_prior_name)
{
  if (distance_prior_name == "Aubourg15_Planck15")
    return move(unique_ptr<cosmobl::modelling::cosmology::Aubourg15_Planck15>(new cosmobl::modelling::cosmology::Aubourg15_Planck15()));
  else if (distance_prior_name == "Aubourg15_WMAP09")
    return move(unique_ptr<cosmobl::modelling::cosmology::Aubourg15_WMAP09>(new cosmobl::modelling::cosmology::Aubourg15_WMAP09()));
  else ErrorCBL("Error in cosmobl::modelling::cosmology::CMB_DistancePrior::Create of Modelling_Cosmology_DistancePrior.h: no such type of CMB_DistancePrior!");

  return NULL;
}


// ============================================================================================


cosmobl::modelling::cosmology::Modelling_Cosmology::Modelling_Cosmology(const shared_ptr<cosmobl::data::Data> dataset, const vector<string> data_type)
{
  m_data = dataset;
  m_data_type = data_type;
}


// ============================================================================================


void cosmobl::modelling::cosmology::Modelling_Cosmology::set_fiducial_cosmology(const cosmobl::cosmology::Cosmology cosmology)
{
  m_cosmology = move(make_shared<cosmobl::cosmology::Cosmology>(cosmology));
}


// ============================================================================================


void cosmobl::modelling::cosmology::Modelling_Cosmology::set_cosmological_parameters(const vector<cosmobl::cosmology::CosmoPar> cosmoPar_name, const vector<cosmobl::statistics::Prior> cosmoPar_prior, const string distance_prior, const vector<string> external_dataset)
{
  (void)distance_prior;
  (void)external_dataset;

  const int nParams = cosmoPar_name.size();
  checkDim(cosmoPar_prior, nParams, "cosmoPar_prior");

  vector<shared_ptr<statistics::Parameter>> ll_parameters;

  for(int i=0; i<nParams; i++){
    m_map_cosmoPar[cosmoPar_name[i]] = move(make_shared<statistics::BaseParameter>(cosmobl::statistics::BaseParameter(cosmoPar_prior[i], cosmobl::cosmology::CosmoPar_name (cosmoPar_name[i])))); 
    ll_parameters.push_back(m_map_cosmoPar[cosmoPar_name[i]]);
  }

  auto z_drag = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("z_drag"));
  ll_parameters.push_back(z_drag);

  auto rs = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("rs"));
  ll_parameters.push_back(rs);

  set_parameters(ll_parameters);

  m_data_model.cosmology = m_cosmology;
  m_data_model.Cpar = cosmoPar_name;
  m_data_model.data_type = m_data_type;

  if(distance_prior != par::defaultString){
    m_data_model.distance_prior = cosmobl::modelling::cosmology::CMB_DistancePrior::Create(distance_prior);
    m_data_fit = cosmobl::data::join_dataset({m_data, m_data_model.distance_prior->dataset()});
    m_fit_range = true;
  }

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_cosmology>(m_data_model);

  // construct the model
  if(distance_prior != par::defaultString)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&cosmobl::modelling::cosmology::cosmological_measurements_model_CMB_DistancePrior, inputs));
  else
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&cosmobl::modelling::cosmology::cosmological_measurements_model, inputs));
}
