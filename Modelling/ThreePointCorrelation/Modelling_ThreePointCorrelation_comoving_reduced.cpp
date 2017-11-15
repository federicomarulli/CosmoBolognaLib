/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation_comoving_reduced.cpp
 *
 *  @brief Methods of the class
 *  Modelling_ThreePointCorrelation_comoving_reduced
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_ThreePointCorrelation_comoving_reduced
 *
 *  @authors Federico Marulli, Michele Moresco
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it
 */


#include "Modelling_ThreePointCorrelation_comoving_reduced.h"

using namespace cosmobl;


// ============================================================================================
	

void cosmobl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced::set_model_nonlinear_localbias (const statistics::Prior bias1_prior, const statistics::Prior bias2_prior)
{
  auto bias1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias1_prior, "b1"));
  auto bias2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias2_prior, "b2"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {bias1, bias2};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Q_nonlinear_localbias, inputs));
}


// ============================================================================================
	

void cosmobl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced::set_model_nonlinear_nonlocalbias (const statistics::Prior bias1_prior, const statistics::Prior bias2_prior, const statistics::Prior g2_prior)
{
  auto bias1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias1_prior, "b1"));
  auto bias2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias2_prior, "b2"));
  auto g2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(g2_prior, "g2"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {bias1, bias2, g2};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Q_nonlinear_nonlocalbias, inputs));
}


// ============================================================================================
  

void cosmobl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced::set_model_nonlinear_nonlocalbias_alpha (const statistics::Prior bias1_prior, const statistics::Prior bias2_prior, const statistics::Prior g2_prior, const statistics::Prior alpha_prior)
{
  auto bias1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias1_prior, "b1"));
  auto bias2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias2_prior, "b2"));
  auto g2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(g2_prior, "g2"));
  auto alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {bias1, bias2, g2, alpha};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Q_nonlinear_nonlocalbias_alpha, inputs));
}
