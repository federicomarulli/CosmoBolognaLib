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
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it
 */


#include "Modelling_ThreePointCorrelation_comoving_reduced.h"

using namespace std;

using namespace cbl;


// ============================================================================================
	

void cbl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced::set_model_nonlinear_localbias (const statistics::PriorDistribution bias1_prior, const statistics::PriorDistribution bias2_prior)
{
  // set the model parameters
  const int nparameters = 2;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "b1";
  parameterName[1] = "b2";

  vector<statistics::PriorDistribution> priors = {bias1_prior, bias2_prior};

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Q_nonlinear_localbias, nparameters, parameterType, parameterName, inputs));
}


// ============================================================================================
	

void cbl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced::set_model_nonlinear_nonlocalbias (const statistics::PriorDistribution bias1_prior, const statistics::PriorDistribution bias2_prior, const statistics::PriorDistribution g2_prior)
{
  // set the model parameters
  const int nparameters = 3;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "b1";
  parameterName[1] = "b2";
  parameterName[2] = "g2";

  vector<statistics::PriorDistribution> priors = {bias1_prior, bias2_prior, g2_prior};

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Q_nonlinear_nonlocalbias, nparameters, parameterType, parameterName, inputs));
}


// ============================================================================================
  

void cbl::modelling::threept::Modelling_ThreePointCorrelation_comoving_reduced::set_model_nonlinear_nonlocalbias_alpha (const statistics::PriorDistribution bias1_prior, const statistics::PriorDistribution bias2_prior, const statistics::PriorDistribution g2_prior, const statistics::PriorDistribution alpha_prior)
{
  // set the model parameters
  const int nparameters = 4;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "b1";
  parameterName[1] = "b2";
  parameterName[2] = "g2";
  parameterName[3] = "alpha";

  vector<statistics::PriorDistribution> priors = {bias1_prior, bias2_prior, g2_prior, alpha_prior};

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Q_nonlinear_nonlocalbias_alpha, nparameters, parameterType, parameterName, inputs));
}
