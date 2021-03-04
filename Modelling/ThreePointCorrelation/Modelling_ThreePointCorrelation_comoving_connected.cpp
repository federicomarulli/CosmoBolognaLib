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
 *  Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation_comoving_connected.cpp
 *
 *  @brief Methods of the class
 *  Modelling_ThreePointCorrelation_comoving_connected
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_ThreePointCorrelation_comoving_connected
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unibo.it
 */


#include "Modelling_ThreePointCorrelation_comoving_connected.h"

using namespace std;

using namespace cbl;


// ============================================================================================
	

void cbl::modelling::threept::Modelling_ThreePointCorrelation_comoving_connected::set_model_RSD (const statistics::PriorDistribution b1_prior, const statistics::PriorDistribution b2_prior, const statistics::PriorDistribution bt_prior, const statistics::PriorDistribution beta_prior)
{
  // set the model parameters
  const int nparameters = 4;

  vector<statistics::ParameterType> parameterType = {statistics::ParameterType::_Base_};

  vector<string> parameterName(nparameters);
  parameterName[0] = "b1";
  parameterName[1] = "b2";
  parameterName[2] = "bt";
  parameterName[3] = "beta";

  vector<statistics::PriorDistribution> priors = {b1_prior, b2_prior, bt_prior, beta_prior};

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model_threept>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&zeta_RSD, nparameters, parameterType, parameterName, inputs));
}
