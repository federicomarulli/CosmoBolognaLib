/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Modelling/Global/Modelling_Distribution.cpp
 *
 *  @brief Methods of the class Modelling_Distribution, used for modelling different
 *  kinds of distributions 
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_Distribution
 *
 *  @author Giorgio Lesci
 *
 *  @author giorgio.lesci2@unibo.it
 */


#include "Modelling_Distribution.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::modelling::distribution::Modelling_Distribution::set_model_Distribution (const statistics::PriorDistribution mean_prior, const statistics::PriorDistribution std_prior, const std::string mean_name, const std::string std_name)
{
  const int nParams = 2;
  vector<statistics::ParameterType> Par_type(nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string(nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the parameter names and priors
  Par_string[0] = mean_name;
  Par_string[1] = std_name;
  param_prior[0] = mean_prior;
  param_prior[1] = std_prior;
  
  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_gaussian, nParams, Par_type, Par_string, NULL));
}


// ============================================================================================


std::vector<double> cbl::modelling::distribution::model_gaussian (const std::vector<double> x, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)inputs;
  
  std::shared_ptr<void> pp;

  std::vector<double> values(x.size());

  for (size_t i=0; i<x.size(); i++)
    values[i] = cbl::gaussian(x[i], pp, {parameter[0], parameter[1]});

  return values;
}


