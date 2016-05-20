/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Statistics/Model.cpp
 *
 *  @brief Methods of the class Model
 *
 *  This file contains the implementation of the methods of the class
 *  Model
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Model.h"
using namespace cosmobl;


// ======================================================================================


cosmobl::statistics::Model::Model (const vector<Parameter> parameters, const shared_ptr<void> model_parameters)
  : m_model_parameters(model_parameters)
{
  for(size_t i=0;i<parameters.size();i++)
    m_parameters.push_back(move(make_shared<Parameter>(parameters[i])));

  m_npar = m_parameters.size();
  m_npar_eff = npar_eff();
}


// ======================================================================================


double cosmobl::statistics::Model::update_parameter (const double new_parameter)
{
  double parameter;

  if (!m_parameters[0]->isFreezed()) {
    parameter = (m_parameters[0]->prior()->isIncluded(new_parameter)) ? new_parameter : closest(new_parameter,m_parameters[0]->prior()->xmin(),m_parameters[0]->prior()->xmax());
    m_parameters[0]->set_value(parameter);
  }
  
  else
    parameter = m_parameters[0]->value();
  
  m_parameters[0]->set_value(parameter);

  return parameter;
}


// ======================================================================================


vector<double> cosmobl::statistics::Model::update_parameters (const vector<double> new_parameters)
{  
  vector<double> parameters;
  int nn = 0;
  
  for (unsigned int i=0; i<m_npar; i++) {
    
    if (!m_parameters[i]->isFreezed()) {
      double value = (m_parameters[i]->prior()->isIncluded(new_parameters[nn])) ? new_parameters[nn] : closest(new_parameters[nn], m_parameters[i]->prior()->xmin(), m_parameters[i]->prior()->xmax());
      m_parameters[i]->set_value(value);
      nn ++;
    }

    parameters.push_back(m_parameters[i]->value());
  }
  
  return parameters;
}

// ======================================================================================


void cosmobl::statistics::Model::set_chains (const int nchains, const int chain_size)
{
  for (unsigned int i=0; i<m_npar; i++) 
    m_parameters[i]->set_chains(nchains, chain_size);
}


// ======================================================================================


void cosmobl::statistics::Model::set_random_number_generator (const default_random_engine generator)
{
  m_generator = generator;
  uniform_real_distribution<double> distribution(0.,1.);
  m_random_numbers = bind(distribution,m_generator);
}
