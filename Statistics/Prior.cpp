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

/** @file Statistics/Prior.cpp
 *
 *  @brief Methods of the class Prior 
 *
 *  This file contains the implementation of the methods of the class
 *  Prior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Prior.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::statistics::Prior::Prior (const Prior_function prior_function, const std::shared_ptr<void> prior_function_inputs)
{
  m_prior_function = prior_function;
  m_prior_function_inputs = prior_function_inputs;
}


// ============================================================================================


double cbl::statistics::Prior::operator() (const vector<double> parameters)
{
  return m_prior_function(parameters, m_prior_function_inputs);
}


// ============================================================================================


double cbl::statistics::Prior::log (const vector<double> parameters)
{
  double prior = this->operator()(parameters);
  return (prior>0) ? std::log(prior) : par::defaultDouble;
}
