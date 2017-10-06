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
 *  @file Statistics/DerivedParameter.cpp
 *
 *  @brief Methods of the  class DerivedParameter 
 *
 *  This file contains the implementation of the methods of the class
 *  DerivedParameter, used to manage derived model parameters in
 *  statistical analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "DerivedParameter.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::statistics::DerivedParameter::set_posterior (const int start, const int thin, const int seed)
{
  m_posterior = m_chain->PosteriorDistribution(start, thin, seed);
}


// ============================================================================================


void cosmobl::statistics::DerivedParameter::set_posterior (const vector<double> chain_values, const int nwalkers, const int seed)
{
  set_chain(chain_values, nwalkers);
  m_posterior = m_chain->PosteriorDistribution(seed);
}


// ============================================================================================


void cosmobl::statistics::DerivedParameter::set_posterior (const vector<vector<double>> chain_values, const int seed)
{
  set_chain(chain_values);
  m_posterior = m_chain->PosteriorDistribution(seed);
}
