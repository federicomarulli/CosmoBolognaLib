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
 *  @file Statistics/Chain.cpp
 *
 *  @brief Methods of the  class Chain 
 *
 *  This file contains the implementation of the methods of the class
 *  Chain, output of the montecarlo
 *  process
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Chain.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::statistics::Chain::Statistics (const int max, const int min)
{
  double Max = (max<=0) ? m_values.size() : max;
  double Min = (min<=0) ? 0 : min;
  vector<double> values;
  
  for (int i=Min; i<Max; i++)
    values.push_back(m_values[i]);

  m_mean = Average(values);
  m_std = Sigma(values);

  m_median = Quartile(values)[1];
}


// ============================================================================================


void cosmobl::statistics::Chain::ComputeDistribution (const int nbin)
{
  vector<double> ww(m_values.size(),1), err;
  distribution(m_var, m_dist, err, m_values, ww, nbin);
}
