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


// ============================================================================================


cosmobl::statistics::Chain::Chain (const vector<double> values, const int nwalkers)
{
  set_values(values, nwalkers);
}


// ============================================================================================


cosmobl::statistics::Chain::Chain (const int size, const int nwalkers)
{
  set_chain(size, nwalkers);
}


// ============================================================================================


cosmobl::statistics::Chain::Chain (vector<vector<double>> values)
{
  set_values(values);
}


// ============================================================================================


void cosmobl::statistics::Chain::set_chain (const int size, const int nwalkers)
{
  m_size = size;
  m_nwalkers = nwalkers;
  reset();
}


// ============================================================================================


void cosmobl::statistics::Chain::reset()
{
  m_values.erase(m_values.begin(), m_values.end());
  m_values.resize(m_size*m_nwalkers, 0);
}


// ============================================================================================


void cosmobl::statistics::Chain::expand(const int append)
{
  vector<double> values = m_values;
  int old_size = m_size;

  m_size += append;

  reset();

  for(int i=0; i<old_size; i++)
    for(int j=0; j<m_nwalkers; j++)
      set_value(i, j, values[i*m_nwalkers+j]);
}


// ============================================================================================


vector<double> cosmobl::statistics::Chain::values(const int start, const int thin) const
{
  vector<double> values;

  for(int i=start; i<size(); i+=thin)
    for(int j=0; j<m_nwalkers; j++)
      values.push_back(value(i, j));

  return values;
}


// ============================================================================================


void cosmobl::statistics::Chain::set_values(const vector<double> values, const int nwalkers)
{
  int size = values.size()/nwalkers;
  
  if(values.size()%nwalkers!=0)
    ErrorCBL("Error in set_values! wrong size of input values of wrong number of walkers!");

  set_chain(size, nwalkers);

  for(int i=0; i< m_size; i++)
    for(int j=0; j< m_nwalkers; j++)
      set_value(i, j, values[i*m_nwalkers+j]);
}


// ============================================================================================


void cosmobl::statistics::Chain::set_values(const vector<vector<double>> values)
{
  vector<double> flatten_val;

  int nwalkers = values[0].size();

  for(size_t i=0; i<values.size(); i++){
    checkDim(values[i], nwalkers, "values["+conv(i, par::fINT)+"]");
    for(size_t j=0; j< values[i].size(); j++)
      flatten_val.push_back(values[i][j]);
  }

  set_values(flatten_val, nwalkers);
}


// ============================================================================================
 

shared_ptr<cosmobl::statistics::Posterior> cosmobl::statistics::Chain::PosteriorDistribution(const int seed) const
{
  const int nbins=50;
  return make_shared<cosmobl::statistics::Posterior> (cosmobl::statistics::Posterior(glob::DistributionType::_DiscreteDistribution_, m_values, {}, nbins, "Spline", seed));
}


// ============================================================================================
 

shared_ptr<cosmobl::statistics::Posterior> cosmobl::statistics::Chain::PosteriorDistribution(const int start, const int thin, const int seed) const
{
  const int nbins=50;
  vector<double> vv = values(start, thin);
  return make_shared<cosmobl::statistics::Posterior> (cosmobl::statistics::Posterior(glob::DistributionType::_DiscreteDistribution_, vv, {}, nbins, "Spline", seed));
}
