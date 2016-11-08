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
 *  @file Statistics/Prior.cpp
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

using namespace cosmobl;
using namespace random;

// ======================================================================================


void cosmobl::statistics::Prior::m_set_prior_normalization()
{
  m_prior_normalization = 1;

  function<double(double)> f = bind(&Prior::operator(), this, std::placeholders::_1);
	   
  glob::STR_generic_integrand pp;
  pp.f = f;

  gsl_function Func;
  Func.function = &generic_integrand;
  Func.params = &pp;

  m_prior_normalization = GSL_integrate_qag(Func,m_xmin, m_xmax);
}

// ======================================================================================


cosmobl::statistics::Prior::Prior (const cosmobl::statistics::PriorType priorType, const double xmin, const double xmax, const int seed) 
{
  if (priorType != statistics::PriorType::_UniformPrior_)
    ErrorCBL("Error in constructor of Prior, this constructor only allows PriorType::_UniformPrior_");

  set_uniform_prior(xmin, xmax, seed);
}


// ======================================================================================


cosmobl::statistics::Prior::Prior (const cosmobl::statistics::PriorType priorType, const vector<double> prior_params, const double xmin, const double xmax , const int seed) 
{
  set_limits(xmin, xmax);

  if (priorType == statistics::PriorType::_GaussianPrior_) {
    if (prior_params.size() != 2)
      ErrorCBL("Error in constructor of Prior, wrong size of prior_params. Gaussian prior needs 2 parameters, the mean and the standard deviation");

    set_gaussian_prior(prior_params[0], prior_params[1], seed);
  }
  
  else if (priorType == statistics::PriorType::_PoissonPrior_) {
    if (prior_params.size() != 1)
      ErrorCBL("Error in constructor of Prior, wrong size of prior_params. Poisson prior needs 1 parameter, the mean");

    set_poisson_prior(prior_params[0], seed);
  }
  
  else 
    ErrorCBL("Error in Prior constructor of Prior.cpp. No such type of prior");
}


// ======================================================================================


cosmobl::statistics::Prior::Prior (const cosmobl::statistics::PriorType priorType, const vector<double> discrete_values, const vector<double> weights, const int seed) 
{
  if (priorType != statistics::PriorType::_DiscretePrior_)
    ErrorCBL("Error in constructor of Prior, this constructor only allows PriorType::_DiscretePrior_");

  set_discrete_values(discrete_values, weights, seed);
}


// ======================================================================================


void cosmobl::statistics::Prior::set_limits(const double xmin, const double xmax)
{
  m_xmin = xmin; 
  m_xmax = xmax;
}


// ======================================================================================


void cosmobl::statistics::Prior::set_uniform_prior(const double xmin, const double xmax, const int seed)
{
  set_limits(xmin, xmax);
  m_prior_func_pars.erase(m_prior_func_pars.begin(), m_prior_func_pars.end());

  m_prior_func_pars.push_back(m_xmax);
  m_prior_func_pars.push_back(m_xmax);

  m_prior_random = make_shared<UniformRandomNumbers> (UniformRandomNumbers(m_xmin, m_xmax, seed));
  m_func = &identity<double>; 

  m_prior_normalization = m_xmax-m_xmin;
}


// ======================================================================================


void cosmobl::statistics::Prior::set_gaussian_prior(const double mean, const double sigma, const int seed)
{
  m_prior_func_pars.erase(m_prior_func_pars.begin(), m_prior_func_pars.end());
  m_prior_func_pars.push_back(mean);
  m_prior_func_pars.push_back(sigma);

  m_prior_random = make_shared<NormalRandomNumbers> (NormalRandomNumbers(mean, sigma, seed, m_xmin, m_xmax));
  m_func = &gaussian<double>; 

  m_prior_normalization = 0.5*(erf((m_xmax-mean)/sigma)-erf((m_xmin-mean)/sigma));
}


// ======================================================================================


void cosmobl::statistics::Prior::set_poisson_prior (const double mean, const int seed) 
{
  m_xmin = nint(m_xmin);
  m_xmax = nint(m_xmax);
  int nbins = m_xmax-m_xmin;

  vector<double> poisson_values = linear_bin_vector(nbins, m_xmin, m_xmax);
  vector<double> weights;
  for(int i=0;i<nbins;i++)
    weights.push_back(poisson(poisson_values[i],NULL,{mean}));

  m_prior_random = make_shared<DiscreteRandomNumbers> (DiscreteRandomNumbers(poisson_values, weights, seed, m_xmin, m_xmax));

  glob::STR_closest_probability parameters;
  parameters.values = poisson_values;
  parameters.weights = weights; 

  m_prior_func_fixed_pars = make_shared<glob::STR_closest_probability>(parameters);
  m_func = &closest_probability; 

  m_prior_normalization = accumulate(weights.begin(), weights.end(), 0);
}


// =====================================================================================


void cosmobl::statistics::Prior::set_discrete_values (const vector<double> discrete_values, const vector<double> weights, const int seed) 
{
  if (discrete_values.size()==0)
    ErrorCBL("Error in set_discrete_values of Prior. Vector of values is empty");

  set_limits(Min(discrete_values), Max(discrete_values));
  m_prior_random = make_shared<DiscreteRandomNumbers> (DiscreteRandomNumbers(discrete_values, weights, seed, m_xmin, m_xmax));

  glob::STR_closest_probability parameters;
  parameters.values = discrete_values;
  parameters.weights = weights; 

  m_prior_func_fixed_pars = make_shared<glob::STR_closest_probability>(parameters);

  m_func = &closest_probability; 
  m_prior_normalization = accumulate(weights.begin(), weights.end(), 0);

}


// =====================================================================================


bool cosmobl::statistics::Prior::isIncluded (const double value) const
{
  if (value > m_xmin && m_xmax > value)
    return true;
  else 
    return false;
}


// =====================================================================================


double cosmobl::statistics::Prior::sample () const
{
  return m_prior_random->operator()();
}


// =====================================================================================


double cosmobl::statistics::Prior::sample (const int seed)
{
  m_prior_random->set_seed(seed);
  return sample();
}


// =====================================================================================


vector<double> cosmobl::statistics::Prior::sample_vector (const int nvalues)
{
  vector<double> values;
  
  for (int i=0; i<nvalues; i++)
    values.push_back(sample());

  return values;

}
