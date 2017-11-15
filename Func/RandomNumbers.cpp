/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Func/RandomNumbers.cpp
 *
 *  @brief Methods of the class RandomNumbers
 *
 *  This file contains the implementation of the methods of the class
 *  RandomNumbers, ad derived classes
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Func.h"

using namespace cosmobl;
using namespace random;


// =====================================================================================


cosmobl::random::RandomNumbers::RandomNumbers (const int seed, const double MinVal, const double MaxVal)
{
  set_seed(seed);
  set_range(MinVal, MaxVal);
}


// =====================================================================================


void cosmobl::random::RandomNumbers::set_seed (const int seed)
{
  m_seed = seed;
  m_generator.seed(m_seed);
}


// =====================================================================================


void cosmobl::random::RandomNumbers::set_range (const double MinVal, const double MaxVal)
{
  m_MinVal = MinVal;
  m_MaxVal = MaxVal;
}


// =====================================================================================


cosmobl::random::UniformRandomNumbers::UniformRandomNumbers (double MinVal, const double MaxVal, const int seed) : RandomNumbers(seed, MinVal, MaxVal)
{
  m_distribution = make_shared<uniform_real_distribution<double>>(uniform_real_distribution<double>(0, 1));
}


// =====================================================================================


double cosmobl::random::UniformRandomNumbers::operator () ()
{
  return (m_MaxVal-m_MinVal)*m_distribution->operator()(m_generator)+m_MinVal;
}


// =====================================================================================


cosmobl::random::PoissonRandomNumbers::PoissonRandomNumbers (const double mean, const int seed, const double MinVal, const double MaxVal) : RandomNumbers(seed, MinVal, MaxVal)
{
  set_mean(mean);
}


// =====================================================================================


void cosmobl::random::PoissonRandomNumbers::set_mean (const double mean)
{
  m_mean = mean;
  m_distribution = make_shared<poisson_distribution<int> >(poisson_distribution<int>(mean));
}


// =====================================================================================


double cosmobl::random::PoissonRandomNumbers::operator () ()
{
  double val = m_distribution->operator()(m_generator);

  while (val>=m_MaxVal || val<=m_MinVal)
    val = m_distribution->operator()(m_generator);

  return val;
}


// =====================================================================================


cosmobl::random::NormalRandomNumbers::NormalRandomNumbers (const double mean, const double sigma, const int seed, const double MinVal, const double MaxVal) : RandomNumbers(seed, MinVal, MaxVal)
{
  set_mean_sigma(mean, sigma);
}


// =====================================================================================


void cosmobl::random::NormalRandomNumbers::set_mean_sigma (const double mean, const double sigma)
{
  m_mean = mean;
  m_sigma = sigma;
  m_distribution = make_shared<normal_distribution<double> >(normal_distribution<double>(mean, sigma));
}


// =====================================================================================


double cosmobl::random::NormalRandomNumbers::operator () ()
{
  double val = m_distribution->operator()(m_generator);

  while (val>=m_MaxVal || val<=m_MinVal)
    val = m_distribution->operator()(m_generator);

  return val;
}


// =====================================================================================


cosmobl::random::DiscreteRandomNumbers::DiscreteRandomNumbers (const vector<double> values, const vector<double> weights, const int seed, const double MinVal, const double MaxVal) : RandomNumbers(seed, MinVal, MaxVal)
{
  set_discrete_values(values, weights);
}


// =====================================================================================


void cosmobl::random::DiscreteRandomNumbers::set_discrete_values (const vector<double> values, const vector<double> weights)
{
  if (weights.size()==0) {
    m_values = values;
    m_weights.erase(m_weights.begin(), m_weights.end());
    m_weights.resize(m_values.size(), 1.);
  }
  else if (weights.size()!=values.size())
    ErrorCBL("Error in set_parameters of DiscreteRandomNumbers.h: value and weight vectors have different sizes!");
  else {
    m_values = values;
    m_weights = weights;
  }

  for (size_t i=0; i<m_values.size(); i++)
    if (m_values[i]>=m_MaxVal || m_values[i]<=m_MinVal)
      m_weights[i]=0;

  m_distribution = make_shared<discrete_distribution<int> >(discrete_distribution<int>(m_weights.begin(), m_weights.end()));
}


// =====================================================================================


double cosmobl::random::DiscreteRandomNumbers::operator () ()
{
  return m_values[m_distribution->operator()(m_generator)];
}


// =====================================================================================


cosmobl::random::DistributionRandomNumbers::DistributionRandomNumbers (const vector<double> xx, const vector<double> distribution_function, const string interpolation_method, const int seed) : RandomNumbers()
{
  set_interpolated_distribution(xx, distribution_function, interpolation_method);
  m_uniform_generator = make_shared<UniformRandomNumbers>(0., 1., seed);
}


// =====================================================================================


void cosmobl::random::DistributionRandomNumbers::set_seed (const int seed)
{
  m_uniform_generator->set_seed(seed);
}


// =====================================================================================


void cosmobl::random::DistributionRandomNumbers::set_interpolated_distribution (const vector<double> xx, const vector<double> distribution_function, const string interpolation_method)
{
  glob::FuncGrid ff(xx,distribution_function,interpolation_method);
  double norm = ff.integrate_qag(Min(xx), Max(xx), 1.e-4);

  vector<double> FX;
  for (size_t i=0; i<xx.size(); i++)
    FX.push_back(ff.integrate_qag(Min(xx), xx[i], 1.e-4)/norm);

  m_distribution = make_shared<glob::FuncGrid>(FX, xx, interpolation_method);
}


// =====================================================================================


double cosmobl::random::DistributionRandomNumbers::operator () ()
{
  return m_distribution->operator()(m_uniform_generator->operator()());
}


// =====================================================================================


cosmobl::random::CustomDistributionRandomNumbers::CustomDistributionRandomNumbers (const distribution_func func, const shared_ptr<void> fixed_pars, const vector<double> pars, const int seed, const double MinVal, const double MaxVal) : RandomNumbers(seed, MinVal, MaxVal)
{
  set_custom_distribution(func, fixed_pars, pars);
  m_uniform_generator = make_shared<UniformRandomNumbers>(0., 1., seed);
}


// =====================================================================================


void cosmobl::random::CustomDistributionRandomNumbers::set_seed (const int seed)
{
  m_uniform_generator->set_seed(seed);
}


// =====================================================================================


void cosmobl::random::CustomDistributionRandomNumbers::set_custom_distribution (const distribution_func func, const shared_ptr<void> fixed_pars, const vector<double> pars)
{
  m_func = func;
  m_func_fixed_pars = fixed_pars;
  m_func_pars = pars;

  m_normalization = gsl::GSL_integrate_qag(m_func, m_func_fixed_pars, m_func_pars, m_MinVal, m_MaxVal);
}


// =====================================================================================


double cosmobl::random::CustomDistributionRandomNumbers::operator () ()
{
  auto f = [this] (double xx) {return gsl::GSL_integrate_qag(m_func, m_func_fixed_pars, m_func_pars, m_MinVal, xx)/m_normalization;};
  return gsl::GSL_root_brent(f, m_uniform_generator->operator()(), m_MinVal, m_MaxVal);
}


