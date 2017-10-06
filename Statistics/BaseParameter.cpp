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
 *  @file Statistics/BaseParameter.cpp
 *
 *  @brief Methods of the  class BaseParameter 
 *
 *  This file contains the implementation of the methods of the class
 *  BaseParameter, used to manage input model parameters in
 *  statistical analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "BaseParameter.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::statistics::BaseParameter::m_set_prior ()
{
  if (fixed())
    m_prior = move(unique_ptr<statistics::Prior>(new statistics::Prior(glob::DistributionType::_ConstantDistribution_, m_value)));
  else
    m_prior = m_prior_free;
}


// ============================================================================================


cosmobl::statistics::BaseParameter::BaseParameter (const double value, const string name)
  : Parameter(cosmobl::statistics::ParameterType::_BaseParameter_, name)
{
  set_value(value);
  set_fixed(true);
}


// ============================================================================================


cosmobl::statistics::BaseParameter::BaseParameter (const cosmobl::statistics::Prior prior, const string name) : Parameter(cosmobl::statistics::ParameterType::_BaseParameter_, name)
{
  set_prior(prior);
  set_fixed(false);

  if (m_prior->distributionType() == glob::DistributionType::_ConstantDistribution_)
    fix (m_prior->sample());
}


// ============================================================================================


cosmobl::statistics::BaseParameter::BaseParameter (const double xmin, const double xmax, const int seed, const string name)
  :  Parameter(cosmobl::statistics::ParameterType::_BaseParameter_, name) 
{  
  set_prior(xmin, xmax, seed);
  set_fixed(false);
}


// ============================================================================================


cosmobl::statistics::BaseParameter::BaseParameter (const glob::DistributionType priorType, const vector<double> prior_params, const double xmin, const double xmax, const int seed, const string name)
  :  Parameter(cosmobl::statistics::ParameterType::_BaseParameter_, name)
{  
  set_prior(priorType, prior_params, xmin, xmax, seed);
  set_fixed(false);
}


// ============================================================================================


double cosmobl::statistics::BaseParameter::value () const
{
  return m_value;
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_value (const double value) 
{
  m_value = value;
}


// ============================================================================================


bool cosmobl::statistics::BaseParameter::fixed () const
{
  return m_fixed;
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_fixed (const bool fix)
{
  m_fixed = fix;
  m_set_prior();
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::free ()
{
  set_fixed(false);
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::fix ()
{
  set_fixed(true);
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::fix (const double value)
{
  set_value(value);
  set_fixed(true);
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::fix_at_bestfit ()
{
  set_value(m_bestfit_value);
  set_fixed(true);
}


// ============================================================================================


shared_ptr<statistics::Prior> cosmobl::statistics::BaseParameter::prior () const
{
  return m_prior;
}


// ============================================================================================


double cosmobl::statistics::BaseParameter::PriorProbability (const double value) const
{
  return m_prior->operator()(value);
}


// ============================================================================================


double cosmobl::statistics::BaseParameter::LogPriorProbability (const double value) const
{
  return m_prior->log_distribution(value);
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_prior (statistics::Prior prior)
{
  auto ptr_prior = make_shared<cosmobl::statistics::Prior>(prior);
  m_prior_free = move(ptr_prior);
  m_set_prior();
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_prior (const double xmin, const double xmax, const int seed) 
{
  cosmobl::statistics::Prior prior(cosmobl::glob::DistributionType::_UniformDistribution_, xmin, xmax, seed);
  set_prior(prior);
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_prior (const glob::DistributionType priorType, const vector<double> prior_params, const double xmin, const double xmax, const int seed)
{
  cosmobl::statistics::Prior prior(priorType, prior_params, xmin, xmax, seed);
  set_prior(prior);
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_prior_seed (const int seed)
{
  m_prior->set_seed(seed);
}


// ============================================================================================


double cosmobl::statistics::BaseParameter::prior_sample ()
{
  return m_prior->sample();
}


// ============================================================================================


vector<double> cosmobl::statistics::BaseParameter::prior_sample (const int sample_size)
{
  return m_prior->sample_vector(sample_size);
}


// ============================================================================================


double cosmobl::statistics::BaseParameter::prior_range (const double epsilon)
{
  return epsilon*(m_prior->distribution_range());
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_posterior (const int start, const int thin, const int seed)
{
  if (fixed()){
    (void)start; (void)thin; (void)seed;
    m_posterior = move(unique_ptr<statistics::Posterior>(new statistics::Posterior(glob::DistributionType::_ConstantDistribution_, m_prior->sample())));
  }
  else{
    m_posterior = m_chain->PosteriorDistribution(start, thin, seed);
  }
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_posterior (const vector<double> chain_values, const int nwalkers, const int seed)
{
  if (fixed()){
    (void)chain_values; (void)nwalkers; (void)seed;
    m_posterior = move(unique_ptr<statistics::Posterior>(new statistics::Posterior(glob::DistributionType::_ConstantDistribution_, m_prior->sample())));
  }
  else{
    set_chain(chain_values, nwalkers);
    m_posterior = m_chain->PosteriorDistribution(seed);
  }
}


// ============================================================================================


void cosmobl::statistics::BaseParameter::set_posterior (const vector<vector<double>> chain_values, const int seed)
{
  if (fixed()){
    (void)chain_values; (void)seed;
    m_posterior = move(unique_ptr<statistics::Posterior>(new statistics::Posterior(glob::DistributionType::_ConstantDistribution_, m_prior->sample())));
  }
  else{
    set_chain(chain_values);
    m_posterior = m_chain->PosteriorDistribution(seed);
  }
}
