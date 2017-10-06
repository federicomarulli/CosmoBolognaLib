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
 *  @file Statistics/Parameter.cpp
 *
 *  @brief Methods of the  class Parameter 
 *
 *  This file contains the implementation of the methods of the class
 *  Parameter, used to manage model parameters in chi2 fitting or monte carlo
 *  analysis
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Parameter.h"

using namespace cosmobl;


// ============================================================================================


string cosmobl::statistics::Parameter::name () const
{
  return m_name;
}


// ============================================================================================


cosmobl::statistics::ParameterType cosmobl::statistics::Parameter::parameterType () const
{
  return m_pType;
}


// ============================================================================================


double cosmobl::statistics::Parameter::bestfit_value () const
{
  if(m_bestfit_set)
    return m_bestfit_value;
  else
    ErrorCBL("Error in bestfit_value of Parameter.cpp. Bestfit value is not set!");

  return 0;
}


// ============================================================================================


shared_ptr<statistics::Posterior> cosmobl::statistics::Parameter::posterior () const
{
  return m_posterior;
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_name (const string name)
{
  m_name = name;
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_pType (const cosmobl::statistics::ParameterType pType)
{
  m_pType = pType;
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_bestfit_value (const double bestfit_value)
{
  m_bestfit_value = bestfit_value;
  m_bestfit_set = true;
}


// ============================================================================================

  
shared_ptr<cosmobl::statistics::Chain> cosmobl::statistics::Parameter::chain () const
{
  return m_chain;
}


// ============================================================================================


int cosmobl::statistics::Parameter::chain_size () const
{
  return m_chain->size();
}


// ============================================================================================


int cosmobl::statistics::Parameter::nwalkers () const
{
  return m_chain->nwalkers();
}


// ============================================================================================


double cosmobl::statistics::Parameter::chain_value (const int pp, const int ww) const
{
  return m_chain->value(pp, ww);
}


// ============================================================================================


vector<double> cosmobl::statistics::Parameter::chain_values (const int start, const int thin) const
{
  return m_chain->values(start, thin);
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chain (shared_ptr<Chain> chain)
{
  m_chain = chain;
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chain (const int chain_size, const int nwalkers)
{
  m_chain = make_shared<cosmobl::statistics::Chain>(cosmobl::statistics::Chain(chain_size, nwalkers));
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chain_value (const int pp, const int ww, const double chain_value)
{
  m_chain->set_value(pp, ww, chain_value);
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chain (const vector<double> chain_values, const int nwalkers)
{
  m_chain = make_shared<cosmobl::statistics::Chain>(cosmobl::statistics::Chain(chain_values, nwalkers));
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chain (const vector<vector<double>> chain_values)
{
  m_chain = make_shared<cosmobl::statistics::Chain>(cosmobl::statistics::Chain(chain_values));
}


// ============================================================================================


void cosmobl::statistics::Parameter::expand_chain (const int append)
{
  m_chain->expand(append);
}


// ============================================================================================


double cosmobl::statistics::Parameter::PosteriorProbability (const double value) const
{
  return m_posterior->operator()(value);
}


// ============================================================================================


double cosmobl::statistics::Parameter::posterior_mean () const
{
  return m_posterior->mean();
}


// ============================================================================================


double cosmobl::statistics::Parameter::posterior_median () const
{
  return m_posterior->median();
}


// ============================================================================================


double cosmobl::statistics::Parameter::posterior_std () const
{
  return m_posterior->std();
}


// ============================================================================================


double cosmobl::statistics::Parameter::posterior_percentile (const unsigned int i) const
{
  return m_posterior->percentile(i);
}


// ============================================================================================


double cosmobl::statistics::Parameter::posterior_sample ()
{
  return m_posterior->sample();
}


// ============================================================================================


vector<double> cosmobl::statistics::Parameter::posterior_sample (const int sample_size)
{
  return m_posterior->sample_vector(sample_size);
}


// ============================================================================================


double cosmobl::statistics::Parameter::posterior_mode () const
{
  return m_posterior->mode();
}
