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


cosmobl::statistics::Parameter::Parameter (const double value, const ParameterType pType, const string name)
{
  m_value = value;
  m_pType = pType;
  m_name = name;
  m_prior = make_shared<Prior>(statistics::Prior(statistics::PriorType::_UniformPrior_, value*0.9, value*1.1));
}


// ============================================================================================


cosmobl::statistics::Parameter::Parameter (const double value, const double xmin, const double xmax, const ParameterType pType, const string name)
{
  m_value = value;
  m_name = name;
  m_pType = pType;
  m_prior = make_shared<statistics::Prior>(statistics::Prior(statistics::PriorType::_UniformPrior_, xmin, xmax));
}


// ============================================================================================


cosmobl::statistics::Parameter::Parameter (const double value, const statistics::PriorType priorType, const vector<double> prior_params, const double xmin, const double xmax, const ParameterType pType, const string name)
{
  m_value = value;
  m_name = name;
  m_pType = pType;
  m_prior = make_shared<statistics::Prior>(statistics::Prior(priorType, prior_params, xmin, xmax));
}


// ============================================================================================


cosmobl::statistics::Parameter::Parameter (const double value, const statistics::PriorType priorType, const vector<double> discrete_values, const vector<double> weights, const ParameterType pType, const string name)
{
  m_value = value;
  m_name = name;
  m_pType = pType;
  m_prior = make_shared<statistics::Prior>(statistics::Prior(priorType, discrete_values, weights));
}

// ============================================================================================


cosmobl::statistics::Parameter::Parameter (const double value, const statistics::Prior prior, const ParameterType pType, const string name)
{
  m_value = value;
  m_name = name;
  m_pType = pType;
  m_prior = make_shared<statistics::Prior>(prior);
}

// ============================================================================================


void cosmobl::statistics::Parameter::set_value (const double value)
{
  if (!isFixed())
    m_value = prior()->isIncluded(value) ? value : closest(value, prior()->xmin(), prior()->xmax());
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chains (const int nchains, const int chain_size)
{
  m_nchains = nchains;
  m_chain_size = chain_size;
  m_chains.resize(nchains);

  for (int i=0;i<m_nchains;i++){
    auto chain = make_shared<Chain>(Chain(m_chain_size));
    m_chains[i] = chain;
  }
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chain_value (const int chain, const int position, const double value)
{
  if (chain >= m_nchains)
    ErrorCBL("Error in set_chain value of Parameter.cpp, chain number >= number of chains");
  if (position >= m_chain_size)
    ErrorCBL("Error in set_chain value of Parameter.cpp, position in chain >= chain size");

  m_chains[chain]->set_chain_value(position, value);
}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chains_values (const int position, const vector<double> values)
{
  if (int(values.size()) != m_nchains)
    ErrorCBL("Error in set_chain value of Parameter.cpp, size of values != number of chains");
  if (position >= m_chain_size)
    ErrorCBL("Error in set_chain value of Parameter.cpp, position in chain >= chain size");

  for (int i=0; i<m_nchains; i++)
    m_chains[i]->set_chain_value(position,values[i]);

}


// ============================================================================================


void cosmobl::statistics::Parameter::set_chains_values_from_prior (const int position)
{
  if (position >= m_chain_size)
    ErrorCBL("Error in set_chain value of Parameter.cpp, position in chain >= chain size");

  vector<double> values = sample_from_prior(m_nchains);
  for (int i=0; i<m_nchains; i++)
    m_chains[i]->set_chain_value(position,values[i]);

}


// ============================================================================================


shared_ptr<statistics::Chain> cosmobl::statistics::Parameter::merge_chains (const int max, const int min, const int thin)
{
  vector<double> values;

  for (auto &&cc : m_chains) {
    int cmin = (min<=0) ? 0 : min;
    int cmax = (max<=0) ? cc->chain_size() : max;
    for (int i=cmin; i<cmax; i+=thin)
      values.push_back(cc->chain_value(i));
  }
  
  shared_ptr<statistics::Chain> chain=make_shared<Chain>(statistics::Chain(values.size()));

  for (size_t i=0; i<values.size(); i++)
    chain->set_chain_value(i,values[i]);

  return chain;
}


// ============================================================================================


double cosmobl::statistics::Parameter::random_value () const
{
  return m_prior->sample();
}


// ============================================================================================


double cosmobl::statistics::Parameter::chains_convergence (const int max, const int min, const int thin)
{
  double RR = 0.;
  
  if (!isFixed()) {

    auto chain = merge_chains(max, min, thin);
    chain->Statistics();
    coutCBL << endl << par::col_green << m_name << par::col_default << ":" << endl;
    coutCBL << setprecision(6) << " mean = "  << chain->mean() << endl << " std = " << chain->std() << endl << " median = " << chain->median() << endl;

    set_value(chain->mean());
    m_std = chain->std();

    vector<double> mean,var;

    for (auto &&cc : chains()) {
      cc->Statistics(max, min);
      mean.push_back(cc->mean());
      var.push_back(cc->std()*cc->std());
    }
    
    double W  = Average(var);
    double B = m_chain_size*Sigma(mean)*Sigma(mean);
    RR = sqrt((double(m_chain_size-1)/m_chain_size*W+B/m_chain_size)/W);
    coutCBL << setprecision(6) << "Convergence parameter (sqrt(R)-1) = " << RR-1 << endl;
  }
  
  //else { coutCBL << endl << m_name  << ":" << endl << " fixed value = " << m_value << endl ; }

  return RR;
}


// ============================================================================================


vector<double> cosmobl::statistics::Parameter::sample_from_prior (const int sample_size)
{
  vector<double> values;
  
  if (!isFixed())  
    values = m_prior->sample_vector(sample_size);
  else 
    values.resize(sample_size, m_value);

  return values;
  
}
