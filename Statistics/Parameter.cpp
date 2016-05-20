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


cosmobl::statistics::Parameter::Parameter(const double value, const bool freeze, const string name)
{
  m_value=value;
  m_freeze = freeze;
  m_name = name;
  m_prior= make_shared<Prior>();
}


// ============================================================================================


cosmobl::statistics::Parameter::Parameter(const double value, const double pmin, const double pmax, const bool freeze, const vector<double> discrete_values, const string name)
{
  m_value = value;
  m_name = name;
  m_freeze = freeze;
  m_prior = make_shared<statistics::Prior>(statistics::Prior(statistics::PriorType::_IdentityPrior_,pmin,pmax,discrete_values));
}


// ============================================================================================


cosmobl::statistics::Parameter::Parameter(const double value, const statistics::PriorType priorType, const vector<double> prior_params, const vector<double> Limits, const bool freeze, const vector<double> discrete_values, const string name)
{
  m_value = value;
  m_name = name;
  m_freeze = freeze;
  m_prior = make_shared<statistics::Prior>(statistics::Prior(priorType,prior_params,Limits,discrete_values));
}


// ============================================================================================


cosmobl::statistics::Parameter::Parameter(const double value, const statistics::prior_func func, const vector<double> prior_params, const vector<double> Limits, const bool freeze, const vector<double> discrete_values, const string name)
{
  m_value = value;
  m_name = name;
  m_freeze = freeze;
  m_prior = make_shared<statistics::Prior>(statistics::Prior(statistics::PriorType::_FunctionPrior_, func,prior_params,Limits,discrete_values));
}

// ============================================================================================


cosmobl::statistics::Parameter::Parameter(const double value, const statistics::Prior prior, const bool freeze, const string name)
{
  m_value = value;
  m_name = name;
  m_freeze = freeze;
  m_prior = make_shared<statistics::Prior>(prior);
}

// ============================================================================================


void cosmobl::statistics::Parameter::set_chains(int nchains, int chain_size)
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


shared_ptr<statistics::Chain> cosmobl::statistics::Parameter::merge_chains(int max, int min, int thin)
{
  vector<double> values;

  for(auto &&cc : m_chains){
    int cmin = (min<=0) ? 0 : min;
    int cmax = (max<=0) ? cc->chain_size() : max;
    for(int i=cmin;i<cmax;i+=thin)
      values.push_back(cc->chain_value(i));
  }
  
  shared_ptr<statistics::Chain> chain=make_shared<Chain>(statistics::Chain(values.size()));

  for(size_t i=0;i<values.size();i++)
    chain->set_chain_value(i,values[i]);

  return chain;
}



// ============================================================================================


double cosmobl::statistics::Parameter::eval_proposed(const double proposed_value){
  m_proposed_value = proposed_value;
  return eval_proposed();
}


// ============================================================================================


double cosmobl::statistics::Parameter::eval_proposed() const
{
  return PriorProbability(m_proposed_value)/PriorProbability();
}


// ============================================================================================


void cosmobl::statistics::Parameter::confirm_proposed_value(){
  m_value = m_proposed_value;
}


// ============================================================================================


double cosmobl::statistics::Parameter::random_value() const
{
  return m_prior->sample();

}
