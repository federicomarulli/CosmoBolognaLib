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


// ======================================================================================


cosmobl::statistics::Prior::Prior () 
{

  m_Discrete = 0;
  m_xmin = -par::defaultDouble; 
  m_xmax = par::defaultDouble;
  m_func = &identity<double>;  

  set_inverse_cumulative_distribution();
}


// ======================================================================================


cosmobl::statistics::Prior::Prior (const cosmobl::statistics::PriorType priorType, const double pmin, const double pmax,  const vector<double> discrete_values) 
{
  if(priorType != statistics::PriorType::_IdentityPrior_)
    ErrorMsg("Error in constructor of Prior, this constructor only allows PriorType::_IdentityPrior_");

  m_Discrete = 0;

  m_func = &identity<double>;
  m_xmin = pmin;
  m_xmax = pmax;

  m_prior_normalization = m_xmax-m_xmin;
  m_prior_max_prob = 1./m_prior_normalization;

  set_discrete_values(discrete_values);
  set_inverse_cumulative_distribution();
}


// ======================================================================================


cosmobl::statistics::Prior::Prior (const cosmobl::statistics::PriorType priorType, const vector<double> prior_params, const vector<double> Limits,  const vector<double> discrete_values) 
{
  m_Discrete = 0;
  
  if (Limits.size()==2){
    m_xmin = Min(Limits);
    m_xmax = Max(Limits);
  }

  if (priorType == statistics::PriorType::_GaussianPrior_) {

   m_func = &gaussian<double>;
   set_parameters(prior_params);
   double x1 = (m_xmax-m_prior_func_pars[0])/(sqrt(2)*m_prior_func_pars[1]);
   double x2 = (m_xmin-m_prior_func_pars[0])/(sqrt(2)*m_prior_func_pars[1]);
   m_prior_normalization = 0.5*(erf(x1)-erf(x2));

   if ( m_prior_func_pars[0] < m_xmax && m_prior_func_pars[0]> m_xmin)
     m_prior_max_prob = (*this)(m_prior_func_pars[0]);
   else if ( m_prior_func_pars[0] > m_xmax)
     m_prior_max_prob = (*this)(m_xmax);
   else if ( m_prior_func_pars[0] < m_xmin)
     m_prior_max_prob = (*this)(m_xmin);
  }
  else if (priorType == statistics::PriorType::_PoissonPrior_)
    set_func_parameters(&poisson<double>, prior_params);
  
  else 
    ErrorMsg("Error in Prior constructor of Prior.cpp. No such type of prior");

  set_inverse_cumulative_distribution();
  set_discrete_values(discrete_values);
}


// ======================================================================================


cosmobl::statistics::Prior::Prior (const cosmobl::statistics::PriorType priorType, const prior_func func, const vector<double> prior_params, const vector<double> Limits,  const vector<double> discrete_values) 
{
  if(priorType != statistics::PriorType::_FunctionPrior_)
    ErrorMsg("Error in constructor of Prior, this constructor only allows PriorType::_FunctionPrior_");

  m_Discrete = 0;
 
  if (Limits.size()==2){
    m_xmin = Min(Limits);
    m_xmax = Max(Limits);
  }

  set_func_parameters(func, prior_params);
  set_discrete_values(discrete_values);

  set_inverse_cumulative_distribution();
}

// =====================================================================================


void cosmobl::statistics::Prior::set_limits (const double pmin, const double pmax) {
  m_xmin = pmin; 
  m_xmax = pmax;
}


// =====================================================================================


void cosmobl::statistics::Prior::set_parameters (const vector<double> pars) 
{
   m_prior_func_pars.erase(m_prior_func_pars.begin(),m_prior_func_pars.end());

   for (unsigned int i=0; i<pars.size(); i++)
      m_prior_func_pars.push_back(pars[i]);
}


// =====================================================================================


void cosmobl::statistics::Prior::set_func_parameters (const prior_func _func, const vector<double> pars) 
{
   m_func = _func;
   
   set_parameters(pars);

   set_distribution_parameters();
}


// =====================================================================================


void cosmobl::statistics::Prior::set_discrete_values (const vector<double> discrete_values) 
{
  if(discrete_values.size()>0){
    m_Discrete = 1;
    for(size_t i=0; i< discrete_values.size();i++)
      m_discrete_values.push_back(discrete_values[i]);
  }
}


// =====================================================================================


void cosmobl::statistics::Prior::set_distribution_parameters () 
{

  int nPt = 10000;
  vector<double> x = linear_bin_vector(nPt, m_xmin, m_xmax);
  vector<double> fx;

  shared_ptr<void> pp = NULL;

  for(int i=0;i<nPt;i++)
    fx.push_back(m_func(x[i], pp, m_prior_func_pars));

  double dx = x[1]-x[0];
  
  double normalization=0;
  for(int i=0;i<nPt;i++){
    normalization += fx[i]*dx;
  }
  
  m_prior_normalization = normalization;
  m_prior_max_prob = Max(fx)/normalization;

}


// =====================================================================================


bool cosmobl::statistics::Prior::isIncluded (const double value) const
{
  if (value > m_xmin && m_xmax > value)
    return 1;
  else 
    return 0;
}


// =====================================================================================


double cosmobl::statistics::Prior::apply_discrete (const double value) const
{
  if (!m_Discrete)
    return value;
  else{
    return closest(value,m_discrete_values);
  }
}


// =====================================================================================


double cosmobl::statistics::Prior::sample ()
{
  time_t seed; time(&seed);
  default_random_engine generator(seed);
  uniform_real_distribution<double> distribution(0.0,1.0);

  double value;
  bool done=0;

  while (!done){
    double random = distribution(generator);
    value = (m_xmax-m_xmin)*random+m_xmin;
    double prob = (*this)(value);
    if(prob > distribution(generator)*m_prior_max_prob)
      done = 1;
  }

  return apply_discrete(value);

}


// =====================================================================================


double cosmobl::statistics::Prior::sample (const int seed)
{
  default_random_engine generator(seed);
  uniform_real_distribution<double> distribution(0.0,1.0);

  double value;
  bool done=0;

  while (!done){
    double random = distribution(generator);
    value = (m_xmax-m_xmin)*random+m_xmin;
    double prob = (*this)(value);
    if(prob > distribution(generator)*m_prior_max_prob)
      done = 1;
  }

  return apply_discrete(value);
}


// =====================================================================================


vector<double> cosmobl::statistics::Prior::sample_vector (const int nvalues)
{
  vector<double> values;

  time_t seed; time(&seed);
  default_random_engine generator(seed);
  uniform_real_distribution<double> distribution(0.0,1.0);

  for(int i=0; i< nvalues;i++){
    double value;
    bool done=0;

    while (!done){
      double random = distribution(generator);
      value = (m_xmax-m_xmin)*random+m_xmin;
      double prob = (*this)(value);
      if(prob > distribution(generator)*m_prior_max_prob)
	done = 1;
    }
    values.push_back(apply_discrete(value));
  }

  return values;

}


// =====================================================================================


void cosmobl::statistics::Prior::set_inverse_cumulative_distribution () 
{
  int nn = (m_xmax-m_xmin)/0.001;

  vector<double> xx = linear_bin_vector(nn, m_xmin, m_xmax);
  vector<double> fx(nn,0), FX(nn,0);

  for(int i=0;i<nn;i++)
    fx[i] = (*this)(xx[i]);

  classfunc::func_grid_GSL ff(xx,fx,"Spline");
  double norm = ff.integrate_qag(m_xmin,m_xmax);

  for(int i=0;i<nn;i++)
    FX[i] = ff.integrate_qag(m_xmin,xx[i])/norm;

  m_prior_inverse_cumulative = make_shared<classfunc::func_grid_GSL>(FX,xx,"Spline");
}


// =====================================================================================


double cosmobl::statistics::Prior::sample (const double prob)
{
  return apply_discrete(m_prior_inverse_cumulative->operator()(prob));
}
