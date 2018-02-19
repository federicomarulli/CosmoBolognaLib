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
 *  @file Func/Distribution.cpp
 *
 *  @brief Methods of the class Distribution
 *
 *  This file contains the implementation of the methods of the class
 *  Distribution
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Distribution.h"

using namespace cosmobl;
using namespace random;


// ======================================================================================


double cosmobl::glob::Distribution::m_percentile_integrator (const double xx)
{
  function<double(double)> f_moment = [this] (double x) {return this->operator()(x);};
  return gsl::GSL_integrate_qag(f_moment, m_xmin, xx);
}


// ======================================================================================


void cosmobl::glob::Distribution::m_set_distribution_normalization ()
{
  function<double(double)> f = bind(m_func, std::placeholders::_1, m_distribution_func_fixed_pars, m_distribution_func_pars);

  m_distribution_normalization = gsl::GSL_integrate_qag(f, m_xmin, m_xmax);
}


// ======================================================================================


cosmobl::glob::Distribution::Distribution (const cosmobl::glob::DistributionType distributionType, const double value) 
{
  if (distributionType != glob::DistributionType::_ConstantDistribution_)
    ErrorCBL("Error in constructor of Distribution, this constructor only allows DistributionType::_ConstaintDistribution_");

  set_constant_distribution(value);
}



// ======================================================================================


cosmobl::glob::Distribution::Distribution (const cosmobl::glob::DistributionType distributionType, const double xmin, const double xmax, const int seed) 
{
  if (distributionType != glob::DistributionType::_UniformDistribution_)
    ErrorCBL("Error in constructor of Distribution, this constructor only allows DistributionType::_UniformDistribution_");

  set_uniform_distribution(xmin, xmax, seed);
}


// ======================================================================================


cosmobl::glob::Distribution::Distribution (const cosmobl::glob::DistributionType distributionType, const vector<double> distribution_params, const double xmin, const double xmax , const int seed) 
{
  set_limits(xmin, xmax);

  if (distributionType == glob::DistributionType::_GaussianDistribution_) {
    if (distribution_params.size() != 2)
      ErrorCBL("Error in constructor of Distribution, wrong size of distribution_params. Gaussian distribution needs 2 parameters, the mean and the standard deviation");

    set_gaussian_distribution(distribution_params[0], distribution_params[1], seed);
  }
  
  else if (distributionType == glob::DistributionType::_PoissonDistribution_) {
    if (distribution_params.size() != 1)
      ErrorCBL("Error in constructor of Distribution, wrong size of distribution_params. Poisson distribution needs 1 parameter, the mean");

    set_poisson_distribution(distribution_params[0], seed);
  }
  
  else 
    ErrorCBL("Error in Distribution constructor of Distribution.cpp. No such type of distribution");
}


// ======================================================================================


 cosmobl::glob::Distribution::Distribution (const DistributionType distributionType, const distribution_func func, const shared_ptr<void> distribution_fixed_pars, const vector<double> distribution_pars, const double xmin, const double xmax, const int seed)
{
  set_limits(xmin, xmax);

  if (distributionType != glob::DistributionType::_CustomDistribution_)
    ErrorCBL("Error in constructor of Distribution, this constructor only allows DistributionType::_CustomDistribution_");

  set_custom_distribution(func, distribution_fixed_pars, distribution_pars, seed);
}


// ======================================================================================


cosmobl::glob::Distribution::Distribution (const cosmobl::glob::DistributionType distributionType, const vector<double> discrete_values, const vector<double> weights, const int seed) 
{
  if (distributionType != glob::DistributionType::_DiscreteDistribution_)
    ErrorCBL("Error in constructor of Distribution, this constructor only allows DistributionType::_DiscreteDistribution_");

  set_discrete_values(discrete_values, weights, seed);
}


// ======================================================================================


cosmobl::glob::Distribution::Distribution (const cosmobl::glob::DistributionType distributionType, const vector<double> var, const vector<double> dist, const int nbin, const string interpolationType, const int seed) 
{
  if (distributionType == glob::DistributionType::_DiscreteDistribution_)
  {
    vector<double> vv, dd, edd;
    distribution(vv, dd, edd, var, dist, nbin);

    set_binned_distribution(vv, dd, interpolationType, seed);
  }
  else if (distributionType == glob::DistributionType::_InterpolatedDistribution_)
  {
    set_binned_distribution(var, dist, interpolationType, seed);
  }
  else
    ErrorCBL("Error in Distribution constructor of Distribution.cpp. No such type of distribution");
}


// ======================================================================================


void cosmobl::glob::Distribution::set_limits(const double xmin, const double xmax)
{
  m_xmin = xmin; 
  m_xmax = xmax;
}


// ======================================================================================


void cosmobl::glob::Distribution::set_constant_distribution (const double value)
{
  m_distributionType = _ConstantDistribution_;

  set_limits(par::defaultDouble, -par::defaultDouble);

  m_distribution_random = make_shared<ConstantRandomNumbers> (ConstantRandomNumbers(value));

  m_func = &identity<double>; 
  m_distribution_normalization = 1.;
}

// ======================================================================================


void cosmobl::glob::Distribution::set_uniform_distribution (const double xmin, const double xmax, const int seed)
{

  m_distributionType = _UniformDistribution_;

  set_limits(xmin, xmax);
  m_distribution_func_pars.erase(m_distribution_func_pars.begin(), m_distribution_func_pars.end());

  m_distribution_func_pars.push_back(m_xmax);
  m_distribution_func_pars.push_back(m_xmax);

  m_distribution_random = make_shared<UniformRandomNumbers> (UniformRandomNumbers(m_xmin, m_xmax, seed));
  m_func = &identity<double>; 

  m_distribution_normalization = m_xmax-m_xmin;
}


// ======================================================================================


void cosmobl::glob::Distribution::set_gaussian_distribution (const double mean, const double sigma, const int seed)
{

  m_distributionType = _GaussianDistribution_;

  m_distribution_func_pars.erase(m_distribution_func_pars.begin(), m_distribution_func_pars.end());
  m_distribution_func_pars.push_back(mean);
  m_distribution_func_pars.push_back(sigma);

  m_distribution_random = make_shared<NormalRandomNumbers> (NormalRandomNumbers(mean, sigma, seed, m_xmin, m_xmax));
  m_func = &gaussian<double>; 

  m_distribution_normalization = 0.5*(erf((m_xmax-mean)/sigma)-erf((m_xmin-mean)/sigma));
}


// ======================================================================================


void cosmobl::glob::Distribution::set_poisson_distribution (const double mean, const int seed) 
{

  m_distributionType = _PoissonDistribution_;

  m_xmin = nint(m_xmin);
  m_xmax = nint(m_xmax);
  int nbins = m_xmax-m_xmin;

  vector<double> poisson_values = linear_bin_vector(nbins, m_xmin, m_xmax);
  vector<double> weights;
  for(int i=0;i<nbins;i++)
    weights.push_back(poisson(poisson_values[i],NULL,{mean}));

  m_distribution_random = make_shared<DiscreteRandomNumbers> (DiscreteRandomNumbers(poisson_values, weights, seed, m_xmin, m_xmax));

  glob::STR_closest_probability parameters;
  parameters.values = poisson_values;
  parameters.weights = weights; 

  m_distribution_func_fixed_pars = make_shared<glob::STR_closest_probability>(parameters);
  m_func = &closest_probability; 

  m_distribution_normalization = accumulate(weights.begin(), weights.end(), 0);
}


// =====================================================================================


void cosmobl::glob::Distribution::set_discrete_values (const vector<double> discrete_values, const vector<double> weights, const int seed) 
{

  m_distributionType = _DiscreteDistribution_;

  if (discrete_values.size()==0)
    ErrorCBL("Error in set_discrete_values of Distribution. Vector of values is empty");

  set_limits(Min(discrete_values), Max(discrete_values));

  m_distribution_random = make_shared<DiscreteRandomNumbers> (DiscreteRandomNumbers(discrete_values, weights, seed, m_xmin, m_xmax));

  vector<double> ww = weights;

  if (ww.size() == 0) {
    ww.resize(discrete_values.size(), 1);
  }


  glob::STR_closest_probability parameters;
  parameters.values = discrete_values;
  parameters.weights = ww; 

  m_distribution_func_fixed_pars = make_shared<glob::STR_closest_probability>(parameters);

  m_func = &closest_probability; 
  m_distribution_normalization = accumulate(ww.begin(), ww.end(), 0);

}


// =====================================================================================


void cosmobl::glob::Distribution::set_custom_distribution (const distribution_func func, const shared_ptr<void> distribution_fixed_pars, const vector<double> distribution_pars, const int seed)
{

  m_distributionType = _CustomDistribution_;

  m_func = func;
  m_distribution_func_fixed_pars = distribution_fixed_pars;
  m_distribution_func_pars = distribution_pars;

  m_distribution_random = make_shared<CustomDistributionRandomNumbers> (CustomDistributionRandomNumbers(m_func, m_distribution_func_fixed_pars, m_distribution_func_pars, seed, m_xmin, m_xmax));

  m_set_distribution_normalization ();

}


// =====================================================================================


void cosmobl::glob::Distribution::set_binned_distribution (const vector<double> var, const vector<double> dist, const string interpolationType, const int seed)
{
  m_distributionType = _InterpolatedDistribution_;

  if (var.size()==0)
    ErrorCBL("Error in set_discrete_values of Distribution. Vector of values is empty");

  set_limits(Min(var), Max(var));
  m_distribution_random = make_shared<DistributionRandomNumbers> (DistributionRandomNumbers(var, dist, interpolationType, seed));

  glob::STR_distribution_probability parameters;
  parameters.func = make_shared<glob::FuncGrid>(glob::FuncGrid(var, dist, interpolationType));

  m_distribution_func_fixed_pars = make_shared<glob::STR_distribution_probability>(parameters);

  m_func = distribution_probability; 
  m_set_distribution_normalization();
}


// =====================================================================================


bool cosmobl::glob::Distribution::isIncluded (const double value) const
{
  if (value > m_xmin && m_xmax > value)
    return true;
  else 
    return false;
}


// =====================================================================================


double cosmobl::glob::Distribution::sample () const
{
  return m_distribution_random->operator()();
}


// =====================================================================================


double cosmobl::glob::Distribution::sample (const int seed)
{
  m_distribution_random->set_seed(seed);
  return sample();
}


// =====================================================================================


vector<double> cosmobl::glob::Distribution::sample_vector (const int nvalues)
{
  vector<double> values;
  
  for (int i=0; i<nvalues; i++)
    values.push_back(sample());

  return values;

}


// =====================================================================================


double cosmobl::glob::Distribution::mean()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_DiscreteDistribution_) {
    shared_ptr<STR_closest_probability> pp  = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    val = Average(pp->values, pp->weights);
  }
  else
  {
    function<double(double)> f = bind(&Distribution::m_moments_integrator, this, std::placeholders::_1, 1);

    val = gsl::GSL_integrate_qag(f, m_xmin, m_xmax);
  }
  m_mean = val;
  return val;
}


// =====================================================================================


double cosmobl::glob::Distribution::variance ()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_DiscreteDistribution_) {
    shared_ptr<STR_closest_probability> pp  = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    val = pow(Sigma(pp->values, pp->weights),2);
  }
  else
  {
    mean();

    function<double(double)> f = [this] (double xx) {return  this->m_central_moments_integrator(xx, 2);};

    val = gsl::GSL_integrate_qag(f, m_xmin, m_xmax);
  }
  m_variance = val;
  return val;
}


// =====================================================================================


double cosmobl::glob::Distribution::std ()
{
  return sqrt(variance());
}


// =====================================================================================


double cosmobl::glob::Distribution::skewness ()
{
  double val = 0.;

  if (m_distributionType == glob::DistributionType::_DiscreteDistribution_)
    ErrorCBL("Work in progress!", glob::ExitCode::_workInProgress_);
  
  else {
    variance();
    
    function<double(double)> f = [this] (double xx) { return  this->m_central_moments_integrator(xx, 3); };
    
    val = sqrt(pow(gsl::GSL_integrate_qag(f, m_xmin, m_xmax),2)*pow(m_variance, -3));
  }
  
  return val;
}


// =====================================================================================


double cosmobl::glob::Distribution::kurtosis ()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_DiscreteDistribution_)
    ErrorCBL("Work in progress!", glob::ExitCode::_workInProgress_);
  
  else {
    variance();

    function<double(double)> f = [this] (double xx) { return  this->m_central_moments_integrator(xx, 4); };

    val= gsl::GSL_integrate_qag(f, m_xmin, m_xmax)*pow(m_variance, -2);
  }
  
  return val;
}


// =====================================================================================


vector<double> cosmobl::glob::Distribution::moments ()
{
  return {mean(), std(), skewness(), kurtosis()};
}


// =====================================================================================


double cosmobl::glob::Distribution::percentile (const unsigned int i)
{
  double Area = double(i)/100.;
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_DiscreteDistribution_) {
    shared_ptr<STR_closest_probability> pp  = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    vector<double> vv = pp->values; sort(vv.begin(), vv.end());
    int perc = nint(double(i)/100.*pp->values.size());
    val = vv[perc];
  }
  
  else {
    function<double(double)> f_integral = [this] (double xx) { return this->m_percentile_integrator(xx); };

    val = gsl::GSL_root_brent(f_integral, Area, m_xmin, m_xmax);
  }

  return val;
}


// =====================================================================================


double cosmobl::glob::Distribution::mode ()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_DiscreteDistribution_) {

    int digits = 4;
    shared_ptr<STR_closest_probability> pp  = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    vector<double> vv = pp->values; 
    for(size_t i=0; i<vv.size(); i++)
      vv[i] = round_to_precision (vv[i], digits);
    
    vector<double> unique_vv = vv; 
    unique_unsorted(unique_vv);
    sort(unique_vv.begin(), unique_vv.end());

    int counts = -1;
    for(size_t i=0; i<unique_vv.size(); i++){
       int cc   = std::count(vv.begin(), vv.end(), unique_vv[i]);
       if(cc>counts){
	 val = unique_vv[i];
	 counts = cc;
       }
    }

  }
  else {
    auto func = [&] (const double param){
      return -this->operator()(param);
    };

    double start = (m_xmax+m_xmin)*0.5;
    val= gsl::GSL_minimize_1D(func, start, m_xmin, m_xmax);
  }
  return val;
}
