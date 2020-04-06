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
 *  @file Distribution/Distribution.cpp
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

using namespace std;

using namespace cbl;
using namespace random;


// ======================================================================================


double cbl::glob::Distribution::m_percentile_integrator (const double xx)
{
  function<double(double)> f_moment = [this] (double x) {return this->operator()(x);};
  return wrapper::gsl::GSL_integrate_cquad(f_moment, m_xmin, xx, 1.e-4);
}


// ======================================================================================


void cbl::glob::Distribution::m_set_distribution_normalization ()
{
  function<double(double)> f = bind(m_func, std::placeholders::_1, m_distribution_func_fixed_pars, m_distribution_func_pars);

  m_distribution_normalization = wrapper::gsl::GSL_integrate_cquad(f, m_xmin, m_xmax, 1.e-4);
  m_log_distribution_normalization = log(m_distribution_normalization);
}


// ======================================================================================


cbl::glob::Distribution::Distribution (const cbl::glob::DistributionType distributionType, const double value) 
{
  if (distributionType!=glob::DistributionType::_Constant_)
    ErrorCBL("this constructor only allows DistributionType::_Constant_ !", "Distribution", "Distribution.cpp");

  set_constant_distribution(value);
}



// ======================================================================================


cbl::glob::Distribution::Distribution (const cbl::glob::DistributionType distributionType, const double xmin, const double xmax, const int seed) 
{
  if (distributionType != glob::DistributionType::_Uniform_)
    ErrorCBL("this constructor only allows DistributionType::_Uniform_ !", "Distribution", "Distribution.cpp");

  set_uniform_distribution(xmin, xmax, seed);
}


// ======================================================================================


cbl::glob::Distribution::Distribution (const cbl::glob::DistributionType distributionType, const std::vector<double> distribution_params, const double xmin, const double xmax , const int seed) 
{
  set_limits(xmin, xmax);

  if (distributionType == glob::DistributionType::_Gaussian_) {
    if (distribution_params.size() != 2)
      ErrorCBL("wrong size of distribution_params: Gaussian distribution needs 2 parameters, the mean and the standard deviation!", "Distribution", "Distribution.cpp");

    set_gaussian_distribution(distribution_params[0], distribution_params[1], seed);
  }
  
  else if (distributionType == glob::DistributionType::_Poisson_) {
    if (distribution_params.size() != 1)
      ErrorCBL("wrong size of distribution_params: poisson distribution needs 1 parameter, the mean", "Distribution", "Distribution.cpp");

    set_poisson_distribution(distribution_params[0], seed);
  }
  
  else 
    ErrorCBL("no such type of distribution", "Distribution", "Distribution.cpp");
}


// ======================================================================================


cbl::glob::Distribution::Distribution (const DistributionType distributionType, const distribution_func func, const std::shared_ptr<void> distribution_fixed_pars, const std::vector<double> distribution_pars, const double xmin, const double xmax, const int seed)
{
  set_limits(xmin, xmax);

  if (distributionType != glob::DistributionType::_Custom_)
    ErrorCBL("this constructor only allows DistributionType::_Custom_ !", "Distribution", "Distribution.cpp");

  set_custom_distribution(func, distribution_fixed_pars, distribution_pars, seed);
}


// ======================================================================================


cbl::glob::Distribution::Distribution (const cbl::glob::DistributionType distributionType, const std::vector<double> discrete_values, const std::vector<double> weights, const int seed) 
{
  if (distributionType != glob::DistributionType::_Discrete_)
    ErrorCBL("this constructor only allows DistributionType::_Discrete_ !", "Distribution", "Distribution.cpp");

  set_discrete_values(discrete_values, weights, seed);
}


// ======================================================================================


cbl::glob::Distribution::Distribution (const cbl::glob::DistributionType distributionType, const std::vector<double> var, const std::vector<double> dist, const int nbin, const std::string interpolationType, const int seed) 
{
  if (distributionType == glob::DistributionType::_Discrete_)
  {
    vector<double> vv, dd, edd;
    get_distribution(vv, dd, edd, var, dist, nbin);

    set_binned_distribution(vv, dd, interpolationType, seed);
  }
  else if (distributionType == glob::DistributionType::_Interpolated_)
  {
    set_binned_distribution(var, dist, interpolationType, seed);
  }
  else
    ErrorCBL("no such type of distribution", "Distribution", "Distribution.cpp");
}


// ======================================================================================


double cbl::glob::Distribution::operator() (const double xx)
{
  if (xx<m_xmin || xx>m_xmax) return 0;
  else return m_func(xx, m_distribution_func_fixed_pars, m_distribution_func_pars)/m_distribution_normalization;
}


// ======================================================================================


double cbl::glob::Distribution::log_distribution (const double xx)
{
  if (xx<m_xmin || xx>m_xmax) return par::defaultDouble;
  else return log(m_func(xx, m_distribution_func_fixed_pars, m_distribution_func_pars))-m_log_distribution_normalization;
}


// ======================================================================================


void cbl::glob::Distribution::set_limits(const double xmin, const double xmax)
{
  m_xmin = xmin; 
  m_xmax = xmax;
}


// ======================================================================================


void cbl::glob::Distribution::set_constant_distribution (const double value)
{
  m_distributionType = glob::DistributionType::_Constant_;

  set_limits(par::defaultDouble, -par::defaultDouble);

  m_distribution_random = make_shared<ConstantRandomNumbers> (ConstantRandomNumbers(value));

  m_func = &identity<double>; 
  m_distribution_normalization = 1.;
  m_log_distribution_normalization = 0.;
}

// ======================================================================================


void cbl::glob::Distribution::set_uniform_distribution (const double xmin, const double xmax, const int seed)
{

  m_distributionType = glob::DistributionType::_Uniform_;

  set_limits(xmin, xmax);
  m_distribution_func_pars.erase(m_distribution_func_pars.begin(), m_distribution_func_pars.end());

  m_distribution_func_pars.push_back(m_xmax);
  m_distribution_func_pars.push_back(m_xmax);

  m_distribution_random = make_shared<UniformRandomNumbers> (UniformRandomNumbers(m_xmin, m_xmax, seed));
  m_func = &identity<double>; 

  m_distribution_normalization = m_xmax-m_xmin;
  m_log_distribution_normalization = log(m_xmax-m_xmin);
}


// ======================================================================================


void cbl::glob::Distribution::set_gaussian_distribution (const double mean, const double sigma, const int seed)
{

  m_distributionType = glob::DistributionType::_Gaussian_;

  m_distribution_func_pars.erase(m_distribution_func_pars.begin(), m_distribution_func_pars.end());
  m_distribution_func_pars.push_back(mean);
  m_distribution_func_pars.push_back(sigma);

  m_distribution_random = make_shared<NormalRandomNumbers> (NormalRandomNumbers(mean, sigma, seed, m_xmin, m_xmax));
  m_func = &gaussian<double>; 

  m_distribution_normalization = 0.5*(erf((m_xmax-mean)/sigma)-erf((m_xmin-mean)/sigma));
  m_log_distribution_normalization = log(m_distribution_normalization);
}


// ======================================================================================


void cbl::glob::Distribution::set_poisson_distribution (const double mean, const int seed) 
{

  m_distributionType = glob::DistributionType::_Poisson_;

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
  m_log_distribution_normalization = log(m_distribution_normalization);
}


// =====================================================================================


void cbl::glob::Distribution::set_discrete_values (const std::vector<double> discrete_values, const std::vector<double> weights, const int seed) 
{
  m_distributionType = glob::DistributionType::_Discrete_;

  if (discrete_values.size()==0)
    ErrorCBL("the input vector of discrete values is empty", "set_discrete_values", "Distribution.cpp");

  set_limits(Min(discrete_values), Max(discrete_values));

  m_distribution_random = make_shared<DiscreteRandomNumbers> (DiscreteRandomNumbers(discrete_values, weights, seed, m_xmin, m_xmax));

  vector<double> ww = weights;

  if (ww.size() == 0) 
    ww.resize(discrete_values.size(), 1);

  glob::STR_closest_probability parameters;
  parameters.values = discrete_values;
  parameters.weights = ww; 

  m_distribution_func_fixed_pars = make_shared<glob::STR_closest_probability>(parameters);

  m_func = &closest_probability; 
  m_distribution_normalization = accumulate(ww.begin(), ww.end(), 0);
  m_log_distribution_normalization = log(m_distribution_normalization);
}


// =====================================================================================


void cbl::glob::Distribution::set_custom_distribution (const distribution_func func, const shared_ptr<void> distribution_fixed_pars, const std::vector<double> distribution_pars, const int seed)
{

  m_distributionType = glob::DistributionType::_Custom_;

  m_func = func;
  m_distribution_func_fixed_pars = distribution_fixed_pars;
  m_distribution_func_pars = distribution_pars;

  m_distribution_random = make_shared<CustomDistributionRandomNumbers> (CustomDistributionRandomNumbers(m_func, m_distribution_func_fixed_pars, m_distribution_func_pars, seed, m_xmin, m_xmax));

  m_set_distribution_normalization();

}


// =====================================================================================


void cbl::glob::Distribution::set_binned_distribution (const std::vector<double> var, const std::vector<double> dist, const std::string interpolationType, const int seed)
{
  m_distributionType = glob::DistributionType::_Interpolated_;

  if (var.size()==0)
    ErrorCBL("the input vector is empty", "set_binned_distribution", "Distribution.cpp");

  set_limits(Min(var), Max(var));
  m_distribution_random = make_shared<DistributionRandomNumbers> (DistributionRandomNumbers(var, dist, interpolationType, seed));

  glob::STR_distribution_probability parameters;
  parameters.func = make_shared<glob::FuncGrid>(glob::FuncGrid(var, dist, interpolationType));

  m_distribution_func_fixed_pars = make_shared<glob::STR_distribution_probability>(parameters);

  m_func = distribution_probability; 
  m_set_distribution_normalization();
}


// =====================================================================================


bool cbl::glob::Distribution::isIncluded (const double value) const
{
  if (value > m_xmin && m_xmax > value)
    return true;
  else 
    return false;
}


// =====================================================================================


double cbl::glob::Distribution::sample () const
{
  return m_distribution_random->operator()();
}


// =====================================================================================


double cbl::glob::Distribution::sample (const int seed)
{
  m_distribution_random->set_seed(seed);
  return sample();
}


// =====================================================================================


vector<double> cbl::glob::Distribution::sample_vector (const int nvalues)
{
  vector<double> values;
  
  for (int i=0; i<nvalues; i++)
    values.push_back(sample());

  return values;

}


// =====================================================================================


double cbl::glob::Distribution::mean()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_Discrete_) {
    shared_ptr<STR_closest_probability> pp = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    val = Average(pp->values, pp->weights);
  }
  else
  {
    function<double(double)> f = bind(&Distribution::m_moments_integrator, this, std::placeholders::_1, 1);

    val = wrapper::gsl::GSL_integrate_cquad(f, m_xmin, m_xmax, 1.e-4);
  }
  m_mean = val;
  return val;
}


// =====================================================================================


double cbl::glob::Distribution::variance ()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_Discrete_) {
    shared_ptr<STR_closest_probability> pp = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    val = pow(Sigma(pp->values, pp->weights),2);
  }
  else
  {
    mean();

    function<double(double)> f = [this] (double xx) {return  this->m_central_moments_integrator(xx, 2);};

    val = wrapper::gsl::GSL_integrate_cquad(f, m_xmin, m_xmax, 1.e-4);
  }
  m_variance = val;
  return val;
}


// =====================================================================================


double cbl::glob::Distribution::std ()
{
  return sqrt(variance());
}


// =====================================================================================


double cbl::glob::Distribution::skewness ()
{
  double val = 0.;

  if (m_distributionType == glob::DistributionType::_Discrete_)
    ErrorCBL("", "skewness", "Distribution.cpp", glob::ExitCode::_workInProgress_);
  
  else {
    variance();
    
    function<double(double)> f = [this] (double xx) { return  this->m_central_moments_integrator(xx, 3); };
    
    val = sqrt(pow(wrapper::gsl::GSL_integrate_qag(f, m_xmin, m_xmax, 1.e-4),2)*pow(m_variance, -3));
  }
  
  return val;
}


// =====================================================================================


double cbl::glob::Distribution::kurtosis ()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_Discrete_)
    ErrorCBL("", "kurtosis", "Distribution.cpp", glob::ExitCode::_workInProgress_);
  
  else {
    variance();

    function<double(double)> f = [this] (double xx) { return  this->m_central_moments_integrator(xx, 4); };

    val= wrapper::gsl::GSL_integrate_qag(f, m_xmin, m_xmax, 1.e-4)*pow(m_variance, -2);
  }
  
  return val;
}


// =====================================================================================


vector<double> cbl::glob::Distribution::moments ()
{
  return {mean(), std(), skewness(), kurtosis()};
}


// =====================================================================================


double cbl::glob::Distribution::percentile (const unsigned int i)
{
  double Area = double(i)/100.;
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_Discrete_) {
    shared_ptr<STR_closest_probability> pp = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    vector<double> vv = pp->values; sort(vv.begin(), vv.end());
    int perc = nint(double(i)/100.*pp->values.size());
    val = vv[perc];
  }
  
  else {
    function<double(double)> f_integral = [this] (double xx) { return this->m_percentile_integrator(xx); };

    val = wrapper::gsl::GSL_root_brent(f_integral, Area, m_xmin, m_xmax);
  }

  return val;
}


// =====================================================================================


double cbl::glob::Distribution::mode ()
{
  double val = 0.;
  
  if (m_distributionType == glob::DistributionType::_Discrete_) {

    int digits = 4;
    shared_ptr<STR_closest_probability> pp = static_pointer_cast<glob::STR_closest_probability>(m_distribution_func_fixed_pars);
    vector<double> vv = pp->values; 
    for (size_t i=0; i<vv.size(); i++)
      vv[i] = round_to_precision (vv[i], digits);
    
    vector<double> unique_vv = vv; 
    unique_unsorted(unique_vv);
    sort(unique_vv.begin(), unique_vv.end());

    int counts = -1;
    for (size_t i=0; i<unique_vv.size(); i++) {
      int cc = std::count(vv.begin(), vv.end(), unique_vv[i]);
      if (cc>counts) {
	val = unique_vv[i];
	counts = cc;
      }
    }
  }
  
  else {
    auto func = [&] (const double param) { return -this->operator()(param); };
    double start = median(); //(m_xmax+m_xmin)*0.5;
    val = wrapper::gsl::GSL_minimize_1D(func, start, m_xmin, m_xmax, 10000, false);
  }
  
  return val;
}


// =====================================================================================


void cbl::glob::Distribution::get_distribution (vector<double> &xx, vector<double> &fx, vector<double> &err, const std::vector<double> FF, const std::vector<double> WW, const int nbin, const bool linear, const std::string file_out, const double fact, const double V1, const double V2, const bool bin_type, const bool conv, const double sigma)
{
  if (xx.size()>0 || fx.size()>0 || FF.size()<=0 || nbin<=0) ErrorCBL("the following conditions have to be satisfied: xx.size()<=0, fx.size()<=0, FF.size()>0 and nbin>0. The values recived are instead: xx.size() = "+cbl::conv(xx.size(), par::fINT)+", fx.size() = "+cbl::conv(fx.size(), par::fINT)+", FF.size() = "+cbl::conv(FF.size(), par::fINT)+" and nbin = "+cbl::conv(nbin, par::fINT)+"!", "get_distribution", "Distribution.cpp");

  double minFF = (V1>cbl::par::defaultDouble) ? V1 : Min(FF)*0.9999;
  double maxFF = (V2>cbl::par::defaultDouble) ? V2 : Max(FF)*1.0001;

  
  // using GSL to create the histogram 

  gsl_histogram *histo = gsl_histogram_alloc(nbin);

  if (linear) gsl_histogram_set_ranges_uniform(histo, minFF, maxFF);

  else {
    vector<double> vv = logarithmic_bin_vector(nbin+1, minFF, maxFF);
    double *vvv = new double[nbin+1]; for (int i=0; i<nbin+1; i++) vvv[i] = vv[i];
    gsl_histogram_set_ranges(histo, vvv, nbin+1);
  }

  vector<double> Weight = WW;
  if (Weight.size()==0) Weight.resize(FF.size(), 1.);
  checkDim(Weight, FF.size(), "WW");
  
  for (size_t i=0; i<FF.size(); i++)
    gsl_histogram_accumulate(histo, FF[i], Weight[i]);
  
  double x1, x2;

  for (int i=0; i<nbin; i++) {

    gsl_histogram_get_range(histo, i, &x1, &x2);
    double val = gsl_histogram_get(histo, i);
    
    if (linear) xx.push_back(0.5*(x1+x2));
    else xx.push_back(pow(10., 0.5*(log10(x1)+log10(x2))));

    if (bin_type) {
      fx.push_back(val/((x2-x1)*fact));
      err.push_back(sqrt(val)/((x2-x1)*fact));
    }
    else {
      fx.push_back(val/((log10(x2)-log10(x1))*fact));
      err.push_back(sqrt(val)/((log10(x2)-log10(x1))*fact));
    }
    
  }

  
  // Gaussian convolution

  if (conv) {
    coutCBL << "The distribution is smoothed with a Gaussian filter" << endl;
    double *func;
    fftw_complex *func_tr;

    if (!linear) ErrorCBL("", "get_distribution", "Distribution.cpp", ExitCode::_workInProgress_);
    int nbinN = 2*nbin;
    int i1 = nbin*0.5, i2 = 1.5*nbin;

    int nbinK = 0.5*nbinN+1;

    func = fftw_alloc_real(nbinN);
    func_tr = fftw_alloc_complex(nbinK);

    for (int i=0; i<nbinN; i++)
      func[i] = 0;
    
    for (int i=i1; i<i2; i++)
      func[i] = fx[i-i1];
    
    for (int i=0; i<nbinK; i++) {
      func_tr[i][0] = 0;
      func_tr[i][1] = 0;
    }

    fftw_plan real2complex;
    real2complex = fftw_plan_dft_r2c_1d(nbinN, func, func_tr, FFTW_ESTIMATE);
    fftw_execute(real2complex);
    fftw_destroy_plan(real2complex);

    double delta = (maxFF-minFF)/nbin;
    double SS = pow(sigma,2);

    double fact = 2*par::pi/(nbinN*delta);
    for (int i=0; i<nbinK; i++) {
      double kk = i*fact;
      func_tr[i][0] = func_tr[i][0]*exp(-0.5*kk*kk*SS);
      func_tr[i][1] = func_tr[i][1]*exp(-0.5*kk*kk*SS);
    }
    
    fftw_plan complex2real;
    complex2real = fftw_plan_dft_c2r_1d(nbinN, func_tr, func, FFTW_ESTIMATE);
    fftw_execute(complex2real);
    fftw_destroy_plan(complex2real);

    for (int i=i1; i<i2; i++)
      fx[i-i1] = func[i]/nbinN;
  }

  
  if (file_out!=par::defaultString && file_out!="") {

    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
    
    for (size_t i=0; i<xx.size(); i++)
      fout << xx[i] << "   " << fx[i] << "   " << err[i] << endl;

    fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file_out << endl;
  }

  gsl_histogram_free(histo);
  fftw_cleanup();

}
