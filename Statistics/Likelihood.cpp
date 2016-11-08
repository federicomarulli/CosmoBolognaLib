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

/** @file Statistics/Likelihood.cpp
 *
 *  @brief Methods of the class Likelihood 
 *
 *  This file contains the implementation of the methods of the class
 *  Likelihood, used for bayesian analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Likelihood.h"

using namespace cosmobl;


// ============================================================================================


cosmobl::statistics::Likelihood::Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const LikelihoodType likelihood_type)
{
  set_likelihood_type(likelihood_type, 0);
  set_data(data);
  set_model(model);
  set_likelihood_function();
}


// ============================================================================================


cosmobl::statistics::Likelihood::Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const LikelihoodType likelihood_type, const LogLikelihood_function loglikelihood_function, const bool cov)
{
  set_likelihood_type(likelihood_type, cov);
  set_data(data);
  set_model(model);
  set_likelihood_function(loglikelihood_function, cov);
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_model (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<statistics::STR_likelihood_parameters> pp  = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);
  pp->model->set_parameter_values(model_parameters);

  vector<double> computed_model(pp->data->ndata(), 0);

#pragma omp parallel num_threads(omp_get_max_threads())
  {

#pragma omp for
    for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
      computed_model[i]=pp->model->operator()(pp->data->xx(i));
  }

  double LogLikelihood=0;
  for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
    LogLikelihood += pow((pp->data->fx(i)-computed_model[i])/computed_model[i],2);

  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_error (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);
  pp->model->set_parameter_values(model_parameters);

  vector<double> computed_model(pp->data->ndata(), 0);

#pragma omp parallel num_threads(omp_get_max_threads())
  {

#pragma omp for
    for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
      computed_model[i] = pp->model->operator()(pp->data->xx(i));
  }

  double LogLikelihood=0;
  for (int i=pp->data->x_down(); i<pp->data->x_up(); i++) 
    LogLikelihood += pow((pp->data->fx(i)-computed_model[i])/pp->data->error_fx(i), 2);
  
  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_covariance (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);
  pp->model->set_parameter_values(model_parameters);
  vector<double> computed_model(pp->data->ndata(),0);

#pragma omp parallel num_threads(omp_get_max_threads())
  {

#pragma omp for
    for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
      computed_model[i]=pp->model->operator()(pp->data->xx(i));
  }


  vector<double> diff(pp->data->ndata(),0);
  for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
    diff[i]=pp->data->fx(i)-computed_model[i];

  double LogLikelihood = 0;
  
  for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
    for (int j=pp->data->x_down(); j<pp->data->x_up(); j++)
      LogLikelihood += diff[i]*pp->data->inverse_covariance(i,j)*diff[j];    

  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_2D_model (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);
  pp->model->set_parameter_values(model_parameters);

  vector<vector<double>> computed_model(pp->data->ndata(),vector<double>(pp->data->ndata(),0));

  int imin = pp->data->x_down();
  int imax = pp->data->x_up();

  int jmin = pp->data->y_down();
  int jmax = pp->data->y_up();

  int deltaj = jmax-jmin;

  int ijmax = (imax-imin)*(jmax-jmin);


#pragma omp parallel num_threads(omp_get_max_threads())

  {

#pragma omp for
    for (int ij=0; ij<ijmax; ++ij)
    {
      int i = ij/deltaj+imin;
      int j = ij%deltaj+jmin;
      computed_model[i][j] = pp->model->operator()(pp->data->xx(i),pp->data->yy(j));
    }
  }

  double LogLikelihood = 0;
  for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
    for (int j=pp->data->y_down(); j<pp->data->y_up(); j++)
      LogLikelihood += pow((pp->data->fxy(i,j)-computed_model[i][j])/computed_model[i][j],2);

 
  return -0.5*LogLikelihood;
}

// ============================================================================================


double statistics::LogLikelihood_Gaussian_2D_error (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);
  pp->model->set_parameter_values(model_parameters);

  vector<vector<double>> computed_model(pp->data->ndata(),vector<double>(pp->data->ndata(),0));

  int imin = pp->data->x_down();
  int imax = pp->data->x_up();

  int jmin = pp->data->y_down();
  int jmax = pp->data->y_up();

  int deltaj = jmax-jmin;

  int ijmax = (imax-imin)*(jmax-jmin);


#pragma omp parallel num_threads(omp_get_max_threads())
  {

#pragma omp for schedule(static, 2)
    for (int ij=0; ij<ijmax; ++ij)
    {
      int i = ij/deltaj+imin;
      int j = ij%deltaj+jmin;
      computed_model[i][j] = pp->model->operator()(pp->data->xx(i),pp->data->yy(j));
    }
  }


  double LogLikelihood = 0;
  for (int i=pp->data->x_down(); i<pp->data->x_up(); i++)
    for (int j=pp->data->y_down(); j<pp->data->y_up(); j++)
      LogLikelihood += pow((pp->data->fxy(i,j)-computed_model[i][j])/computed_model[i][j],2);

  return -0.5*LogLikelihood;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_likelihood_type (const LikelihoodType likelihood_type, const bool cov)
{
  m_likelihood_type = likelihood_type;

  switch(m_likelihood_type)
    {
    case(LikelihoodType::_GaussianLikelihood_Model_):
      m_cov = 0;
      break;

    case(LikelihoodType::_GaussianLikelihood_Error_):
      m_cov = 0;
      break; 

    case(LikelihoodType::_GaussianLikelihood_Covariance_):
      m_cov = 1;
      break;

    case(LikelihoodType::_UserDefinedLikelihood_):
      m_cov = cov;
      break;

    default:
      ErrorCBL("Error in set_likelihood_type of Likelihood.cpp, Type of likelihood not recognized");
      break;
    }
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_data (const shared_ptr<data::Data> data)
{
  if (m_likelihood_type == cosmobl::statistics::LikelihoodType::_Likelihood_NotSet_)
    ErrorCBL("Error in set_data of Likelihood.cpp, Likelihood type not setted, please set it first before setting the dataset");

  m_data = data;

  if (m_cov)
    data->invert_covariance();
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_model (const shared_ptr<Model> model)
{
  m_model = model;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_likelihood_function()
{
  if (m_data->dataType() == data::DataType::_1D_data_) {
    switch(m_likelihood_type)
      {
      case(LikelihoodType::_GaussianLikelihood_Model_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_model;
	break;

      case(LikelihoodType::_GaussianLikelihood_Error_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_error;
	break; 

      case(LikelihoodType::_GaussianLikelihood_Covariance_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_covariance;
	break;

      default:
	ErrorCBL("Error in set_likelihood_type of Likelihood.cpp, Type of likelihood not recognized or not yet implemented");
	break;
      }
  }
  else if (m_data->dataType() == data::DataType::_2D_data_) {
    switch(m_likelihood_type)
      {
      case(LikelihoodType::_GaussianLikelihood_Model_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_2D_model;
	break;

      case(LikelihoodType::_GaussianLikelihood_Error_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_2D_error;
	break; 

      default:
	ErrorCBL("Error in set_likelihood_type of Likelihood.cpp, Type of likelihood not recognized or not yet implemented");
	break;
      }
  }

}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_likelihood_function(const LogLikelihood_function loglikelihood_function, const bool cov)
{
  if (m_likelihood_type != LikelihoodType::_UserDefinedLikelihood_)
    ErrorCBL("Error in set_likelihood_function of Likelihood.cpp, This function only works with LikelihoodType::_UserDefinedLikelihood_");

  m_log_likelihood_function = loglikelihood_function;
  m_cov = cov;
}

// ======================================================================================


double cosmobl::statistics::Likelihood::LogPriorProbability ()
{
  double log_prior = 0;

  for (unsigned int i=0; i<m_model->npar(); i++)
    log_prior += m_model->parameter(i)->isFixed() ? 0 : Log(m_model->parameter(i)->PriorProbability());

  return log_prior;
}

// ======================================================================================


double cosmobl::statistics::Likelihood::LogPriorProbability (const vector<double> parameter_values)
{
  double log_prior = 0;

  if (parameter_values.size() != m_model->npar())
    ErrorCBL("Error in LogPriorProbability of Likelihood.cpp, wrong number of parameters");

  for (unsigned int i=0; i<m_model->npar(); i++)
    {
      log_prior += m_model->parameter(i)->isFixed() ? 0 : Log(m_model->parameter(i)->PriorProbability(parameter_values[i]));
    }

  return log_prior;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::minimize_LogLikelihood (const bool usePrior, const unsigned int max_iter, const double tol)
{
  if (m_likelihood_type != LikelihoodType::_UserDefinedLikelihood_) 
    set_likelihood_function();

  unsigned int npar = m_model->npar();
  if (npar==0)
    ErrorCBL("Error in minimize of minimize of Likelihood.cpp, there is no parameter to vary");

  auto fixed_parameters = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data,m_model));

  auto likelihood_prior = [&](vector<double> xx, shared_ptr<void> pp) { return -m_log_likelihood_function(xx, pp) + ((usePrior) ? -LogPriorProbability(xx) : 0); };

  vector<double> step_size(npar, 1);
  int nn = 0;
  for (unsigned int i=0; i<npar; i++) {
    if (!m_model->parameter(i)->isFixed()) { 
      if ( m_model->parameter(i)->interval_size()<-par::defaultDouble) 
	step_size[nn]=0.05*m_model->parameter(i)->interval_size();
      nn ++;
    }
  }

  // Starting the minimization...
  
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  x = gsl_vector_alloc (npar);
  ss = gsl_vector_alloc (npar);
  for (size_t i=0; i<npar; i++) {
    gsl_vector_set(x, i, m_model->parameter_values()[i]);
    gsl_vector_set(ss, i, step_size[i]);
  }

  typedef function<double(vector<double>)> likelihood;
  likelihood func = bind(&Likelihood::operator(), this, std::placeholders::_1);

  minex_func.n = npar;
  minex_func.f = [](const gsl_vector *gsl_x, void *p) {
    auto f = (likelihood *)p;
    vector<double> xx;
    for (unsigned int i=0; i<gsl_x->size; i++)
      xx.push_back(gsl_vector_get(gsl_x,i));

    return -f->operator()(xx); //perator()(xx); //likelihood(xx);
  };
  minex_func.params = &func;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do {
    iter++;

    status = gsl_multimin_fminimizer_iterate(s);

    if (status) break;

    size = gsl_multimin_fminimizer_size (s);
    
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
      printf ("-----> converged to minimum \n");
  }
  
  while (status == GSL_CONTINUE && iter <max_iter);

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  for (unsigned int i=0; i<m_model->npar(); i++) {
    if (m_model->parameter(i)->isFixed())
      coutCBL << par::col_green << m_model->parameter(i)->name() << par::col_default << ": fixed value = " << m_model->parameter(i)->value() << endl;
    else
      coutCBL << par::col_green << m_model->parameter(i)->name() << par::col_default << ": best-fit value from the chi2 minimization =" << m_model->parameter(i)->value() << endl;
  }


  coutCBL << "-2*LogLikelihood = " << likelihood_prior(m_model->parameter_values(), fixed_parameters) << endl << endl;

}


// ============================================================================================


double cosmobl::statistics::Likelihood::sample_stretch_move (const int nchains, const int chain_size, const int seed, const bool do_write_chain, const string output_dir, const string output_file)
{
  coutCBL << "Sampling the likelihood" << endl;
  
  // set the likelihood type, according to data and parameters, if not yet setted
  if (m_likelihood_type != LikelihoodType::_UserDefinedLikelihood_) 
    set_likelihood_function();

  // set seed for priors
  random::UniformRandomNumbers prior_seeds(0, 23412432, seed);

  for (unsigned int i=0; i<m_model->npar(); i++)
    m_model->parameter(i)->set_prior_seed(int(prior_seeds()));

  // initialize chains
  m_nchains=nchains;
  m_chain_size=chain_size;
  vector<double> chains_index = linear_bin_vector(m_nchains, 0., m_nchains-1.);
  vector<double> chains_weights(m_nchains,1);
  random::DiscreteRandomNumbers chains(chains_index, chains_weights, int(prior_seeds()), 0, m_nchains-1);
  
  for (size_t i=0; i<m_model->npar(); i++) {
    m_model->parameter(i)->set_chains(m_nchains, m_chain_size);
    m_model->parameter(i)->set_chains_values_from_prior(0);
  }
 
  // initialize LogLikelihood
  auto fixed_parameters = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data, m_model));

  vector<double> log_likelihood(m_nchains, 0);
  vector<double> log_prior(m_nchains, 0);
  
  for (int i=0; i<m_nchains; i++) {
    log_likelihood[i] = m_log_likelihood_function(m_model->parameter_values_from_chain(i,0), fixed_parameters);
    log_prior[i] = LogPriorProbability(m_model->parameter_values_from_chain(i,0));
  }

  
  // stretch-move

  // sample the distribuion function g(Z), to be moved somewhere else
  
  double gzpar = 2; // User defined ?

  double zmin = 1./gzpar;
  double zmax = gzpar;
  vector<double> zz = linear_bin_vector(1000,zmin,zmax), gzz;
  for (auto &&zzz : zz)
    gzz.push_back(1./sqrt(zzz));
 
  random::UniformRandomNumbers MH_random(0., 1., int(prior_seeds()));
  random::DistributionRandomNumbers GZ(zz, gzz, "Spline", int(prior_seeds()));

  int dp = cout.precision();
  
  for (int n=1; n<m_chain_size; n++) {
    for (int i=0; i<m_nchains; i++)
    {
      int kk = i;
      while (kk==i)
	kk = int(chains());
        
      vector<double> parameters_i;
      vector<double> parameters_k;

      bool isIncluded = 0;
      double gen_z;
      
      while (!isIncluded) {
	gen_z = GZ();
	parameters_i = m_model->parameter_values_from_chain(i, n-1);
	parameters_k = m_model->parameter_values_from_chain(kk, n-1);
	vector<bool> included(m_model->npar(), true);
	for (unsigned int p=0; p<m_model->npar(); p++) {
	  if (!m_model->parameter(p)->isFixed()) {
	    parameters_i[p] = parameters_k[p] + gen_z*(parameters_i[p]-parameters_k[p]);
	    included[p] =  m_model->parameter(p)->prior()->isIncluded(parameters_i[p]);
	  }
	}
	isIncluded = accumulate(included.begin(), included.end(), 1, std::multiplies<bool>());
      }

      m_model->set_parameter_values(parameters_i);

      double proposed_loglikelihood = m_log_likelihood_function(parameters_i, fixed_parameters);
      double proposed_prior = LogPriorProbability(parameters_i);

      double lnprob = min(1.,pow(gen_z,(m_model->npar_eff()-1))*exp(proposed_loglikelihood-proposed_prior-log_likelihood[i]+log_prior[i]));

      if (MH_random() <lnprob) {
	log_prior[i] = proposed_prior;
	log_likelihood[i] = proposed_loglikelihood;
	for (unsigned int p=0; p<m_model->npar(); p++)
	  m_model->parameter(p)->set_chain_value (i, n, parameters_i[p]);
      }
      else {
	parameters_i = m_model->parameter_values_from_chain(i, n-1);
	for (unsigned int p=0; p<m_model->npar(); p++)
	  m_model->parameter(p)->set_chain_value (i, n, parameters_i[p]);
      }

      double progress = double((n+1)*m_nchains)/(m_nchains*m_chain_size)*100;
      coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();
    }

    if (do_write_chain)
      m_write_chain (output_dir, output_file, 0, n, 1);
  }

  coutCBL << endl;
  cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
  
  for (size_t i = 0; i<m_model->npar(); i++)
    m_model->parameter(i)->chains_convergence(m_chain_size, m_chain_size*0.5);

  return 0;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::write_chain (const string output_dir, const string output_file, const double start, const double stop, const int thin)
{
  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {}
  
  string file = output_dir+output_file; 
  ofstream fout(file.c_str()); checkIO(fout, file);
  
  double new_start = int(m_chain_size*start);
  double new_stop = int(m_chain_size*stop);
  int nn = 0;

  if (m_model->npar()) {
    vector< shared_ptr<Parameter> > pp = m_model->parameters();
    for (int i=0; i<m_nchains; i++) {
      for (int j=new_start; j<new_stop; j+=thin) {
	fout << nn << " ";
	vector<double> pars;
	for (size_t k = 0; k<pp.size(); k++) {
	  pars.push_back(pp[k]->chain(i)->chain_value(j));
	  fout << " " << pp[k]->chain(i)->chain_value(j) << " ";
	}
	fout << endl;
	nn ++;
      }
    }
  }

  fout.clear(); fout.close();
  coutCBL << "I wrote the file: " << file << endl;
}


// ============================================================================================


double cosmobl::statistics::Likelihood::sample_tabulated_likelihood (const int nstep_p1, const int nstep_p2, const string interpolation_method, const int nchains, const int chain_size, const int seed, const bool do_write_chain, const string output_dir, const string output_file)
{
  if (m_model->npar()!=2)
    ErrorCBL("Error in sample tabulated_likelihood of Likelihood.cpp, it only works with 2 parameters");

  coutCBL << "Sampling tabulated likelihood" << endl;
  
  // set the likelihood type, according to data and parameters, if not yet setted
  if (m_likelihood_type != LikelihoodType::_UserDefinedLikelihood_) 
    set_likelihood_function();

  vector<double> parameter1(nstep_p1);
  vector<double> parameter2(nstep_p2);
  vector<vector<double> > tabulated_likelihood(nstep_p1, vector<double>(nstep_p2,0));

  for (int i=0; i<nstep_p1; i++) {
    parameter1[i] = m_model->parameter(0)->prior()->xmin()+(i)*(m_model->parameter(0)->prior()->xmax()-m_model->parameter(0)->prior()->xmin())/(nstep_p1-1);
    for (int j=0; j<nstep_p2; j++) {
      parameter2[j] = m_model->parameter(1)->prior()->xmin()+(j)*(m_model->parameter(1)->prior()->xmax()-m_model->parameter(1)->prior()->xmin())/(nstep_p2-1);
      tabulated_likelihood[i][j] = this->operator()({parameter1[i], parameter2[j]});
    }
  }


  // Set seed for priors
  random::UniformRandomNumbers prior_seeds(0, 23412432, seed);

  for (unsigned int i=0; i< m_model->npar(); i++)
    m_model->parameter(i)->set_prior_seed(int(prior_seeds()));
  

  // Initialize chains
  m_nchains=nchains;
  m_chain_size=chain_size;
  vector<double> chains_index = linear_bin_vector(m_nchains, 0., m_nchains-1.);
  vector<double> chains_weights(m_nchains,1);
  random::DiscreteRandomNumbers chains(chains_index, chains_weights, int(prior_seeds()), 0, m_nchains-1);

  for (size_t i=0; i<m_model->npar(); i++) {
    m_model->parameter(i)->set_chains(m_nchains, m_chain_size);
    m_model->parameter(i)->set_chains_values_from_prior(0);
  }


  // initialize LogLikelihood;
  auto fixed_parameters = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data,m_model));

  vector<double> log_likelihood(m_nchains, 0);
  vector<double> log_prior(m_nchains, 0);
  
  for (int i=0; i<m_nchains; i++) {
    vector<double> params = m_model->parameter_values_from_chain(i,0);
    log_likelihood[i] = interpolated_2D(params[0], params[1], parameter1, parameter2, tabulated_likelihood, interpolation_method);
    log_prior[i] = LogPriorProbability(m_model->parameter_values_from_chain(i,0));
  }


  // ----- stretch-move -----

  // sample the distribuion function g(Z), to be moved somewhere else
  
  double gzpar = 2; // (user defined ?)

  double zmin = 1./gzpar;
  double zmax = gzpar;
  vector<double> zz = linear_bin_vector(1000,zmin,zmax), gzz;
  for (auto &&zzz : zz)
    gzz.push_back(1./sqrt(zzz));

  random::UniformRandomNumbers MH_random(0., 1., int(prior_seeds()));
  random::DistributionRandomNumbers GZ(zz, gzz, "Spline", int(prior_seeds()));

  for (int n=1; n<m_chain_size; n++) {
    for (int i=0; i< m_nchains; i++) {
      int kk = i;
      while(kk==i)
	kk = int(chains());
        
      vector<double> parameters_i;
      vector<double> parameters_k;

      bool isIncluded = 0;
      double gen_z;

      while (!isIncluded) {
	gen_z = GZ();
	parameters_i = m_model->parameter_values_from_chain(i, n-1);
	parameters_k = m_model->parameter_values_from_chain(kk, n-1);
	vector<bool> included(m_model->npar());
	for (unsigned int p=0; p<m_model->npar(); p++) {
	  parameters_i[p] = parameters_k[p] + gen_z*(parameters_i[p]-parameters_k[p]);
	  included[p] =  m_model->parameter(p)->prior()->isIncluded(parameters_i[p]);
	}
	isIncluded = accumulate(included.begin(), included.end(), 1, std::multiplies<bool>());
      }

      m_model->set_parameter_values(parameters_i);


      double proposed_loglikelihood = interpolated_2D(parameters_i[0], parameters_i[1], parameter1, parameter2, tabulated_likelihood, interpolation_method);
//m_log_likelihood_function(parameters_i, fixed_parameters);
      double proposed_prior = LogPriorProbability(parameters_i);

      double lnprob = min(1.,pow(gen_z,(m_model->npar_eff()-1))*exp(proposed_loglikelihood+proposed_prior-log_likelihood[i]-log_prior[i]));

      if (MH_random() < lnprob) {
	log_prior[i] = proposed_prior;
	log_likelihood[i] = proposed_loglikelihood;
	for (unsigned int p=0; p<m_model->npar(); p++)
	  m_model->parameter(p)->set_chain_value (i, n, parameters_i[p]);
      }
      else {
	parameters_i = m_model->parameter_values_from_chain(i, n-1);
	for (unsigned int p=0; p<m_model->npar(); p++)
	  m_model->parameter(p)->set_chain_value (i, n, parameters_i[p]);
      }

      double progress = double((n+1)*m_nchains)/(m_nchains*m_chain_size)*100;
      coutCBL << setprecision(2) << setiosflags(ios::fixed) << setw(8) << progress << "% \r"; cout.flush();
    }

    coutCBL << endl;
    
    if (do_write_chain)
      m_write_chain(output_dir, output_file, 0, n, 1);
  }
  
  for (size_t i=0; i<m_model->npar(); i++)
    m_model->parameter(i)->chains_convergence(m_chain_size, m_chain_size*0.5);

  return 0;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::m_write_chain (const string output_dir, const string output_file, const int start, const int stop, const int thin)
{
  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {}
  
  string file = output_dir+output_file; 
  ofstream fout(file.c_str()); checkIO(fout, file);
  
  int new_start = (start<0) ?  0 : start;
  int new_stop = (start<0) ?  m_chain_size : stop;

  int nn = 0;

  if (m_model->npar()) {
    vector< shared_ptr<Parameter> > pp = m_model->parameters();
    for (int j=new_start; j<new_stop; j+=thin) {
      for (int i=0; i<m_nchains; i++) {
	fout << nn << " ";
	vector<double> pars;
	for (size_t k = 0;k<pp.size();k++) {
	  pars.push_back(pp[k]->chain(i)->chain_value(j));
	  fout << " " << pp[k]->chain(i)->chain_value(j) << " ";
	}
	fout << i << " " << j << endl;
	nn ++;
      }
    }
  }

  fout.clear(); fout.close();
  coutCBL << "I wrote the file: " << file << endl;
}
