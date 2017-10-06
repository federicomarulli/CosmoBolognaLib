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

#include "Sampler.h"
#include "Likelihood.h"

using namespace cosmobl;


// ============================================================================================


cosmobl::statistics::Likelihood::Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const vector<shared_ptr<Parameter>> parameters, const LikelihoodType likelihood_type)
{
  set_data(data);
  set_model(model);
  set_parameters(parameters);
  set_likelihood_function(likelihood_type);
}


// ============================================================================================


cosmobl::statistics::Likelihood::Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const vector<shared_ptr<Parameter>> parameters, const LogLikelihood_function loglikelihood_function)
{
  set_data(data);
  set_model(model);
  set_parameters(parameters);
  set_likelihood_function(loglikelihood_function);
}


// ============================================================================================


cosmobl::statistics::Likelihood::Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const shared_ptr<LikelihoodParameters> parameters, const LikelihoodType likelihood_type)
{
  set_data(data);
  set_model(model);
  m_parameters=parameters;
  set_likelihood_function(likelihood_type);
}


// ============================================================================================


cosmobl::statistics::Likelihood::Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const shared_ptr<LikelihoodParameters> parameters, const LogLikelihood_function loglikelihood_function)
{
  set_data(data);
  set_model(model);
  m_parameters = parameters;
  set_likelihood_function(loglikelihood_function);
}
										
// ============================================================================================


double cosmobl::statistics::Likelihood::prior (vector<double> &pp) const
{
  double prior_value = 1;
  vector<double> param_prior = m_parameters->PriorProbability(pp);

  for (size_t i=0;i<param_prior.size(); i++) {
    prior_value *= param_prior[i];
  }

  return prior_value;
}
										
// ============================================================================================


double cosmobl::statistics::Likelihood::log_prior (vector<double> &pp) const
{
  double prior_value = 0;
  vector<double> param_prior = m_parameters->LogPriorProbability(pp);

  for (size_t i=0;i<param_prior.size(); i++)
    prior_value += param_prior[i];

  return prior_value;
}


// ============================================================================================


double cosmobl::statistics::Likelihood::likelihood (vector<double> &pp) const
{
  return exp(log_likelihood(pp));
}


// ============================================================================================


double cosmobl::statistics::Likelihood::log_likelihood (vector<double> &pp) const
{
  return m_log_likelihood_function(pp, m_likelihood_parameters);
}


// ============================================================================================


double cosmobl::statistics::Likelihood::likelihood_and_priors (vector<double> &pp) const
{
  double pr=prior(pp);
  return ((pr==0) ? 0 : likelihood(pp)*pr);
}


// ============================================================================================


double cosmobl::statistics::Likelihood::log_likelihood_and_priors (vector<double> &pp) const
{
  double pr=prior(pp);
  return ((pr==0) ? par::defaultDouble : log_likelihood(pp)+log(pr));
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_data (const shared_ptr<data::Data> data)
{
  m_data = data;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_model (const shared_ptr<Model> model)
{
  m_model = model;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_parameters (vector<shared_ptr<cosmobl::statistics::Parameter>> parameters)
{
  m_parameters = make_shared<cosmobl::statistics::LikelihoodParameters>(cosmobl::statistics::LikelihoodParameters(parameters));
}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_likelihood_function (const LikelihoodType likelihood_type)
{
  m_likelihood_type = likelihood_type;

  if (m_data->dataType()==data::DataType::_1D_data_ || m_data->dataType()==data::DataType::_1D_data_extra_) {

    switch (m_likelihood_type)
    {
      case (LikelihoodType::_GaussianLikelihood_Error_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_error;
	break; 

      case (LikelihoodType::_GaussianLikelihood_Covariance_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_1D_covariance;
	break;

      default:
	ErrorCBL("Error in set_likelihood_function() of Likelihood.cpp: type of likelihood not recognized or not yet implemented!");
	break;
    }
  }

  else if (m_data->dataType() == data::DataType::_2D_data_ || m_data->dataType() == data::DataType::_2D_data_extra_) {
    switch (m_likelihood_type)
    {
      case (LikelihoodType::_GaussianLikelihood_Error_):
	m_log_likelihood_function = &LogLikelihood_Gaussian_2D_error;
	break; 

      default:
	ErrorCBL("Error in set_likelihood_function() of Likelihood.cpp: type of likelihood not recognized or not yet implemented!");
	break;
    }
  }

  else ErrorCBL("Error in set_likelihood_function() of Likelihood.cpp: data type not recognized or not yet implemented!");

}


// ============================================================================================


void cosmobl::statistics::Likelihood::set_likelihood_function (const LogLikelihood_function loglikelihood_function)
{
  m_likelihood_type = LikelihoodType::_UserDefinedLikelihood_;
  m_log_likelihood_function = loglikelihood_function;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::maximize (vector<double> &guess, const int ntry, const int prior_seed, const bool usePrior, const unsigned int max_iter, const double tol, const double epsilon)
{
  unsigned int npar = m_parameters->nparameters_free();
  if(guess.size() != npar)
    ErrorCBL("Error in maximize of Likelihood.cpp, guess has wrong size!");

  vector<double> guess_complete = m_parameters->full_parameters(guess);

  m_likelihood_parameters = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data,m_model));

  function<double(vector<double> &)> ll;

  if (usePrior) 
    ll = [this](vector<double> & pp) { return -this->log_likelihood_and_priors(pp); };
  else
    ll = [this](vector<double> & pp) { return -this->log_likelihood(pp); };

  double min_loglik = -par::defaultDouble;

  coutCBL << "input values from prior" << endl;

  m_parameters->set_prior_seed(prior_seed);

  for (int i=0; i<ntry; i++) {
    vector<double> pars = m_parameters->prior_sample();
    double temp_loglik = ll(pars);

    if (temp_loglik< min_loglik) {
      min_loglik=temp_loglik;
      guess_complete = pars;
    }
  }

  vector<double> step_size = m_parameters->prior_range(epsilon);


  coutCBL << "Maximizing the likelihood..." << endl;
  guess_complete = cosmobl::gsl::GSL_minimize_nD(ll, guess_complete, step_size, max_iter, tol);
  coutCBL << "Done!" << endl << endl;

  m_parameters->write_bestfit_info(guess_complete);

  if (usePrior)
    coutCBL << "LogLikelihood = " << log_likelihood(guess_complete) << endl << endl;
  else
    coutCBL << "LogLikelihood = " << log_likelihood(guess_complete) << endl << endl;

  guess = guess_complete;
}



// ============================================================================================


void cosmobl::statistics::Likelihood::maximize (vector<double> &guess, const bool usePrior, const unsigned int max_iter, const double tol, const double epsilon)
{
  unsigned int npar = m_parameters->nparameters_free();

  vector<double> guess_complete = m_parameters->full_parameters(guess);

  if (npar==0)
    ErrorCBL("Error in maximize of minimize of Likelihood.cpp, there is no parameter free to vary");

  m_likelihood_parameters = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data, m_model));

  function<double(vector<double> &)> ll;

  if (usePrior) 
    ll = [this](vector<double> & pp) { return -this->log_likelihood_and_priors(pp); };
  else
    ll = [this](vector<double> & pp) { return -this->log_likelihood(pp); };

  vector<double> step_size = m_parameters->prior_range(epsilon);
  ll(guess_complete);

  coutCBL << "Maximizing the likelihood..." << endl;
  guess_complete = cosmobl::gsl::GSL_minimize_nD(ll, guess_complete, step_size, max_iter, tol);
  coutCBL << "Done!" << endl << endl;

  m_parameters->write_bestfit_info(guess_complete);

  if (usePrior)
    coutCBL << "LogLikelihood = " << log_likelihood(guess_complete) << endl << endl;
  else
    coutCBL << "LogLikelihood = " << log_likelihood(guess_complete) << endl << endl;

  guess = guess_complete;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::sample_stretch_move (const int seed, const double aa)
{
  if (m_likelihood_type==statistics::LikelihoodType::_Likelihood_NotSet_)
    ErrorCBL("Error in cosmobl::statistics::Likelihood::sample_stretch_move of Likelihood.h: the Likelihood function is not set!");

  coutCBL << "Sampling the likelihood..." << endl;

  vector<vector<double> > start = m_parameters->chain_values(0);

  m_likelihood_parameters = make_shared<cosmobl::statistics::STR_likelihood_parameters>(cosmobl::statistics::STR_likelihood_parameters(m_data, m_model));

  auto likelihood_function = [this] (vector<double> &pp) { return this->log_likelihood_and_priors(pp); }; 

  cosmobl::statistics::Sampler sampler(m_parameters->nparameters(), m_parameters->nparameters_free(), likelihood_function); 
  sampler.sample_stretch_move (m_parameters->chain_size(), m_parameters->nwalkers(), start, seed, aa);
  
  vector<vector<double>> chain_values;
  sampler.get_chain_function_acceptance(chain_values, m_likelihood_values, m_acceptance);

  m_parameters->set_chains(chain_values, m_parameters->nwalkers());

  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::sample_stretch_move_parallel (const int seed, const double aa)
{
  if (m_likelihood_type==statistics::LikelihoodType::_Likelihood_NotSet_)
    ErrorCBL("Error in cosmobl::statistics::Likelihood::sample_stretch_move_parallel of Likelihood.h: the Likelihood function is not set!");

  coutCBL << "Sampling the likelihood..." << endl;
  
  vector<vector<double> > start = m_parameters->chain_values(0);

  m_likelihood_parameters = make_shared<cosmobl::statistics::STR_likelihood_parameters>(cosmobl::statistics::STR_likelihood_parameters(m_data, m_model));

  auto likelihood_function = [this] (vector<double> &pp) { return this->log_likelihood_and_priors(pp); }; 

  cosmobl::statistics::Sampler sampler(m_parameters->nparameters(), m_parameters->nparameters_free(), likelihood_function); 
  
  sampler.sample_stretch_move_parallel(m_parameters->chain_size(), m_parameters->nwalkers(), start, seed, aa);
  
  vector<vector<double>> chain_values;
  sampler.get_chain_function_acceptance(chain_values, m_likelihood_values, m_acceptance);

  m_parameters->set_chains(chain_values, m_parameters->nwalkers());

  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::write_chain (const string output_dir, const string output_file, const int start, const int thin)
{
  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {}
  
  string file = output_dir+output_file; 
  ofstream fout(file.c_str()); checkIO(fout, file);
  fout.precision(10);
  
  int nn = 0;

  for (int j=start; j<m_parameters->chain_size(); j+=thin) {
    for (int i=0; i<m_parameters->nwalkers(); i++) {
      int index = i+j*m_parameters->nwalkers();
      fout << nn << " ";
      vector<double> pp;
      for (int k = 0; k<m_parameters->nparameters(); k++) {
	pp.push_back(m_parameters->chain_value(j, i, k));
	fout << " " << m_parameters->chain_value(j, i, k) << " ";
      }

      double pr=prior(pp);
      fout << m_likelihood_values[index] -log(pr) << " " << pr << endl;
      nn ++;
    }
  }

  fout.clear(); fout.close();
  coutCBL << "I wrote the file: " << file << endl;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::read_chain (const string input_dir, const string input_file, const int nwalkers, const int skip_header)
{

  string file = input_dir+input_file;
  coutCBL << "Reading the chain file " << file << endl;

  ifstream fin(file.c_str()); checkIO(fin, file);

  string line;
  for(int i=0; i<skip_header; i++)
    getline(fin, line);

  vector<vector<double>> chain_values;
  vector<double> likelihood_values;

  while(getline(fin, line))
  {
    stringstream ss(line);
    double NUM;
    vector<double> ll, params;

    while(ss>>NUM) ll.push_back(NUM);
    for(size_t i=1;i<ll.size()-1; i++){
      params.push_back(ll[i]);
    }

    likelihood_values.push_back(ll[ll.size()-1]);
    chain_values.push_back(params);
  }
  int chain_size = chain_values.size()/nwalkers;

  checkDim(chain_values, nwalkers*chain_size, m_parameters->nparameters(), "chain_from_file");
  checkDim(likelihood_values, nwalkers*chain_size, "likelihood_from_file");

  chain_values = cosmobl::transpose(chain_values);

  m_parameters->set_chains(chain_values, nwalkers);
  m_likelihood_values = likelihood_values;

  fin.clear(); fin.close();
  
  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::initialize_chains (const int chain_size, const int nwalkers, const int seed)
{
  m_parameters->set_chains(chain_size, nwalkers);

  m_parameters->initialize_chains_from_prior(seed);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::initialize_chains (const int chain_size, const int nwalkers, const int seed, const double radius, const int ntry, const int prior_seed, const unsigned int max_iter, const double tol, const double epsilon)
{
  m_parameters->set_chains(chain_size, nwalkers);

  vector<double> guess;
  maximize(guess, ntry, prior_seed, true, max_iter, tol, epsilon);

  m_parameters->initialize_chains_around_bestfit_values(radius, seed);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::initialize_chains (const int chain_size, const int nwalkers, const int seed, const double radius, vector<double> &guess, const unsigned int max_iter, const double tol, const double epsilon)
{
  vector<double> guess_complete = m_parameters->full_parameters(guess);

  m_parameters->set_chains(chain_size, nwalkers);

  maximize(guess_complete, true, max_iter, tol, epsilon);

  m_parameters->initialize_chains_around_bestfit_values(radius, seed);

  guess = guess_complete;
}


// ============================================================================================


void cosmobl::statistics::Likelihood::initialize_chains (const int chain_size, const int nwalkers, const int seed, vector<double> &values, const double radius)
{
  m_parameters->set_chains(chain_size, nwalkers);

  m_parameters->initialize_chains_around_values(values, radius, seed);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::initialize_chains (const int chain_size, const vector<vector<double>> chain_values)
{
  const int nwalkers = chain_values[0].size();
  m_parameters->set_chains(chain_size, nwalkers);

  for(int pp=0; pp<m_parameters->nparameters(); pp++)
    m_parameters->set_chain_values(0, pp, chain_values[pp]);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::initialize_chains (const int chain_size, const int nwalkers, const string input_dir, const string input_file)
{
  string last_step_file = input_dir+input_file+"_LastStep";
  string get_last_step = "tail -n "+conv(nwalkers, par::fINT)+" "+input_dir+input_file+" > "+last_step_file;
  if (system(get_last_step.c_str())) {}

  ifstream fin(last_step_file);
  string line;

  vector<vector<double>> chain_values;

  while(getline(fin, line))
  {
    stringstream ss(line);
    double NUM;
    vector<double> ll, params;

    while(ss>>NUM) ll.push_back(NUM);
    for(size_t i=1;i<ll.size()-2; i++){
      params.push_back(ll[i]);
    }

    chain_values.push_back(params);
  }
  fin.clear(); fin.close();

  string rm_last_step = "rm -r "+last_step_file;
  if (system(rm_last_step.c_str())) {}

  checkDim(chain_values, nwalkers, m_parameters->nparameters(), "chain_from_LastStep_file");
  chain_values = cosmobl::transpose(chain_values);

  initialize_chains (chain_size, chain_values);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::show_results (const int start, const int thin, const int seed)
{
  m_parameters->show_results(start, thin, seed);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::write_results (const string dir, const string file, const int start, const int thin, const int seed)
{
  m_parameters->write_results(dir, file, start, thin, seed);
  write_chain(dir, file+"_chain.dat", start, thin);
}


// ============================================================================================


void cosmobl::statistics::Likelihood::sample_tabulated_likelihood (const int nstep_p1, const int nstep_p2, const string interpolation_method, const int nchains, const int chain_size, const int seed, const bool do_write_chain, const string output_dir, const string output_file)
{
  (void)nstep_p1; (void)nstep_p2; (void)interpolation_method; (void)nchains; (void)chain_size; 
  (void)seed; (void)do_write_chain; (void)output_dir; (void)output_file;

  ErrorCBL("Error in sample_tabulated_likelihood, work in progress");
}
/*
{
  if (m_model->npar()!=2)
    ErrorCBL("Error in sample tabulated_likelihood of Likelihood.cpp, it only works with 2 parameters");

  coutCBL << "Sampling tabulated likelihood..." << endl;
  
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

      double lnprob = min(1.,pow(gen_z,(m_model->npar_free()-1))*exp(proposed_loglikelihood+proposed_prior-log_likelihood[i]-log_prior[i]));

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
*/
