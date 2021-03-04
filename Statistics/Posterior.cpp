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

/** @file Statistics/Posterior.cpp
 *
 *  @brief Methods of the class Posterior 
 *
 *  This file contains the implementation of the methods of the class
 *  Posterior, used for bayesian analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "FITSwrapper.h"
#include "PosteriorParameters.h"
#include "Sampler.h"
#include "Posterior.h"
#include "ChainMesh.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::statistics::Posterior::m_set_seed (const int seed)
{
  m_seed = seed;
  m_seed_generator = make_shared<random::UniformRandomNumbers_Int>(random::UniformRandomNumbers_Int(0, std::numeric_limits<int>::max(), m_seed));
  m_model_parameters->set_prior_distribution_seed(m_seed_generator);
}


// ============================================================================================


cbl::statistics::Posterior::Posterior (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const int seed)
{
  set(prior_distributions, data, model, likelihood_type, x_index, w_index, seed);
}


// ============================================================================================


cbl::statistics::Posterior::Posterior (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const Likelihood &likelihood, const int seed) : Likelihood(likelihood)
{
  m_model_parameters = make_shared<PosteriorParameters>(PosteriorParameters(m_model->parameters()->nparameters(), prior_distributions, m_model->parameters()->type(), m_model->parameters()->name()));

  m_model->set_parameters(m_model_parameters);

  m_prior = m_model_parameters->prior();

  m_likelihood_inputs = make_shared<STR_likelihood_inputs>(STR_likelihood_inputs(m_data, m_model, m_x_index, m_w_index));

  m_set_seed(seed);
}


// ============================================================================================


cbl::statistics::Posterior::Posterior (const std::string file_name, const std::string path_file, const std::vector<int> usecols, const std::vector<std::string> parameter_names, const int skip_header)
{
  std::vector<std::vector<double>> chain = cbl::read_file(file_name, path_file, usecols, skip_header);

  m_Nparameters = usecols.size()-1;
  m_log_posterior = chain[m_Nparameters];
  m_parameters.resize(m_Nparameters);

  for (int N=0; N<m_Nparameters; N++)
    m_parameters[N] = chain[N];

  m_seed = -1;

  m_model_parameters = make_shared<ModelParameters> (ModelParameters(m_Nparameters, parameter_names));
}


// ============================================================================================



double cbl::statistics::Posterior::operator () (std::vector<double> &pp) const
{
  pp = m_model_parameters->full_parameter(pp);
  double prior = m_prior->operator()(pp);

  if (prior>0)
    return (m_use_grid) ? m_likelihood_function_grid(pp, m_likelihood_inputs)*prior : m_likelihood_function(pp, m_likelihood_inputs)*prior;

  return 0.;
}


// ============================================================================================


vector<double> cbl::statistics::Posterior::weight (const int start, const int thin) const
{
  vector<double> ww;
  const int chain_size = m_model_parameters->chain_size();
  const int n_walkers = m_model_parameters->chain_nwalkers();

  for (int j=start; j<chain_size; j+=thin)
    for (int i=0; i<n_walkers; i++)
      ww.push_back(m_weight[j*n_walkers+i]);

  return ww;
}


// ============================================================================================


double cbl::statistics::Posterior::log (std::vector<double> &pp) const
{
  pp = m_model_parameters->full_parameter(pp);
  const double logprior = m_prior->log(pp);
 
  double val;
  if (logprior>par::defaultDouble)
    val = (m_use_grid) ? m_log_likelihood_function_grid(pp, m_likelihood_inputs)+logprior : m_log_likelihood_function(pp, m_likelihood_inputs)+logprior;
  else
    val = par::defaultDouble;

  return val;
}


// ============================================================================================


double cbl::statistics::Posterior::chi2 (const std::vector<double> parameter) const
{
  vector<double> pp = (parameter.size()>0) ? parameter : m_model_parameters->bestfit_value();
  
  return (m_use_grid) ? -2.*m_log_likelihood_function_grid(pp, m_likelihood_inputs) : -2.*m_log_likelihood_function(pp, m_likelihood_inputs);
}


// ============================================================================================


void cbl::statistics::Posterior::set_model (const std::shared_ptr<Model> model, const std::shared_ptr<ModelParameters> model_parameters)
{
  switch (model->dimension()) {
  case Dim::_1D_:
    m_model = make_shared<Model1D>(*static_pointer_cast<Model1D>(model));
    break;
  case Dim::_2D_:
    m_model = make_shared<Model2D>(*static_pointer_cast<Model2D>(model));
    break;
  default:
    ErrorCBL("the model dimension must be Dim::_1D_ or Dim::_2D_ !", "set_model", "Posterior.cpp");
  }

  if (model_parameters!=NULL)
    m_model->set_parameters(model_parameters);
  else
    m_model->set_parameters(m_model_parameters);

}


// ============================================================================================


void cbl::statistics::Posterior::set (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const int seed)
{
  set_data(data);
  set_model(model);

  m_model_parameters = make_shared<PosteriorParameters>(PosteriorParameters(m_model->parameters()->nparameters(), prior_distributions, m_model->parameters()->type(), m_model->parameters()->name()));

  m_model->set_parameters(m_model_parameters);
  m_model_parameters->set_prior_distribution(prior_distributions);

  m_prior = m_model_parameters->prior();

  set_function(likelihood_type, x_index, w_index);
  m_likelihood_inputs = make_shared<STR_likelihood_inputs>(STR_likelihood_inputs(m_data, m_model, m_x_index, m_w_index));

  m_set_seed(seed);
}


// ============================================================================================


void cbl::statistics::Posterior::maximize (const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon)
{
  vector<double> starting_par;

  const unsigned int npar_free = m_model->parameters()->nparameters_free();
  const unsigned int npar = m_model->parameters()->nparameters();
  
  if (start.size()==npar_free) 
    starting_par = start;
  else if (start.size()==npar)
    for (size_t i=0; i<npar_free; i++)
      starting_par.push_back(start[m_model->parameters()->free_parameter()[i]]);
  else
    ErrorCBL("check your inputs: start.size()="+conv(start.size(), par::fINT)+" must be equal to either npar_free="+conv(npar_free, par::fINT)+" or npar="+conv(npar, par::fINT)+"!", "maximize", "Posterior.cpp");
  
  function<double(vector<double> &)> post = [this] (vector<double> &pp) { return -this->log(pp); };


  // extra check on epsilon
  
  function<bool(vector<double> &)> checkWrong = [&] (vector<double> &pp)
						{
						  bool ch = true;
						  if (post(pp)<-par::defaultDouble)
						    ch = false;
						  return ch;
						};

  vector<double> par = starting_par;
  if (checkWrong(par))
    ErrorCBL("The starting position is outside the prior range: the first input parameter must be changed!", "maximize", "Posterior.cpp");

  // loop on simplex side
  for (size_t i=0; i<npar_free; i++) {
    par = starting_par;
    par[i] += epsilon;
    if (checkWrong(par))
      ErrorCBL("The simplex side is outside prior range: the epsilon parameter or the starting position must be changed.", "maximize", "Posterior.cpp");
  }

  // everything is fine up to here... let's go
  coutCBL << "Maximizing the posterior..." << endl;
  vector<double> result = cbl::wrapper::gsl::GSL_minimize_nD(post, starting_par, {}, max_iter, tol, epsilon);
  
  // check if the result is inside the prior ranges
  if (m_prior->log(result)>par::defaultDouble) {
    coutCBL << "Done!" << endl << endl;
    m_model_parameters->set_bestfit_values(result);
    m_model_parameters->write_bestfit_info();
    coutCBL << "-log(posterior) = " << post(result) << endl << endl;
  }
  else
    ErrorCBL("the maximization ended with parameter values out of the priors: check your inputs or change the epsilon value!", "maximize", "Posterior.cpp");

}


// ============================================================================================


void cbl::statistics::Posterior::sample_stretch_move (const double aa, const bool parallel, const string outputFile, const int start, const int thin, const int nbins)
{
  if (parallel && outputFile!=cbl::par::defaultString)
    WarningMsgCBL("no run-time output available for parallel stretch-move algorithm: this option will be ignored!", "sample_stretch_move", "Posterior.cpp");

  coutCBL << "Sampling the posterior..." << endl;
  
  const int seed = m_generate_seed();
  const int nparameters = m_model_parameters->nparameters();
  const int nparameters_free = m_model_parameters->nparameters_free();
  const int chain_size = m_model_parameters->chain_size();
  const int n_walkers = m_model_parameters->chain_nwalkers();

  vector<vector<double>> Start(n_walkers, vector<double>(nparameters, 0));

  for (int i=0; i<nparameters; i++)
    for (int j=0; j<n_walkers; j++)
      Start[j][i] = m_model_parameters->chain_value(i, 0, j);

  auto posterior = [this] (vector<double> &pp) { return log(pp); }; 

  cbl::statistics::Sampler sampler(nparameters, nparameters_free, posterior); 
  if (parallel)
    sampler.sample_stretch_move_parallel(chain_size, n_walkers, Start, seed, aa);
  else
    sampler.sample_stretch_move(chain_size, n_walkers, Start, seed, aa, outputFile);
  
  vector<vector<double>> chain_values;

  sampler.get_chain_function_acceptance(chain_values, m_log_posterior, m_acceptance);

  m_model_parameters->set_chain_values(chain_values, n_walkers);

  // set the best-fit parameters to the meadian values of the MCMC chain
  m_model_parameters->set_bestfit_values(start, thin, nbins, m_generate_seed());

  m_weight.erase(m_weight.begin(), m_weight.end());
  m_weight.resize(m_log_posterior.size(), 1.);
}


// ============================================================================================


void cbl::statistics::Posterior::importance_sampling (const std::string input_dir, const std::string input_file, const vector<size_t> column, const int header_lines_to_skip, const bool is_FITS_format, const bool apply_to_likelihood, const int n_walkers)
{ 
  // import the chain
  read_chain(input_dir, input_file, n_walkers, column, header_lines_to_skip, is_FITS_format);
 
  const int n_parameters = m_model_parameters->nparameters();
  const int chain_size = m_model_parameters->chain_size();

  vector<double> log_dist(m_log_likelihood.size(), 1.), input_log_dist(m_log_likelihood.size(), 1.);
 
  // loop on the chain
  for (int chain_position=0; chain_position<chain_size; ++chain_position) {
    
    // loop to parallelize the computation
#pragma omp parallel for schedule(dynamic)
    for (int walker_index=0; walker_index<n_walkers; ++walker_index) {
      
      vector<double> pp(n_parameters);
      for (int parameter=0; parameter<n_parameters; ++parameter) 
	pp[parameter] = m_model_parameters->chain_value(parameter, chain_position, walker_index);
      
      // compute the logarithm of the posterior, of likelihood, at the parameters of the input chain
      log_dist[chain_position*n_walkers+walker_index] = (apply_to_likelihood) ? Likelihood::log(pp) : log(pp);
      
      // value of the logarithm of the posterior, or likelihood, of the input chain
      input_log_dist[chain_position*n_walkers+walker_index] =
	(apply_to_likelihood) ? m_log_likelihood[chain_position*n_walkers+walker_index]
	: m_log_posterior[chain_position*n_walkers+walker_index];
    }
    
    coutCBL << int(double(chain_position)/(chain_size)*100.) << "% \r"; cout.flush();
  }
  
  
  // compute the weights as posterior, or likelihood, normalized ratios
  
  const double shift = Max(input_log_dist)-Max(log_dist);

  for (int chain_position=0; chain_position<chain_size; ++chain_position) 
#pragma omp parallel for schedule(dynamic)
    for (int walker_index=0; walker_index<n_walkers; ++walker_index) 
      m_weight[chain_position*n_walkers+walker_index] = exp(log_dist[chain_position*n_walkers+walker_index]-input_log_dist[chain_position*n_walkers+walker_index]+shift);
  
}


// ============================================================================================


void cbl::statistics::Posterior::read_chain (const string input_dir, const string input_file, const int n_walkers, const vector<size_t> column, const int header_lines_to_skip, const bool is_FITS_format)
{
  if (is_FITS_format) 
    read_chain_fits(input_dir, input_file, n_walkers, column);
  else 
    read_chain_ascii(input_dir, input_file, n_walkers, column, header_lines_to_skip);
}


// ============================================================================================


void cbl::statistics::Posterior::read_chain_ascii (const string input_dir, const string input_file, const int n_walkers, const vector<size_t> column, const int skip_header)
{
  const string file = input_dir+input_file;
  coutCBL << "Reading the chain file " << file << endl;

  const int nparameters = m_model_parameters->nparameters();

  std::function<vector<double>(const vector<double> params)> assign;

  if (column.size()==0) {
    assign = [&] (const vector<double> params) {
      vector<double> pp(nparameters);
      for (int i=0; i<nparameters; i++)
	pp[i] = params[i+1];
      return pp;
    };
  }
  else {
    checkDim(column, nparameters, "columns");
    assign = [&] (const vector<double> params) {
      vector<double> pp(nparameters);
      for (int i=0; i<nparameters; i++)
	pp[i] = params[column[i]+1];
      return pp;
    };
  }

  ifstream fin(file.c_str()); checkIO(fin, file);

  string line;
  for (int i=0; i<skip_header; i++)
    getline(fin, line);

  vector<vector<double>> chain_value;

  m_log_likelihood.erase(m_log_likelihood.begin(), m_log_likelihood.end());
  m_log_posterior.erase(m_log_posterior.begin(), m_log_posterior.end());
  m_weight.erase(m_weight.begin(), m_weight.end());

  while (getline(fin, line)) {
    stringstream ss(line);
    double NUM;
    vector<double> ll;

    while (ss>>NUM) ll.push_back(NUM);

    vector<double> params = assign(ll); 

    m_log_likelihood.push_back(ll[ll.size()-4]);
    m_log_posterior.push_back(ll[ll.size()-2]);
    m_weight.push_back(ll[ll.size()-1]);

    chain_value.push_back(params);
  }
  const int chain_size = chain_value.size()/n_walkers;
  
  fin.clear(); fin.close();

  checkDim(chain_value, n_walkers*chain_size, m_model_parameters->nparameters(), "chain_value");
  checkDim(m_log_posterior, n_walkers*chain_size, "m_log_posterior");
  checkDim(m_weight, n_walkers*chain_size, "m_weight");

  chain_value = cbl::transpose(chain_value);

  m_model_parameters->set_chain_values(chain_value, n_walkers);

  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cbl::statistics::Posterior::read_chain_fits (const string input_dir, const string input_file, const int n_walkers, const vector<size_t> column)
{
  (void)column;
  
  const int nparameters = m_model_parameters->nparameters();

  string file = input_dir+input_file;
  coutCBL << "Reading the chain file " << file << endl;

  vector<string> names = m_model_parameters->name();
  names.push_back("Log(Likelihood)");
  names.push_back("Log(Posterior)");
  names.push_back("Weight");

  vector<vector<double>> chain_value(nparameters);
  m_log_likelihood.erase(m_log_likelihood.begin(), m_log_likelihood.end());
  m_log_posterior.erase(m_log_posterior.begin(), m_log_posterior.end());
  m_weight.erase(m_weight.begin(), m_weight.end());

  vector<vector<double>> values = cbl::wrapper::ccfits::read_table_fits(file, names);

  for (int i=0; i<nparameters; i++)
    chain_value[i] = values[i];

  m_log_likelihood = values[nparameters-2];
  m_log_posterior = values[nparameters];
  m_weight = values[nparameters+1];

  int chain_size = m_log_posterior.size()/n_walkers;
 
  checkDim(chain_value, nparameters, n_walkers*chain_size, "chain_value");
  checkDim(m_log_likelihood, n_walkers*chain_size, "m_log_likelihood"); 
  checkDim(m_log_posterior, n_walkers*chain_size, "m_log_posterior"); 
  checkDim(m_weight, n_walkers*chain_size, "m_weight"); 

  m_model_parameters->set_chain_values(chain_value, n_walkers);

  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cbl::statistics::Posterior::write_chain (const string output_dir, const string output_file, const int start, const int thin, const bool is_FITS_format, const int prec, const int ww)
{
  if (is_FITS_format) 
    write_chain_fits(output_dir, output_file, start, thin);
  else 
    write_chain_ascii(output_dir, output_file, start, thin, prec, ww);
}


// ============================================================================================


void cbl::statistics::Posterior::write_chain_ascii (const string output_dir, const string output_file, const int start, const int thin, const int prec, const int ww)
{
  const int nparameters = m_model_parameters->nparameters();
  const int chain_size = m_model_parameters->chain_size();
  const int n_walkers = m_model_parameters->chain_nwalkers();
  
  checkDim(m_log_posterior, n_walkers*chain_size, "m_log_posterior", false);
  
  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {}
  
  string file = output_dir+output_file; 

  ofstream fout(file.c_str()); checkIO(fout, file);
  fout.precision(10);

  fout << "# step ";
  for (int k=0; k<nparameters; k++) 
    fout << m_model_parameters->name(k) << "  ";
  fout << " log(Likelihood)  log(Prior)  log(Posterior)  Weight" << endl;

  int nn = 0;
  
  for (int j=start; j<chain_size; j+=thin) {
    for (int i=0; i<n_walkers; i++) {
      fout << setw(6) << nn++ << "  ";
      
      vector<double> pp(nparameters);
      
      for (int k=0; k<nparameters; k++) {
	pp[k] = m_model_parameters->chain_value(k, j, i);
	cbl::Print(m_model_parameters->chain_value(k, j, i), prec, ww, "", "  ", false, fout);
      }

      if (m_log_likelihood.size()>size_t(j*n_walkers+i)) {
	cbl::Print(m_log_likelihood[j*n_walkers+i], prec, ww, "", "  ", false, fout);
	cbl::Print(m_log_posterior[j*n_walkers+i]-m_log_likelihood[j*n_walkers+i], prec, ww, "", "  ", false, fout);
      }
      else {
	const double pr = m_prior->log(pp);
	cbl::Print(m_log_posterior[j*n_walkers+i]-pr, prec, ww, "", "  ", false, fout);
	cbl::Print(pr, prec, ww, "", "  ", false, fout);
      }

      cbl::Print(m_log_posterior[j*n_walkers+i], prec, ww, "", "  ", false, fout);
      cbl::Print(m_weight[j*n_walkers+i], prec, ww, "", "\n", false, fout);
    }
  }

  fout.clear(); fout.close();

  coutCBL << "I wrote the file: " << file << endl;
}


// ============================================================================================


void cbl::statistics::Posterior::write_chain_fits (const string output_dir, const string output_file, const int start, const int thin)
{  
  const int nparameters = m_model_parameters->nparameters();
  const int chain_size = m_model_parameters->chain_size();
  const int n_walkers = m_model_parameters->chain_nwalkers();

  checkDim(m_log_posterior, n_walkers*chain_size, "m_log_posterior", false);
  
  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {}

  vector<string> names;
  names.push_back("Step");
  for (int k=0; k<nparameters; k++) 
    names.emplace_back(m_model_parameters->name(k));

  names.push_back("Log(Likelihood)");
  names.push_back("Log(Prior)");
  names.push_back("Log(Posterior)");
  names.push_back("Weight");
  
  vector<vector<double>> value(nparameters+5);

  int n = 0;
  for (int j=start; j<chain_size; j+=thin) 
    for (int i=0; i<n_walkers; i++) {
      value[0].emplace_back(n);
      vector<double> pp(nparameters);
      for (int k=0; k<nparameters; k++) {
	value[k+1].emplace_back(m_model_parameters->chain_value(k, j, i));
	pp[k] = m_model_parameters->chain_value(k, j, i);
      }

      if (m_log_likelihood.size()>size_t(j*n_walkers+i)) {
	value[nparameters+1].emplace_back(m_log_likelihood[j*n_walkers+i]);
	value[nparameters+2].emplace_back(m_log_posterior[j*n_walkers+i]-m_log_likelihood[j*n_walkers+i]);
      }
      else {
	const double lpr = m_prior->log(pp);
	value[nparameters+1].emplace_back(m_log_posterior[j*n_walkers+i]-lpr);
	value[nparameters+2].emplace_back(lpr);
      }
      
      value[nparameters+3].emplace_back(m_log_posterior[j*n_walkers+i]);
      value[nparameters+4].emplace_back(m_weight[j*n_walkers+i]);
      n ++;
    }
  
  cbl::wrapper::ccfits::write_table_fits(output_dir, output_file, names, value);
  
}


// ============================================================================================


void cbl::statistics::Posterior::initialize_chains (const int chain_size, const int n_walkers)
{
  m_model_parameters->set_chain(chain_size, n_walkers);
  m_model_parameters->initialize_chain_from_prior();
}


// ============================================================================================


void cbl::statistics::Posterior::initialize_chains (const int chain_size, const int n_walkers, const double radius, const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon)
{
  maximize(start, max_iter, tol, epsilon);

  m_model_parameters->set_chain(chain_size, n_walkers);
  m_model_parameters->initialize_chain_ball_bestfit(radius, m_generate_seed());
}


// ============================================================================================


void cbl::statistics::Posterior::initialize_chains (const int chain_size, const int n_walkers, std::vector<double> &value, const double radius)
{
  m_model_parameters->set_chain(chain_size, n_walkers);
  m_model_parameters->initialize_chain_ball(value, radius, m_generate_seed());
}


// ============================================================================================


void cbl::statistics::Posterior::initialize_chains (const int chain_size, const std::vector<std::vector<double>> chain_value)
{
  const int n_walkers = chain_value[0].size();
  m_model_parameters->set_chain(chain_size, n_walkers);

  for (size_t pp=0; pp<m_model_parameters->nparameters(); pp++)
    for (int ww=0; ww<n_walkers; ww++)
      m_model_parameters->set_chain_value(pp, 0, ww, chain_value[pp][ww]);
}


// ============================================================================================


void cbl::statistics::Posterior::initialize_chains (const int chain_size, const int n_walkers, const std::string input_dir, const std::string input_file)
{
  string last_step_file = input_dir+input_file+"_LastStep";
  string get_last_step = "tail -n "+conv(n_walkers, par::fINT)+" "+input_dir+input_file+" > "+last_step_file;
  if (system(get_last_step.c_str())) {}

  ifstream fin(last_step_file);
  string line;

  vector<vector<double>> chain_value;

  while (getline(fin, line))
    {
      stringstream ss(line);
      double NUM;
      vector<double> ll, params;

      while (ss>>NUM) ll.push_back(NUM);
      for (size_t i=1; i<ll.size()-3; i++) 
	params.push_back(ll[i]);

      chain_value.push_back(params);
    }
  fin.clear(); fin.close();

  string rm_last_step = "rm -r "+last_step_file;
  if (system(rm_last_step.c_str())) {}

  checkDim(chain_value, n_walkers, m_model_parameters->nparameters(), "chain_from_LastStep_file");
  chain_value = cbl::transpose(chain_value);

  initialize_chains(chain_size, chain_value);
}


// ============================================================================================


void cbl::statistics::Posterior::write_maximization_results (const string dir_output, const string file)
{
  coutCBL << "Writing results of posterior maximization on " << dir_output+file << endl;
  vector<double> bestFitValues = m_model_parameters->bestfit_value();
  string name = LikelihoodTypeNames ()[static_cast<int>(m_likelihood_type)];
  double posteriorValue = this->log(bestFitValues);

  string mkdir = "mkdir -p "+dir_output;
  if (system(mkdir.c_str())) {}

  ofstream fout(dir_output+file);

  fout << "#Parameters information" << endl;
  fout << "number of parameters = " << bestFitValues.size() << endl;

  for (size_t i=0; i<bestFitValues.size(); i++) {
    fout << "par" << i+1 << "_name = " << m_model_parameters->name(i) << endl;
    fout << "par" << i+1 << "_status = " << m_model_parameters->status(i) << endl;
    fout << "par" << i+1 << "_bestfit_value = " << bestFitValues[i] << endl;
  }

  fout << "#Likelihood information" << endl;
  fout << "likelihoodType = " << name << endl;
  fout << "logPosteriorValue = " << posteriorValue << endl;

  fout.clear(); fout.close();
  coutCBL << "I wrote the file " << dir_output+file << endl;
}


// ============================================================================================


void cbl::statistics::Posterior::show_results (const int start, const int thin, const int nbins, const bool show_mode, const int ns, const int nb)
{
  m_model_parameters->show_results(start, thin, nbins, m_generate_seed(), show_mode, ns, nb, weight(start, thin));
}


// ============================================================================================


void cbl::statistics::Posterior::write_results (const string output_dir, const string root_file, const int start, const int thin, const int nbins, const bool fits, const bool compute_mode, const int ns, const int nb)
{
  const string extension = (fits) ? "_chain.fits" : "_chain.dat";
  write_chain(output_dir, root_file+extension, start, thin, fits);

  m_model_parameters->write_results(output_dir, root_file, start, thin, nbins, m_generate_seed(), compute_mode, ns, nb, weight(start, thin));
}


// ============================================================================================


void cbl::statistics::Posterior::write_model_from_chain (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const int start, const int thin)
{
  switch (m_model->dimension()) {

  case Dim::_1D_: 
    {
      vector<double> xvec = xx;
      if (xx.size()==0)
	xvec = m_data->xx();

      m_model->write_from_chains(output_dir, output_file, xvec, start, thin);
    }

    break;
  case Dim::_2D_: 
    {
      vector<double> xvec = xx, yvec = yy;
      if (xx.size()==0)
	xvec = m_data->xx();
      if (yy.size()==0)
	yvec = m_data->yy();

      m_model->write_from_chains(output_dir, output_file, xvec, yvec, start, thin);
    }

    break;
  default:
    ErrorCBL("the input dimension must be Dim::_1D_ or Dim::_2D_ !", "write_model_from_chain", "Posterior.cpp");
  }
}
