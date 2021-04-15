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

/** @file Statistics/CombinedPosterior.cpp
 *
 *  @brief Methods of the class CombinedPosterior
 *
 *  This file contains the implementation of the methods of the class
 *  CombinedPosterior, used for bayesian analyses
 *
 *  @authors Davide Pelliciari, Sofia Contarini, Giorgio Lesci
 *
 *  @authors davide.pelliciari@studio.unibo.it, sofia.contarini3@unibo.it, giorgio.lesci2@unibo.it
 */

#include "FITSwrapper.h"
#include "PosteriorParameters.h"
#include "Sampler.h"
#include "CombinedPosterior.h"
#include "ChainMesh.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::statistics::CombinedPosterior::CombinedPosterior (const std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par)
{
  m_posteriors = posteriors;
  m_Nposteriors = m_posteriors.size();
  std::vector<bool> is_from_chain(m_Nposteriors);

  // eventualmente dentro una funzione che controlla e setta is_from_chain
  for(int N=0; N<m_Nposteriors; N++)
    if(m_posteriors[N]->m_get_seed() != -1) is_from_chain[N] = true;
    else is_from_chain[N] = false;

  if(std::find(is_from_chain.begin(), is_from_chain.end(), false) == is_from_chain.end()) {    // All false
    m_set_parameters_priors(m_posteriors, repeated_par);
    m_set_independent_probes();
  }

  else if (std::find(is_from_chain.begin(), is_from_chain.end(), true) == is_from_chain.end()) // All true
    {
      impsampling = true;
      for(int N=1; N<m_Nposteriors; N++)
	if(m_posteriors[N]->get_Nparameters()!=m_posteriors[0]->get_Nparameters())
	  ErrorCBL("Different number of parameters for the combination", "CombinedPosterior", "CombinedPosterior.cpp");
      m_Nparameters = m_posteriors[0]->get_Nparameters();
      m_model_parameters = m_posteriors[0]->get_model_parameters();
    }

  // at least one is false
  else
    ErrorCBL("Different kind of Posterior objects given as input!", "CombinedPosterior", "CombinedPosterior.cpp");
}


// ============================================================================================


cbl::statistics::CombinedPosterior::CombinedPosterior (const std::vector<std::vector<std::shared_ptr<Posterior>>> posteriors, const std::vector<std::shared_ptr<data::CovarianceMatrix>> covariance, const std::vector<std::shared_ptr<Posterior>> independent_posteriors, std::vector<std::string> repeated_par)
{
  // Check the covariance matrices
  if (posteriors.size() != covariance.size())
    ErrorCBL("The number of sets of dependent posteriors must match the number of covariance matrices!", "CombinedPosterior", "CombinedPosterior.cpp");

  for (size_t i=0; i<posteriors.size(); i++) {
    int dataset_length = 0;
    for (size_t j=0; j<posteriors[i].size(); j++)
      dataset_length += posteriors[i][j]->get_m_data()->xx().size();
    std::vector<std::vector<double>> cov = covariance[i]->operator()();
    if (dataset_length != (int)(cov.size()))
      ErrorCBL("Dataset "+cbl::conv(i+1, cbl::par::fINT)+": The dimension of the covariance matrix does not match the length of the data vector!", "CombinedPosterior", "CombinedPosterior.cpp");
  }

  // Set the core variables
  m_Nposteriors = posteriors.size() + independent_posteriors.size();

  m_use_grid.resize(m_Nposteriors);
  m_likelihood_inputs.resize(m_Nposteriors);
  m_log_likelihood_functions.resize(m_Nposteriors);
  m_likelihood_functions.resize(m_Nposteriors);
  m_log_likelihood_functions_grid.resize(m_Nposteriors);
  m_likelihood_functions_grid.resize(m_Nposteriors);

  // Set parameters and priors
  std::vector<std::shared_ptr<Posterior>> dummy_posteriors;
  for (size_t i=0; i<posteriors.size(); i++)
    for (size_t j=0; j<posteriors[i].size(); j++)
      dummy_posteriors.emplace_back(posteriors[i][j]);
  for (size_t i=0; i<independent_posteriors.size(); i++)
    dummy_posteriors.emplace_back(independent_posteriors[i]);

  m_set_parameters_priors(dummy_posteriors, repeated_par);

  m_parameter_indexes2.resize(dummy_posteriors.size());
  
  // Set the models and the likelihood functions for the dependent probes
  int name_idx = 0;
  for (size_t i=0; i<posteriors.size(); i++) {
    auto inputs = make_shared<STR_DependentProbes_data_model>(m_data_model);
    inputs->models.resize(posteriors[i].size(), NULL);
    inputs->xx.resize(posteriors[i].size());
    inputs->par_indexes.resize(posteriors[i].size());
    for (size_t j=0; j<posteriors[i].size(); j++) {
      if (posteriors[i][j]->get_m_model()->dimension() != Dim::_1D_)
	ErrorCBL("The model dimension of the statistically dependent probes must be 1!", "CombinedPosterior", "CombinedPosterior.cpp");
      inputs->models[j] = posteriors[i][j]->get_m_model();
      m_models.emplace_back(posteriors[i][j]->get_m_model());  // Used only for writing the models at percentile values
      m_datasets.emplace_back(posteriors[i][j]->get_m_data()); // Used only for writing the models at percentile values
      inputs->xx[j] = posteriors[i][j]->get_m_data()->xx();
      for (size_t k=0; k<posteriors[i][j]->get_model_parameters()->name().size(); k++)
	for(int rr=0; rr<m_Nparameters; rr++)
	  if (m_parameter_names[name_idx][k] == m_model_parameters->name(rr)) {
	    inputs->par_indexes[j].emplace_back(rr);
	    m_parameter_indexes2[name_idx].emplace_back(rr);
	  }
      for (size_t k=0; k<posteriors[i][j]->get_m_data()->xx().size(); k++)
	inputs->flat_data.emplace_back(posteriors[i][j]->get_m_data()->data(k));
      name_idx ++;
    }
    inputs->covariance = covariance[i];

    m_likelihood_inputs[i] = inputs;
    m_log_likelihood_functions[i] = &LogLikelihood_Gaussian_combined;
    m_likelihood_functions[i] = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_functions[i](par, input)); };

    m_use_grid[i] = false;
    m_log_likelihood_functions_grid[i] = NULL;
    m_likelihood_functions_grid[i] = NULL;
  }

  // Set the models and the likelihood functions for the independent probes
  int idx = 0;
  if (independent_posteriors.size() > 0) {
    for(int i=posteriors.size(); i<m_Nposteriors; i++) {
      m_use_grid[i] = independent_posteriors[idx]->get_m_use_grid();
      m_likelihood_inputs[i] = independent_posteriors[idx]->get_m_likelihood_inputs();
      m_log_likelihood_functions[i] = independent_posteriors[idx]->get_m_log_likelihood_function();
      m_likelihood_functions[i] = independent_posteriors[idx]->get_m_likelihood_function();

      m_models.emplace_back(independent_posteriors[idx]->get_m_model());
      m_datasets.emplace_back(independent_posteriors[idx]->get_m_data());
      
      idx ++;
    }
  }

  // Set the parameter indexes for the dependent posteriors
  m_parameter_indexes.resize(m_Nposteriors);

  for (size_t N=0; N<posteriors.size(); N++)
    for (size_t j=0; j<m_model_parameters->nparameters(); j++)
      m_parameter_indexes[N].emplace_back(j);

  // Set the parameter indexes for the independent posteriors
  for (int N=(int)(posteriors.size()); N<m_Nposteriors; N++) {
    for (size_t k=0; k<m_parameter_names[name_idx].size(); k++) {
      for(int j=0; j<m_Nparameters; j++) {
	if (m_parameter_names[name_idx][k] == m_model_parameters->name(j)) {
	  m_parameter_indexes[N].emplace_back(j);
	  m_parameter_indexes2[name_idx].emplace_back(j);
	}
      }
    }
    name_idx ++;
  }

  m_set_seed(321);
  
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_set_parameters_priors (std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par)
{
  const int dummy_Nposteriors = (int)(posteriors.size());
  
  // Set the parameter names and the priors
  std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions = posteriors[0]->get_model_parameters()->prior_distribution();
  std::vector<ParameterType> parameter_types = posteriors[0]->get_model_parameters()->type();
  std::vector<std::string> parameter_names = posteriors[0]->get_model_parameters()->name();
  
  m_parameter_names.resize(posteriors.size());
  m_parameter_names[0] = parameter_names;
  
  if (dummy_Nposteriors > 1) {
    for (int N=1; N<dummy_Nposteriors; N++) {
      for (size_t k=0; k<posteriors[N]->get_model_parameters()->name().size(); k++) {
	const bool is_in_parnames = std::count(parameter_names.begin(), parameter_names.end(), posteriors[N]->get_model_parameters()->name(k));
	const bool is_repeated = std::count(repeated_par.begin(), repeated_par.end(), posteriors[N]->get_model_parameters()->name(k));
	if (is_in_parnames && is_repeated == false) {
	  m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k));
	} else {
	  if ( (is_in_parnames && is_repeated) || (is_in_parnames == false && is_repeated) ) {
	    m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior"+cbl::conv(N+1, cbl::par::fINT));
	    parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior"+cbl::conv(N+1, cbl::par::fINT));
	  } else if (is_in_parnames == false && is_repeated == false) {
	    m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k));
	    parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
	  }
	  parameter_types.emplace_back(posteriors[N]->get_model_parameters()->type(k));
	  prior_distributions.emplace_back(posteriors[N]->get_model_parameters()->prior_distribution(k));
	}
      }
    }
  }

  m_Nparameters = (int)(prior_distributions.size());
  
  cbl::statistics::PosteriorParameters posterior_par (m_Nparameters, prior_distributions, parameter_types, parameter_names);

  // Check if the strings in repeated_par correspond to the parameter names defined internally
  for (size_t i=0; i<repeated_par.size(); i++)
    if (std::count(parameter_names.begin(), parameter_names.end(), repeated_par[i]))
      continue;
    else
      ErrorCBL("Wrong parameter name declaration in repeated_par!", "m_set_parameter_priors", "CombinedPosterior.cpp");
  
  // Set m_model_parameters
  m_model_parameters = std::make_shared<PosteriorParameters>(posterior_par);  
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_set_independent_probes ()
{
  m_use_grid.resize(m_Nposteriors);
  m_models.resize(m_Nposteriors);
  m_datasets.resize(m_Nposteriors);
  m_likelihood_inputs.resize(m_Nposteriors);
  m_log_likelihood_functions_grid.resize(m_Nposteriors);
  m_log_likelihood_functions.resize(m_Nposteriors);
  m_likelihood_functions.resize(m_Nposteriors);
  m_likelihood_functions_grid.resize(m_Nposteriors);

  for(int N=0; N<m_Nposteriors; N++) {
    m_use_grid[N] = m_posteriors[N]->get_m_use_grid();
    m_models[N] = m_posteriors[N]->get_m_model();
    m_datasets[N] = m_posteriors[N]->get_m_data();
    m_likelihood_inputs[N] = m_posteriors[N]->get_m_likelihood_inputs();
    m_log_likelihood_functions_grid[N] = m_posteriors[N]->get_m_log_likelihood_function_grid();
    m_log_likelihood_functions[N] = m_posteriors[N]->get_m_log_likelihood_function();
    m_likelihood_functions[N] = m_posteriors[N]->get_m_likelihood_function();
    m_likelihood_functions_grid[N] = m_posteriors[N]->get_m_likelihood_function_grid();
  }

  m_seed = m_posteriors[0]->m_get_seed();
  m_set_seed(m_seed);

  // Set the parameter indexes
  m_parameter_indexes.resize(m_Nposteriors);
  m_parameter_indexes2.resize(m_Nposteriors);

  for(int N=0; N<m_Nposteriors; N++)
    for (size_t k=0; k<m_parameter_names[N].size(); k++)
      for(int j=0; j<m_Nparameters; j++)
	if (m_parameter_names[N][k] == m_model_parameters->name(j)) {
	  m_parameter_indexes[N].emplace_back(j);
	  m_parameter_indexes2[N].emplace_back(j);
	}
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::set_parameters (const std::vector<std::vector<double>> parametersA, const std::vector<std::vector<double>> parametersB)
{
  std::vector<std::vector<double>> parameters (m_Nparameters);
  for(int N=0; N<m_Nparameters; N++){
    parameters[N].insert(parameters[N].end(), parametersA[N].begin(), parametersA[N].end());
    parameters[N].insert(parameters[N].end(), parametersB[N].begin(), parametersB[N].end());
  }
  m_parameters = parameters;
  parameters.clear();
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::set_log_posterior (const std::vector<double> logpostA, const std::vector<double> logpostB)
{
  m_log_posterior.reserve(logpostA.size() + logpostB.size()); // preallocate memory
  m_log_posterior.insert(m_log_posterior.end(), logpostA.begin(), logpostA.end());
  m_log_posterior.insert(m_log_posterior.end(), logpostB.begin(), logpostB.end());
}


// ============================================================================================

void cbl::statistics::CombinedPosterior::set_weight (const std::vector<double> weightsA, const std::vector<double> weightsB)
{

  m_weight.reserve(weightsA.size() + weightsB.size()); // preallocate memory
  m_weight.insert(m_weight.end(), weightsA.begin(), weightsA.end());
  m_weight.insert(m_weight.end(), weightsB.begin(), weightsB.end());

}

// ============================================================================================

void cbl::statistics::CombinedPosterior::importance_sampling (const int distNum, const double cell_size, const double rMAX, const double cut_sigma)
{
  std::vector<std::vector<double>> parametersA = m_posteriors[0]->get_parameters();
  std::vector<std::vector<double>> parametersB = m_posteriors[1]->get_parameters();
  std::vector<double> logpostA = m_posteriors[0]->get_log_posterior();
  std::vector<double> logpostB = m_posteriors[1]->get_log_posterior();

  set_parameters(parametersA, parametersB);
  set_log_posterior(logpostA, logpostB);

  // initialize the Chain Mesh for interpolation
  cbl::chainmesh::ChainMesh chmesh(cell_size, m_Nparameters);

  // interpolate using chainmesh
  std::vector<double> logpostA_interpolated = chmesh.interpolate(parametersA, logpostA, parametersB, distNum, rMAX);
  std::vector<double> logpostB_interpolated = chmesh.interpolate(parametersB, logpostB, parametersA, distNum, rMAX);

  // compute shifts

  double shift_A = *max_element(logpostA.begin(), logpostA.end()) - *max_element(logpostB_interpolated.begin(), logpostB_interpolated.end());
  double shift_B = *max_element(logpostB.begin(), logpostB.end()) - *max_element(logpostA_interpolated.begin(), logpostA_interpolated.end());

  // compute m_weight

  std::vector<double> weights_A(logpostA.size());
  std::vector<double> weights_B(logpostB.size());

  for(size_t ii=0; ii<logpostA.size(); ii++) weights_A[ii] = exp(logpostB_interpolated[ii] - logpostA[ii] + shift_A);
  for(size_t ii=0; ii<logpostB.size(); ii++) weights_B[ii] = exp(logpostA_interpolated[ii] - logpostB[ii] + shift_B);

  // cut weights distribution if cut_sigma!=-1
  if(cut_sigma>0)
    {
      const double mean_A = cbl::Average(weights_A);
      const double mean_B = cbl::Average(weights_B);
      const double sigma_A = cbl::Sigma(weights_A);
      const double sigma_B = cbl::Sigma(weights_B);

      for(size_t ii=0; ii<weights_A.size(); ii++)
	if(weights_A[ii]>mean_A+cut_sigma*sigma_A) weights_A[ii] = mean_A;
      for(size_t ii=0; ii<weights_B.size(); ii++)
	if(weights_B[ii]>mean_B+cut_sigma*sigma_B) weights_B[ii] = mean_B;
    }

  set_weight(weights_A, weights_B);
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::importance_sampling (const std::string output_path, const std::string model_nameA, const std::string model_nameB, const std::vector<double> start, const int chain_size, const int nwalkers, const int burn_in, const int thin)
{
  if(m_Nposteriors != 2) ErrorCBL("You can't do importance sampling for a number of probes > 2", "importance_sampling", "CombinedPosterior.cpp");

  impsampling = true;
  std::vector<std::string> modelnames(m_Nposteriors);
  modelnames[0] = model_nameA;
  modelnames[1] = model_nameB;

  std::string file_AB = model_nameA+"_post_"+model_nameB;
  std::string file_BA = model_nameB+"_post_"+model_nameA;

  for(int N=0; N<m_Nposteriors; N++){
    m_posteriors[N]->initialize_chains(chain_size, nwalkers, 1.e-5, start);
    coutCBL << "Sampling the posterior distribution for " << modelnames[N] << endl << endl;
    m_posteriors[N]->sample_stretch_move(2);
    m_posteriors[N]->write_results(output_path, modelnames[N], burn_in, thin);
    cout << endl;
  }

  coutCBL << "Doing importance sampling for " << modelnames[0] << "_post_" << modelnames[1] << ".." << endl;
  m_posteriors[0]->importance_sampling(output_path, modelnames[1]+"_chain.dat");
  m_posteriors[0]->write_results(output_path, file_AB);
  cout << endl;
  coutCBL << "Doing importance sampling for " << modelnames[1] << "_post_" << modelnames[0] << ".." << endl;
  m_posteriors[1]->importance_sampling(output_path, modelnames[0]+"_chain.dat");
  m_posteriors[1]->write_results(output_path, file_BA);

  std::vector<std::vector<double>> chainA;
  std::vector<std::vector<double>> chainB;

  chainA = read(output_path, file_AB+"_chain.dat");
  chainB = read(output_path, file_BA+"_chain.dat");

  std::vector<std::vector<double>> parametersA (m_Nparameters, std::vector<double> (chainA.size()));
  std::vector<std::vector<double>> parametersB (m_Nparameters, std::vector<double> (chainB.size()));
  std::vector<double> weightsAB = m_posteriors[0]->weight();
  std::vector<double> weightsBA = m_posteriors[1]->weight();

  for(int N=1; N<m_Nparameters+1; N++){
    for(size_t ii=0; ii<chainA.size(); ii++){
      parametersA[N-1][ii] = chainA[ii][N];
      parametersB[N-1][ii] = chainB[ii][N];
    }
  }

  set_log_posterior(m_posteriors[0]->get_log_posterior(), m_posteriors[1]->get_log_posterior());
  set_parameters(parametersA, parametersB);
  set_weight(weightsAB, weightsBA);
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::initialize_chains (const int chain_size, const int n_walkers, const double radius, const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon)
{
  maximize(start, max_iter, tol, epsilon);
  m_model_parameters->set_chain(chain_size, n_walkers);
  m_model_parameters->initialize_chain_ball_bestfit(radius, m_generate_seed());
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::initialize_chains (const int chain_size, const int n_walkers)
{
  m_model_parameters->set_chain(chain_size, n_walkers);
  m_model_parameters->initialize_chain_from_prior();
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::sample_stretch_move (const double aa, const bool parallel, const std::string outputFile, const int start, const int thin, const int nbins)
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
  m_model_parameters->set_bestfit_values(start, thin, nbins, m_generate_seed());
  m_weight.erase(m_weight.begin(), m_weight.end());
  m_weight.resize(m_log_posterior.size(), 1.);
    
  // Set the chain values also for each Model
  for (size_t i=0; i<m_models.size(); i++) {
    vector<vector<double>> probe_chain_values;
    probe_chain_values.resize(m_parameter_indexes2[i].size());
    for (size_t j=0; j<m_parameter_indexes2[i].size(); j++)
      probe_chain_values[j] = chain_values[m_parameter_indexes2[i][j]];
    m_models[i]->parameters()->set_chain_values(probe_chain_values, n_walkers);
    m_models[i]->parameters()->set_bestfit_values(start, thin, nbins, m_generate_seed());
  }

}


// ============================================================================================


double cbl::statistics::CombinedPosterior::operator () (std::vector<double> &pp) const
{
  pp = m_model_parameters->full_parameter(pp);
  double prior = m_model_parameters->prior()->operator()(pp);

  double val = 1.;

  for(int N=0; N<m_Nposteriors; N++)
    {
      if(prior<=0)
	{
	  val = 0.;
	  break;
	}
      std::vector<double> pp_single (m_parameter_indexes.size(), 0);
      for (size_t ii=0; ii<pp_single.size(); ii++)
	pp_single[ii] = pp[m_parameter_indexes[N][ii]];
      val *= (m_use_grid[N]) ? m_likelihood_functions_grid[N](pp_single, m_likelihood_inputs[N]) : m_likelihood_functions[N](pp_single, m_likelihood_inputs[N]);
    }
  
  return val*prior;
}


// ============================================================================================


double cbl::statistics::CombinedPosterior::log (std::vector<double> &pp) const
{
  pp = m_model_parameters->full_parameter(pp);  
  const double logprior = m_model_parameters->prior()->log(pp);  

  double val = 0.;
  
  for(int N=0; N<m_Nposteriors; N++)
    {
      if (logprior>par::defaultDouble) {
	std::vector<double> pp_single (m_parameter_indexes[N].size(), 0);
	for (size_t ii=0; ii<pp_single.size(); ii++)
	  pp_single[ii] = pp[m_parameter_indexes[N][ii]];
	val += (m_use_grid[N]) ? m_log_likelihood_functions_grid[N](pp_single, m_likelihood_inputs[N]) : m_log_likelihood_functions[N](pp_single, m_likelihood_inputs[N]);
      }
      else
	{
	  val = par::defaultDouble;
	  break;
	}
    }

  return val+logprior;
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::maximize (const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon)
{
  vector<double> starting_par;

  unsigned int npar = m_model_parameters->nparameters();
  unsigned int npar_free = m_model_parameters->nparameters_free();
  
  if (start.size()==npar_free)
    starting_par = start;
  else if (start.size()==npar)
    for (size_t i=0; i<npar_free; i++)
      starting_par.push_back(start[m_model_parameters->free_parameter()[i]]);
  else
    ErrorCBL("check your inputs: start.size()="+conv(start.size(), par::fINT)+" must be equal to either npar_free="+conv(npar_free, par::fINT)+" or npar="+conv(npar, par::fINT)+"!", "maximize", "CombinedPosterior.cpp");

  function<double(vector<double> &)> post = [this](vector<double> &pp) { return -this->log(pp); };


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
    ErrorCBL("The starting position is outside the prior range: the first input parameter must be changed!", "maximize", "CombinedPosterior.cpp");  
  
  // loop on simplex side
  for (size_t i=0; i<npar_free; i++) {
    par = starting_par;
    par[i] += epsilon;
    if (checkWrong(par))
      ErrorCBL("The simplex side is outside prior range: the epsilon parameter or the starting position must be changed.", "maximize", "CombinedPosterior.cpp");
  }  

  // everything is fine up to here... let's go
  coutCBL << "Maximizing the posterior..." << endl;
  vector<double> result = cbl::wrapper::gsl::GSL_minimize_nD(post, starting_par, {}, max_iter, tol, epsilon);
  // check if the result is inside the prior ranges

  if(m_model_parameters->prior()->log(result)<=par::defaultDouble)
    ErrorCBL("the maximization ended with parameter values out of the priors: check your inputs or change the epsilon value!", "maximize", "CombinedPosterior.cpp");

  coutCBL << "Done!" << endl << endl;

  m_model_parameters->set_bestfit_values(result);
  m_model_parameters->write_bestfit_info();

  coutCBL << "log(posterior) = " << post(result) << endl << endl;

}


// ============================================================================================


void cbl::statistics::CombinedPosterior::show_results (const int start, const int thin, const int nbins, const bool show_mode, const int ns, const int nb)
{
  if(!impsampling) Posterior::show_results(start, thin, nbins, show_mode, ns, nb);
  else{
    for (int i=0; i<m_Nparameters; i++){
      coutCBL << "Parameter: " << par::col_yellow << this->parameters()->name(i) << par::col_default << endl;
      coutCBL << "Weighted Average: " << cbl::Average(m_parameters[i], m_weight) << endl;
      coutCBL << "Standard Deviation: " << cbl::Sigma(m_parameters[i], m_weight) << endl << endl;
    }
  }
}
// ============================================================================================


void cbl::statistics::CombinedPosterior::write_results (const string output_dir, const string root_file, const int start, const int thin, const int nbins, const bool fits, const bool compute_mode, const int ns, const int nb)
{
  if(!impsampling){
    const string extension = (fits) ? "_chain.fits" : "_chain.dat";
    write_chain(output_dir, root_file+extension, start, thin, fits);
    m_model_parameters->write_results(output_dir, root_file, start, thin, nbins, m_generate_seed(), compute_mode, ns, nb, weight(start, thin));
  }
  else {
    string file = output_dir+root_file+"_chain.dat";
    ofstream fout(file.c_str()); checkIO(fout, file);

    fout << "# step" << setw(25);
    for (int k=0; k<m_Nparameters; k++){
      fout << this->parameters()->name(k) << setw(25);
    }
    fout << "log(Posterior)" << setw(25) << "Weight" << endl;

    fout << std::fixed;
    fout << setprecision(7);

    for(size_t ii=0; ii<m_weight.size(); ii++)
      {
	fout << ii;
	for(int N=0; N<m_Nparameters; N++){
	  fout << setw(25) << std::fixed << m_parameters[N][ii];
	}
	fout << setw(25) << std::fixed << m_log_posterior[ii] << setw(25) << std::scientific << m_weight[ii] << endl;
      }
    fout.clear(); fout.close();

    coutCBL << "I wrote the file: " << file << endl;
  }
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::write_chain (const string output_dir, const string output_file, const int start, const int thin, const bool is_FITS_format, const int prec, const int ww)
{
  if (is_FITS_format)
    write_chain_fits(output_dir, output_file, start, thin);
  else
    write_chain_ascii(output_dir, output_file, start, thin, prec, ww);
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::write_chain_ascii (const string output_dir, const string output_file, const int start, const int thin, const int prec, const int ww)
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

      double pr = m_model_parameters->prior()->log(pp);
      cbl::Print(m_log_posterior[j*n_walkers+i]-pr, prec, ww, "", "  ", false, fout);
      cbl::Print(pr, prec, ww, "", "  ", false, fout);

      cbl::Print(m_log_posterior[j*n_walkers+i], prec, ww, "", "  ", false, fout);
      cbl::Print(m_weight[j*n_walkers+i], prec, ww, "", "\n", false, fout);
    }
  }

  fout.clear(); fout.close();

  coutCBL << "I wrote the file: " << file << endl;
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::write_chain_fits (const string output_dir, const string output_file, const int start, const int thin)
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


      double lpr = m_model_parameters->prior()->log(pp);
      value[nparameters+1].emplace_back(m_log_posterior[j*n_walkers+i]-lpr);
      value[nparameters+2].emplace_back(lpr);

      value[nparameters+3].emplace_back(m_log_posterior[j*n_walkers+i]);
      value[nparameters+4].emplace_back(m_weight[j*n_walkers+i]);
      n ++;
    }

  cbl::wrapper::ccfits::write_table_fits(output_dir, output_file, names, value);

}


// ============================================================================================


void cbl::statistics::CombinedPosterior::write_model_from_chain (const std::string output_dir, const std::string output_file, const int start, const int thin)
{
  for(size_t N=0; N<m_models.size(); N++) {
  
    switch (m_models[N]->dimension()) {

    case Dim::_1D_:
      {
	vector<double> xvec = m_datasets[N]->xx();
	m_models[N]->write_from_chains(output_dir, std::to_string(N)+"_"+output_file, xvec, start, thin);
      }
      break;
    
    case Dim::_2D_:
      {
	vector<double> xvec = m_datasets[N]->xx();
	vector<double> yvec = m_datasets[N]->yy();
	m_models[N]->write_from_chains(output_dir, std::to_string(N)+"_"+output_file, xvec, yvec, start, thin);
      }
      break;
    
    default:
      ErrorCBL("the input dimension must be Dim::_1D_ or Dim::_2D_ !", "write_model_from_chain", "CombinedPosterior.cpp");
      
    }
  }
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::write_maximization_results (const std::string dir_output, const std::string file)
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


std::vector<std::vector<double>> cbl::statistics::CombinedPosterior::read (const std::string path, const std::string filename)
{
  std::vector< std::vector<double>> table;
  std::fstream ifs;

  ifs.open(path+filename);

  while (true)
    {
      std::string line;
      double buf;
      getline(ifs, line);

      std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

      if (!ifs)
	// mainly catch EOF
	break;

      if (line[0] == '#' || line.empty())
	// catch empty lines or comment lines
	continue;

      std::vector<double> row;

      while (ss >> buf)
	row.push_back(buf);

      table.push_back(row);
    }

  ifs.close();

  return table;
}


// ============================================================================================


double cbl::statistics::LogLikelihood_Gaussian_combined (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters -----    
  shared_ptr<statistics::STR_DependentProbes_data_model> pp = static_pointer_cast<statistics::STR_DependentProbes_data_model>(fixed_parameter);

  // ----- compute the model values -----   
  vector<double> computed_model;
  for (size_t i=0; i<pp->models.size(); i++) {
    std::vector<double> single_par ((int)(pp->par_indexes[i].size()));
    for (size_t jj=0; jj<single_par.size(); jj++)
      single_par[jj] = likelihood_parameter[pp->par_indexes[i][jj]];
    std::vector<double> single_probe_model = pp->models[i]->operator()(pp->xx[i], single_par);
    for (size_t kk=0; kk<single_probe_model.size(); kk++)
      computed_model.emplace_back(single_probe_model[kk]);
  }
   
  // ----- compute the difference between model and data at each bin -----    
  vector<double> diff(pp->flat_data.size(), 0);
  for (size_t i=0; i<pp->flat_data.size(); i++)
    diff[i] = pp->flat_data[i]-computed_model[i];

  // ----- define the inverse covariance matrix -----
  data::CovarianceMatrix cov = *pp->covariance;
  std::vector<std::vector<double>> inverse_covariance = cov.precision();
  
  // ----- estimate the Gaussian log-likelihood -----    
  double LogLikelihood = 0.;
  for (size_t i=0; i<pp->flat_data.size(); i++)
    for (size_t j=0; j<pp->flat_data.size(); j++)
      LogLikelihood += diff[i]*inverse_covariance[i][j]*diff[j];

  return -0.5*LogLikelihood - sqrt(pow(2*cbl::par::pi, computed_model.size()) * cov.determinant());
}
