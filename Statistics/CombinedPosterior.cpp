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


cbl::statistics::CombinedPosterior::CombinedPosterior (const std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par)
{      
  m_posteriors = posteriors;
  m_Nposteriors = m_posteriors.size();
  std::vector<bool> is_from_chain(m_Nposteriors);

  // eventualmente dentro una funzione che controlla e setta is_from_chain
  for(int N=0; N<m_Nposteriors; N++)
    if(m_posteriors[N]->m_get_seed() != -1) is_from_chain[N] = true;
    else is_from_chain[N] = false;

  if(std::find(is_from_chain.begin(), is_from_chain.end(), false) == is_from_chain.end()) {    // All false
    m_set_parameters_priors(m_posteriors, repeated_par, common_repeated_par);
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


cbl::statistics::CombinedPosterior::CombinedPosterior (const std::vector<std::vector<std::shared_ptr<Posterior>>> posteriors, const std::vector<std::shared_ptr<data::CovarianceMatrix>> covariance, const std::vector<cbl::statistics::LikelihoodType> likelihood_types, const std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par)
{
  // Check the likelihood types and the covariance matrices
  int N_userdefined = 0;
  std::vector<bool> userdefined;
  for (size_t i=0; i<posteriors.size(); i++) {
    int dataset_length = 0;
    
    switch (likelihood_types[i])
      {
	
      case (LikelihoodType::_Poissonian_): case (LikelihoodType::_Gaussian_Covariance_):
	for (size_t j=0; j<posteriors[i].size(); j++)
	  dataset_length += posteriors[i][j]->get_m_data()->xx().size();
	if (dataset_length != (int)(covariance[i]->operator()().size()))
	  ErrorCBL("Dataset with index "+cbl::conv(i, cbl::par::fINT)+": The dimension of the covariance matrix does not match the length of the data vector!", "CombinedPosterior", "CombinedPosterior.cpp");
	userdefined.emplace_back(false);
	break;

      case (LikelihoodType::_UserDefined_):
	if (posteriors[i].size() != 1)
	  ErrorCBL("If the likelihood type is _UserDefined_, you must provide one Posterior object. You provided "+cbl::conv((int)(posteriors.size()), cbl::par::fINT)+" Posterior objects in posteriors["+cbl::conv(i, cbl::par::fINT)+"].", "CombinedPosterior", "CombinedPosterior.cpp");

	userdefined.emplace_back(true);
	N_userdefined ++;
	break;

      default:
	ErrorCBL("Wrong likelihood type declaration for the set of variables with index "+cbl::conv(i, cbl::par::fINT)+". Choose among: _Gaussian_Covariance_, _Poissonian_, _UserDefined_", "CombinedPosterior", "CombinedPosterior.cpp");
	break;
	
      }
  }
  
  if (N_userdefined > 0) {
    bool userdefined_started = false;
    for (size_t i=0; i<userdefined.size(); i++) {
      if (userdefined[i] == true)
	userdefined_started = true;
      else if (userdefined[i] == false && userdefined_started == true)
	ErrorCBL("You must sort the posteriors elements so that the sets of probes described by UserDefined likelihoods are the last elements in the vector!", "CombinedPosterior", "CombinedPosterior.cpp");
      else
	continue;
    }    
    if (covariance.size() != likelihood_types.size()-N_userdefined)
      ErrorCBL("If a set of probes is described by a UserDefined likelihood, you can't give in input a CovarianceMatrix object, in this constructor, related to that set of probes!\nThe number of required CovarianceMatrix objects is "+cbl::conv((int)(likelihood_types.size()-N_userdefined), cbl::par::fINT)+", i.e. the number of non-UserDefined likelihoods.", "CombinedPosterior", "CombinedPosterior.cpp");
  }
  if (posteriors.size() != likelihood_types.size() && N_userdefined == 0)
    ErrorCBL("The number of sets of dependent posteriors must match the number of likelihood types!", "CombinedPosterior", "CombinedPosterior.cpp");

  if (posteriors.size() != covariance.size()+N_userdefined)
    ErrorCBL("The number of sets of dependent posteriors must match the number of covariance matrices!", "CombinedPosterior", "CombinedPosterior.cpp");  
  
  
  // Set the core variables
  m_Nposteriors = posteriors.size();

  m_parameter_indexes.resize(m_Nposteriors);
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

  m_set_parameters_priors(dummy_posteriors, repeated_par, common_repeated_par);

  m_parameter_indexes2.resize(dummy_posteriors.size());
  
  // Set the models and the likelihood functions
  int name_idx = 0;
  for (size_t i=0; i<posteriors.size(); i++) {
    auto inputs = make_shared<STR_DependentProbes_data_model>(m_data_model);
    
    switch (likelihood_types[i])
      {
	
      case (LikelihoodType::_Gaussian_Covariance_): case (LikelihoodType::_Poissonian_):

	if (covariance[i]->isSet_SSC())
	  if (covariance[i]->SSC()->Sij_dimension() != (int)(posteriors[i].size()))
	    ErrorCBL("The dimension of the Sij matrix in the super-sample covariance is "+cbl::conv(covariance[i]->SSC()->Sij_dimension(), cbl::par::fINT)+", but posteriors["+cbl::conv((int)(i), cbl::par::fINT)+" contains "+cbl::conv((int)(posteriors.size()), cbl::par::fINT)+" Posterior objects! These dimensions must be the same.", "CombinedPosterior", "CombinedPosterior.cpp");
	
	inputs->models.resize(posteriors[i].size(), NULL);
	inputs->xx.resize(posteriors[i].size());
	inputs->par_indexes.resize(posteriors[i].size());
	for (size_t j=0; j<posteriors[i].size(); j++) {
	  if (posteriors[i][j]->get_m_model()->dimension() != Dim::_1D_)
	    ErrorCBL("The model dimension of the statistically dependent probes must be 1!", "CombinedPosterior", "CombinedPosterior.cpp");
	  inputs->models[j] = posteriors[i][j]->get_m_model();
	  m_models.emplace_back(posteriors[i][j]->get_m_model());  // Used only for writing the models at percentiles
	  m_datasets.emplace_back(posteriors[i][j]->get_m_data()); // Used only for writing the models at percentiles
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
	inputs->cosmoPar_indexes = m_cosmoPar_indexes;

	for (size_t j=0; j<m_model_parameters->nparameters(); j++)
	  m_parameter_indexes[i].emplace_back(j);
	
	break;

      case (LikelihoodType::_UserDefined_):

	for (size_t k=0; k<m_parameter_names[name_idx].size(); k++)
	  for(int j=0; j<m_Nparameters; j++)
	    if (m_parameter_names[name_idx][k] == m_model_parameters->name(j)) {
	      m_parameter_indexes[i].emplace_back(j);
	      m_parameter_indexes2[name_idx].emplace_back(j);
	    }
	name_idx ++;
	
	break;

      default:
	break;
      }
    
    switch (likelihood_types[i])
      {      
      case (LikelihoodType::_Gaussian_Covariance_):
	if (covariance[i]->isSet_SSC())
	  ErrorCBL("The super-sample covariance in the Gaussian likelihood case is not available yet!", "CombinedPosterior", "CombinedPosterior.cpp");
	m_log_likelihood_functions[i] = &LogLikelihood_Gaussian_combined;
	m_likelihood_inputs[i] = inputs;
	break;
	
      case (LikelihoodType::_Poissonian_):
	if (covariance[i]->isSet_SSC())
	  m_log_likelihood_functions[i] = &LogLikelihood_Poissonian_SSC_combined;
	else
	  m_log_likelihood_functions[i] = &LogLikelihood_Poissonian_combined;
	m_likelihood_inputs[i] = inputs;
	break;

      case (LikelihoodType::_UserDefined_):
        WarningMsgCBL("The likelihood type is UserDefined. Note that the user must set the natural logarithm of the likelihood.", "CombinedPosterior", "CombinedPosterior.cpp");
	m_log_likelihood_functions[i] = posteriors[i][0]->get_m_log_likelihood_function();
	m_likelihood_inputs[i] = posteriors[i][0]->get_m_likelihood_inputs();
	break;

      default:
	break;
      }

    m_likelihood_functions[i] = [&] (vector<double> &par, const shared_ptr<void> input) { return exp(m_log_likelihood_functions[i](par, input)); };
    m_use_grid[i] = false;
    m_log_likelihood_functions_grid[i] = NULL;
    m_likelihood_functions_grid[i] = NULL;
  }

  m_set_seed(321);
  
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_check_repeated_par (int dummy_Nposteriors, std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par)
{
  std::vector<bool> at_least_one (repeated_par.size(), false);
  for (int i=0; i<dummy_Nposteriors; i++) {
    std::vector<std::string> names = posteriors[i]->get_model_parameters()->name();
    for (size_t j=0; j<repeated_par.size(); j++) {
      if (std::count(names.begin(), names.end(), repeated_par[j]))
	at_least_one[j] = true;
      else
	continue;
    }
  }
  for (size_t j=0; j<repeated_par.size(); j++)
    if (at_least_one[j] == false)
      ErrorCBL("Wrong parameter name declaration in repeated_par! The parameter \""+repeated_par[j]+"\" does not exist or is not set.", "m_check_repeated_par", "CombinedPosterior.cpp");
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_check_common_repeated_par (int dummy_Nposteriors, std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par, std::vector<std::vector<std::vector<int>>> common_repeated_par)
{
  if (common_repeated_par.size() != repeated_par.size())
    ErrorCBL("The size of common_repeated_par must match the size of repeated_par!", "m_check_common_repeated_par", "CombinedPosterior.cpp");
  for (size_t cc=0; cc<common_repeated_par.size(); cc++) {
    std::vector<int> indexes;
    for (size_t dd=0; dd<common_repeated_par[cc].size(); dd++) {
      if (common_repeated_par[cc][dd].size() < 2)
	ErrorCBL("Error in common_repeated_par["+cbl::conv(cc, cbl::par::fINT)+"]["+cbl::conv(dd, cbl::par::fINT)+"]: the vectors must have dimension > 1!", "m_set_parameters_priors", "CombinedPosterior.cpp");
      for (size_t ee=0; ee<common_repeated_par[cc][dd].size(); ee++) {
	if (common_repeated_par[cc][dd][ee] > (int)(dummy_Nposteriors-1) || common_repeated_par[cc][dd][ee] < 0)
	  ErrorCBL("Error in common_repeated_par["+cbl::conv(cc, cbl::par::fINT)+"]["+cbl::conv(dd, cbl::par::fINT)+"]: the index "+cbl::conv(common_repeated_par[cc][dd][ee], cbl::par::fINT)+" does not correspond to any probe given in input!", "m_check_common_repeated_par", "CombinedPosterior.cpp");
	std::vector<std::string> names = posteriors[common_repeated_par[cc][dd][ee]]->get_model_parameters()->name();
	if (std::count(names.begin(), names.end(), repeated_par[cc]) == false)
	  ErrorCBL("Error in common_repeated_par: the parameter "+repeated_par[cc]+" is not set for the probe with index "+cbl::conv(common_repeated_par[cc][dd][ee], cbl::par::fINT)+"!", "m_check_common_repeated_par", "CombinedPosterior.cpp");
	else
	  indexes.emplace_back(common_repeated_par[cc][dd][ee]);
      }
    }
    if (cbl::different_elements(indexes).size() != indexes.size())
      ErrorCBL("Posterior index declared multiple times in common_repeated_par["+cbl::conv(cc, cbl::par::fINT)+"]!", "m_check_common_repeated_par", "CombinedPosterior.cpp");
  }
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_add_prior (bool par_is_repeated, std::vector<std::shared_ptr<Posterior>> posteriors, const int N, const int k, std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> &prior_distributions)
{
  switch (posteriors[N]->get_model_parameters()->type(k)) {
  case (cbl::statistics::ParameterType::_Base_):
    prior_distributions.emplace_back(posteriors[N]->get_model_parameters()->prior_distribution(k));
    break;
  case (cbl::statistics::ParameterType::_Derived_):
    break;
  case (cbl::statistics::ParameterType::_Correlated_):
    if (par_is_repeated)
      ErrorCBL("For the moment, the code can manage only Base and Derived parameters (in the case of repeated parameters)!", "m_add_prior", "CombinedPosterior.cpp");
    else
      prior_distributions.emplace_back(posteriors[N]->get_model_parameters()->prior_distribution(k));
    break;
  }
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_set_common_repeated_par (std::vector<std::shared_ptr<Posterior>> posteriors, const bool is_in_parnames, const int N, const int k, const std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par, std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> &prior_distributions, std::vector<std::string> &parameter_names, std::vector<std::string> &original_names, std::vector<ParameterType> &parameter_types)
{
  auto itr = std::find(repeated_par.begin(), repeated_par.end(), posteriors[N]->get_model_parameters()->name(k));
  int the_i = std::distance(repeated_par.begin(), itr);
  if (common_repeated_par[the_i].size() == 0) {
    if (is_in_parnames == false)
      original_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
    m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior_"+cbl::conv(N+1, cbl::par::fINT));
    parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior_"+cbl::conv(N+1, cbl::par::fINT));
    parameter_types.emplace_back(posteriors[N]->get_model_parameters()->type(k));
    m_add_prior(true, posteriors, N, k, prior_distributions);
  } else {
    for (size_t dd=0; dd<common_repeated_par[the_i].size(); dd++) {
      std::string posterior_name = "_Posterior_";
      if (std::count(common_repeated_par[the_i][dd].begin(), common_repeated_par[the_i][dd].end(), N)) {
	for (size_t ee=0; ee<common_repeated_par[the_i][dd].size(); ee++) {
	  if (ee == common_repeated_par[the_i][dd].size()-1)
	    posterior_name += cbl::conv(common_repeated_par[the_i][dd][ee]+1, cbl::par::fINT);
	  else
	    posterior_name += cbl::conv(common_repeated_par[the_i][dd][ee]+1, cbl::par::fINT)+"_";
	}
	m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k)+posterior_name);
	if (N == common_repeated_par[the_i][dd][0] && std::count(common_repeated_par[the_i][dd].begin(), common_repeated_par[the_i][dd].end(), 0) == false) {
	  if (is_in_parnames == false)
	    original_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
	  parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k)+posterior_name);
	  parameter_types.emplace_back(posteriors[N]->get_model_parameters()->type(k));
	  m_add_prior(true, posteriors, N, k, prior_distributions);
	}
      } else if ( std::count(common_repeated_par[the_i][dd].begin(), common_repeated_par[the_i][dd].end(), N) == false ) {
	bool N_in_one = false;
	for (size_t new_dd=0; new_dd<common_repeated_par[the_i].size(); new_dd++)
	  if (std::count(common_repeated_par[the_i][new_dd].begin(), common_repeated_par[the_i][new_dd].end(), N))
	    N_in_one = true;
	if (N_in_one == false && dd == 0) {
	  if (is_in_parnames == false)
	    original_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
	  m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior_"+cbl::conv(N+1, cbl::par::fINT));
	  parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior_"+cbl::conv(N+1, cbl::par::fINT));
	  parameter_types.emplace_back(posteriors[N]->get_model_parameters()->type(k));
	  m_add_prior(true, posteriors, N, k, prior_distributions);
	} else
	  continue;
      }
    } 
  }
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_set_repeated_par (std::vector<std::shared_ptr<Posterior>> posteriors, const bool is_in_parnames, const int N, const int k, const std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par, std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> &prior_distributions, std::vector<std::string> &parameter_names, std::vector<std::string> &original_names, std::vector<ParameterType> &parameter_types)
{
  // First, check if there are common repeated parameters
  if (common_repeated_par.size() > 0) {
    m_set_common_repeated_par(posteriors, is_in_parnames, N, k, repeated_par, common_repeated_par, prior_distributions, parameter_names, original_names, parameter_types);
  }
  else { 
    if (is_in_parnames == false)
      original_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
    m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior_"+cbl::conv(N+1, cbl::par::fINT));
    parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k)+"_Posterior_"+cbl::conv(N+1, cbl::par::fINT));
    parameter_types.emplace_back(posteriors[N]->get_model_parameters()->type(k));
    m_add_prior(true, posteriors, N, k, prior_distributions);
  }
}


// ============================================================================================


void cbl::statistics::CombinedPosterior::m_set_parameters_priors (std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par)
{
  const int dummy_Nposteriors = (int)(posteriors.size());

  // Check repeated_par
  if (repeated_par.size() > 0)
    m_check_repeated_par(dummy_Nposteriors, posteriors, repeated_par);

  // Check common_repeated_par
  if (common_repeated_par.size() > 0)
    m_check_common_repeated_par(dummy_Nposteriors, posteriors, repeated_par, common_repeated_par);
  
  // Set the parameter types, names and priors, starting from the first Posterior object
  std::vector<ParameterType> parameter_types = posteriors[0]->get_model_parameters()->type();
  std::vector<std::string> parameter_names = posteriors[0]->get_model_parameters()->name();

  std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions;
  for (size_t k=0; k<parameter_types.size(); k++)
    m_add_prior(false, posteriors, 0, k, prior_distributions);
  
  std::vector<std::string> original_names = parameter_names;
  
  m_parameter_names.resize(posteriors.size());
  m_parameter_names[0] = parameter_names;

  // If the number Posterior objects is >1, set the other parameters
  if (dummy_Nposteriors > 1) {
    for (int N=1; N<dummy_Nposteriors; N++) {
      for (size_t k=0; k<posteriors[N]->get_model_parameters()->name().size(); k++) {
	
	const bool is_in_parnames = std::count(original_names.begin(), original_names.end(), posteriors[N]->get_model_parameters()->name(k));
	const bool is_repeated = std::count(repeated_par.begin(), repeated_par.end(), posteriors[N]->get_model_parameters()->name(k));

	if (is_in_parnames && is_repeated == false) {
	  m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k));
	}
	else if (is_in_parnames == false && is_repeated == false) {
	  original_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
	  
	  m_parameter_names[N].emplace_back(posteriors[N]->get_model_parameters()->name(k));
	  parameter_names.emplace_back(posteriors[N]->get_model_parameters()->name(k));
	  parameter_types.emplace_back(posteriors[N]->get_model_parameters()->type(k));
	  m_add_prior(false, posteriors, N, k, prior_distributions);
	}
	else if (is_repeated) {
	  m_set_repeated_par(posteriors, is_in_parnames, N, k, repeated_par, common_repeated_par, prior_distributions, parameter_names, original_names, parameter_types);
	}
      }
    }
  }

  m_Nparameters = (int)(parameter_types.size());
  

  // Rename the parameters of the first probe, if necessary
  if (repeated_par.size() > 0) {
    if (common_repeated_par.size() > 0) {
      for (size_t i=0; i<posteriors[0]->get_model_parameters()->name().size(); i++) {
	if (std::count(repeated_par.begin(), repeated_par.end(), posteriors[0]->get_model_parameters()->name(i))) {
	  auto itr = std::find(repeated_par.begin(), repeated_par.end(), posteriors[0]->get_model_parameters()->name(i));
	  int the_i = std::distance(repeated_par.begin(), itr);
	  if (common_repeated_par[the_i].size() > 0) {
	    for (size_t jj=0; jj<common_repeated_par[the_i].size(); jj++) {
	      if (common_repeated_par[the_i][jj].size() > 0) {
		if (std::count(common_repeated_par[the_i][jj].begin(), common_repeated_par[the_i][jj].end(), 0)) {
		  parameter_names[i] += "_Posterior_";
		  m_parameter_names[0][i] += "_Posterior_"; 
		  for (size_t dd=0; dd<common_repeated_par[the_i][jj].size(); dd++) {
		    if (dd == common_repeated_par[the_i][jj].size()-1) {
		      parameter_names[i] += cbl::conv(common_repeated_par[the_i][jj][dd]+1, cbl::par::fINT);
		      m_parameter_names[0][i] += cbl::conv(common_repeated_par[the_i][jj][dd]+1, cbl::par::fINT);
		    } else {
		      parameter_names[i] += cbl::conv(common_repeated_par[the_i][jj][dd]+1, cbl::par::fINT)+"_";
		      m_parameter_names[0][i] += cbl::conv(common_repeated_par[the_i][jj][dd]+1, cbl::par::fINT)+"_";
		    }
		  }
		} else if (std::count(common_repeated_par[the_i][jj].begin(), common_repeated_par[the_i][jj].end(), 0) == false) {
		  bool zero_in_one = false;
		  for (size_t new_jj=0; new_jj<common_repeated_par[the_i].size(); new_jj++)
		    if (std::count(common_repeated_par[the_i][new_jj].begin(), common_repeated_par[the_i][new_jj].end(), 0))
		      zero_in_one = true;
		  if (zero_in_one == false && jj == 0) {
		    parameter_names[i] += "_Posterior_1";
		    m_parameter_names[0][i] += "_Posterior_1";
		  } else
		    continue;
		}
	      } else if (common_repeated_par[the_i][jj].size() == 0) {
		parameter_names[i] += "_Posterior_1";
		m_parameter_names[0][i] += "_Posterior_1";
	      }
	    }
	  } else {
	    parameter_names[i] += "_Posterior_1";
	    m_parameter_names[0][i] += "_Posterior_1";
	  }
	}
      }
    } else if (common_repeated_par.size() == 0) {
      for (size_t i=0; i<posteriors[0]->get_model_parameters()->name().size(); i++)
	if (std::count(repeated_par.begin(), repeated_par.end(), posteriors[0]->get_model_parameters()->name(i))) {
	  parameter_names[i] += "_Posterior_1";
	  m_parameter_names[0][i] += "_Posterior_1";
	}
    }
  }
  
  // Set m_model_parameters
  cbl::statistics::PosteriorParameters posterior_par (m_Nparameters, prior_distributions, parameter_types, parameter_names);
  m_model_parameters = std::make_shared<PosteriorParameters>(posterior_par);

  // Print the parameters for each probe on screen
  std::cout<<std::endl;
  for (size_t i=0; i<m_parameter_names.size(); i++) {
    coutCBL<<"For the probe "<<i+1<<" (with internal index "<<i<<"), the following parameters are set:\n";
    for (size_t j=0; j<m_parameter_names[i].size(); j++) {
      if (j == m_parameter_names[i].size()-1)
	std::cout<<m_parameter_names[i][j]<<std::endl<<std::endl;
      else if (j==0)
	coutCBL<<m_parameter_names[i][j]<<", ";
      else
	std::cout<<m_parameter_names[i][j]<<", ";
    }
  }
  coutCBL<<"Thus the guess vector should be sorted in this way:\n";
  for (size_t j=0; j<parameter_names.size(); j++)
    if (j == parameter_names.size()-1)
      std::cout<<parameter_names[j]<<std::endl<<std::endl;
    else if (j==0)
      coutCBL<<parameter_names[j]<<", ";
    else
      std::cout<<parameter_names[j]<<", ";

  
  // Find the indexes of the cosmological parameters, if any.
  // This is useful for the super-sample covariance.
  std::vector<std::string> cosmoNames = cbl::cosmology::CosmologicalParameterNames();
  for (size_t i=0; i<parameter_names.size(); i++)
    if (std::count(cosmoNames.begin(), cosmoNames.end(), parameter_names[i]))
      m_cosmoPar_indexes.emplace_back(i);

  // Check if any cosmological parameter is repeated
  for (size_t i=0; i<repeated_par.size(); i++)
    if (std::count(cosmoNames.begin(), cosmoNames.end(), repeated_par[i]))
      ErrorCBL("You cannot have more than one posterior for a cosmological parameter ("+repeated_par[i]+" in this case)!", "m_set_parameters_priors", "CombinedPosterior.cpp");
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
  for(int N=0; N<m_Nposteriors; N++) {
    for (size_t k=0; k<m_parameter_names[N].size(); k++) {
      for(int j=0; j<m_Nparameters; j++) {
	if (m_parameter_names[N][k] == m_model_parameters->name(j)) {
	  m_parameter_indexes[N].emplace_back(j);
	  m_parameter_indexes2[N].emplace_back(j);
	}
      }
    }
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


void cbl::statistics::CombinedPosterior::initialize_chains (const int chain_size, const int n_walkers, const std::string input_dir, const std::string input_file, const int seed)
{
  m_set_seed(seed);
  
  m_model_parameters->set_chain(chain_size, n_walkers);

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
      for (size_t i=1; i<ll.size()-4; i++) 
	params.push_back(ll[i]);
  
      chain_value.push_back(params);
    }
  fin.clear(); fin.close();
  
  string rm_last_step = "rm -r "+last_step_file;
  if (system(rm_last_step.c_str())) {}
  
  checkDim(chain_value, n_walkers, m_model_parameters->nparameters(), "chain_from_LastStep_file");
  chain_value = cbl::transpose(chain_value);
    
  for (size_t pp=0; pp<m_model_parameters->nparameters(); pp++)
    for (int ww=0; ww<n_walkers; ww++)
      m_model_parameters->set_chain_value(pp, 0, ww, chain_value[pp][ww]);
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
      for (size_t ii=0; ii<pp_single.size(); ii++)
	  pp[m_parameter_indexes[N][ii]] = pp_single[ii]; // this is necessary in presence of derived parameters
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
	for (size_t ii=0; ii<pp_single.size(); ii++)
	  pp[m_parameter_indexes[N][ii]] = pp_single[ii]; // this is necessary in presence of derived parameters
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

  return -0.5*LogLikelihood - log(sqrt(pow(2*cbl::par::pi, computed_model.size()) * cov.determinant()));
}


// ============================================================================================


double cbl::statistics::LogLikelihood_Poissonian_combined (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
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
  
  // ----- estimate the Poissonian log-likelihood -----    
  double LogLikelihood = 0.;
  for (size_t i=0; i<pp->flat_data.size(); i++)
    LogLikelihood += pp->flat_data[i]*cbl::Ln(computed_model[i],1.e-50)-computed_model[i]-gsl_sf_lnfact(int(pp->flat_data[i]));

  return LogLikelihood;
}


// ============================================================================================


double cbl::statistics::LogLikelihood_Poissonian_SSC_combined (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters -----    
  shared_ptr<statistics::STR_DependentProbes_data_model> pp = static_pointer_cast<statistics::STR_DependentProbes_data_model>(fixed_parameter);

  // ----- define the covariance matrix -----
  data::CovarianceMatrix cov = *pp->covariance;

  // ----- define the Sij matrix for the super-sample covariance -----
  std::vector<double> cosmoPar (pp->cosmoPar_indexes.size(), 0.);
  for (size_t i=0; i<cosmoPar.size(); i++)
    cosmoPar[i] = likelihood_parameter[pp->cosmoPar_indexes[i]];  
  
  std::vector<std::vector<double>> Sij = cov.SSC()->operator()(cosmoPar);
  vector<vector<double>> invSij; cbl::invert_matrix(Sij, invSij);
  const double detSij = cbl::determinant_matrix(Sij);

  // ----- for the super-sample covariance, extract the values of delta_b, i.e. the background matter contrast variation -----
  srand(time(0));
  std::vector<double> delta_b (pp->models.size(), 0.);
  for (size_t i=0; i<pp->models.size(); i++) {
    cbl::random::NormalRandomNumbers extraction (0, sqrt(Sij[i][i]), rand(), -sqrt(Sij[i][i]), sqrt(Sij[i][i]));
    delta_b[i] = extraction();
  }

  // ----- compute the model values -----   
  vector<double> computed_model, computed_response;
  vector<double> delta_b_index;
  for (size_t i=0; i<pp->models.size(); i++) {
    std::vector<double> single_par ((int)(pp->par_indexes[i].size()));
    for (size_t jj=0; jj<single_par.size(); jj++)
      single_par[jj] = likelihood_parameter[pp->par_indexes[i][jj]];
    std::vector<double> single_probe_model = pp->models[i]->operator()(pp->xx[i], single_par);
    std::vector<double> single_probe_response = cov.SSC()->get_response(i, pp->xx[i], single_par);
    for (size_t kk=0; kk<single_probe_model.size(); kk++) {
      computed_model.emplace_back(single_probe_model[kk]);
      computed_response.emplace_back(single_probe_response[kk]);
      delta_b_index.emplace_back(i);
    }
  }
  
  // ----- estimate the log-likelihood -----    
  double LogLikelihood = 0.;
  
  for (size_t i=0; i<pp->flat_data.size(); i++)
    LogLikelihood += pp->flat_data[i]*cbl::Ln(computed_model[i]+delta_b[delta_b_index[i]]*computed_response[i],1.e-50)-(computed_model[i]+delta_b[delta_b_index[i]]*computed_response[i])-gsl_sf_lnfact(int(pp->flat_data[i]));

  for (size_t i=0; i<delta_b.size(); i++)
    for (size_t j=0; j<delta_b.size(); j++)
      LogLikelihood += -0.5 * delta_b[i]*invSij[i][j]*delta_b[j];
  
  LogLikelihood = LogLikelihood - log( sqrt(detSij*pow(2.*cbl::par::pi, delta_b.size())) );

  return LogLikelihood;
}
