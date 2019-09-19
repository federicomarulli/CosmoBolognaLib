/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Modelling/Global/Modelling.cpp
 *
 *  @brief Methods of the class Modelling, used for modelling any kind
 *  of measurements
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::modelling::Modelling::m_set_prior (vector<statistics::PriorDistribution> prior_distribution)
{
  // set the priors
  m_parameter_priors.erase(m_parameter_priors.begin(), m_parameter_priors.end());
  for (size_t i=0; i<prior_distribution.size(); i++)
    m_parameter_priors.emplace_back(make_shared<statistics::PriorDistribution>(prior_distribution[i]));
}


// ============================================================================================


void cbl::modelling::Modelling::m_set_posterior (const int seed)
{
  if (m_likelihood!=NULL && m_parameter_priors.size()==m_model->parameters()->nparameters_base())
    m_posterior = make_shared<statistics::Posterior>(statistics::Posterior(m_parameter_priors, *m_likelihood, seed));
  else
    ErrorCBL("either the posterior is not defined or a wrong number of prior distributions has been provided!", "m_set_posterior", "Modelling.cpp");
}

// ============================================================================================


shared_ptr<statistics::Likelihood> cbl::modelling::Modelling::likelihood ()
{
  if (m_likelihood!=NULL)
    return m_likelihood;
  else
    ErrorCBL("the likelihood is not defined!", "likelihood", "Modelling.cpp");
  return NULL;
}


// ============================================================================================


shared_ptr<cbl::statistics::Posterior> cbl::modelling::Modelling::posterior ()
{
  if (m_posterior!=NULL)
    return move(m_posterior);
  else
    ErrorCBL("the posterior is not defined!", "posterior", "Modelling.cpp");
  return NULL;
}


// ============================================================================================


shared_ptr<cbl::statistics::ModelParameters> cbl::modelling::Modelling::likelihood_parameters ()
{
  if (m_likelihood!=NULL)
    return m_likelihood->parameters();
  else
    ErrorCBL("the likelihood is not defined!", "likelihood_parameters", "Modelling.cpp");
  return NULL;
}


// ============================================================================================


shared_ptr<cbl::statistics::ModelParameters> cbl::modelling::Modelling::posterior_parameters ()
{
  if (m_posterior!=NULL)
    return m_posterior->parameters();
  else
    ErrorCBL("the posterior is not defined!", "posterior_parameters", "Modelling.cpp");
  return NULL;
}


// ============================================================================================


void cbl::modelling::Modelling::set_likelihood (const statistics::LikelihoodType likelihood_type, const vector<size_t> x_index, const int w_index)
{
  if (m_model==NULL)
    ErrorCBL("undefined  model!", "set_likelihood", "Modelling.cpp");

  if (m_fit_range) {
    if (m_data_fit==NULL)
      ErrorCBL("undefined fit range!", "set_likelihood", "Modelling.cpp");
    m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data_fit, m_model, likelihood_type, x_index, w_index));
  }
  else  {
    if (m_data == NULL)
      ErrorCBL("Error in set_likelihood of Modelling.cpp. Undefined dataset!", "set_likelihood", "Modelling.cpp");
    m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type, x_index, w_index));
  }
}


// ============================================================================================


void cbl::modelling::Modelling::maximize_likelihood (const vector<double> start, const vector<vector<double>> parameter_ranges, const unsigned int max_iter, const double tol, const double epsilon)
{
  m_likelihood->maximize(start, parameter_ranges, max_iter, tol, epsilon);
}


// ============================================================================================


void cbl::modelling::Modelling::maximize_posterior (const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon, const int seed)
{
  m_set_posterior(seed);
  m_posterior->maximize(start, max_iter, tol, epsilon);
}


// ============================================================================================


void cbl::modelling::Modelling::sample_posterior (const int chain_size, const int nwalkers, const int seed, const double aa, const bool parallel)
{
  m_set_posterior(seed);
  m_posterior->initialize_chains(chain_size, nwalkers);
  m_posterior->sample_stretch_move(aa, parallel);
}

// ============================================================================================


void cbl::modelling::Modelling::sample_posterior (const int chain_size, const int nwalkers, const double radius, const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon, const int seed, const double aa, const bool parallel)
{ 
  m_set_posterior(seed);
  m_posterior->initialize_chains(chain_size, nwalkers, radius, start, max_iter, tol, epsilon);
  m_posterior->sample_stretch_move(aa, parallel);
}

// ============================================================================================


void cbl::modelling::Modelling::sample_posterior (const int chain_size, const int nwalkers, std::vector<double> &value, const double radius, const int seed, const double aa, const bool parallel)
{
  m_set_posterior(seed);
  m_posterior->initialize_chains(chain_size, nwalkers, value, radius);
  m_posterior->sample_stretch_move(aa, parallel);
}


// ============================================================================================


void cbl::modelling::Modelling::sample_posterior (const int chain_size, const std::vector<std::vector<double>> chain_value, const int seed, const double aa, const bool parallel)
{
  m_set_posterior(seed);
  m_posterior->initialize_chains(chain_size, chain_value);
  m_posterior->sample_stretch_move(aa, parallel);
}

      
// ============================================================================================


void cbl::modelling::Modelling::sample_posterior (const int chain_size, const int nwalkers, const std::string input_dir, const std::string input_file, const int seed, const double aa, const bool parallel)
{
  m_set_posterior(seed);
  m_posterior->initialize_chains(chain_size, nwalkers, input_dir, input_file);
  m_posterior->sample_stretch_move(aa, parallel);
}


// ============================================================================================
      

void cbl::modelling::Modelling::write_chain (const string output_dir, const string output_file, const int start, const int thin, const bool fits)
{
  m_posterior->write_chain(output_dir, output_file, start, thin, fits);
}


// ============================================================================================


void cbl::modelling::Modelling::read_chain (const string input_dir, const string input_file, const int nwalkers, const int skip_header, const bool fits)
{
  m_posterior->read_chain(input_dir, input_file, nwalkers, skip_header, fits);
}


// ============================================================================================


void cbl::modelling::Modelling::show_results (const int start, const int thin, const int nbins, const bool show_mode)
{
  m_posterior->show_results(start, thin, nbins, show_mode);
}

// ============================================================================================


void cbl::modelling::Modelling::write_results (const string dir, const string file, const int start, const int thin, const int nbins, const bool fits, const bool compute_mode)
{
  m_posterior->write_results(dir, file, start, thin, nbins, fits, compute_mode);
}

