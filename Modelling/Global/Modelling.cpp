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

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::Modelling::set_parameters (const vector<shared_ptr<cosmobl::statistics::Parameter>> parameters)
{
  m_parameters = make_shared<cosmobl::statistics::LikelihoodParameters>(cosmobl::statistics::LikelihoodParameters(parameters));
}


// ============================================================================================


void cosmobl::modelling::Modelling::set_likelihood (const statistics::LikelihoodType likelihood_type)
{
  if(m_fit_range)
    m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data_fit, m_model, m_parameters, likelihood_type));
  else {
    m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, m_parameters, likelihood_type));
  }
}


// ============================================================================================


void cosmobl::modelling::Modelling::set_likelihood (const statistics::LogLikelihood_function likelihood_func)
{
  if (m_fit_range)
    m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data_fit, m_model, m_parameters, likelihood_func));
  else
    m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, m_parameters, likelihood_func));
}


// ============================================================================================


void cosmobl::modelling::Modelling::maximize_likelihood (vector<double> &guess, const int ntry, const int prior_seed, const bool usePrior, const unsigned int max_iter, const double tol, const double epsilon)
{
  m_likelihood->maximize(guess, ntry, prior_seed, usePrior, max_iter, tol, epsilon);
}


// ============================================================================================


void cosmobl::modelling::Modelling::maximize_likelihood (vector<double> &guess, const bool usePrior, const unsigned int max_iter, const double tol, const double epsilon)
{
  m_likelihood->maximize(guess, usePrior, max_iter, tol, epsilon);
}


// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const int chain_size, const int nwalkers, const int seed, const double aa)
{
  m_likelihood->initialize_chains(chain_size, nwalkers, seed);

  m_likelihood->sample_stretch_move_parallel(seed, aa);
}


// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const int chain_size, const int nwalkers, const int seed, vector<double> &start, const double radius, const double aa)
{
  m_likelihood->initialize_chains(chain_size, nwalkers, seed, start, radius);

  m_likelihood->sample_stretch_move_parallel(seed, aa);
}

// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const int chain_size, const vector<vector<double>> chain_values, const int seed, const double aa)
{
  m_likelihood->initialize_chains(chain_size, chain_values);

  m_likelihood->sample_stretch_move_parallel(seed, aa);
}

// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const int chain_size, const int nwalkers, const int seed, const string input_dir, const string input_file, const double aa)
{
  m_likelihood->initialize_chains(chain_size, nwalkers, input_dir, input_file);

  m_likelihood->sample_stretch_move_parallel(seed, aa);
}


// ============================================================================================


void cosmobl::modelling::Modelling::show_results (const int start, const int thin, const int seed)
{
  m_parameters->show_results(start, thin, seed);
}

// ============================================================================================


void cosmobl::modelling::Modelling::write_results (const string dir, const string file, const int start, const int thin, const int seed)
{
  m_likelihood->write_results(dir, file, start, thin, seed);
}

// ============================================================================================


void cosmobl::modelling::Modelling::read_chain (const string dir, const string file, const int nwalkers, const int skip_header)
{
  m_likelihood->read_chain(dir, file, nwalkers, skip_header);
}
