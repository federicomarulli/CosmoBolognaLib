/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Modelling/Modelling.cpp
 *
 *  @brief Methods of the class Modelling, used for modelling any kind
 *  of measurements
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling.h"

using namespace cosmobl;


// ============================================================================================


void modelling::Modelling::sample_likelihood(const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_data->set_limits(xmin, xmax);

  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type));

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);

}


// ============================================================================================


void modelling::Modelling::sample_likelihood(const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type));

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);

}
// ============================================================================================



void modelling::Modelling::sample_likelihood(const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_data->set_limits(xmin, xmax);

  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type, loglikelihood_function, cov));

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);
}


// ============================================================================================


void modelling::Modelling::sample_likelihood(const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type, loglikelihood_function, cov));

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);
}
