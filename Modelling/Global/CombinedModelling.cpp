/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Modelling/Global/CombinedModelling.cpp
 *
 *  @brief Methods of the class CombinedModelling, used for modelling any kind
 *  of combined modelling 
 *
 *  This file contains the implementation of the methods of the class
 *  CombinedModelling
 *
 *  @author Giorgio Lesci, Sofia Contarini (and Federico Marulli)
 *
 *  @author giorgio.lesci2@unibo.it, sofia.contarini3@unibo.it (and federico.marulli3@unibo.it)
 */


#include "CombinedModelling.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::CombinedModelling::CombinedModelling (std::vector<std::shared_ptr<modelling::Modelling>> modelling, std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par)
{
  std::vector<std::shared_ptr<statistics::Posterior>> posteriors(modelling.size());
  for (size_t i=0; i<modelling.size(); i++) {
    modelling[i]->m_set_posterior(321);
    auto posterior_ptr = modelling[i]->posterior();
    posteriors[i] = std::move(posterior_ptr);
  }  
  m_combined_posterior = make_shared<statistics::CombinedPosterior>(statistics::CombinedPosterior(posteriors, repeated_par, common_repeated_par));
}


// ============================================================================================


cbl::modelling::CombinedModelling::CombinedModelling (std::vector<std::vector<std::shared_ptr<modelling::Modelling>>> modelling, std::vector<std::shared_ptr<data::CovarianceMatrix>> covariance, const std::vector<cbl::statistics::LikelihoodType> likelihood_types, const std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par)
{
  std::vector<std::vector<std::shared_ptr<statistics::Posterior>>> posteriors(modelling.size());
  for (size_t i=0; i<modelling.size(); i++) {
    posteriors[i].resize(modelling[i].size());
    for (size_t j=0; j<modelling[i].size(); j++) {
      modelling[i][j]->m_set_posterior(321);
      auto posterior_ptr = modelling[i][j]->posterior();
      posteriors[i][j] = std::move(posterior_ptr);
    }
  }
  m_combined_posterior = make_shared<statistics::CombinedPosterior>(statistics::CombinedPosterior(posteriors, covariance, likelihood_types, repeated_par, common_repeated_par));
}


// ============================================================================================


void cbl::modelling::CombinedModelling::maximize_combined_posterior (const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon)
{
  m_combined_posterior->maximize(start, max_iter, tol, epsilon);
}


// ============================================================================================


void cbl::modelling::CombinedModelling::sample_combined_posterior (const int chain_size, const int nwalkers, const double aa, const bool parallel)
{
  m_combined_posterior->initialize_chains(chain_size, nwalkers);
  m_combined_posterior->sample_stretch_move(aa, parallel);
}


// ============================================================================================


void cbl::modelling::CombinedModelling::sample_combined_posterior (const int chain_size, const int n_walkers, const double radius, const std::vector<double> start, const unsigned int max_iter, const double tol, const double epsilon, const double aa, const bool parallel)
{
  m_combined_posterior->initialize_chains(chain_size, n_walkers, radius, start, max_iter, tol, epsilon);
  m_combined_posterior->sample_stretch_move(aa, parallel);
}


// ============================================================================================


void cbl::modelling::CombinedModelling::sample_combined_posterior (const int chain_size, const int nwalkers, const std::string input_dir, const std::string input_file, const int seed, const double aa, const bool parallel)
{
  m_combined_posterior->initialize_chains(chain_size, nwalkers, input_dir, input_file, seed);
  m_combined_posterior->sample_stretch_move(aa, parallel);
}


// ============================================================================================


void cbl::modelling::CombinedModelling::write_combined_results (const std::string output_dir, const std::string root_file, const int start, const int thin, const int nbins, const bool fits, const bool compute_mode, const int ns)
{
  m_combined_posterior->write_results(output_dir, root_file, start, thin, nbins, fits, compute_mode, ns);
}


// ============================================================================================


void cbl::modelling::CombinedModelling::write_model_from_combined_chain (const std::string output_dir, const std::string output_file, const int start, const int thin)
{
  m_combined_posterior->write_model_from_chain(output_dir, output_file, start, thin);
}
