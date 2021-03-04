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
 *  @author Giorgio Lesci (and Federico Marulli)
 *
 *  @author giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "CombinedModelling.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::CombinedModelling::CombinedModelling (std::vector<std::shared_ptr<modelling::Modelling>> modelling)
{
  m_Nmodellings = modelling.size();
  m_modelling = modelling;
  m_set_combined_posterior();
}


// ============================================================================================


void cbl::modelling::CombinedModelling::m_set_combined_posterior ()
{
  std::vector<std::shared_ptr<statistics::Posterior>> posteriors;
  for (int i=0; i<m_Nmodellings; i++){
    m_modelling[i]->m_set_posterior(321);
    auto posterior_ptr = m_modelling[i]->posterior();
    posteriors.push_back(std::move(posterior_ptr));
  }  
  m_combined_posterior = make_shared<statistics::CombinedPosterior>(statistics::CombinedPosterior(posteriors));
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


void cbl::modelling::CombinedModelling::write_combined_results (const std::string output_dir, const std::string root_file, const int start, const int thin, const int nbins, const bool fits, const bool compute_mode, const int ns)
{
  m_combined_posterior->write_results(output_dir, root_file, start, thin, nbins, fits, compute_mode, ns);
}
