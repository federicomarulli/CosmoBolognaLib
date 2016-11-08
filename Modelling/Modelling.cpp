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
 *  @file Modelling/Modelling.cpp
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


void cosmobl::modelling::Modelling::set_likelihood (const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type)
{
  m_data->set_limits(xmin, xmax);

  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type));
}


// ============================================================================================


void cosmobl::modelling::Modelling::set_likelihood (const double xmin, const double xmax, const double ymin, const double ymax, const statistics::LikelihoodType likelihood_type)
{
  m_data->set_limits(xmin, xmax, ymin, ymax);

  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type));
}


// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  set_likelihood(xmin, xmax, likelihood_type);

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);
}


// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const double xmin, const double xmax, const double ymin, const double ymax, const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  set_likelihood(xmin, xmax, ymin, ymax, likelihood_type);

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);

}

// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type));

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);
}


// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_likelihood = make_shared<statistics::Likelihood> (statistics::Likelihood(m_data, m_model, likelihood_type, loglikelihood_function, cov));

  m_likelihood->sample_stretch_move(n_chains, chain_size, seed, do_write_chain, dir_output, chain_file);

  m_likelihood->write_chain(dir_output, chain_file, start, stop, thin);
}


// ============================================================================================


void cosmobl::modelling::Modelling::sample_likelihood (const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_data->set_limits(xmin, xmax);

  sample_likelihood(likelihood_type, loglikelihood_function, cov, n_chains, chain_size, seed, dir_output, chain_file, do_write_chain, start, stop, thin);
}


// ============================================================================================



void cosmobl::modelling::Modelling::sample_likelihood (const double xmin, const double xmax, const double ymin, const double ymax, const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain, const double start, const double stop, const int thin)
{
  m_data->set_limits(xmin, xmax, ymin, ymax);

  sample_likelihood(likelihood_type, loglikelihood_function, cov, n_chains, chain_size, seed, dir_output, chain_file, do_write_chain, start, stop, thin);
}


// ============================================================================================


void cosmobl::modelling::Modelling::write_parameters (const string dir, const string file) const
{
  string output = dir+file; 
  ofstream fout(output); checkIO(fout, output);

  for (unsigned int i=0; i<m_model->npar(); i++)
    fout << m_model->parameter(i)->value() << "   " << m_model->parameter(i)->std() << endl;
  
  fout.clear(); fout.close();
  coutCBL << "I wrote the file: " << output << endl;
}
