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

/**
 *  @file Statistics/LikelihoodParameters.cpp
 *
 *  @brief Methods of the class LikelihoodParameters 
 *
 *  This file contains the implementation of the methods of the class
 *  LikelihoodParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "LikelihoodParameters.h"

using namespace cosmobl;


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::full_parameters (const vector<double> free_parameter_values) const
{

  if((int)free_parameter_values.size() == (int)m_nparameters_free){
    vector<double> all_parameters(m_nparameters, 0);

    for(int i=0; i<m_nparameters_free; i++)
      all_parameters[m_free_parameters[i]] = free_parameter_values[i];

    for(int i=0; i<m_nparameters_fixed; i++)
      all_parameters[m_fixed_parameters[i]] = parameter(m_fixed_parameters[i])->value();

    for(int i=0; i<m_nparameters_output; i++)
      all_parameters[m_output_parameters[i]] = 0.;

    return all_parameters;
  }
  else if ((int)free_parameter_values.size() == (int)m_nparameters)
    return free_parameter_values;
  else
    ErrorCBL("Error in parameters of LikelihoodParameters.cpp, the vector of values of free parameters has the wrong size!");

  vector<double> vv;
  return vv;
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::m_set_parameter_type (vector<shared_ptr<Parameter>> parameters)
{
  m_parameters = parameters;

  m_nparameters = m_parameters.size();
  for (int i=0; i<m_nparameters; i++) {
    switch (m_parameters[i]->parameterType()) {
      case statistics::ParameterType::_BaseParameter_:
	if (m_parameters[i]->fixed()) {
	  m_nparameters_fixed +=1;
	  m_fixed_parameters.push_back(i);
	}
	else {
	  m_nparameters_free +=1;
	  m_free_parameters.push_back(i);
	}
	break;

      case statistics::ParameterType::_DerivedParameter_:
	m_nparameters_output += 1;
	m_output_parameters.push_back(i);
	break;

      default:
	ErrorCBL("Error in cosmobl::statistics::LikelihoodParameters of LikelihoodParameters.cpp: no such kind of parameter!");
    }

	
  }
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::m_set_parameter_posterior (const int start, const int thin, const int seed)
{
  for (int i=0; i<m_nparameters; i++) 
    m_parameters[i]->set_posterior(start, thin, seed);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::m_set_parameter_covariance(const int start, const int thin)
{
  vector<vector<double>> chains;

  for (int i=0; i<m_nparameters; i++)
    chains.push_back(m_parameters[i]->chain_values(start, thin));

  cosmobl::covariance_matrix(transpose(chains), m_parameter_covariance, false);
}


// ============================================================================================


cosmobl::statistics::LikelihoodParameters::LikelihoodParameters (const vector<shared_ptr<Parameter>> parameters)
{
  set_parameters(parameters);
}


// ============================================================================================

  
int cosmobl::statistics::LikelihoodParameters::nparameters () const
{
  return m_nparameters;
}


// ============================================================================================


int cosmobl::statistics::LikelihoodParameters::nparameters_free () const
{
  return m_nparameters_free;
}


// ============================================================================================


int cosmobl::statistics::LikelihoodParameters::nparameters_fixed () const
{
  return m_nparameters_fixed;
}


// ============================================================================================


int cosmobl::statistics::LikelihoodParameters::nparameters_output () const
{
  return m_nparameters_output;
}


// ============================================================================================


shared_ptr<cosmobl::statistics::Parameter> cosmobl::statistics::LikelihoodParameters::parameter (const unsigned int i) const
{
  return m_parameters[i];
}



// ============================================================================================


vector<shared_ptr<cosmobl::statistics::Parameter>> cosmobl::statistics::LikelihoodParameters::parameters () const
{
  return m_parameters;
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_parameters (const vector<shared_ptr<cosmobl::statistics::Parameter>> parameters)
{
  m_set_parameter_type(parameters);
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::parameter_covariance (const int i, const int j) const
{
  return m_parameter_covariance[i][j];
}


// ============================================================================================


vector<vector<double>> cosmobl::statistics::LikelihoodParameters::parameter_covariance () const 
{
  return m_parameter_covariance;
}


// ============================================================================================


int cosmobl::statistics::LikelihoodParameters::chain_size () const
{
  return m_chain_size;
}


// ============================================================================================


int cosmobl::statistics::LikelihoodParameters::nwalkers () const
{
  return m_nwalkers;
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::free (const int p)
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      m_parameters[p]->free();
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::fix (const int p)
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      m_parameters[p]->fix();
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::fix (const int p, const double value)
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      m_parameters[p]->fix(value);
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::fix_at_bestfit (const int p)
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      m_parameters[p]->fix_at_bestfit();
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }

}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::PriorProbability (const int p, const double value) const
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      m_parameters[p]->PriorProbability(value);
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
  return 0;
}



// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::PriorProbability (const vector<double> value) const
{
  vector<double> prior(m_nparameters, 1);

  for (int i=0; i<m_nparameters_free; i++) {
    int pp = m_free_parameters[i];
    prior[pp] = m_parameters[pp]->PriorProbability(value[pp]);
  }

  return prior;
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::LogPriorProbability (const int p, const double value) const
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      m_parameters[p]->LogPriorProbability(value);
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
  return 0;
}



// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::LogPriorProbability (const vector<double> value) const
{
  vector<double> prior(m_nparameters, 1);

  for (int i=0; i<m_nparameters_free; i++) {
    int pp = m_free_parameters[i];
    prior[pp] = m_parameters[pp]->LogPriorProbability(value[pp]);
  }

  return prior;
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_prior_seed (const int seed)
{
  cosmobl::random::UniformRandomNumbers ran(0, 321434523, seed);
  
  for (int i=0; i<m_nparameters_free; i++) {
    int pp = m_free_parameters[i];
    m_parameters[pp]->set_prior_seed(int(ran()));
  }
  
  for (int i=0; i<m_nparameters_fixed; i++) {
    int pp = m_fixed_parameters[i];
    m_parameters[pp]->set_prior_seed(int(ran()));
  }


}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::prior_sample (const int p)
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      return m_parameters[p]->prior_sample();
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }

  return 0;
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::prior_sample (const int p, const int sample_size)
{
  switch (m_parameters[p]->parameterType()) {

    case statistics::ParameterType::_BaseParameter_:
      return m_parameters[p]->prior_sample(sample_size);
      break;

    case statistics::ParameterType::_DerivedParameter_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameters[p]->name()+" is an output parameter");
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }

  vector<double> vv;
  return vv;
}

// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::prior_sample ()
{
  vector<double> prior_values(nparameters(), 0);
  
  for (int i=0; i<m_nparameters_free; i++) {
    int pp = m_free_parameters[i];
    prior_values[pp] = m_parameters[pp]->prior_sample();
  }
  
  for (int i=0; i<m_nparameters_fixed; i++) {
    int pp = m_fixed_parameters[i];
    prior_values[pp] = m_parameters[pp]->prior_sample();
  }

  return prior_values;
}


// ============================================================================================
	

double cosmobl::statistics::LikelihoodParameters::prior_range (const int p, const double epsilon)
{
  switch (m_parameters[p]->parameterType()) {
    case statistics::ParameterType::_BaseParameter_:
      if (m_parameters[p]->fixed()) {
	return 0.;
      }
      else {
	return m_parameters[p]->prior_range(epsilon); 
      }
      break;

    case statistics::ParameterType::_DerivedParameter_:
      return 0;
      break;

    default:
      ErrorCBL("Error in cosmobl::statistics::LikelihoodParameters of LikelihoodParameters.cpp: no such kind of parameter!");
  }
  return 0;
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::prior_range (const double epsilon)
{
  vector<double> ss(m_nparameters, 0);
  for (int i=0; i<m_nparameters_free; i++)
    ss[m_free_parameters[i]] = prior_range(m_free_parameters[i], epsilon);

  return ss;
}


// ============================================================================================

  
double cosmobl::statistics::LikelihoodParameters::bestfit_value (const int p) const
{
  return m_parameters[p]->bestfit_value();
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::bestfit_values () const
{
  vector<double> bestfit(nparameters(),0);

  for (int i=0; i<m_nparameters; i++)
    bestfit[i] = bestfit_value(i);

  return bestfit;
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_bestfit_value (const int p, const double bestfit_value)
{
  m_parameters[p]->set_bestfit_value(bestfit_value);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_bestfit_value (const vector<double> bestfit_value)
{
  for (int i=0; i<m_nparameters; i++)
    m_parameters[i]->set_bestfit_value(bestfit_value[i]);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::write_bestfit_info (const vector<double> bestfit_value)
{
  set_bestfit_value(bestfit_value);
  for (int i=0; i<m_nparameters; i++) {

    switch (m_parameters[i]->parameterType()) {

      case statistics::ParameterType::_BaseParameter_:
	if (m_parameters[i]->fixed()) {
	  coutCBL << "Parameter: " << par::col_yellow << m_parameters[i]->name() << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	  coutCBL << "value = " << m_parameters[i]->bestfit_value() << endl;
	  cout << endl;
	}
	else {
	  coutCBL << "Parameter: " << par::col_yellow << m_parameters[i]->name() << par::col_default << " --> status: " << par::col_green << "FREE" << endl;
	  coutCBL << "value = " << m_parameters[i]->bestfit_value() << endl;
	  cout << endl;
	}
	
	break;

      case statistics::ParameterType::_DerivedParameter_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameters[i]->name() << par::col_default << " --> status: " << par::col_bred << "OUTPUT" << endl;
	coutCBL << "value = " << m_parameters[i]->bestfit_value() << endl;
	cout << endl;
	break;

      default:
	ErrorCBL("Error in cosmobl::statistics::LikelihoodParameters of LikelihoodParameters.cpp: no such kind of parameter!");
    }
  }

}


// ============================================================================================

  
double cosmobl::statistics::LikelihoodParameters::PosteriorProbability (const int p, const double value) const
{
  return m_parameters[p]->PosteriorProbability(value);
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::PosteriorProbability (const vector<double> value) const
{
  vector<double> vv;

  for (int i=0; i<m_nparameters; i++)
    vv.push_back(PosteriorProbability(i, value[i]));

  return vv;
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::posterior_mean (const int p) const
{
  return m_parameters[p]->posterior_mean();
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::posterior_mean () const
{
  vector<double> vv;

  for (int i=0; i<m_nparameters; i++)
    vv.push_back(posterior_mean(i));

  return vv;
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::posterior_median (const int p) const
{
  return m_parameters[p]->posterior_median();
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::posterior_median () const
{ 
  vector<double> vv;

  for (int i=0; i<m_nparameters; i++)
    vv.push_back(posterior_median(i));

  return vv;
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::posterior_std (const int p) const
{
  return m_parameters[p]->posterior_std();
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::posterior_std () const
{
  vector<double> vv;

  for (int i=0; i<m_nparameters; i++)
    vv.push_back(posterior_std(i));

  return vv;
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::posterior_percentile (const int p, const unsigned int i) const
{
  return m_parameters[p]->posterior_percentile(i);
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::posterior_percentile (const unsigned int i) const
{
  vector<double> vv;

  for (int p=0; p<m_nparameters; p++)
    vv.push_back(posterior_percentile(p, i));

  return vv;
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::posterior_sample (const int p)
{
  return m_parameters[p]->posterior_sample();
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::posterior_sample (const int p, const int sample_size)
{
  return m_parameters[p]->posterior_sample(sample_size);
}


// ============================================================================================


double cosmobl::statistics::LikelihoodParameters::chain_value (const int pos, const int ww, const int par)
{
  return m_parameters[par]->chain_value(pos, ww);
}


// ============================================================================================


vector<double> cosmobl::statistics::LikelihoodParameters::chain_values (const int pos, const int par)
{
  vector<double> vv(nwalkers(),0);
  for (int i=0; i<nwalkers(); i++)
    vv[i] = m_parameters[par]->chain_value(pos, i);
  return vv;
}


// ============================================================================================


vector<vector<double>> cosmobl::statistics::LikelihoodParameters::chain_values (const int pos)
{
  vector<vector<double>> vv;

  for (int i=0;i<m_nparameters;i++)
    vv.push_back(chain_values(pos, i));

  return cosmobl::transpose(vv);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_chain_value (const int pos, const int ww, const int par, const double value)
{
  m_parameters[par]->set_chain_value(pos, ww, value);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_chain_values (const int pos, const int par, const vector<double> values)
{
  for (size_t i=0;i<values.size();i++)
    m_parameters[par]->set_chain_value(pos, i, values[i]);
}



// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_chain_values (const int pos, const vector<vector<double>> values)
{
  for (int p=0;p<m_nparameters;p++)
      set_chain_values(pos, p, values[p]);

}
  

// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_chains (const int chain_size, const int nwalkers)
{
  for (int p=0;p<m_nparameters;p++)
    m_parameters[p]->set_chain(chain_size, nwalkers);

  m_nwalkers = nwalkers;
  m_chain_size = chain_size;
}
  

// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_chains (const vector<vector<double>> chain_values, const int nwalkers)
{
  for (int p=0;p<m_nparameters;p++)
    m_parameters[p]->set_chain(chain_values[p], nwalkers);

  m_nwalkers = nwalkers;
  m_chain_size = m_parameters[0]->chain_size();
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::set_chains (const vector<vector<vector<double>>> chain_values)
{
  for (int p=0;p<m_nparameters;p++)
    m_parameters[p]->set_chain(chain_values[p]);

  m_nwalkers = m_parameters[0]->nwalkers();
  m_chain_size = m_parameters[0]->chain_size();
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::initialize_chains_from_prior (const int seed)
{
  set_prior_seed(seed);

  for (int p=0;p<nparameters_free(); p++) {
    int ii = m_free_parameters[p];
    vector<double> vv(nwalkers(), 0);

    if (m_parameters[ii]->parameterType()==statistics::ParameterType::_BaseParameter_) {
      vv = m_parameters[ii]->prior_sample(nwalkers());
    }
    set_chain_values(0, ii, vv);
  }

  for (int p=0;p<nparameters_fixed(); p++) {
    int ii = m_fixed_parameters[p];
    vector<double> vv(nwalkers(), m_parameters[ii]->value());
    set_chain_values(0, ii, vv);
  }
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::initialize_chains_around_bestfit_values (const double radius, const int seed)
{
  initialize_chains_around_values(bestfit_values(), radius, seed);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::initialize_chains_around_values (const vector<double> values, const double radius, const int seed)
{
  if ((int)values.size()!=nparameters())
    ErrorCBL("Error in initialize_chains_around_values of LikelihoodParameeres, starting values must have the same dimension of the number of parameters");

   cosmobl::random::NormalRandomNumbers rn(0, radius, seed);

  for (int p=0;p<nparameters_free(); p++) {
    int ii = m_free_parameters[p];
    vector<double> vv(nwalkers(), 0);

    if (m_parameters[ii]->parameterType()==statistics::ParameterType::_BaseParameter_) {
      for (int i=0;i<nwalkers(); i++){
	bool ok = false;
	while(!ok){
	  vv[i] = rn()+values[ii];
	  ok = m_parameters[ii]->prior()->isIncluded (vv[i]);
	}
      }
    }
    set_chain_values(0, ii, vv);
  }

  for (int p=0;p<nparameters_fixed(); p++) {
    int ii = m_fixed_parameters[p];
    vector<double> vv(nwalkers(), m_parameters[ii]->value());
    set_chain_values(0, ii, vv);
  }
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::show_results (const int start, const int thin, const int seed)
{
  m_set_parameter_posterior(start, thin, seed);

  m_set_parameter_covariance(start, thin);

  const int dp = cout.precision(); cout.precision(4);
  cout << endl;
  
  for (int i=0; i<m_nparameters; i++) {

    switch (m_parameters[i]->parameterType()) {

      case statistics::ParameterType::_BaseParameter_:
	if (m_parameters[i]->fixed()) {
	  coutCBL << "Parameter: " << par::col_yellow << m_parameters[i]->name() << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	  coutCBL << "value = " << m_parameters[i]->value() << endl;
	  cout << endl;
	}
	else {
	  coutCBL << "Parameter: " << par::col_yellow << m_parameters[i]->name() << par::col_default << " --> status: " << par::col_green << "FREE" << endl;
	  coutCBL << "Posterior mean = " << m_parameters[i]->posterior_mean() << endl;
	  coutCBL << "Posterior standard deviation = " << m_parameters[i]->posterior_std() << endl;
	  coutCBL << "Posterior median = " << m_parameters[i]->posterior_median() << endl;
	  coutCBL << "Posterior 18th percentile = " << m_parameters[i]->posterior_median()-m_parameters[i]->posterior_percentile(18) << endl;
	  coutCBL << "Posterior 82th percentile = " << m_parameters[i]->posterior_percentile(82)-m_parameters[i]->posterior_median() << endl;
	  coutCBL << "Posterior mode = " << m_parameters[i]->posterior_mode() << endl;
	  cout << endl;
	}
	break;

      case statistics::ParameterType::_DerivedParameter_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameters[i]->name() << par::col_default << " --> status: " << par::col_bred << "OUTPUT" << endl;
	coutCBL << "Posterior mean = " << m_parameters[i]->posterior_mean() << endl;
	coutCBL << "Posterior standard deviation = " << m_parameters[i]->posterior_std() << endl;
	coutCBL << "Posterior median = " << m_parameters[i]->posterior_median() << endl;
	coutCBL << "Posterior 18th percentile = " << m_parameters[i]->posterior_median()-m_parameters[i]->posterior_percentile(18) << endl;
	coutCBL << "Posterior 82th percentile = " << m_parameters[i]->posterior_percentile(82)-m_parameters[i]->posterior_median() << endl;
	coutCBL << "Posterior mode = " << m_parameters[i]->posterior_mode() << endl;
	coutCBL << endl;
	break;

      default:
	ErrorCBL("Error in cosmobl::statistics::LikelihoodParameters of LikelihoodParameters.cpp: no such kind of parameter!");
    }
  }
  
  cout.precision(dp);
}


// ============================================================================================


void cosmobl::statistics::LikelihoodParameters::write_results (const string dir, const string file, const int start, const int thin, const int seed)
{
  m_set_parameter_posterior(start, thin, seed);
  m_set_parameter_covariance(start, thin);
  
  string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}

  string file_parameters = dir+file+"_parameters.dat";
  string file_covariance = dir+file+"_covariance.dat";

  ofstream fout(file_parameters.c_str());

  fout << "### Parameter # status ###" << endl << "### Posterior mean # Posterior standard deviation # Posterior median # Posterior 18th percentile # Posterior 82th percentile # Posterior mode ###" << endl << endl;
  
  for (int i=0; i<m_nparameters; i++) 
    if (m_parameters[i]->parameterType()==statistics::ParameterType::_BaseParameter_ && m_parameters[i]->fixed())
      fout << "# " << m_parameters[i]->name() << " FIXED #" << endl;
    else 
      fout << "# " << m_parameters[i]->name() << " FREE #" << endl;
   
  for (int i=0; i<m_nparameters; i++) 
    if (m_parameters[i]->parameterType()==statistics::ParameterType::_BaseParameter_ && m_parameters[i]->fixed())
      fout << m_parameters[i]->value() << " 0 0 0 0 0 0" << endl;
    else 
      fout << m_parameters[i]->posterior_mean() << " " << m_parameters[i]->posterior_std() << " " << m_parameters[i]->posterior_median() << " " << m_parameters[i]->posterior_median()-m_parameters[i]->posterior_percentile(18) << " " << m_parameters[i]->posterior_percentile(82)-m_parameters[i]->posterior_median() << " " << m_parameters[i]->posterior_mode() << endl;
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file_parameters << endl;


  fout.open(file_covariance.c_str());
  for (int i=0; i<m_nparameters; i++) {
    for (int j=0; j<m_nparameters; j++)
      fout << i << " " << j << " " << m_parameter_covariance[i][j] << endl;
    fout << endl;
  }
  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file_covariance << endl;

}
