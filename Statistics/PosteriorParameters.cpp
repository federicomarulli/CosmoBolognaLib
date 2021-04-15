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
 *  @file Statistics/PosteriorParameters.cpp
 *
 *  @brief Methods of the class PosteriorParameters
 *
 *  This file contains the implementation of the methods of the class
 *  PosteriorParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "PosteriorParameters.h"

using namespace std;

using namespace cbl;
using namespace statistics;


// ============================================================================================


string cbl::statistics::PosteriorParameters::status (const int p) const
{
  string stat;

  switch (m_parameter_type[p]) {
  case statistics::ParameterType::_Base_:
    if (m_parameter_prior[p]->distributionType()==cbl::glob::DistributionType::_Constant_)
      stat = "FIXED";
    else
      stat = "FREE";
    break;

  case statistics::ParameterType::_Correlated_:
    stat = "FREE";
    break;

  case statistics::ParameterType::_Derived_:
    stat = "OUTPUT";
    break;

  default:
    ErrorCBL("no such kind of parameter!", "status", "PosteriorParameters.cpp");
  }

  return stat;
}


// ============================================================================================


vector<string> cbl::statistics::PosteriorParameters::status () const
{
  vector<string> stat;
  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
    case statistics::ParameterType::_Base_:
      if (m_parameter_prior[i]->distributionType()==cbl::glob::DistributionType::_Constant_)
	stat.push_back("FIXED");
      else
	stat.push_back("FREE");
      break;

    case statistics::ParameterType::_Correlated_:
      stat.push_back("FREE");
      break;

    case statistics::ParameterType::_Derived_:
      stat.push_back("OUTPUT");
      break;

    default:
      ErrorCBL("no such kind of parameter!", "status", "PosteriorParameters.cpp");
    }
  }
  return stat;
}
// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::full_parameter (const std::vector<double> parameter_value) const
{
  if (parameter_value.size()==m_nparameters_free) {
    vector<double> all_parameters(m_nparameters, 0);

    for (size_t i=0; i<m_nparameters_free; i++)
      all_parameters[m_free_parameter[i]] = parameter_value[i];

    for (size_t i=0; i<m_nparameters_fixed; i++)
      all_parameters[m_fixed_parameter[i]] = m_parameter_prior[m_fixed_parameter[i]]->sample();

    for (size_t i=0; i<m_nparameters_derived; i++)
      all_parameters[m_derived_parameter[i]] = 0.;

    return all_parameters;
  }
  else if (parameter_value.size() == m_nparameters) {
    vector<double> all_parameters = parameter_value;

    for (size_t i=0; i<m_nparameters_fixed; i++)
      all_parameters[m_fixed_parameter[i]] = m_parameter_prior[m_fixed_parameter[i]]->sample();

    return all_parameters;
  }
  else
    ErrorCBL("the vector of free parameters has the wrong size!", "full_parameter", "PosteriorParameters.cpp");

  vector<double> vv;
  return vv;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::m_set_parameter_type ()
{
  m_nparameters_free = 0;
  m_nparameters_fixed = 0;
  m_nparameters_base = 0;
  m_nparameters_correlated = 0;
  m_nparameters_derived = 0;

  m_base_parameter.erase(m_base_parameter.begin(), m_base_parameter.end());
  m_fixed_parameter.erase(m_fixed_parameter.begin(), m_fixed_parameter.end());
  m_free_parameter.erase(m_free_parameter.begin(), m_free_parameter.end());
  m_derived_parameter.erase(m_derived_parameter.begin(), m_derived_parameter.end());

  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
    case statistics::ParameterType::_Base_:
      if (m_parameter_prior[i] != NULL) {
	if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_) {
	  m_nparameters_fixed ++;
	  m_fixed_parameter.push_back(i);
	}
	else {
	  m_nparameters_free ++;
	  m_free_parameter.push_back(i);
	}
	m_nparameters_base ++;
	m_base_parameter.push_back(i);
      }
      break;

    case statistics::ParameterType::_Correlated_:
      m_nparameters_free ++;
      m_free_parameter.push_back(i);
      m_nparameters_correlated ++;
      m_nparameters_base ++;
      m_base_parameter.push_back(i);
      break;

    case statistics::ParameterType::_Derived_:
      m_nparameters_derived ++;
      m_derived_parameter.push_back(i);
      break;

    default:
      ErrorCBL("no such kind of parameter!", "m_set_parameter_type", "PosteriorParameters.cpp");
    }
  }
}


// ============================================================================================


cbl::statistics::PosteriorParameters::PosteriorParameters (const size_t nparameters, const std::vector<std::shared_ptr<PriorDistribution>> priorDistribution, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
{
  set_parameters(nparameters, priorDistribution, parameterTypes, parameterNames);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_parameters (const size_t nparameters, const std::vector<std::shared_ptr<PriorDistribution>> priorDistribution, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
{
  // check parameterTypes size

  if (nparameters==0)
    ErrorCBL("nparameters must be >0!", "set_parameters", "PosteriorParameters.cpp");

  if ((parameterTypes.size()!=nparameters) && (parameterTypes.size()!=0))
    ErrorCBL("the size of parameterTypes is incorrect!", "set_parameters", "PosteriorParameters.cpp");

  if ((parameterNames.size()!=nparameters) && (parameterNames.size()!=0))
    ErrorCBL("the size of parameterNames is incorrect!", "set_parameters", "PosteriorParameters.cpp");


  if ((parameterTypes.size()==nparameters) && (parameterNames.size()==nparameters)) {
    m_nparameters = nparameters;
    m_parameter_type = parameterTypes;
    m_parameter_name = parameterNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)) {
    m_nparameters = nparameters;
    vector<ParameterType> pTypes(m_nparameters);
    vector<string> pNames(m_nparameters);
    for (size_t i=0; i<m_nparameters; i++) {
      pTypes[i] = ParameterType::_Base_;
      pNames[i] = "par_"+conv(i+1, par::fINT);
    }
    m_parameter_type = pTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==nparameters) && (parameterNames.size()==0)) {
    m_nparameters = nparameters;
    vector<string> pNames(m_nparameters);
    for (size_t i=0; i<m_nparameters; i++)
      pNames[i] = "par_"+conv(i+1, par::fINT);

    m_parameter_type = parameterTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)) {
    m_nparameters = nparameters;
    vector<ParameterType> pTypes(m_nparameters);
    for (size_t i=0; i<m_nparameters; i++)
      pTypes[i] = ParameterType::_Base_;

    m_parameter_type = pTypes;
    m_parameter_name = parameterNames;
  }

  ModelParameters::m_set_parameter_type();
  m_parameter_bestfit_value.erase(m_parameter_bestfit_value.begin(), m_parameter_bestfit_value.end());

  set_prior_distribution(priorDistribution);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_parameter_covariance (const int start, const int thin)
{
  vector<vector<double>> values;

  for (size_t i=0; i<m_nparameters; i++)
    values.push_back(parameter_chain_values(i, start, thin));

  covariance_matrix (transpose(values), m_parameter_covariance);
}


// ============================================================================================


double cbl::statistics::PosteriorParameters::parameter_covariance (const int i, const int j) const
{
  return m_parameter_covariance[i][j];
}


// ============================================================================================


std::vector<std::vector<double>> cbl::statistics::PosteriorParameters::parameter_covariance () const
{
  return m_parameter_covariance;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_prior_distribution (const int p, const std::shared_ptr<PriorDistribution> priorDistribution)
{
  switch (m_parameter_type[p]) {

  case statistics::ParameterType::_Base_:
    m_parameter_prior[p] = priorDistribution;
    break;

  case statistics::ParameterType::_Correlated_:
    m_parameter_prior[p] = priorDistribution;
    break;

  case statistics::ParameterType::_Derived_:
    m_parameter_prior[p] = NULL;
    WarningMsgCBL(m_parameter_name[p]+" is a derived parameter!", "set_prior_distribution", "PosteriorParameters.cpp");
    break;

  default:
    ErrorCBL("no such kind of parameter!", "set_prior_distribution", "PosteriorParameters.cpp");
  }
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_prior_distribution (const std::vector<std::shared_ptr<PriorDistribution>> priorDistribution)
{
  size_t priorSize = priorDistribution.size();
  size_t diff = priorSize-m_nparameters_base+m_nparameters_correlated;

  m_parameter_prior.erase(m_parameter_prior.begin(), m_parameter_prior.end());
  m_parameter_prior.resize(m_nparameters, NULL);
  m_correlated_parameter_prior.erase(m_correlated_parameter_prior.begin(), m_correlated_parameter_prior.end());
  m_correlated_parameter_prior.resize(diff, NULL);
  m_multipriors.erase(m_multipriors.begin(), m_multipriors.end());
  m_multipriors.resize(diff, 0);

  if (m_nparameters_correlated==0 and diff>0)
    ErrorCBL ("The size of the prior vector is incorrect!", "set_prior_distribution", "PosteriorParameters.cpp");

  else if (m_nparameters_correlated>0) {
    WarningMsgCBL("I'm assuming that the correlated parameters have been included as last, with "+conv(diff, par::fINT)+" multidimensional prior(s)!", "set_prior_distribution", "PosteriorParameters.cpp");

    for (size_t i=diff; i --> 0;) {
      m_correlated_parameter_prior[i] = priorDistribution[i+priorSize-diff];
      m_multipriors[i] = priorDistribution[i+priorSize-diff]->get_size_distribution();
    }
  }

  size_t index_corr = 0;
  size_t nset_prior = 0;

  for (size_t p=0; p<m_nparameters_base; p++) {

    if (m_parameter_type[m_base_parameter[p]]==statistics::ParameterType::_Correlated_) {

      if (index_corr>(m_multipriors[nset_prior]-1)) {
	nset_prior++;
	index_corr = 0;
      }

      std::shared_ptr<cbl::glob::Distribution> Distr = priorDistribution[nset_prior+m_nparameters_base-m_nparameters_correlated]->get_distribution(index_corr);
      PriorDistribution priorDistr;
      if (Distr->distributionType()==cbl::glob::DistributionType::_Gaussian_)
	priorDistr = PriorDistribution(cbl::glob::DistributionType::_Gaussian_, {Distr->get_mean(), Distr->get_sigma()}, Distr->xmin(), Distr->xmax(), Distr->get_seed());
      else priorDistr = PriorDistribution(cbl::glob::DistributionType::_Uniform_, Distr->xmin(), Distr->xmax(), Distr->get_seed());

      set_prior_distribution(m_base_parameter[p], std::make_shared<PriorDistribution>(priorDistr));
      index_corr++;
    }

    else set_prior_distribution(m_base_parameter[p], priorDistribution[p]);
  }

  m_set_parameter_type();
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_prior_distribution_seed (const std::shared_ptr<random::UniformRandomNumbers_Int> ran_generator)
{
  for (size_t i=0; i<m_nparameters_base; i++)
    m_parameter_prior[m_base_parameter[i]]->set_seed(ran_generator->operator()());
}


// ============================================================================================


std::shared_ptr<cbl::statistics::Prior> cbl::statistics::PosteriorParameters::prior () const
{
  auto prior_function = [&] (const vector<double> parameters, const shared_ptr<void> prior_inputs)
    {
      (void)prior_inputs;
      double prior_value = 1;

      vector<double> new_pars = parameters;
      for (size_t i=0; i<m_nparameters_derived; i++)
	new_pars.erase(new_pars.begin()+m_derived_parameter[i]);

      for (size_t i=0; i<m_nparameters_base-m_nparameters_correlated; i++)
	prior_value *= m_parameter_prior[m_base_parameter[i]]->operator()(new_pars[m_base_parameter[i]]);

      int add_previous = 0;
      for (size_t i=0; i<m_multipriors.size(); i++) {
	int start = m_nparameters_base - m_nparameters_correlated + add_previous;
	add_previous += m_multipriors[i];
	vector<double> correlated_parameters(new_pars.begin()+start, new_pars.begin()+start+m_multipriors[i]);
	prior_value *= m_correlated_parameter_prior[i]->operator[](correlated_parameters);
      }

      return prior_value;
    };

  shared_ptr<void> prior_function_inputs = NULL;

  return make_shared<Prior>(prior_function, prior_function_inputs);

}


// ============================================================================================


double cbl::statistics::PosteriorParameters::bestfit_value (const int p) const
{
  if (m_parameter_bestfit_value.size()==0)
    ErrorCBL("the best-fit values have not been computed!", "bestfit_value", "PosteriorParameters.cpp");

  return m_parameter_bestfit_value[p];
}

// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::bestfit_value () const
{
  if (m_parameter_bestfit_value.size()==0)
    ErrorCBL("the best-fit values have not been computed!", "bestfit_value", "PosteriorParameters.cpp");

  return m_parameter_bestfit_value;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_bestfit_values (const std::vector<double> bestfit_value)
{
  if (bestfit_value.size()!=m_nparameters)
    ErrorCBL("the size of the input vector is incorrect!", "set_bestfit_values", "PosteriorParameters.cpp");

  m_parameter_bestfit_value.erase(m_parameter_bestfit_value.begin(), m_parameter_bestfit_value.end());

  for (size_t i=0; i<m_nparameters; i++)
    m_parameter_bestfit_value.push_back(bestfit_value[i]);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_bestfit_values (const int start, const int thin, const int nbins, const int seed)
{
  set_posterior_distribution(start, thin, nbins, seed);
  set_parameter_covariance(start, thin);

  vector<double> bestfit_parameter;

  for (size_t i=0; i<m_nparameters; i++) {

    auto posterior = m_posterior_distribution[i];

    switch (m_parameter_type[i]) {

    case statistics::ParameterType::_Base_:
      if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_)
	bestfit_parameter.emplace_back(m_parameter_prior[i]->sample());
      else
	bestfit_parameter.emplace_back(posterior->median());
      break;

    case statistics::ParameterType::_Correlated_:
      bestfit_parameter.emplace_back(posterior->median());
      break;

    case statistics::ParameterType::_Derived_:
      bestfit_parameter.emplace_back(posterior->median());
      break;

    default:
      ErrorCBL("no such kind of parameter!", "set_bestfit_values", "PosteriorParameters.cpp");
    }
  }

  set_bestfit_values(bestfit_parameter);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::write_bestfit_info ()
{
  if (m_parameter_bestfit_value.size()==m_nparameters) {
    for (size_t i=0; i<m_nparameters; i++) {

      switch (m_parameter_type[i]) {
      case statistics::ParameterType::_Base_:
	if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_)
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	else
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " <<  par::col_green << "FREE" << endl;
	break;

      case statistics::ParameterType::_Correlated_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " <<  par::col_green << "FREE" << endl;
	break;

      case statistics::ParameterType::_Derived_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: "  << par::col_bred << "OUTPUT" << endl;
	break;

      default:
	ErrorCBL("no such kind of parameter!", "write_bestfit_info", "PosteriorParameters.cpp");
      }

      Print(m_parameter_bestfit_value[i], 5, 10, "value = ", "\n", true, std::cout);
    }
  }

  else
    ErrorCBL("the best-fit values have not been computed!", "write_bestfit_info", "PosteriorParameters.cpp");
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_chain (const size_t size, const size_t nwalkers)
{
  m_chain_size = size;
  m_chain_nwalkers = nwalkers;
  reset_chain();
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::reset_chain ()
{
  m_chain_value.erase(m_chain_value.begin(), m_chain_value.end());
  m_chain_value.resize(m_nparameters, vector<double>(m_chain_size*m_chain_nwalkers, 0));
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::expand_chain (const int append)
{
  vector<vector<double>> values = m_chain_value;
  size_t old_size = m_chain_size;

  m_chain_size += append;

  reset();

  for (size_t i=0; i<old_size; i++)
    for (size_t j=0; j<m_chain_nwalkers; j++)
      for (size_t p=0; p<m_nparameters; p++)
	set_chain_value(p, i, j, values[p][i*m_chain_nwalkers+j]);
}


// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::chain_value_parameter (const int pos, const int ww) const
{
  vector<double> vv(m_nparameters, 0);
  for (size_t i=0; i<m_nparameters; i++)
    vv[i] = chain_value(i, pos, ww);
  return vv;
}

// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::parameter_chain_values (const int param, const int start, const int thin) const
{
  vector<double> values;

  for (size_t i=start; i<m_chain_size; i+=thin)
    for (size_t j=0; j<m_chain_nwalkers; j++)
      values.push_back(chain_value(param, i, j));

  return values;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_chain_values (const std::vector<std::vector<double>> values, const int nwalkers)
{
  int size = values[0].size()/nwalkers;

  if (values[0].size()%nwalkers!=0)
    ErrorCBL("the size of the input values or the number of walkers is incorrect!", "set_chain_values", "PosteriorParameters.cpp");

  set_chain(size, nwalkers);

  for (size_t p=0; p<m_nparameters; p++)
    for (size_t i=0; i<m_chain_size; i++)
      for (size_t j=0; j<m_chain_nwalkers; j++)
	set_chain_value(p, i, j, values[p][i*m_chain_nwalkers+j]);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_chain_values (const std::vector<std::vector<std::vector<double>>> values)
{
  vector<vector<double>> flatten_val;

  int nwalkers = values[0][0].size();

  for (size_t i=0; i<values.size(); i++) {
    vector<double> vv;
    for (size_t j=0; j<values[0].size(); j++) {
      checkDim(values[i][j], nwalkers, "values["+conv(i, par::fINT)+"]"+"["+conv(j, par::fINT)+"]");
      for (size_t k=0; k<values[i][j].size(); k++)
	vv.push_back(values[i][j][k]);
    }
    flatten_val.push_back(vv);
  }

  set_chain_values(flatten_val, nwalkers);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain (const std::vector<std::vector<double>> values)
{
  for (size_t j=0; j<values[0].size(); j++) {
    vector<double> val;
    for (size_t i=0; i<values.size(); i++)
      val.emplace_back(values[i][j]);

    val = full_parameter(val);
    for (size_t i=0; i<val.size(); i++)
      set_chain_value(i, 0, j, val[i]);
  }
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain_from_prior ()
{
  std::vector<std::vector<double>> values(m_nparameters, std::vector<double>(m_chain_nwalkers, 0));

  for (size_t i=0; i<m_nparameters_base; i++) {
    const int index = m_base_parameter[i];
    for (size_t j=0; j<m_chain_nwalkers; j++) {
      values[index][j] = m_parameter_prior[index]->sample();
    }
  }
  initialize_chain(values);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain_ball (const std::vector<double> center, const double radius, const double seed)
{
  vector<vector<double>> value(m_nparameters, vector<double>(m_chain_nwalkers, 0));

  // verify if the number of given parameters is fine and, in case
  // only the free parameters are given in inputs, fill the empy
  // elements
  vector<double> cen = full_parameter(center);

  cbl::random::UniformRandomNumbers ran(-radius, radius, seed);

  for (size_t i=0; i<m_nparameters_free; i++) {
    int index = m_free_parameter[i];

    if (cen[index]+radius<m_parameter_prior[index]->xmin() or cen[index]-radius>m_parameter_prior[index]->xmax())
      ErrorCBL("Wrong parameter range!", "initialize_chain_ball", "PosteriorParameters.cpp");

    for (size_t j=0; j<m_chain_nwalkers; j++) {
      double val = ran()+cen[index];
      while (!m_parameter_prior[index]->isIncluded(val))
	val = ran()+cen[index];
      value[index][j] = val;
    }
  }
  initialize_chain(value);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain_ball_bestfit (const double radius, const double seed)
{
  initialize_chain_ball(bestfit_value(), radius, seed);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_posterior_distribution (const int start, const int thin, const int nbins, const int seed, const vector<double> weight)
{
  m_posterior_distribution.erase(m_posterior_distribution.begin(), m_posterior_distribution.end());
  m_posterior_distribution.resize(m_nparameters);

  for (size_t i=0; i<m_nparameters_base; i++) {
    const int index = m_base_parameter[i];
    if (m_parameter_prior[index]->distributionType()==glob::DistributionType::_Constant_)
      m_posterior_distribution[index] = make_shared<PosteriorDistribution>(PosteriorDistribution(glob::DistributionType::_Constant_, m_parameter_prior[index]->sample()));
    else
      m_posterior_distribution[index] = make_shared<PosteriorDistribution>(PosteriorDistribution(glob::DistributionType::_Discrete_, parameter_chain_values(index, start, thin), weight, nbins, "Spline", seed));
  }

  for (size_t i=0; i<m_nparameters_derived; i++) {
    const int index = m_derived_parameter[i];
    m_posterior_distribution[index] = make_shared<PosteriorDistribution>(PosteriorDistribution(glob::DistributionType::_Discrete_, parameter_chain_values(index, start, thin), weight, nbins, "Spline", seed));
  }
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::show_results (const int start, const int thin, const int nbins, const int seed, const bool show_mode, const int ns, const int nb, const vector<double> weight)
{
  set_posterior_distribution(start, thin, nbins, seed, weight);
  set_parameter_covariance(start, thin);

  const int dp = cout.precision(); cout.precision(4);
  cout << endl;

  // correction for covariance matrix uncertainties (Percival et al. 2014)
  double corr = 1.;
  if (ns>0 && nb>0) {
    const double AA = 2./double((ns-nb-1.)*(ns-nb-4.));
    const double BB = (ns-nb-2.)/double((ns-nb-1.)*(ns-nb-4.));
    corr = sqrt((1.+BB*(nb-m_nparameters))/(1.+AA+BB*(m_nparameters-1.)));
  }

  for (size_t i=0; i<m_nparameters; i++) {

    switch (m_parameter_type[i]) {

      case statistics::ParameterType::_Base_:
      if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_) {
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	Print(m_parameter_prior[i]->sample(), 5, 10, "value =", "\n", true, std::cout);
      }
      else {
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_green << "FREE" << endl;

	auto posterior = m_posterior_distribution[i];
	const double std = posterior->std()*corr;
	const double std_diff = (posterior->std()*(corr-1.))*0.5;

	Print(posterior->mean(), 5, 10, "Posterior mean = ", "\n", true, std::cout);
	Print(std, 5, 10, "Posterior standard deviation = ", "\n", true, std::cout);
	Print(posterior->median(), 5, 10, "Posterior median = ", "\n", true, std::cout);
	Print((posterior->median()-posterior->percentile(16))-std_diff, 5, 10, "Posterior 16th percentile = ", "\n", true, std::cout);
	Print((posterior->percentile(84)-posterior->median())+std_diff, 5, 10, "Posterior 84th percentile = ", "\n", true, std::cout);
	if (show_mode) Print(posterior->mode(), 5, 10, "Posterior mode = ", "\n", true, std::cout);
	cout << endl;
      }
      break;

    case statistics::ParameterType::_Correlated_:
      {
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_green << "FREE" << endl;

	auto posterior = m_posterior_distribution[i];
	const double std = posterior->std()*corr;
	const double std_diff = (posterior->std()*(corr-1.))*0.5;

	Print(posterior->mean(), 5, 10, "Posterior mean = ", "\n", true, std::cout);
	Print(std, 5, 10, "Posterior standard deviation = ", "\n", true, std::cout);
	Print(posterior->median(), 5, 10, "Posterior median = ", "\n", true, std::cout);
	Print((posterior->median()-posterior->percentile(16))-std_diff, 5, 10, "Posterior 16th percentile = ", "\n", true, std::cout);
	Print((posterior->percentile(84)-posterior->median())+std_diff, 5, 10, "Posterior 84th percentile = ", "\n", true, std::cout);
	if (show_mode) Print(posterior->mode(), 5, 10, "Posterior mode = ", "\n", true, std::cout);
	cout << endl;
	break;
      }

    case statistics::ParameterType::_Derived_:
      {
	auto posterior = m_posterior_distribution[i];
	const double std = posterior->std()*corr;
	const double std_diff = (posterior->std()*(corr-1.))*0.5;
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_bred << "OUTPUT" << endl;
	Print(posterior->mean(), 5, 10, "Posterior mean = ", "\n", true, std::cout);
	Print(std, 5, 10, "Posterior standard deviation = ", "\n", true, std::cout);
	Print(posterior->median(), 5, 10, "Posterior median = ", "\n", true, std::cout);
	Print((posterior->median()-posterior->percentile(16))-std_diff, 5, 10, "Posterior 16th percentile = ", "\n", true, std::cout);
	Print((posterior->percentile(84)-posterior->median())+std_diff, 5, 10, "Posterior 84th percentile = ", "\n", true, std::cout);
	if (show_mode) Print(posterior->mode(), 5, 10, "Posterior mode = ", "\n", true, std::cout);
	cout << endl;
	break;
      }

      default:
      ErrorCBL("no such kind of parameter!", "show_results", "PosteriorParameters.cpp");
    }
  }

  cout.precision(dp);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::write_results (const string dir, const string file, const int start, const int thin, const int nbins, const int seed, const bool compute_mode, const int ns, const int nb, const vector<double> weight)
{
  set_posterior_distribution(start, thin, nbins, seed, weight);
  set_parameter_covariance(start, thin);

  const string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}

  const string file_parameters = dir+file+"_parameters.dat";
  const string file_covariance = dir+file+"_covariance.dat";


  // correction for covariance matrix uncertainties (Percival et al. 2014)

  double corr = 1.;
  if (ns>0 && nb>0) {
    const double AA = 2./double((ns-nb-1.)*(ns-nb-4.));
    const double BB = (ns-nb-2.)/double((ns-nb-1.)*(ns-nb-4.));
    corr = sqrt((1.+BB*(nb-m_nparameters))/(1.+AA+BB*(m_nparameters-1.)));
  }

  ofstream fout(file_parameters.c_str());

  if (compute_mode)
    fout << "### Parameter # status # Posterior mean # Posterior standard deviation # Posterior median # Posterior 16th percentile # Posterior 84th percentile # Posterior mode ###" << endl << endl;
  else
    fout << "### Parameter # status # Posterior mean # Posterior standard deviation # Posterior median # Posterior 16th percentile # Posterior 84th percentile ###" << endl;

  for (size_t i=0; i<m_nparameters; i++) {
    auto posterior = m_posterior_distribution[i];

    if (m_parameter_type[i]==statistics::ParameterType::_Base_ && m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_) {
      fout << m_parameter_name[i] << " FIXED " << m_parameter_prior[i]->sample() << " 0 0 0 0 0";
      if (compute_mode)
	fout << "0";
    }
    else {

      const double std = posterior->std()*corr;
      const double std_diff = (posterior->std()*(corr-1.))*0.5;
      fout << m_parameter_name[i] << " FREE " << posterior->mean() << " " << std << " " << posterior->median() << " " << (posterior->median()-posterior->percentile(16))-std_diff << " " << (posterior->percentile(84)-posterior->median())+std_diff;
      if (compute_mode)
	fout << " " << posterior->mode() << endl;
    }
    fout << endl;
  }

  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file_parameters << endl;


  fout.open(file_covariance.c_str());
  for (size_t i=0; i<m_nparameters; i++) {
    for (size_t j=0; j<m_nparameters; j++)
      fout << i << " " << j << " " << m_parameter_covariance[i][j] << endl;
    fout << endl;
  }

  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file_covariance << endl;

}
