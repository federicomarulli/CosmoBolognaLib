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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "PosteriorParameters.h"

using namespace std;

using namespace cbl;
using namespace statistics;


// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::full_parameters (const std::vector<double> parameter_values) const
{
  if (parameter_values.size() == m_nparameters_free){
    vector<double> all_parameters(m_nparameters, 0);

    for (size_t i=0; i<m_nparameters_free; i++)
      all_parameters[m_free_parameters[i]] = parameter_values[i];

    for (size_t i=0; i<m_nparameters_fixed; i++)
      all_parameters[m_fixed_parameters[i]] = m_parameter_prior[m_fixed_parameters[i]]->sample();

    for (size_t i=0; i<m_nparameters_derived; i++)
      all_parameters[m_derived_parameters[i]] = 0.;

    return all_parameters;
  }
  else if (parameter_values.size() == m_nparameters) {
    vector<double> all_parameters = parameter_values;

    for (size_t i=0; i<m_nparameters_fixed; i++)
      all_parameters[m_fixed_parameters[i]] = m_parameter_prior[m_fixed_parameters[i]]->sample();
    
    return all_parameters;
  }
  else
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::full_parameters of PosteriorParameters.cpp, the vector of values of free parameters has the wrong size!");

  vector<double> vv;
  return vv;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::m_set_parameter_type ()
{
  m_nparameters_free = 0;
  m_nparameters_fixed = 0;
  m_nparameters_base = 0;
  m_nparameters_derived = 0;

  m_base_parameters.erase(m_base_parameters.begin(), m_base_parameters.end());
  m_fixed_parameters.erase(m_fixed_parameters.begin(), m_fixed_parameters.end());
  m_free_parameters.erase(m_free_parameters.begin(), m_free_parameters.end());
  m_derived_parameters.erase(m_derived_parameters.begin(), m_derived_parameters.end());

  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
      case statistics::ParameterType::_Base_:
	if (m_parameter_prior[i] != NULL){
	  if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_) {
	    m_nparameters_fixed +=1;
	    m_fixed_parameters.push_back(i);
	  }
	  else {
	    m_nparameters_free +=1;
	    m_free_parameters.push_back(i);
	  }
	  m_nparameters_base +=1;
	  m_base_parameters.push_back(i);
	}
	break;

      case statistics::ParameterType::_Derived_:
	m_nparameters_derived += 1;
	m_derived_parameters.push_back(i);
	break;

      default:
	ErrorCBL("Error in cbl::statistics::PosteriorParameters::m_set_parameter_type of PosteriorParameters.cpp: no such kind of parameter!");
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
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::set_parameters() of ModelParameters.cpp: nparameters must be >0!");

  if ((parameterTypes.size()!=nparameters) && (parameterTypes.size()!=0))
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::set_parameters() of ModelParameters.cpp: wrong size of the vector parameterTypes.!");

  if ((parameterNames.size()!=nparameters) && (parameterNames.size()!=0))
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::set_parameters() of ModelParameters.cpp: wrong size of the vector parameterNames.!");


  if ((parameterTypes.size()==nparameters) && (parameterNames.size()==nparameters)){
    m_nparameters=nparameters;
    m_parameter_type = parameterTypes;
    m_parameter_name = parameterNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)){
    m_nparameters=nparameters;
    vector<ParameterType> pTypes(m_nparameters);
    vector<string> pNames(m_nparameters);
    for (size_t i=0; i<m_nparameters; i++){
      pTypes[i] = ParameterType::_Base_;
      pNames[i] = "par_"+conv(i+1, par::fINT);
    }
    m_parameter_type = pTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==nparameters) && (parameterNames.size()==0)){
    m_nparameters=nparameters;
    vector<string> pNames(m_nparameters);
    for (size_t i=0; i<m_nparameters; i++)
      pNames[i] = "par_"+conv(i+1, par::fINT);
    
    m_parameter_type = parameterTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)){
    m_nparameters=nparameters;
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


size_t cbl::statistics::PosteriorParameters::nparameters_free () const
{
  return m_nparameters_free;
}


// ============================================================================================


size_t cbl::statistics::PosteriorParameters::nparameters_fixed () const
{
  return m_nparameters_fixed;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_parameter_covariance (const int start, const int thin)
{
  vector<vector<double>> values(m_nparameters, vector<double>(m_nparameters, 0));

  for (size_t i=0; i<m_nparameters; i++)
    values.push_back(parameter_chain_values(i, start, thin));

  covariance_matrix (values, m_parameter_covariance);
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

    case statistics::ParameterType::_Derived_:
      m_parameter_prior[p] = NULL;
      WarningMsg("Warning in set_prior_distribution of PosteriorParameters, "+m_parameter_name[p]+" is a derived parameter");
      break;

    default:
      ErrorCBL("Error in cbl::statistics::set_prior_distribution() of PosteriorParameters.cpp: no such kind of parameter!");
  }
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_prior_distribution (const std::vector<std::shared_ptr<PriorDistribution>> priorDistribution)
{
  if (m_nparameters_base != priorDistribution.size())
    ErrorCBL ("Error in cbl::statistics::PosteriorParameters::set_prior_distribution() of PosteriorParameters.cpp: wrong size of the prior vector!");

  m_parameter_prior.erase(m_parameter_prior.begin(), m_parameter_prior.end());
  m_parameter_prior.resize(m_nparameters, NULL);

  for (size_t p=0; p< m_nparameters_base; p++) 
    set_prior_distribution(m_base_parameters[p], priorDistribution[p]);

  m_set_parameter_type ();
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_prior_distribution_seed (const std::shared_ptr<random::UniformRandomNumbers_Int> ran_generator)
{
  for (size_t i=0; i< m_nparameters_base; i++) 
    m_parameter_prior[m_base_parameters[i]]->set_seed(ran_generator->operator()());
}


// ============================================================================================


std::shared_ptr<cbl::statistics::Prior> cbl::statistics::PosteriorParameters::prior () const
{
  auto prior_function = [&] (const vector<double> parameters, const shared_ptr<void> prior_inputs)
  {
    (void)prior_inputs;
    double prior_value = 1;
    for (size_t i=0; i<m_nparameters_base; i++) 
      prior_value *= m_parameter_prior[m_base_parameters[i]]->operator()(parameters[m_base_parameters[i]]);
    
    return prior_value;
  };

  shared_ptr<void> prior_function_inputs = NULL;

  return make_shared<Prior>(prior_function, prior_function_inputs);
}


// ============================================================================================


double cbl::statistics::PosteriorParameters::bestfit_value (const int p) const
{
  if (m_parameter_bestfit_value.size() == 0) 
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::bestfit_value() of PosteriorParameters: can't found best fit values!"); 

  return m_parameter_bestfit_value[p];
}

// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::bestfit_values () const
{
  if (m_parameter_bestfit_value.size() == 0) 
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::bestfit_values() of PosteriorParameters: the best-fit values couldn't be found!"); 

  return m_parameter_bestfit_value;
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_bestfit_value (const std::vector<double> bestfit_value)
{
  if (bestfit_value.size() != m_nparameters)
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::set_bestfit_value() of PosteriorParameters: wrong size for the input vector!"); 

  m_parameter_bestfit_value.erase(m_parameter_bestfit_value.begin(), m_parameter_bestfit_value.end());

  for (size_t i=0; i<m_nparameters; i++) 
    m_parameter_bestfit_value.push_back(bestfit_value[i]);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::write_bestfit_info ()
{
  if (m_parameter_bestfit_value.size() == m_nparameters) {
    for (size_t i=0; i<m_nparameters; i++) {

      switch (m_parameter_type[i]) {
	case statistics::ParameterType::_Base_:
	if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_) 
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	else 
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " <<  par::col_green << "FREE" << endl;
	break;

	case statistics::ParameterType::_Derived_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: "  << par::col_bred << "OUTPUT" << endl;
	break;

	default:
	ErrorCBL("Error in cbl::statistics::PosteriorParameters::write_bestfit_info() of PosteriorParameters.cpp: no such kind of parameter!");
      }

      coutCBL << "value = " << m_parameter_bestfit_value[i] << endl << endl;
    }
  }
  
  else
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::write_bestfit_info() of PosteriorParameters: can't found best fit values!"); 
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_chain (const size_t size, const size_t nwalkers)
{
  m_chain_size = size;
  m_chain_nwalkers = nwalkers;
  reset_chain();
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::reset_chain()
{
  m_chain_values.erase(m_chain_values.begin(), m_chain_values.end());
  m_chain_values.resize(m_nparameters, vector<double>(m_chain_size*m_chain_nwalkers, 0));
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::expand_chain (const int append)
{
  vector<vector<double>> values = m_chain_values;
  size_t old_size = m_chain_size;

  m_chain_size += append;

  reset();

  for (size_t i=0; i<old_size; i++)
    for (size_t j=0; j<m_chain_nwalkers; j++)
      for (size_t p=0; p<m_nparameters; p++)
         set_chain_value(p, i, j, values[p][i*m_chain_nwalkers+j]);
}


// ============================================================================================


std::vector<double> cbl::statistics::PosteriorParameters::chain_value_parameters (const int pos, const int ww) const
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
    ErrorCBL("Error in cbl::statistics::PosteriorParameters::set_chain_values() of PosteriorParameters.cpp: wrong size of input values of wrong number of walkers!");

  set_chain(size, nwalkers);

  for (size_t p=0; p< m_nparameters; p++)
    for (size_t i=0; i< m_chain_size; i++)
      for (size_t j=0; j< m_chain_nwalkers; j++)
	set_chain_value(p, i, j, values[p][i*m_chain_nwalkers+j]);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_chain_values(const std::vector<std::vector<std::vector<double>>> values)
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

    val = full_parameters(val);
    for (size_t i=0; i<val.size(); i++)
      set_chain_value(i, 0, j, val[i]);
  }
} 


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain_from_prior ()
{
  std::vector<std::vector<double>> values(m_nparameters, std::vector<double>(m_chain_nwalkers, 0));

  for (size_t i=0; i<m_nparameters_base; i++) {
    const int index = m_base_parameters[i];
    for (size_t j=0; j<m_chain_nwalkers; j++) {
      values[index][j] = m_parameter_prior[index]->sample();
    }
  }
  initialize_chain(values);
} 


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain_ball (const std::vector<double> center, const double radius, const double seed)
{
  vector<vector<double>> values(m_nparameters, vector<double>(m_chain_nwalkers, 0));
  vector<double> cen = full_parameters(center);

  cbl::random::UniformRandomNumbers ran(-radius, radius, seed);
  for (size_t i=0; i<m_nparameters_free; i++) {
    const int index = m_free_parameters[i];
    for (size_t j=0; j<m_chain_nwalkers; j++)      
      values[index][j] = ran()+cen[index];
  }

  initialize_chain(values);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::initialize_chain_ball_bestfit (const double radius, const double seed)
{
  initialize_chain_ball(bestfit_values(), radius, seed);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::set_posterior_distribution (const int start, const int thin, const int nbins, const int seed)
{
  m_parameter_posterior.erase(m_parameter_posterior.begin(), m_parameter_posterior.end());
  m_parameter_posterior.resize(m_nparameters);

  for (size_t i=0; i<m_nparameters_base; i++) {
    const int index = m_base_parameters[i];
    if (m_parameter_prior[index]->distributionType()==glob::DistributionType::_Constant_)  
      m_parameter_posterior[index] = make_shared<PosteriorDistribution>(PosteriorDistribution(glob::DistributionType::_Constant_, m_parameter_prior[index]->sample()));
    else 
      m_parameter_posterior[index] = make_shared<PosteriorDistribution>(PosteriorDistribution(glob::DistributionType::_Discrete_, parameter_chain_values(index, start, thin), {}, nbins, "Spline", seed));
  }

  for (size_t i=0; i<m_nparameters_derived; i++) {
    const int index = m_derived_parameters[i];
    m_parameter_posterior[index] = make_shared<PosteriorDistribution>(PosteriorDistribution(glob::DistributionType::_Discrete_, parameter_chain_values(index, start, thin), {}, nbins, "Spline", seed));
  }
}

// ============================================================================================


void cbl::statistics::PosteriorParameters::show_results (const int start, const int thin, const int nbins, const int seed, const bool show_mode)
{
  set_posterior_distribution(start, thin, nbins, seed);
  set_parameter_covariance(start, thin);

  const int dp = cout.precision(); cout.precision(4);
  cout << endl;

  for (size_t i=0; i<m_nparameters; i++) {

    auto posterior = m_parameter_posterior[i];

    switch (m_parameter_type[i]) {

      case statistics::ParameterType::_Base_:
	if (m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_) {
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	  coutCBL << "value = " << m_parameter_prior[i]->sample() << endl;
	  cout << endl;
	}
	else {
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_green << "FREE" << endl;
	  coutCBL << "Posterior mean = " << posterior->mean() << endl;
	  coutCBL << "Posterior standard deviation = " << posterior->std() << endl;
	  coutCBL << "Posterior median = " << posterior->median() << endl;
	  coutCBL << "Posterior 18th percentile = " << posterior->median()-posterior->percentile(18) << endl;
	  coutCBL << "Posterior 82th percentile = " << posterior->percentile(82)-posterior->median() << endl;
	  if (show_mode) coutCBL << "Posterior mode = " << posterior->mode() << endl;
	  cout << endl;
	}
	break;

      case statistics::ParameterType::_Derived_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_bred << "OUTPUT" << endl;
	coutCBL << "Posterior mean = " << posterior->mean() << endl;
	coutCBL << "Posterior standard deviation = " << posterior->std() << endl;
	coutCBL << "Posterior median = " << posterior->median() << endl;
	coutCBL << "Posterior 18th percentile = " << posterior->median()-posterior->percentile(18) << endl;
	coutCBL << "Posterior 82th percentile = " << posterior->percentile(82)-posterior->median() << endl;
	if (show_mode) coutCBL << "Posterior mode = " << posterior->mode() << endl;
	coutCBL << endl;
	break;

      default:
	ErrorCBL("Error in cbl::statistics::PosteriorParameters::show_results() of PosteriorParameters.cpp: no such kind of parameter!");
    }
  }
  
  cout.precision(dp);
}


// ============================================================================================


void cbl::statistics::PosteriorParameters::write_results (const string dir, const string file, const int start, const int thin, const int nbins, const int seed, const bool compute_mode)
{
  set_posterior_distribution (start, thin, nbins, seed);
  set_parameter_covariance(start, thin);
  
  string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}

  string file_parameters = dir+file+"_parameters.dat";
  string file_covariance = dir+file+"_covariance.dat";

  ofstream fout(file_parameters.c_str());

  if (compute_mode)
    fout << "### Parameter # status ###" << endl << "### Posterior mean # Posterior standard deviation # Posterior median # Posterior 18th percentile # Posterior 82th percentile # Posterior mode ###" << endl << endl;
  else
    fout << "### Parameter # status ###" << endl << "### Posterior mean # Posterior standard deviation # Posterior median # Posterior 18th percentile # Posterior 82th percentile ###" << endl << endl;
  
  for (size_t i=0; i<m_nparameters; i++) 
    if (m_parameter_type[i]==statistics::ParameterType::_Base_ && m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_)
      fout << "# " << m_parameter_name[i] << " FIXED #" << endl;
    else 
      fout << "# " << m_parameter_name[i] << " FREE #" << endl;

  for (size_t i=0; i<m_nparameters; i++) {
    auto posterior = m_parameter_posterior[i];
    if (m_parameter_type[i]==statistics::ParameterType::_Base_ && m_parameter_prior[i]->distributionType()==glob::DistributionType::_Constant_)
      fout << posterior->sample() << " 0 0 0 0 0 0" << endl;
    else 
      fout << posterior->mean() << " " << posterior->std() << " " << posterior->median() << " " << posterior->median()-posterior->percentile(18) << " " << posterior->percentile(82)-posterior->median();
    if (compute_mode)
      fout << " " << posterior->mode() << endl;
    else
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

