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
 *  @file
 *  Modelling/MassObservableRelation/Modelling_MassObservableRelation.cpp
 *
 *  @brief Methods of the class Modelling_MassObservableRelation
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_MassObservableRelation, i.e. the common functions to model
 *  the galaxy cluster mass - mass proxy relation
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "Modelling_MassObservableRelation.h"

using namespace std;

using namespace cbl;


// ===========================================================================================

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const std::vector<double> redshift, const double redshift_pivot, const double proxy_or_mass_pivot, const double log_base)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units
  m_data_model.cluster = make_shared<catalogue::Cluster>(cluster);
  
  m_data_model.redshift = redshift;
  m_data_model.redshift_pivot = redshift_pivot;
  
  m_data_model.proxy_or_mass_pivot = proxy_or_mass_pivot;
  m_data_model.log_base = log_base;
}


// ===========================================================================================

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_model_MassObservableRelation_cosmology (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior)
{
  const size_t nParams = 8;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors for the mass-observable relation parameters
  Par_string[0] = "alpha";
  param_prior[0] = alpha_prior;
  Par_string[1] = "beta";
  param_prior[1] = beta_prior;
  Par_string[2] = "gamma";
  param_prior[2] = gamma_prior;
  Par_string[3] = "scatter0";
  param_prior[3] = scatter0_prior;
  Par_string[4] = "scatterM";
  param_prior[4] = scatterM_prior;
  Par_string[5] = "scatterM_exponent";
  param_prior[5] = scatterM_exponent_prior;
  Par_string[6] = "scatterz";
  param_prior[6] = scatterz_prior;
  Par_string[7] = "scatterz_exponent";
  param_prior[7] = scatterz_exponent_prior;

  // input data used to construct the model
  auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_direct_z, nParams, Par_type, Par_string, inputs));
}



// ===========================================================================================

void cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_model_MassObservableRelation_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior)
{
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+8; // The total number of parameters is given by the cosmological ones + 8

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the mass-observable relation parameters
  Par_string[cosmo_param.size()] = "alpha";
  param_prior[cosmo_param.size()] = alpha_prior;
  Par_string[cosmo_param.size()+1] = "beta";
  param_prior[cosmo_param.size()+1] = beta_prior;
  Par_string[cosmo_param.size()+2] = "gamma";
  param_prior[cosmo_param.size()+2] = gamma_prior;
  Par_string[cosmo_param.size()+3] = "scatter0";
  param_prior[cosmo_param.size()+3] = scatter0_prior;
  Par_string[cosmo_param.size()+4] = "scatterM";
  param_prior[cosmo_param.size()+4] = scatterM_prior;
  Par_string[cosmo_param.size()+5] = "scatterM_exponent";
  param_prior[cosmo_param.size()+5] = scatterM_exponent_prior;
  Par_string[cosmo_param.size()+6] = "scatterz";
  param_prior[cosmo_param.size()+6] = scatterz_prior;
  Par_string[cosmo_param.size()+7] = "scatterz_exponent";
  param_prior[cosmo_param.size()+7] = scatterz_exponent_prior;

  // input data used to construct the model
  auto inputs = make_shared<STR_MOrelation_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&model_E_z, nParams, Par_type, Par_string, inputs));
}


// ===========================================================================================

double cbl::modelling::massobsrel::scaling_relation (cbl::catalogue::Cluster cluster, const double div_proxy_or_mass, const double f_z)
{  
  return cluster.alpha_scaling_rel() + cluster.beta_scaling_rel() * div_proxy_or_mass + cluster.gamma_scaling_rel() * f_z;
}


// ===========================================================================================

std::vector<double> cbl::modelling::massobsrel::model_E_z (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_alpha_scaling_rel(parameter[pp->Cpar.size()]);
  cluster.set_beta_scaling_rel(parameter[pp->Cpar.size()+1]);
  cluster.set_gamma_scaling_rel(parameter[pp->Cpar.size()+2]);
  cluster.set_scatter0_scaling_rel(parameter[pp->Cpar.size()+3]);
  cluster.set_scatterM_scaling_rel(parameter[pp->Cpar.size()+4]);
  cluster.set_scatterM_exponent_scaling_rel(parameter[pp->Cpar.size()+5]);
  cluster.set_scatterz_scaling_rel(parameter[pp->Cpar.size()+6]);
  cluster.set_scatterz_exponent_scaling_rel(parameter[pp->Cpar.size()+7]);

  std::vector<double> res(proxy_or_mass.size());
  for (size_t j=0; j<proxy_or_mass.size(); j++) {
    double log1 = log(proxy_or_mass[j]/pp->proxy_or_mass_pivot)/log(pp->log_base);
    double log2 = log(cosmo.EE(pp->redshift[j])/cosmo.EE(pp->redshift_pivot))/log(pp->log_base);
    res[j] = scaling_relation(cluster, log1, log2);
  }

  return res;
}


// ===========================================================================================

std::vector<double> cbl::modelling::massobsrel::model_direct_z (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_MOrelation_data_model> pp = static_pointer_cast<STR_MOrelation_data_model>(inputs);

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;
  
  // set the cluster parameters
  cluster.set_alpha_scaling_rel(parameter[0]);
  cluster.set_beta_scaling_rel(parameter[1]);
  cluster.set_gamma_scaling_rel(parameter[2]);
  cluster.set_scatter0_scaling_rel(parameter[3]);
  cluster.set_scatterM_scaling_rel(parameter[4]);
  cluster.set_scatterM_exponent_scaling_rel(parameter[5]);
  cluster.set_scatterz_scaling_rel(parameter[6]);
  cluster.set_scatterz_exponent_scaling_rel(parameter[7]);

  std::vector<double> res(proxy_or_mass.size());
  for (size_t j=0; j<proxy_or_mass.size(); j++) {
    double log1 = log(proxy_or_mass[j]/pp->proxy_or_mass_pivot)/log(pp->log_base);
    double log2 = log((1+pp->redshift[j])/(1+pp->redshift_pivot))/log(pp->log_base);
    res[j] = scaling_relation(cluster, log1, log2);
  }

  return res;
}
