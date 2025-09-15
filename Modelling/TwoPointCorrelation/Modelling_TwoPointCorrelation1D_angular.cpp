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
 *  @file
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D_angular.cpp
 *
 *  @brief Methods of the class
 *  Modelling_TwoPointCorrelation1D_angular
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation1D_angular, used to model the angular
 *  of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, massimilia.romanell2@unibo.it
 */


#include "Modelling_TwoPointCorrelation1D_angular.h"

using namespace cbl;
using namespace std;

  
// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_data_model (const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double z_min, const double z_max, const cbl::CoordinateUnits theta_units, const std::string method_Pk, const bool NL, const int norm, const double k_min, const double k_max, std::vector<double> dN_par, const std::string distribution_type, const bool redshift_space, const std::string interpType)
{
  m_data_model = make_shared<modelling::twopt::STR_data_model>(modelling::twopt::STR_data_model());

  m_data_model->z_abs_err = z_abs_err;
  m_data_model->cosmology = move(cosmology);
  m_data_model->z_min=z_min;
  m_data_model->z_max=z_max;
  m_data_model->theta_units=theta_units;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->norm = norm;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->dN_par=dN_par;
  m_data_model->distribution_type=distribution_type;
  m_data_model->redshift_space=redshift_space;
  m_data_model->interpType=interpType;

  if (norm==1) {
    cosmology::PkXi PX(m_data_model->cosmology);
    const double s8 = PX.sigma8_Pk(method_Pk, 0., false, "test", NL, k_min, k_max);
    m_data_model->cosmology->set_parameter("sigma8", s8);
  }
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_data_model_MF_theory (const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double z_min, const double z_max, const cbl::CoordinateUnits theta_units, const std::string method_Pk, const bool NL, const int norm, const double k_min, const double k_max, const double area_degrees, const bool redshift_space, const string interpType, const double Delta, const bool isDelta_critical, const std::string model_MF, const double Mass_min, const double Mass_max, const int mass_step)
{
  m_data_model = make_shared<STR_data_model>(STR_data_model());
  m_data_model->z_abs_err = z_abs_err;
  m_data_model->cosmology = move(cosmology);
  m_data_model->z_min = z_min;
  m_data_model->z_max = z_max;
  m_data_model->theta_units = theta_units;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->norm = norm;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->area_degrees = area_degrees;
  m_data_model->interpType = interpType;
  m_data_model->redshift_space = redshift_space;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter((z_max+z_min)*0.5) : Delta;
  m_data_model->model_MF = model_MF;
  m_data_model->Mass_min = Mass_min;
  m_data_model->Mass_max = Mass_max;
  m_data_model->mass_step = mass_step;

  if(norm!=0){
    cosmology::PkXi PX(m_data_model->cosmology);
    double s8 = PX.sigma8_Pk(method_Pk, 0., false, "test", NL, k_min, k_max);
    m_data_model->cosmology->set_parameter("sigma8", s8);
  }

}


//===================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_data_model_bias_MF_theory (const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> cluster_mass, const std::vector<double> cluster_redshift, const double z_min, const double z_max, const cbl::CoordinateUnits theta_units, const std::string method_Pk, const bool NL, const int norm, const double k_min, const double k_max, const std::string output_dir, const std::string output_root, const double area_degrees, const bool redshift_space, const string interpType, const double Delta, const bool isDelta_critical, const std::string model_MF, const std::string model_bias, const int mass_step, const bool bias_from_sample, const bool include_quadrupole)
{
  m_data_model = make_shared<STR_data_model>(STR_data_model());
  m_data_model->z_abs_err = z_abs_err;
  double redshift=cbl::Average(cluster_redshift); 
  m_data_model->cosmology = move(cosmology);
  m_data_model->z_min = z_min;
  m_data_model->z_max = z_max;
  m_data_model->theta_units = theta_units;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->norm = norm;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->output_dir = output_dir;
  m_data_model->output_root = output_root;
  m_data_model->area_degrees = area_degrees;
  m_data_model->interpType = interpType;
  m_data_model->redshift_space = redshift_space;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;
  m_data_model->model_MF = model_MF;
  m_data_model->model_bias = model_bias;
  m_data_model->mass_step = mass_step;
  m_data_model->bias_from_sample = bias_from_sample;
  m_data_model->include_quadrupole = include_quadrupole;

  if (norm!=0) {
    cosmology::PkXi PX(m_data_model->cosmology);
    const double s8 = PX.sigma8_Pk(method_Pk, 0., false, output_root, NL, k_min, k_max);
    m_data_model->cosmology->set_parameter("sigma8", s8);
  }

  if (cluster_mass.size()==0 || cluster_redshift.size()==0 || cluster_mass.size()!=cluster_redshift.size())
    ErrorCBL("The mass and redshift vectors must contain at least one object each, and must have the same size.", "set_data_model_bias_MF_theory", "Modelling_TwoPointCorrelation1D_angular.cpp");
  m_data_model->mass=cluster_mass;
  m_data_model->cluster_redshift=cluster_redshift;

  m_data_model->Mass_min = *min_element(cluster_mass.begin(), cluster_mass.end());
  m_data_model->Mass_max = *max_element(cluster_mass.begin(), cluster_mass.end());


  // Build a dummy dataset for the angular power spectrum object (useful only for wtheta covariance).
  std::vector<double> l = linear_bin_vector(200, 1., 200.);
  std::shared_ptr<cbl::data::Data> pow_dataset = std::make_shared<cbl::data::Data1D>(l,l,l); 
  cbl::modelling::angularpk::Modelling_PowerSpectrum_angular model_pow(pow_dataset);
  m_data_model->angularpk = make_shared<modelling::angularpk::Modelling_PowerSpectrum_angular>(model_pow);
  
  (m_data_model->angularpk)->set_data_model_bias_MF_theory(z_abs_err, cosmology, cluster_mass, cluster_redshift, z_min, z_max, method_Pk, NL, norm, k_min, k_max, output_dir, output_root, area_degrees, true, 50., {}, {}, redshift_space, interpType, Delta, isDelta_critical, model_MF, model_bias, mass_step, false, 1.e-5);
  
}


// ===================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_data_model_scaling_relation_cosmology (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk, const bool NL, const double k_min, const double k_max, const int norm, const bool redshift_space, std::vector<double> dN_par, const std::string distribution_type, const double z_min, const double z_max, const cbl::CoordinateUnits theta_units, const std::string output_dir, const std::string output_root, const double prec, const std::string file_par, const double Delta, const bool isDelta_critical, const std::string model_bias, const std::string interpType)
{
  if (file_par!=par::defaultString)
    WarningMsgCBL("check the consistency between the parameters of the object cosmology, provided in input, and the ones in the parameter file", "set_data_model", "Modelling_TwoPointCorrelation1D_angular.cpp");
  
  m_data_model = make_shared<modelling::twopt::STR_data_model>(modelling::twopt::STR_data_model());
  
  if (z_abs_err < 0 || proxy_rel_err < 0)
    ErrorCBL("The errors on z and proxy must be > 0.", "set_data_model", "Modelling_TwoPointCorrelation1D_angular.cpp");
  m_data_model->z_abs_err = z_abs_err;
  m_data_model->proxy_rel_err = proxy_rel_err;
  m_data_model->cosmology = std::move(cosmology);
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->redshift_space=redshift_space;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->z_min = z_min;
  m_data_model->z_max = z_max;
  m_data_model->theta_units=theta_units;
  m_data_model->output_dir = output_dir;
  m_data_model->norm = norm;
  m_data_model->dN_par=dN_par;
  m_data_model->distribution_type=distribution_type;
  m_data_model->prec = prec;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;
  m_data_model->model_bias = model_bias;
  m_data_model->interpType=interpType;

  if (cosmology->sigma8()>0) 
    m_data_model->sigma8_z = cosmology->sigma8(redshift);  
  else { 
    coutCBL << "sigma8 is not set, it will be computed from the power spectrum with " << method_Pk << endl;
    cosmology::PkXi PX(cosmology);
    m_data_model->sigma8_z = PX.sigma8_Pk(method_Pk, redshift, false, output_root, NL, k_min, k_max, prec, file_par);    
    coutCBL << "--> sigma8(z=" << redshift << ") = " << m_data_model->sigma8_z << endl << endl;
  }
  
  if (cluster_mass_proxy.size() == 0 || cluster_redshift.size() == 0 || cluster_mass_proxy.size() != cluster_redshift.size())
    ErrorCBL("The mass proxy and redshift vectors must contain at least one object each, and must have the same size.", "set_data_model", "Modelling_TwoPointCorrelation1D_angular.cpp");
  
  // Build a dummy dataset for the scaling relation Modelling object.
  // Only the x values of the dataset (the mass proxy values) are used in the computation
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(cluster_mass_proxy, cluster_mass_proxy, cluster_mass_proxy));
  
  // Build the scaling relation Modelling object
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation (dataset);
  m_data_model->scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);
  
  (m_data_model->scaling_relation)->set_data_model(cosmology, cluster_redshift, redshift_pivot, proxy_pivot, log_base);

  (m_data_model->scaling_relation)->set_mass_pivot(mass_pivot);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_data_model_scaling_relation_cosmology_MF_theory (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk, const bool NL, const double k_min, const double k_max, const int norm, const bool redshift_space, const double area_degrees, const double z_min, const double z_max, const double proxy_min, const double proxy_max, const cbl::CoordinateUnits theta_units, const std::string output_dir, const std::string output_root, const double prec, const std::string file_par, const double Delta, const bool isDelta_critical, const std::string model_bias, const string interpType, const std::string model_MF, const std::vector<double> SF_weights)
{
  if (file_par!=par::defaultString)
    WarningMsgCBL("check the consistency between the parameters of the object cosmology, provided in input, and the ones in the parameter file", "set_data_model", "Modelling_TwoPointCorrelation1D_angular.cpp");

  if (z_abs_err < 0 || proxy_rel_err < 0)
    ErrorCBL("The errors on z and proxy must be > 0.", "set_data_model", "Modelling_TwoPointCorrelation1D_angular.cpp");

  m_data_model = make_shared<modelling::twopt::STR_data_model>(modelling::twopt::STR_data_model());
  m_data_model->z_abs_err = z_abs_err;
  m_data_model->proxy_rel_err = proxy_rel_err;
  m_data_model->cosmology = move(cosmology);
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->norm = norm;
  m_data_model->theta_units=theta_units;
  m_data_model->z_min = z_min;
  m_data_model->z_max = z_max;
  m_data_model->proxy_min = proxy_min;
  m_data_model->proxy_max = proxy_max;
  m_data_model->output_dir = output_dir;
  m_data_model->prec = prec;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;
  m_data_model->model_bias = model_bias;
  m_data_model->redshift_space = redshift_space;
  m_data_model->interpType = interpType;
  m_data_model->area_degrees = area_degrees;
  m_data_model->model_MF = model_MF;
  m_data_model->SF_weights = SF_weights;

  if (cosmology->sigma8()>0) 
    m_data_model->sigma8_z = cosmology->sigma8(redshift);  
  else { 
    coutCBL << "sigma8 is not set, it will be computed from the power spectrum with " << method_Pk << endl;
    cosmology::PkXi PX(cosmology);
    m_data_model->sigma8_z = PX.sigma8_Pk(method_Pk, redshift, false, output_root, NL, k_min, k_max, prec, file_par);    
    coutCBL << "--> sigma8(z=" << redshift << ") = " << m_data_model->sigma8_z << endl << endl;
  }
  
  if (cluster_mass_proxy.size() == 0 || cluster_redshift.size() == 0 || cluster_mass_proxy.size() != cluster_redshift.size())
    ErrorCBL("The mass proxy and redshift vectors must contain at least one object each, and must have the same size.", "set_data_model", "Modelling_TwoPointCorrelation1D_angular.cpp");
  
  // Build a dummy dataset for the scaling relation Modelling object.
  // Only the x values of the dataset (the mass proxy values) are used in the computation
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(cluster_mass_proxy, cluster_mass_proxy, cluster_mass_proxy));
  
  // Build the scaling relation Modelling object
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation(dataset);
  m_data_model->scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);
  
  (m_data_model->scaling_relation)->set_data_model(cosmology, cluster_redshift, redshift_pivot, proxy_pivot, log_base);

  (m_data_model->scaling_relation)->set_mass_pivot(mass_pivot);
}


// =================================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior)
{
  m_data_model->Cpar = cosmo_param;

  const int nparameters = cosmo_param.size()+3; 
  
  std::vector<statistics::ParameterType> parameterType (nparameters-2, statistics::ParameterType::_Base_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  std::vector<std::string> parameterName(nparameters);
  std::vector<statistics::PriorDistribution> priors(nparameters-2);
  
  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    parameterName[i] = cosmo_param[i];
    priors[i] = cosmo_prior[i];
  }
  
  priors[nparameters-3] = bias_prior;
  parameterName[nparameters-3] = "bias";
  parameterName[nparameters-2] = "S8";
  parameterName[nparameters-1] = "Baryonic fraction";

  // set the priors
  m_set_prior(priors);

  // construct the model
  m_model = std::make_shared<statistics::Model1D>(statistics::Model1D(&wtheta, nparameters, parameterType, parameterName, m_data_model));
}


// =================================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution bias_prior)
{
  m_data_model->Cpar = cosmo_param;

  const int nparameters = cosmo_param.size()+3; 
  
  std::vector<statistics::ParameterType> parameterType (nparameters-2, statistics::ParameterType::_Base_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  std::vector<std::string> parameterName(nparameters);
  std::vector<statistics::PriorDistribution> priors(nparameters-2);
  
  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    parameterName[i] = cosmo_param[i];
    priors[i] = cosmo_prior[i];
  }
  
  priors[nparameters-3] = bias_prior;
  parameterName[nparameters-3] = "bias";
  parameterName[nparameters-2] = "S8";
  parameterName[nparameters-1] = "Baryonic fraction";

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = std::make_shared<statistics::Model1D>(statistics::Model1D(&wtheta_MF_theory, nparameters, parameterType, parameterName, m_data_model));

}


// =================================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_bias_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior)
{
  m_data_model->Cpar = cosmo_param;

  const int nparameters = cosmo_param.size()+2; 
  
  std::vector<statistics::ParameterType> parameterType (nparameters-2, statistics::ParameterType::_Base_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  std::vector<std::string> parameterName(nparameters);
  std::vector<statistics::PriorDistribution> priors(nparameters-2);
  
  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    parameterName[i] = cosmo_param[i];
    priors[i] = cosmo_prior[i];
  }
  
  parameterName[nparameters-2] = "S8";
  parameterName[nparameters-1] = "bias";

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = std::make_shared<statistics::Model1D>(statistics::Model1D(&wtheta_bias_MF_theory, nparameters, parameterType, parameterName, m_data_model));
}


// ========================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_bias_MF_theory_BAO (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution SigmaNL_prior)
{

  m_data_model->Cpar = cosmo_param;

  const int nparameters = cosmo_param.size()+2; 
  
  std::vector<statistics::ParameterType> parameterType (nparameters-1, statistics::ParameterType::_Base_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  std::vector<std::string> parameterName(nparameters);
  std::vector<statistics::PriorDistribution> priors(nparameters-1);
  
  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    parameterName[i] = cosmo_param[i];
    priors[i] = cosmo_prior[i];
  }
  parameterName[cosmo_param.size()] = "SigmaNL";
  parameterName[cosmo_param.size()+1] = "bias";

  priors[cosmo_param.size()] = SigmaNL_prior;

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = std::make_shared<statistics::Model1D>(statistics::Model1D(&wtheta_bias_MF_theory_BAO, nparameters, parameterType, parameterName, m_data_model));

}


// ========================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_bias_MF_theory_BAO_poly (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution SigmaNL_prior, const statistics::PriorDistribution A0_prior, const statistics::PriorDistribution A1_prior, const statistics::PriorDistribution A2_prior)
{

  m_data_model->Cpar = cosmo_param;

  const int nparameters = cosmo_param.size()+6; 
  
  std::vector<statistics::ParameterType> parameterType (nparameters-2, statistics::ParameterType::_Base_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Derived_);
  std::vector<std::string> parameterName(nparameters);
  std::vector<statistics::PriorDistribution> priors(nparameters-2);
  
  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    parameterName[i] = cosmo_param[i];
    priors[i] = cosmo_prior[i];
  }
  parameterName[cosmo_param.size()] = "SigmaNL";
  parameterName[cosmo_param.size()+1] = "A0";
  parameterName[cosmo_param.size()+2] = "A1";
  parameterName[cosmo_param.size()+3] = "A2";
  parameterName[cosmo_param.size()+4] = "S8";
  parameterName[cosmo_param.size()+5] = "bias";

  priors[cosmo_param.size()] = SigmaNL_prior;
  priors[cosmo_param.size()+1] = A0_prior;
  priors[cosmo_param.size()+2] = A1_prior;
  priors[cosmo_param.size()+3] = A2_prior;
  
  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = std::make_shared<statistics::Model1D>(statistics::Model1D(&wtheta_bias_MF_theory_BAO_poly, nparameters, parameterType, parameterName, m_data_model));

}


// ========================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_scaling_relation_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::string z_evo)
{

  // Set the scaling relation Modelling object
  if ((m_data_model->scaling_relation)->data()->xx().size() == 0)
    ErrorCBL("The mass-observable relation is not set! Use the correct set_data_model().", "set_model_scaling_relation_cosmology", "Modelling_TwoPointCorrelation1D_angular.cpp");
   
  (m_data_model->scaling_relation)->set_model_MassObservableRelation_cosmology(z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);
  (m_data_model->scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {}); // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
     
  // Set the parameter names and priors
  m_data_model->Cpar = cosmo_param;
  const size_t nParams_scaling_relation = (m_data_model->scaling_relation)->likelihood()->parameters()->nparameters();
  const size_t nParams_base = nParams_scaling_relation; // scaling relation parameters (which include the cosmological parameters)
  const size_t nParams_derived = 1; // the effective bias is derived from the scaling relation
   
  const size_t nParams = nParams_base + nParams_derived;
  
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[nParams-1] = statistics::ParameterType::_Derived_; //change the last parameter, derived
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams_base);
  // Cosmological and scaling relation parameters
  for (size_t i=0; i<nParams_base; i++) {
    Par_string[i] = (m_data_model->scaling_relation)->likelihood()->parameters()->name(i);
    param_prior[i] = *(m_data_model->scaling_relation)->get_prior(i);
  }
  // bias
  Par_string[nParams-1] = "bias";
  
  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&wtheta_scaling_relation_cosmology, nParams, Par_type, Par_string, m_data_model));
}


// =================================================================================================

  
void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_scaling_relation_cosmology_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<statistics::PriorDistribution> Plambda_prior, const std::string z_evo)
{
  if (Plambda_prior.size()!=3)
    ErrorCBL("Plambda_prior must have size = 3.","set_model_scaling_relation_cosmology_MF_theory","Modelling_TwoPointCorrelation1D_angular.cpp");
   
  // Set the scaling relation Modelling object
  if ((m_data_model->scaling_relation)->data()->xx().size() == 0)
    ErrorCBL("The mass-observable relation is not set! Use the correct set_data_model().", "set_model_scaling_relation_cosmology_MF_theory", "Modelling_TwoPointCorrelation1D_angular.cpp");
  (m_data_model->scaling_relation)->set_model_MassObservableRelation_cosmology(z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  (m_data_model->scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {}); // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
   
  // Set the parameter names and priors
  m_data_model->Cpar = cosmo_param;
  
  //const size_t nParams_scaling_relation = (m_data_model->scaling_relation)->likelihood()->parameters()->nparameters();
  const size_t nParams_base = cosmo_param.size()+11; // scaling relation parameters (which include the cosmological parameters)
  const size_t nParams_derived = 1; // the effective bias is derived from the scaling relation

  const size_t nParams = nParams_base + nParams_derived;
   
  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[nParams-1] = statistics::ParameterType::_Derived_; // change the last parameter, derived
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams_base);
  
  // Cosmological and scaling relation parameters
  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++) {
    Par_string[i] = cosmo_param[i];
    param_prior[i] = cosmo_prior[i];
  }

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
  Par_string[cosmo_param.size()+8] = "Plambda_a";
  param_prior[cosmo_param.size()+8] = Plambda_prior[0];
  Par_string[cosmo_param.size()+9] = "Plambda_b";
  param_prior[cosmo_param.size()+9] = Plambda_prior[1];
  Par_string[cosmo_param.size()+10] = "Plambda_c";
  param_prior[cosmo_param.size()+10] = Plambda_prior[2];
  
  // Set the functional form for the redshift evolution in the scaling relation
  
  if (z_evo == "E_z")
    m_data_model->fz = [] (const double z, const double z_piv, const std::shared_ptr<void> cosmo)
      {
	auto cosmology = std::static_pointer_cast<cbl::cosmology::Cosmology>(cosmo);
	return cosmology->Hubble(z)/cosmology->Hubble(z_piv);
      };

  else if (z_evo == "direct")
    m_data_model->fz = [] (const double z, const double z_piv, const std::shared_ptr<void> cosmo)
      {
	(void)cosmo; return (1+z)/(1+z_piv);
      };
  
  else
    cbl::ErrorCBL("Error in the input parameter z_evo: no such a possibility for f(z)!", "set_model_scaling_relation_cosmology_MF_theory", "Modelling_TwoPointCorrelation1D_angular.cpp");
  
  // bias
  Par_string[nParams-1] = "bias";
  
  // set prior
  m_set_prior(param_prior);
  
  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&wtheta_scaling_relation_cosmology_MF_theory, nParams, Par_type, Par_string, m_data_model));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_angular::set_model_covariance_bias_MF_theory (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution alpha_prior)
{
  // Set the angularpk Modelling object
  if ((m_data_model->angularpk)->data()->xx().size() == 0)
    ErrorCBL("The angularpk is not set! Use the correct set_data_model_bias_MF_theory().", "set_model_covariance_bias_MF_theory", "Modelling_TwoPointCorrelation1D_angular.cpp");

  //add fake prior for the modelling of angularpk
  cbl::statistics::PriorDistribution shot_noise_prior = cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Constant_, 0.);
  (m_data_model->angularpk)->set_model_bias_MF_theory(cosmo_param, cosmo_prior, shot_noise_prior);
  
  (m_data_model->angularpk)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {}); // Set the likelihood for the angularpk (only to avoid internal errors, of course it is not used)
  
  m_data_model->Cpar = cosmo_param;
  
  const int nparameters = cosmo_param.size()+2; 
  
  vector<statistics::ParameterType> parameterType (nparameters, statistics::ParameterType::_Base_);
  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    parameterName[i] = cosmo_param[i];
    priors[i] = cosmo_prior[i];
  }
  
  //beta prior
  priors[nparameters-2] = beta_prior;
  parameterName[nparameters-2] = "beta";

  //alpha prior
  priors[nparameters-1] = alpha_prior;
  parameterName[nparameters-1] = "alpha";
   
  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&covariance_wtheta_bias_MF_theory, nparameters, parameterType, parameterName, m_data_model));

}
