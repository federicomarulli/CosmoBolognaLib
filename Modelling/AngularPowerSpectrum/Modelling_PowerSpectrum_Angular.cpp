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
 *  Modelling/AngularPowerSpectrum/Modelling_PowerSpectrum_Angular.cpp
 *
 *  @brief Methods of the class
 *  Modelling_PowerSpectrum_angular
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_PowerSpectrum_angular, used to model the
 *  angular power spectrum
 *
 *  @authors Federico Marulli, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */


#include "Modelling_PowerSpectrum_Angular.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::modelling::angularpk::Modelling_PowerSpectrum_angular::set_data_model (const cbl::cosmology::Cosmology cosmology, const double z_min, const double z_max, const std::string method_Pk, const bool NL, const int norm, const double k_min, const double k_max, std::vector<double> dN_par, const double fsky, std::vector<double> ll, std::vector<std::vector<double>> mixing_matrix)
{
  m_data_model = make_shared<STR_data_model>(STR_data_model());
  
  m_data_model->cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model->z_min=z_min;
  m_data_model->z_max=z_max;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->norm = norm;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->dN_par=dN_par;
  m_data_model->fsky=fsky;
  m_data_model->ll=ll;
  m_data_model->mixing_matrix=mixing_matrix;

  if(norm==1){
    double s8=m_data_model->cosmology->sigma8_Pk(method_Pk, 0., false, "test", NL, k_min, k_max);
    m_data_model->cosmology->set_sigma8(s8);
  }
}


// ============================================================================================


void cbl::modelling::angularpk::Modelling_PowerSpectrum_angular::set_model_limber (const statistics::PriorDistribution bias_prior, cbl::statistics::PriorDistribution offset_slope_prior)
{
  // set the model parameters
  const int nparameters = 3;
  vector<statistics::ParameterType> parameterType(nparameters-2, statistics::ParameterType::_Base_);
  parameterType.emplace_back(cbl::statistics::ParameterType::_Correlated_);  
  parameterType.emplace_back(cbl::statistics::ParameterType::_Correlated_);   //offset and slope are correlated parameters

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters-1);
  priors[0] = bias_prior;
  priors[1] = offset_slope_prior;      //offset_prior and slope_prior;
  
  parameterName[0] = "bias";
  parameterName[1] = "offset";
  parameterName[2] = "slope";

  //set the priors
  m_set_prior(priors);
  
  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Cl_limber, nparameters, parameterType, parameterName, m_data_model));

}


// =================================================================================================


void cbl::modelling::angularpk::Modelling_PowerSpectrum_angular::set_model_limber (const statistics::PriorDistribution bias_prior)
{
  // set the model parameters
  const int nparameters = 1;
  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters);
  priors[0] = bias_prior;
  
  parameterName[0] = "bias";

  //set the priors
  m_set_prior(priors);
  
  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&Cl_limber, nparameters, parameterType, parameterName, m_data_model));

}


// =================================================================================================



