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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_projected.cpp
 *
 *  @brief Methods of the class
 *  Modelling_TwoPointCorrelation_projected
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_projected, used to model the
 *  projected two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation_projected.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_projected::set_model_linearBias (const statistics::PriorDistribution bsigma8_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_wpDM();

  // set the model parameters
  const int nparameters = 3;

  vector<statistics::ParameterType> parameterType = {statistics::ParameterType::_Base_, statistics::ParameterType::_Base_, statistics::ParameterType::_Base_};

  vector<string> parameterName = {"alpha", "f*sigma8", "b*sigma8"};

  statistics::PriorDistribution alpha_prior (glob::DistributionType::_Constant_, 1);
  statistics::PriorDistribution fsigma8_prior (glob::DistributionType::_Constant_, 0);

  vector<statistics::PriorDistribution> priors = {alpha_prior, fsigma8_prior, bsigma8_prior};

  // input data used to construct the model
  m_data_model->poly_order = 0;

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_projected::set_fiducial_wpDM ()
{
  coutCBL << "Setting up the fiducial model for the projected correlation function of the dark matter" << endl;

  const vector<double> rad = logarithmic_bin_vector(m_data_model->step, max(m_data_model->r_min, 1.e-4), min(m_data_model->r_max, 100.));
  vector<double> wpDM(m_data_model->step);
  
  for (size_t i=0; i<(size_t)m_data_model->step; i++) 
    wpDM[i] = m_data_model->cosmology->wp_DM(rad[i], m_data_model->method_Pk, m_data_model->NL, m_data_model->redshift, m_data_model->pi_max, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->r_min, m_data_model->r_max, m_data_model->k_min, m_data_model->k_max, m_data_model->aa, m_data_model->GSL, m_data_model->prec, m_data_model->file_par);
  
  m_data_model->func_xi = make_shared<glob::FuncGrid>(glob::FuncGrid(rad, wpDM, "Spline"));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_projected::set_model_HOD (const statistics::PriorDistribution Mmin_prior, const statistics::PriorDistribution sigmalgM_prior, const statistics::PriorDistribution M0_prior, const statistics::PriorDistribution M1_prior, const statistics::PriorDistribution alpha_prior)
{
  // compute the fiducial dark matter power spectrum
  set_fiducial_PkDM();

  // compute the fiducial mass variance and its logarithmic derivative
  set_fiducial_sigma();

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType = {statistics::ParameterType::_Base_};

  vector<string> parameterName(nparameters);
  parameterName[0] = "Mmin";
  parameterName[1] = "sigmalgM";
  parameterName[2] = "M0";
  parameterName[3] = "M1";
  parameterName[4] = "alpha";

  vector<statistics::PriorDistribution> priors = {Mmin_prior, sigmalgM_prior, M0_prior, M1_prior, alpha_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&wp_HOD, nparameters, parameterType, parameterName, m_data_HOD));
  
}
