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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation_projected.h"

using namespace cosmobl;

// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_projected::set_model_linearBias (const statistics::Prior bsigma8_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_wpDM();

  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(statistics::Prior(glob::DistributionType::_ConstantDistribution_, 1), "alpha"));

  m_fsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(statistics::Prior(glob::DistributionType::_ConstantDistribution_, 0.), "f*sigma8"));
  
  m_bsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bsigma8_prior, "b*sigma8"));
  
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_alpha, m_fsigma8, m_bsigma8};
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear, inputs));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_projected::set_fiducial_wpDM ()
{
  coutCBL << "Setting up the fiducial model for the projected correlation function of the dark matter" << endl;

  const vector<double> rad = logarithmic_bin_vector(m_data_model.step, max(m_data_model.r_min, 1.e-4), min(m_data_model.r_max, 100.));
  vector<double> wpDM(m_data_model.step);
  
  for (size_t i=0; i<(size_t)m_data_model.step; i++) 
    wpDM[i] = m_data_model.cosmology->wp_DM(rad[i], m_data_model.method_Pk, m_data_model.redshift, m_data_model.pi_max, m_data_model.output_root, m_data_model.NL, m_data_model.norm, m_data_model.r_min, m_data_model.r_max, m_data_model.k_min, m_data_model.k_max, m_data_model.aa, m_data_model.GSL, m_data_model.prec, m_data_model.file_par);

  m_data_model.func_xi = make_shared<glob::FuncGrid>(glob::FuncGrid(rad, wpDM, "Spline"));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_projected::set_model_HOD (const statistics::Prior Mmin_prior, const statistics::Prior sigmalgM_prior, const statistics::Prior M0_prior, const statistics::Prior M1_prior, const statistics::Prior alpha_prior)
{
  // compute the fiducial dark matter power spectrum
  set_fiducial_PkDM();

  // compute the fiducial mass variance and its logarithmic derivative
  set_fiducial_sigma();
  
  // set the parameter of the HOD model

  m_Mmin = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Mmin_prior, "Mmin"));

  m_sigmalgM = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigmalgM_prior, "sigmalgM"));

  m_M0 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(M0_prior, "M0"));
  
  m_M1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(M1_prior, "M1"));

  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));

  set_parameters({m_Mmin, m_sigmalgM, m_M0, m_M1, m_alpha});

  
  // input data used to construct the model
  auto inputs = make_shared<STR_data_HOD>(m_data_HOD);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&wp_HOD, inputs));
  
}
