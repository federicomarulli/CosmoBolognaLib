/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Modelling/Modelling_TwoPointCorrelation_monopole.cpp
 *
 *  @brief Methods of the class Modelling, used for modelling the monopole
 *  of the 2pcf
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_monopole
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling_TwoPointCorrelation_monopole.h"

using namespace cosmobl;


// ============================================================================================


cosmobl::modelling::Modelling_TwoPointCorrelation_monopole::Modelling_TwoPointCorrelation_monopole (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
{
  m_data = twop->dataset();
}


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_monopole::set_fiducial_twop ()
{
  cout << "Setting up the fiducial two point correlation function model" << endl;
  
  m_twop_parameters.fiducial_twop.erase(m_twop_parameters.fiducial_twop.begin(), m_twop_parameters.fiducial_twop.end());

  if (m_twop_parameters.sigmaNL ==0) {
    for(size_t i=0; i<m_twop_parameters.model_scales.size(); i++) {
      double xi = m_twop_parameters.cosmology->xi_DM(m_twop_parameters.model_scales[i], m_twop_parameters.method, m_twop_parameters.redshift, m_twop_parameters.output_root, m_twop_parameters.NL, m_twop_parameters.norm, m_twop_parameters.k_min, m_twop_parameters.k_max, m_twop_parameters.aa, m_twop_parameters.GSL, m_twop_parameters.prec, m_twop_parameters.file_par); 
      m_twop_parameters.fiducial_twop.push_back(xi);
    }
  }
  
  else {
    vector<double> kk = logarithmic_bin_vector(500, m_twop_parameters.k_min+1.e-4, m_twop_parameters.k_max), Pk;

    for(size_t i=0; i<kk.size(); i++)
      Pk.push_back(m_twop_parameters.cosmology->Pk_DeWiggle (kk[i], m_twop_parameters.redshift, m_twop_parameters.sigmaNL, m_twop_parameters.output_root, m_twop_parameters.norm, m_twop_parameters.k_min, m_twop_parameters.k_max, m_twop_parameters.aa, m_twop_parameters.prec));

    m_twop_parameters.fiducial_twop = Xi0(m_twop_parameters.model_scales, kk, Pk);
  }

  m_twop_parameters.func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(m_twop_parameters.model_scales, m_twop_parameters.fiducial_twop, "Spline"));
  
  cout << "Done!" << endl;
}


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_monopole::set_model_AP_isotropic (const statistics::Prior alpha_prior, const statistics::Prior B_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior, const statistics::ParameterType pT_alpha, const statistics::ParameterType pT_B, const statistics::ParameterType pT_A0, const statistics::ParameterType pT_A1, const statistics::ParameterType pT_A2)
{
  set_fiducial_twop();

  vector<statistics::Parameter> model_parameters;

  statistics::Parameter alpha(alpha_prior.sample(), alpha_prior, pT_alpha, "alpha");
  statistics::Parameter B(B_prior.sample(), B_prior, pT_B, "B");
  statistics::Parameter A0(A0_prior.sample(), A0_prior, pT_A0, "A0");
  statistics::Parameter A1(A1_prior.sample(), A1_prior, pT_A1, "A1");
  statistics::Parameter A2(A2_prior.sample(), A2_prior, pT_A1, "A2");

  model_parameters.push_back(B);
  model_parameters.push_back(alpha);
  model_parameters.push_back(A0);
  model_parameters.push_back(A1);
  model_parameters.push_back(A2);

  auto fixed_parameters = make_shared<glob::STR_twop_model>(m_twop_parameters);

  m_model = make_shared<statistics::Model1D>(statistics::Model1D(model_parameters, fixed_parameters, &glob::xi_alpha_B_poly)); 
}
