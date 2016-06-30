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
 *  @file Modelling/Modelling_TwoPointCorrelation_projected.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation_projected, 
 *  used for modelling projected 2pcf
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_projected
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling_TwoPointCorrelation_projected.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_projected::set_fiducial_twop ()
{
  cout << "Setting up the fiducial two-point correlation function model" << endl;
  m_twop_parameters.fiducial_twop.erase(m_twop_parameters.fiducial_twop.begin(), m_twop_parameters.fiducial_twop.end());

  for(size_t i=0; i<m_twop_parameters.model_scales.size(); i++){
    double wp = m_twop_parameters.cosmology->wp_DM(m_twop_parameters.model_scales[i], m_twop_parameters.method, m_twop_parameters.redshift, m_twop_parameters.output_root, m_twop_parameters.NL, m_twop_parameters.norm, m_twop_parameters.r_min, m_twop_parameters.r_max, m_twop_parameters.k_min, m_twop_parameters.k_max, m_twop_parameters.aa, m_twop_parameters.GSL, m_twop_parameters.prec, m_twop_parameters.file_par); 
    m_twop_parameters.fiducial_twop.push_back(wp);
  }
  
  m_twop_parameters.func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(m_twop_parameters.model_scales, m_twop_parameters.fiducial_twop, "Spline"));
  cout << "Done!" << endl;

}

// ============================================================================================


cosmobl::modelling::Modelling_TwoPointCorrelation_projected::Modelling_TwoPointCorrelation_projected(const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
{
  m_data = twop->dataset();
}


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation_projected::set_model_bias(const statistics::Prior bias_prior, const statistics::ParameterType pT_bias)
{
  set_fiducial_twop();
  statistics::Parameter bias(bias_prior.sample(), bias_prior, pT_bias, "bias");

  vector<statistics::Parameter> model_parameters;
  model_parameters.push_back(bias);

  auto fixed_parameters = make_shared<glob::STR_twop_model>(m_twop_parameters);

  m_model = make_shared<statistics::Model1D>(statistics::Model1D(model_parameters, fixed_parameters, &glob::wp_bias)); 
}
