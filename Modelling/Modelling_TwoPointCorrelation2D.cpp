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
 *  @file Modelling/Modelling_TwoPointCorrelation2D.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation2D
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation2D
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling_TwoPointCorrelation2D.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation2D::set_model_dispersionModel (const double alpha_perp_guess, const statistics::Prior alpha_perp_prior, const statistics::ParameterType alpha_perp_type, const double alpha_par_guess, const statistics::Prior alpha_par_prior, const statistics::ParameterType alpha_par_type, const double fsigma8_guess, const statistics::Prior fsigma8_prior, const statistics::ParameterType fsigma8_type, const double bsigma8_guess, const statistics::Prior bsigma8_prior, const statistics::ParameterType bsigma8_type, const double sigma12_guess, const statistics::Prior sigma12_prior, const statistics::ParameterType sigma12_type)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  
  // ----- set the parameter of the model xi0_linear (implemented in ModelFunction.cpp) -----

  statistics::Parameter alpha_perp((alpha_perp_guess>par::defaultDouble) ? alpha_perp_guess : alpha_perp_prior.sample(), alpha_perp_prior, alpha_perp_type, "alpha_perp");
  statistics::Parameter alpha_par((alpha_par_guess>par::defaultDouble) ? alpha_par_guess : alpha_par_prior.sample(), alpha_par_prior, alpha_par_type, "alpha_par");
  statistics::Parameter fsigma8((fsigma8_guess>par::defaultDouble) ? fsigma8_guess : fsigma8_prior.sample(), fsigma8_prior, fsigma8_type, "f*sigma8");
  statistics::Parameter bsigma8((bsigma8_guess>par::defaultDouble) ? bsigma8_guess : bsigma8_prior.sample(), bsigma8_prior, bsigma8_type, "bias*sigma8");
  statistics::Parameter sigma12((sigma12_guess>par::defaultDouble) ? sigma12_guess : sigma12_prior.sample(), sigma12_prior, sigma12_type, "bias*sigma8");

  vector<statistics::Parameter> model_parameters = {alpha_perp, alpha_par, fsigma8, bsigma8, sigma12};
  

  // input data used to construct the model
  auto inputs = make_shared<STR_twop_model>(m_twop_parameters);

  // construct the model
  m_model = make_shared<statistics::Model2D>(statistics::Model2D(model_parameters, inputs, &xi2D_dispersionModel)); 
}
