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
 *  @file Modelling/Modelling_TwoPointCorrelation1D.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation1D
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation1D.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation1D::set_model_monopole (const double alpha_guess, const statistics::Prior alpha_prior, const statistics::ParameterType alpha_type, const double fsigma8_guess, const statistics::Prior fsigma8_prior, const statistics::ParameterType fsigma8_type, const double bsigma8_guess, const statistics::Prior bsigma8_prior, const statistics::ParameterType bsigma8_type, const double A0_guess, const statistics::Prior A0_prior, const statistics::ParameterType A0_type, const double A1_guess, const statistics::Prior A1_prior, const statistics::ParameterType A1_type, const double A2_guess, const statistics::Prior A2_prior, const statistics::ParameterType A2_type)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  
  // ----- set the parameter of the model xi0_linear (implemented in ModelFunction.cpp) -----

  statistics::Parameter alpha((alpha_guess>par::defaultDouble) ? alpha_guess : alpha_prior.sample(), alpha_prior, alpha_type, "alpha");
  statistics::Parameter fsigma8((fsigma8_guess>par::defaultDouble) ? fsigma8_guess : fsigma8_prior.sample(), fsigma8_prior, fsigma8_type, "f*sigma8");
  statistics::Parameter bsigma8((bsigma8_guess>par::defaultDouble) ? bsigma8_guess : bsigma8_prior.sample(), bsigma8_prior, bsigma8_type, "bias*sigma8");
  statistics::Parameter A0((A0_guess>par::defaultDouble) ? A0_guess : A0_prior.sample(), A0_prior, A0_type, "A0");
  statistics::Parameter A1((A1_guess>par::defaultDouble) ? A1_guess : A1_prior.sample(), A1_prior, A1_type, "A1");
  statistics::Parameter A2((A2_guess>par::defaultDouble) ? A2_guess : A2_prior.sample(), A2_prior, A2_type, "A2");

  vector<statistics::Parameter> model_parameters = {alpha, fsigma8, bsigma8, A0, A1, A2};
  

  // input data used to construct the model
  auto inputs = make_shared<STR_twop_model>(m_twop_parameters);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(model_parameters, inputs, &xi0_linear));

}
