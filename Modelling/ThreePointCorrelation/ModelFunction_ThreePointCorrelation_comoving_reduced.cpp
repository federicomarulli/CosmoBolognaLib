/*******************************************************************
 *  Copyright (C) 2017 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/** 
 *  @file
 *  Modelling/ThreePointCorrelation/ModelFunction_ThreePointCorrelation_comoving_reduced.cpp
 *
 *  @brief Functions to model the reduced three-point correlation
 *  function in comoving coordinates
 *
 *  This file contains the implementation of the functions used to
 *  model the reduced three-point correlation function in comoving
 *  coordinates
 *
 *
 *  @authors Federico Marulli, Michele Moresco
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it
 */

#include "ModelFunction_ThreePointCorrelation_comoving_reduced.h"

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<double> cbl::modelling::threept::Q_nonlinear_localbias (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)theta;

  // structure contaning the required input data
  shared_ptr<STR_data_model_threept> pp = static_pointer_cast<STR_data_model_threept>(inputs);

  // input likelihood parameters

  // bias
  double bias1 = parameter[0];
  double bias2 = parameter[1];
  
  vector<double> Q_nl_lb(pp->Q_DM.size());
  for (size_t i=0; i<Q_nl_lb.size(); ++i) 
  	Q_nl_lb[i] = 1./bias1*(pp->Q_DM[i]+bias2/bias1);

  return Q_nl_lb;
}


// ============================================================================================

std::vector<double> cbl::modelling::threept::Q_nonlinear_nonlocalbias (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)theta;

  // structure contaning the required input data
  shared_ptr<STR_data_model_threept> pp = static_pointer_cast<STR_data_model_threept>(inputs);

  // input likelihood parameters

  // bias
  double bias1 = parameter[0];
  double bias2 = parameter[1];
  double g2 = parameter[2];

  vector<double> Q_nl_nlb(pp->theta.size());
  Q_nl_nlb = pp->cosmology->Q_halo (pp->r1, pp->r2, pp->theta, bias1, bias2, g2, pp->model, pp->kk, pp->Pk_DM);

  return Q_nl_nlb;
}


// ============================================================================================


std::vector<double> cbl::modelling::threept::Q_nonlinear_nonlocalbias_alpha (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)theta;

  // structure contaning the required input data
  shared_ptr<STR_data_model_threept> pp = static_pointer_cast<STR_data_model_threept>(inputs);

  // input likelihood parameters

  // bias
  double bias1 = parameter[0];
  double bias2 = parameter[1];
  double g2 = parameter[2];
  double alpha = parameter[3];

  vector<double> Q_nl_nlb_alpha(pp->theta.size());

  Q_nl_nlb_alpha = pp->cosmology->Q_halo (pp->r1*alpha, pp->r2*alpha, pp->theta, bias1, bias2, g2, pp->model, pp->kk, pp->Pk_DM);

  return Q_nl_nlb_alpha;
}
