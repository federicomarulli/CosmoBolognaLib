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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation2D_cartesian.cpp
 *
 *  @brief Functions to model the 2D two-point correlation
 *  function in Cartesian coordinates
 *
 *  This file contains the implementation of the functions used to
 *  model the 2D two-point correlation function in Cartesian
 *  coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_TwoPointCorrelation2D_cartesian.h"

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::xi2D_dispersionModel (const std::vector<double> rp, const std::vector<double> pi, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);


  // ----- input parameters -----

  // AP parameter: D_A,1(z)/D_A,2(z)
  double AP1 = parameter[0];

  // AP parameter: H_2(z)/H_1(z)
  double AP2 = parameter[1];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[2];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[3];

  // the dispersion in the pairwise random peculiar velocities
  double sigma12 = parameter[4];


  // ----- derived parameters -----

  // the linear bias
  double bias = bsigma8/pp->sigma8_z;
  
  // the distortion parameter
  double beta = fsigma8/bsigma8;
  
  vector<vector<double>> model(rp.size(), vector<double>(pi.size(), 0.));

  // return the 2D correlation function in Cartesian coordinates modelled with the dispersion model
  if (sigma12 == 0)
    for (size_t i=0; i<rp.size(); i++)
      for (size_t j=0; j<pi.size(); j++) 
	model[i][j] = xi2D_lin_model(AP1*rp[i], AP2*pi[j], beta, bias, pp->func_xi, pp->func_xi_, pp->func_xi__, pp->bias_nl, pp->bA);
  else
    for (size_t i=0; i<rp.size(); i++)
      for (size_t j=0; j<pi.size(); j++) 
	model[i][j] = xi2D_model(AP1*rp[i], AP2*pi[j], beta, bias, sigma12, pp->func_xi, pp->func_xi_, pp->func_xi__, pp->var, pp->FV, pp->bias_nl, pp->bA, pp->v_min, pp->v_max, pp->step_v);

  return model;
}
