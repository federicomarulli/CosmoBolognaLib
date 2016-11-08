/*******************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file Modelling/ModelFunction.cpp
 *
 *  @brief implementation of the function for data modelling
 *
 *  This file contains the implementation of the functions used to
 *  model any kind of data 
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "ModelFunction.h"


// ============================================================================================


double cosmobl::modelling::xi0_linear (const double rad, const shared_ptr<void> inputs, vector<double> parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(inputs);

  
  // ----- input parameters -----

  // AP parameter that contains the distance information
  double alpha = parameter[0];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[1];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[2];

  // polynomial parameter A0 (used for BAO analyses)
  double A0 = parameter[3];

  // polynomial parameter A1 (used for BAO analyses)
  double A1 = parameter[4];

  // polynomial parameter A2 (used for BAO analyses)
  double A2 = parameter[5];

  
  // ----- derived parameters -----

  // rescaled radius
  double new_rad = alpha*rad;

  // polynomial used to marginalize over signals caused by systematics
  // not fully taken into account (see e.g. Anderson et al. 2012, and
  // reference therein)
  double poly = A0+A1/rad+A2/(rad*rad);
  

  // return the redshift-space monopole of the two-point correlation function
  return xi_ratio(fsigma8, bsigma8)*pow(bsigma8/pp->sigma8_z, 2)*pp->func_xi->operator()(new_rad)+poly;
  
}


// ============================================================================================


double cosmobl::modelling::xi0_linear_cosmology (const double rad, const shared_ptr<void> inputs, vector<double> parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(inputs);

  
  // ----- input parameters -----

  // AP parameter that contains the distance information
  double alpha = parameter[0];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[1];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[2];

  // polynomial parameter A0 (used for BAO analyses)
  double A0 = parameter[3];

  // polynomial parameter A1 (used for BAO analyses)
  double A1 = parameter[4];

  // polynomial parameter A2 (used for BAO analyses)
  double A2 = parameter[5];

  
  // ----- derived parameters -----

  // rescaled radius
  double new_rad = alpha*rad;

  // polynomial used to marginalize over signals caused by systematics
  // not fully taken into account (see e.g. Anderson et al. 2012, and
  // reference therein)
  double poly = A0+A1/rad+A2/(rad*rad);
  
  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<parameter.size(); ++i)
    pp->cosmology->set_parameter(pp->Cpar[i], parameter[i]);

  
  // return the redshift-space monopole of the two-point correlation function
  return xi_ratio(fsigma8, bsigma8)*pp->cosmology->xi_DM(new_rad, pp->method_Pk, pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par)/pow(pp->sigma8_z, 2)+poly;
  
}


// ============================================================================================


double cosmobl::modelling::xi2D_dispersionModel (const double rp, const double pi, const shared_ptr<void> inputs, vector<double> parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(inputs);

  
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

  
  // return the 2D correlation function in Cartesian coordinates modelled with the dispersion model
  return xi2D_model(AP1*rp, AP2*pi, beta, bias, sigma12, pp->funcXiR, pp->funcXiR_, pp->funcXiR__, pp->var, pp->FV, pp->bias_nl, pp->bA, pp->v_min, pp->v_max, pp->step_v);
  
}


