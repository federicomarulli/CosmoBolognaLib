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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_projected.cpp
 *
 *  @brief Functions to model the projected two-point correlation
 *  function
 *
 *  This file contains the implementation of the functions used to
 *  model the projected two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_TwoPointCorrelation_projected.h"

using namespace cosmobl;


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_from_xi_approx (function<vector<double>(const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter)> func, const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  vector<double> wp(rp.size(), 0);

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<wp.size(); i++) {
      
      auto integrand = [&] (double rad) { return func({rad}, inputs, parameter)[0]/sqrt(rad*rad-rp[i]*rp[i]); };
      
      double r_out = sqrt(pow(rp[i], 2)+pow(pp->r_max_int, 2));
      wp[i] = 2.*gsl::GSL_integrate_qag(integrand, rp[i], r_out);
      
    }

  }

  return wp;
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_1halo_approx (const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  return wp_from_xi_approx(modelling::twopt::xi_1halo, rp, inputs, parameter); 
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_2halo_approx (const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  return wp_from_xi_approx(modelling::twopt::xi_2halo, rp, inputs, parameter); 
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_HOD_approx (const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  return wp_from_xi_approx(modelling::twopt::xi_HOD, rp, inputs, parameter); 
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_from_xi (function<double(const double rp, const double pi, const shared_ptr<void> inputs, vector<double> &parameter)> func, const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  // structure contaning the required input data
   shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);
  
  vector<double> wp(rp.size(), 0);

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<wp.size(); i++) {
      
      auto integrand = [&] (double pi) { return func(rp[i], pi, inputs, parameter); };
      
      wp[i] = 2.*gsl::GSL_integrate_qag(integrand, 0., pp->pi_max);
      
    }

  }

  return wp;
}
						       

// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_1halo (const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  return wp_from_xi(modelling::twopt::xi_1halo_zspace, rp, inputs, parameter); 
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_2halo (const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  return wp_from_xi(modelling::twopt::xi_2halo_zspace, rp, inputs, parameter); 
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_HOD (const vector<double> rp, const shared_ptr<void> inputs, vector<double> &parameter)
{
  return wp_from_xi(modelling::twopt::xi_HOD_zspace, rp, inputs, parameter); 
}
