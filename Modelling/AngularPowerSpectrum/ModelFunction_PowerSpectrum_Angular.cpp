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
 *  Modelling/AngularPowerSpectrum/ModelFunction_PowerSpectrum_Angular.cpp
 *
 *  @brief Functions to model the angular power spectrum
 *
 *  This file contains the implementation of the functions used to
 *  model the angular power spectrum
 *
 *  @authors Federico Marulli, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */


#include "ModelFunction_PowerSpectrum_Angular.h"

using namespace std;

using namespace cbl;


// ============================================================================================


double cbl::modelling::angularpk::integrand_limber(double redshift, std::shared_ptr<void> inputs, std::vector<double> par){
  
  // structure contaning the required input data

  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;
  //std::vector<double> dN_par= pp->dN_par;
  double distribution=0;
  for (size_t i=0; i<pp->dN_par.size(); ++i)
    distribution += pp->dN_par[i]*pow(redshift, i);  
  distribution+=par[1]+par[2]*redshift;
  //(void)pp;
  double l=par[0];
  double _cc=1/cbl::par::cc;
  double Pk=cosmo.Pk_matter({double((l+0.5)/cosmo.D_C(redshift))}, pp->method_Pk, pp->NL, redshift, false, "test", pp->norm, pp->k_min, pp->k_max)[0];
  double inv_d2= 1/(cosmo.D_C(redshift)*cosmo.D_C(redshift));
  return distribution*distribution*Pk*inv_d2*cosmo.HH(redshift)*_cc;
  
}


// =====================================================================


double cbl::modelling::angularpk::integral_limber (const double l, double offset, double slope, std::shared_ptr<void> inputs){

  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  double lower_limit = pp->z_min;
  double upper_limit = pp->z_max;

  // relative measurement precision (it raises a warning message if not satisfied)
  const double prec = 1.e-5;
  
  // default values (refer to the GSL documentation for different choiches)
  const int limit_size = 1000;
  const int rule = 6;
  
  return cbl::wrapper::gsl::GSL_integrate_qag(&integrand_limber,inputs,{l, offset, slope}, lower_limit, upper_limit, prec, limit_size, rule);
   
}


// ===========================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_limber (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
 
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  vector<double> Cl(l.size(), 0);

  // input parameters
  double bias = parameter[0];
  double offset = 0;
  double slope = 0;
  
  if(parameter.size()>1){
    offset = parameter[1];
    slope = parameter[2];
  }
  for (size_t i =0; i<l.size(); i++) Cl[i] = cbl::modelling::angularpk::integral_limber(l[i], offset, slope, inputs);
  for (size_t i =0; i<Cl.size(); i++) Cl[i] *= bias*bias;
  return Cl;
   
}


// ============================================================================================




