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
 *  Modelling/NumberCounts/ModelFunction_NumberCounts1D_Redshift.cpp
 *
 *  @brief Functions to model the redshift number counts
 *
 *  This file contains the implementation of the functions used to
 *  model redshift number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts1D_Redshift.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_density_redshift (const std::vector<double> redshift, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_NC_data_model> pp = static_pointer_cast<STR_NC_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // input likelihood parameters

  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  // compute the power spectrum
  std::vector<double> Pk = cosmo.Pk(pp->kk, pp->method_Pk, false, 0., pp->output_dir, pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  std::vector<std::vector<double>> mass_function = cbl::modelling::numbercounts::mass_function(redshift, pp->Mass_vector, cosmo, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_Vir, pp->kk, Pk, "Spline", pp->k_max);

  std::vector<double> number_density(redshift.size());

  for (size_t i=0; i<redshift.size(); i++) {
    glob::FuncGrid interpMF(pp->Mass_vector, mass_function[i], "Spline");
    number_density[i] = interpMF.integrate_qag(pp->Mass_min, pp->Mass_max);
  }

  return number_density;
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_counts_redshift (const std::vector<double> redshift, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_NC_data_model> pp = static_pointer_cast<STR_NC_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // input likelihood parameters

  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<pp->Cpar.size(); ++i) {
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  }

  // compute the power spectrum
  std::vector<double> Pk = cosmo.Pk(pp->kk, pp->method_Pk, false, 0., pp->output_dir, pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  std::vector<std::vector<double>> mass_function = cbl::modelling::numbercounts::mass_function (redshift, pp->Mass_vector, cosmo, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_Vir, pp->kk, Pk, "Spline", pp->k_max);

  std::vector<double> number_counts(redshift.size());
  double deltaZ = redshift[1]-redshift[0];

  for (size_t i=0; i<redshift.size(); i++) {
    glob::FuncGrid interpMF(pp->Mass_vector, mass_function[i], "Spline");
    number_counts[i] = pp->area_rad*interpMF.integrate_qag(pp->Mass_min, pp->Mass_max)*cosmo.dV_dZdOmega(redshift[i], true)*deltaZ;
  }

  /*
  std::vector<std::vector<double>> mass_function = cbl::modelling::numbercounts::mass_function (pp->z_std::vector, pp->Mass_vector, cosmo, pp->model_MF, pp->Delta, pp->isDelta_Vir, pp->kk, Pk, "Spline", pp->k_max);

  std::vector<double> dV_dzdOmega (pp->z_vector.size());

  for(size_t i=0; i<pp->z_vector.size(); i++)
    dV_dzdOmega[i] = pp->area_rad*cosmo.dV_dZdOmega(pp->z_vector[i], true);

  std::vector<double> number_density(redshift.size());

  for (size_t i=0; i<pp->z_vector.size(); i++) {
    glob::FuncGrid interpMF(pp->Mass_vector, mass_function[i], "Linear");
    number_density[i] = interpMF.integrate_qag(pp->Mass_min, pp->Mass_max);
  }

  glob::FuncGrid interp_dVdZdOmega(pp->z_vector, dV_dzdOmega, "Linear");
  glob::FuncGrid interp_nz(pp->z_vector, number_density, "Linear");

  std::vector<double> number_counts(redshift.size());
  double deltaZ = redshift[1]-redshift[0];

  for (size_t i=0; i<redshift.size(); i++) {

    auto integrand = [&] (const double z) 
    {
      return interp_nz(z)*interp_dVdZdOmega(z);
    };

    number_counts[i] = cbl::gsl::GSL_integrate_qag(integrand, redshift[0]+i*deltaZ, redshift[0]+(i+1)*deltaZ, 1.e-4);
  }
  */

  return number_counts;
}
