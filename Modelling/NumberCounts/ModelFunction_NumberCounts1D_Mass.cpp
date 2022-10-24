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
 *  Modelling/NumberCounts/ModelFunction_NumberCounts1D_Mass.cpp
 *
 *  @brief Functions to model the mass number counts
 *
 *  This file contains the implementation of the functions used to
 *  model mass number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts1D_Mass.h"

using namespace std;

using namespace cbl;

// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::mass_function_mass (const std::vector<double> mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
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
  std::vector<double> Pk = cosmo.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  return cbl::modelling::numbercounts::mass_function (mass, cosmo, pp->redshift, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_critical, pp->kk, Pk, "Spline", pp->k_max);
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_density_mass (const std::vector<double> mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
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
  std::vector<double> Pk = cosmo.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  std::vector<double> mass_function = cbl::modelling::numbercounts::mass_function (pp->Mass_vector, cosmo, pp->redshift, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_critical, pp->kk, Pk, "Spline", pp->k_max);

  glob::FuncGrid interp_MF (pp->Mass_vector, mass_function, "Spline");

  std::vector<double> number_density(mass.size());

  for (size_t i=0; i<mass.size(); i++) 
    number_density[i] = interp_MF.integrate_qag(pp->edges_x[i], pp->edges_x[i+1]);

  return number_density;
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_counts_mass (const std::vector<double> mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)mass;
  
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
  std::vector<double> Pk = cosmo.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (pp->Mass_vector, cosmo, pp->kk, Pk, "Spline", pp->k_max);

  std::vector<double> number_counts(pp->edges_x.size()-1);

  for (size_t j=0; j<pp->edges_x.size()-1; j++) {
    number_counts[j] = cbl::modelling::numbercounts::number_counts(pp->z_min, pp->z_max, pp->edges_x[j], pp->edges_x[j+1], cosmo, pp->area_rad, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_critical, interp[0], interp[1]);
  }

  return number_counts;
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_counts_mass_snapshot (const std::vector<double> mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)mass;
  
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
  std::vector<double> Pk = cosmo.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  std::vector<double> mass_function = cbl::modelling::numbercounts::mass_function (pp->Mass_vector, cosmo, pp->redshift, pp->model_MF, pp->store_output, pp->Delta, pp->isDelta_critical, pp->kk, Pk, "Spline", pp->k_max);

  glob::FuncGrid interp_MF (pp->Mass_vector, mass_function, "Spline");

  std::vector<double> number_counts(pp->edges_x.size()-1);

  for (size_t i=0; i<pp->edges_x.size()-1; i++) {
    number_counts[i] = pp->Volume*interp_MF.integrate_qag(pp->edges_x[i], pp->edges_x[i+1]); //*(maxM-minM);
  }

  return number_counts;
}

