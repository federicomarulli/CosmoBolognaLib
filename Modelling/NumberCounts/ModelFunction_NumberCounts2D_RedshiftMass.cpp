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
 *  Modelling/NumberCounts/ModelFunction_NumberCounts2D_RedshiftMass.cpp
 *
 *  @brief Functions to model the redshift-mass number counts
 *
 *  This file contains the implementation of the functions used to
 *  model redshift-mass number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts2D_RedshiftMass.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::mass_function_mass_redshift (const std::vector<double> redshift, const std::vector<double> mass, const shared_ptr<void> inputs, std::vector<double> &parameter)
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
  std::vector<double> Pk = cosmo.Pk(pp->kk, pp->method_Pk, false, 0., pp->output_dir, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par);

  return cbl::modelling::numbercounts::mass_function (redshift, mass, cosmo, pp->model_MF, pp->Delta, pp->isDelta_Vir, pp->kk, Pk, "Spline", pp->k_max);
}


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::number_density_mass_redshift (const std::vector<double> redshift, const std::vector<double> mass, const shared_ptr<void> inputs, std::vector<double> &parameter)
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
  std::vector<double> Pk = cosmo.Pk(pp->kk, pp->method_Pk, false, 0., pp->output_dir, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par);

  std::vector<std::vector<double>> mass_function = cbl::modelling::numbercounts::mass_function (pp->z_vector, pp->Mass_vector, cosmo, pp->model_MF, pp->Delta, pp->isDelta_Vir, pp->kk, Pk, "Spline", pp->k_max);

  glob::FuncGrid2D interp_MF (pp->z_vector, pp->Mass_vector, mass_function, "Linear");

  std::vector<std::vector<double>> number_density(redshift.size(), std::vector<double>(mass.size()));

  double deltaz = redshift[1]-redshift[0];
  double logMinMass = log10(mass[0]);
  double deltaLogM = log10(mass[1])-logMinMass;

  for (size_t i=0; i<redshift.size(); i++) {
    std::vector<double> zlimits = {redshift[0]+i*deltaz, redshift[0]+(i+1)*deltaz};
    for (size_t j=0; j<mass.size(); j++) {
      std::vector<double> mlimits = {pow(10., logMinMass+i*deltaLogM), pow(10, logMinMass+(i+1)*deltaLogM)};
      number_density[i][j] = interp_MF.IntegrateVegas(zlimits[0], zlimits[1], mlimits[0], mlimits[1]);
    }
  }

  return number_density;
}


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::number_counts_mass_redshift (const std::vector<double> redshift, const std::vector<double> mass, const shared_ptr<void> inputs, std::vector<double> &parameter)
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

  std::vector<std::vector<double>> number_counts(redshift.size(), std::vector<double>(mass.size()));

  // compute the power spectrum
  const std::vector<double> Pk = cosmo.Pk(pp->kk, pp->method_Pk, false, 0., pp->output_dir, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (pp->Mass_vector, cosmo, pp->kk, Pk, "Spline", pp->k_max);
  
  double deltaz = 0.5*(redshift[1]-redshift[0]);
  double deltaLogM = 0.5*(log10(mass[1])-log10(mass[0]));

  for (size_t i=0; i<redshift.size(); i++){ 
    std::vector<double> zlim = {redshift[i]-deltaz, redshift[i]+deltaz};
    for (size_t j=0; j<mass.size(); j++) {
      std::vector<double> Mlim = {pow(10, log10(mass[j])-deltaLogM), pow(10, log10(mass[j])+deltaLogM)};
      number_counts[i][j] = cbl::modelling::numbercounts::number_counts(zlim[0], zlim[1], Mlim[0], Mlim[1], cosmo, pp->area_rad, pp->model_MF, pp->Delta, pp->isDelta_Vir, interp[0], interp[1], 10, 10);

    }
  }

  return number_counts;
}
