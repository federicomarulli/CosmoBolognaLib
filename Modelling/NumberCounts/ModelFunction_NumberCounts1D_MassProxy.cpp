/********************************************************************
 *  Copyright (C) 2021 by Giorgio Lesci and Federico Marulli        *
 *  giorgio.lesci2@unibo.it, federico.marulli3@unibo.it             *
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
 *  Modelling/NumberCounts/ModelFunction_NumberCounts1D_MassProxy.cpp
 *
 *  @brief Functions to model the mass number counts as a function of a mass proxy
 *
 *  This file contains the implementation of the functions used to
 *  model mass number counts as a function of mass proxy
 *
 *  @authors Giorgio Lesci, Federico Marulli
 *
 *  @authors giorgio.lesci2@unibo.it, federico.marulli3@unibo.it
 */


#include "ModelFunction_NumberCounts1D_MassProxy.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_counts_proxy (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)proxy;
  
  // structure contaning the required input data
  shared_ptr<STR_NC_data_model> pp = static_pointer_cast<STR_NC_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  // compute the power spectrum
  std::vector<double> Pk = cosmo.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  // interpolate sigmaM and its derivative
  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (pp->Mass_vector, cosmo, pp->kk, Pk, "Spline", pp->k_max);

  // interpolate the growth factor
  const std::vector<double> z_for_DN = cbl::linear_bin_vector(200, 0.0001, 3.);
  std::vector<double> DN (z_for_DN.size(), 0.);
  for (size_t i=0; i<z_for_DN.size(); i++)
    DN[i] = cosmo.DN(z_for_DN[i]);
  cbl::glob::FuncGrid DN_interp (z_for_DN, DN, "Spline");

  // Compute the counts
  std::vector<double> number_counts(pp->edges_x.size()-1);

  for (size_t j=0; j<pp->edges_x.size()-1; j++) {
    number_counts[j] = cbl::modelling::numbercounts::counts_proxy(parameter[pp->Cpar.size()], parameter[pp->Cpar.size()+1], parameter[pp->Cpar.size()+2], parameter[pp->Cpar.size()+3], parameter[pp->Cpar.size()+4], parameter[pp->Cpar.size()+5], parameter[pp->Cpar.size()+6], parameter[pp->Cpar.size()+7], parameter[pp->Cpar.size()+8], parameter[pp->Cpar.size()+9], parameter[pp->Cpar.size()+10], parameter[pp->Cpar.size()+11], parameter[pp->Cpar.size()+12], parameter[pp->Cpar.size()+13], parameter[pp->Cpar.size()+14], pp->fz, pp->z_error, pp->proxy_error, pp->response_fact, pp->z_min, pp->z_max, pp->edges_x[j], pp->edges_x[j+1], cosmo, pp->area_rad, pp->model_MF, pp->model_bias, pp->store_output, pp->Delta, pp->isDelta_critical, interp[0], interp[1], DN_interp, pp->proxy_pivot, pp->z_pivot, pp->mass_pivot, pp->log_base, pp->SF_weights[j]);
  }

  return number_counts;
}
