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
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  // compute the power spectrum
  cosmology::PkXi PX(cosmo);
  std::vector<double> Pk = PX.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par, true);

  // interpolate sigmaM and its derivative
  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (pp->Mass_vector, cosmo, pp->kk, Pk, "Spline", pp->k_max);

  // interpolate the growth factor
  const std::vector<double> z_for_DN = cbl::linear_bin_vector(200, 0.0001, 3.);
  std::vector<double> DN (z_for_DN.size(), 0.);
  for (size_t i=0; i<z_for_DN.size(); i++)
    DN[i] = cosmo->DN(z_for_DN[i]);
  cbl::glob::FuncGrid DN_interp (z_for_DN, DN, "Spline");

  // Compute the counts
  std::vector<double> number_counts(pp->edges_x.size()-1);

  for (size_t j=0; j<pp->edges_x.size()-1; j++)
    number_counts[j] = cbl::modelling::numbercounts::counts_proxy(parameter[pp->Cpar.size()], parameter[pp->Cpar.size()+1], parameter[pp->Cpar.size()+2], parameter[pp->Cpar.size()+3], parameter[pp->Cpar.size()+4], parameter[pp->Cpar.size()+5], parameter[pp->Cpar.size()+6], parameter[pp->Cpar.size()+7], parameter[pp->Cpar.size()+8], parameter[pp->Cpar.size()+9], parameter[pp->Cpar.size()+10], parameter[pp->Cpar.size()+11], parameter[pp->Cpar.size()+12], parameter[pp->Cpar.size()+13], parameter[pp->Cpar.size()+14], pp->fz, pp->z_error, pp->proxy_error, pp->response_fact, pp->z_min, pp->z_max, pp->edges_x[j], pp->edges_x[j+1], cosmo, pp->area_rad, pp->model_MF, pp->model_bias, pp->store_output, pp->Delta, pp->isDelta_critical, interp[0], interp[1], DN_interp, pp->proxy_pivot, pp->z_pivot, pp->mass_pivot, pp->log_base, pp->SF_weights[j]);

  return number_counts;
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::number_counts_proxy_classic (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  (void)proxy;
  
  // structure contaning the required input data
  shared_ptr<STR_NC_data_model> pp = static_pointer_cast<STR_NC_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  // scaling relation parameters  
  const double alpha = parameter[pp->Cpar.size()];
  const double beta = parameter[pp->Cpar.size()+1];
  const double gamma = parameter[pp->Cpar.size()+2];
  const double scatter0 = parameter[pp->Cpar.size()+3];
  const double scatterM = parameter[pp->Cpar.size()+4];
  const double scatterM_exp = parameter[pp->Cpar.size()+5];
  const double scatterz = parameter[pp->Cpar.size()+6];
  const double scatterz_exp = parameter[pp->Cpar.size()+7];

  // compute the power spectrum
  cosmology::PkXi PX(cosmo);
  std::vector<double> Pk = PX.Pk_matter(pp->kk, pp->method_Pk, false, 0., pp->store_output, pp->output_root, pp->norm, cbl::Min(pp->kk), cbl::Max(pp->kk), pp->prec, cbl::par::defaultString, true);

  // interpolate sigmaM and its derivative
  const std::vector<cbl::glob::FuncGrid> interp = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (cbl::logarithmic_bin_vector(250, 1.e8, 1.e16), cosmo, pp->kk, Pk, "Spline", cbl::Max(pp->kk));

  // interpolate the growth factor
  const std::vector<double> redshift_vector = cbl::linear_bin_vector(20, 0.0001, 2.);
  std::vector<double> DN (redshift_vector.size(), 0.);
  for (size_t i=0; i<redshift_vector.size(); i++)
    DN[i] = cosmo->DN(redshift_vector[i]);
  cbl::glob::FuncGrid DN_interp (redshift_vector, DN, "Spline");

  // Compute the counts
  std::vector<double> number_counts(pp->edges_x.size()-1);

  for (size_t j=0; j<pp->edges_x.size()-1; j++)
    number_counts[j] = cbl::modelling::numbercounts::counts_proxy_classic(cosmo, pp->area_rad, pp->z_min, pp->z_max, pp->edges_x[j], pp->edges_x[j+1], pp->z_error_normalised, pp->Plambda_mean_fc, pp->Plambda_A_mu, pp->Plambda_B_mu, pp->Plambda_C_mu, pp->Plambda_std_fc, pp->Plambda_A_sigma, pp->Plambda_B_sigma, pp->Plambda_C_sigma, interp[0], interp[1], DN_interp, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, pp->purity[j], pp->completeness_interp, pp->scaling_relation, pp->z_pivot, pp->proxy_pivot, pp->mass_pivot, pp->log_base, pp->isDelta_critical, pp->Delta, pp->model_MF, pp->response_fact, pp->model_bias);

  return number_counts;
}
