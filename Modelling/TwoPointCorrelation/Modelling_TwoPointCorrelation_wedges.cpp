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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_wedges.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation_wedges
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_wedges, used to model the wedges of
 *  the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "Data1D.h"
#include "ModelFunction_TwoPointCorrelation_wedges.h"
#include "Modelling_TwoPointCorrelation_wedges.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop, const int nWedges, const std::vector<std::vector<double>> mu_integral_limits)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nWedges(nWedges), m_mu_integral_limits(mu_integral_limits)
{
  m_ModelIsSet = false;
  
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  for (int j=0; j<m_nWedges; j++)
    for (int i=0; i<m_data->ndata()/m_nWedges; i++)
      m_wedges_order.push_back(j);

  m_use_wedge.resize(m_nWedges, true);
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<data::Data> twop_dataset, const int nWedges, const std::vector<std::vector<double>> mu_integral_limits)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nWedges(nWedges), m_mu_integral_limits(mu_integral_limits)
{
  m_ModelIsSet = false;
  
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());
  m_use_wedge.resize(m_nWedges, false);
  
  for (int j=0; j<m_nWedges; j++) {
    m_use_wedge[j] = true;
    for (int i=0; i<m_data->ndata()/m_nWedges; i++)
      m_wedges_order.push_back(j);
  }
  
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const double xmin, const double xmax, const int nWedges)
{
  vector<vector<double>> fit_range(m_nWedges, vector<double>(2, -1.));

  int mp = (0<nWedges && nWedges<m_nWedges) ? nWedges : m_nWedges;

  for (int i=0; i<mp; i++) {
    fit_range[i][0] = xmin;
    fit_range[i][1] = xmax;
  }

  set_fit_range(fit_range);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const std::vector<std::vector<double>> fit_range)
{
  if ((int)fit_range.size()!=m_nWedges)
    ErrorCBL("wrong number of wedges provided, "+conv(fit_range.size(), cbl::par::fINT)+" instead of "+conv(m_nWedges, cbl::par::fINT)+"!", "set_fit_range", "Modelling_TwoPointCorrelation_wedges.cpp");

  m_use_wedge = {false, false};
  
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  const int size = m_data->ndata()/m_nWedges;
  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;
  
  for (int j=0; j<m_nWedges; j++) {
    for (int i=0; i<size; i++) {
      if (fit_range[j][0]<m_data->xx(i+j*size) && m_data->xx(i+j*size)<fit_range[j][1]) {
	m_wedges_order.push_back(j);
	xx.push_back(m_data->xx(i+j*size));
	m_use_wedge[j] = true;
	mask[i+j*size] = true;
      }
    }
  }
  m_data_fit = m_data->cut(mask);
  
  m_fit_range = true; 

  if (m_ModelIsSet)
    m_data_model->dataset_order = m_wedges_order;
  
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_PkDM ()
{
  m_data_model->nmultipoles = 3;
  m_data_model->nWedges = 2;

  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));
  vector<double> Pk = m_data_model->cosmology->Pk_matter(m_data_model->kk, m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
  
  if (m_data_model->Pk_mu_model=="dispersion_dewiggled") {
    vector<double> PkNW = m_data_model->cosmology->Pk_matter(m_data_model->kk, "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_NW);
  }

  else if (m_data_model->Pk_mu_model=="dispersion_modecoupling") {
    vector<double> kk_1loop, Pk_1loop;
    for (size_t i=0; i<(size_t)m_data_model->step; i++) {
      if (m_data_model->kk[i]<par::pi) {
	kk_1loop.push_back(m_data_model->kk[i]);
	Pk_1loop.push_back(m_data_model->cosmology->Pk_1loop(m_data_model->kk[i], m_data_model->func_Pk, 0,  m_data_model->k_min, 5., m_data_model->prec));
      }
    }
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk1loop = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk_1loop, Pk_1loop, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk1loop);
  }

  else if (m_data_model->Pk_mu_model=="dispersion_Gauss" || m_data_model->Pk_mu_model=="dispersion_Lorentz") {
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
  }

  else if (m_data_model->Pk_mu_model=="Scoccimarro_Pezzotta_Gauss" || m_data_model->Pk_mu_model=="Scoccimarro_Pezzotta_Lorentz" || m_data_model->Pk_mu_model=="Scoccimarro_Bel_Gauss" || m_data_model->Pk_mu_model=="Scoccimarro_Bel_Lorentz") {
    vector<double> Pknonlin = m_data_model->cosmology->Pk_matter(m_data_model->kk, m_data_model->method_Pk, true, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_nonlin = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pknonlin, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_nonlin);
  }

  else if (m_data_model->Pk_mu_model=="Scoccimarro_Gauss" || m_data_model->Pk_mu_model=="Scoccimarro_Lorentz") {
    vector<vector<double>> Pk_terms = m_data_model->cosmology->Pk_TNS_dd_dt_tt(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[2], "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
  }

  else if (m_data_model->Pk_mu_model=="TNS_Gauss" || m_data_model->Pk_mu_model=="TNS_Lorentz") {
    vector<vector<double>> Pk_terms = m_data_model->cosmology->Pk_TNS_dd_dt_tt(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
    vector<vector<double>> Pk_AB = m_data_model->cosmology->Pk_TNS_AB_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
    
    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[2], "Spline"));

    m_data_model->func_Pk_A11 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[0], "Spline"));
    m_data_model->func_Pk_A12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[1], "Spline"));
    m_data_model->func_Pk_A22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[2], "Spline"));
    m_data_model->func_Pk_A23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[3], "Spline"));
    m_data_model->func_Pk_A33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[4], "Spline"));
    m_data_model->func_Pk_B12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[5], "Spline"));
    m_data_model->func_Pk_B13 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[6], "Spline"));
    m_data_model->func_Pk_B14 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[7], "Spline"));
    m_data_model->func_Pk_B22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[8], "Spline"));
    m_data_model->func_Pk_B23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[9], "Spline"));
    m_data_model->func_Pk_B24 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[10], "Spline"));
    m_data_model->func_Pk_B33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[11], "Spline"));
    m_data_model->func_Pk_B34 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[12], "Spline"));
    m_data_model->func_Pk_B44 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[13], "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A11);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B13);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B14);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B24);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B34);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B44);
  }

  else if (m_data_model->Pk_mu_model=="eTNS_Gauss" || m_data_model->Pk_mu_model=="eTNS_Lorentz") {
    vector<vector<double>> Pk_AB = m_data_model->cosmology->Pk_TNS_AB_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
    vector<vector<double>> Pk_eTNS_terms = m_data_model->cosmology->Pk_eTNS_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[2], "Spline"));

    m_data_model->func_Pk_A11 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[0], "Spline"));
    m_data_model->func_Pk_A12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[1], "Spline"));
    m_data_model->func_Pk_A22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[2], "Spline"));
    m_data_model->func_Pk_A23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[3], "Spline"));
    m_data_model->func_Pk_A33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[4], "Spline"));
    m_data_model->func_Pk_B12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[5], "Spline"));
    m_data_model->func_Pk_B13 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[6], "Spline"));
    m_data_model->func_Pk_B14 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[7], "Spline"));
    m_data_model->func_Pk_B22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[8], "Spline"));
    m_data_model->func_Pk_B23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[9], "Spline"));
    m_data_model->func_Pk_B24 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[10], "Spline"));
    m_data_model->func_Pk_B33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[11], "Spline"));
    m_data_model->func_Pk_B34 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[12], "Spline"));
    m_data_model->func_Pk_B44 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[13], "Spline"));

    m_data_model->func_Pk_b2d = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[3], "Spline"));
    m_data_model->func_Pk_b2v = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[4], "Spline"));
    m_data_model->func_Pk_b22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[5], "Spline"));
    m_data_model->func_Pk_bs2d = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[6], "Spline"));
    m_data_model->func_Pk_bs2v = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[7], "Spline"));
    m_data_model->func_Pk_b2s2 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[8], "Spline"));
    m_data_model->func_Pk_bs22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[9], "Spline"));
    m_data_model->func_sigma32Pklin = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[10], "Spline"));

    //-------------------------

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A11);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B13);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B14);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B24);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B34);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B44);

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b2d);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b2v);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_bs2d);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_bs2v);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b2s2);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_bs22);
    m_data_model->funcs_pk.push_back(m_data_model->func_sigma32Pklin);
  }
  
  else ErrorCBL("the chosen model ("+m_data_model->Pk_mu_model+") is not currently implemented!", "set_fiducial_PkDM", "Modelling_TwoPointCorrelation_wedges.cpp");
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial dark matter two-point correlation function model" << endl;

  m_data_model->nmultipoles = 3;
  m_data_model->nWedges = 2;
  
  const vector<double> rad = linear_bin_vector(m_data_model->step, m_data_model->r_min, m_data_model->r_max);

  m_data_model->rr = rad;
  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));

  vector<double> Pk =  m_data_model->cosmology->Pk_matter(m_data_model->kk, m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
  vector<double> PkNW =  m_data_model->cosmology->Pk_matter(m_data_model->kk, "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);


  m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
  m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

  std::vector<double> template_parameter = {m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->linear_growth_rate_z, m_data_model->bias, 0.};
  cbl::Print(template_parameter);

  vector<vector<double>> xil = Xi_l(rad, m_data_model->nmultipoles, "dispersion_dewiggled", {m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->linear_growth_rate_z, m_data_model->bias, 0.}, {m_data_model->func_Pk, m_data_model->func_Pk_NW}, m_data_model->prec, 1, 1);

  m_data_model->func_multipoles.erase(m_data_model->func_multipoles.begin(), m_data_model->func_multipoles.end());

  for (int i=0; i< m_data_model->nmultipoles; i++)
    m_data_model->func_multipoles.push_back(make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(rad, xil[i], "Spline")));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameter)
{
  if (m_likelihood==NULL) ErrorCBL("this function requires the likelihood to be defined (with the function set_likelihood)!", "write_model", "Modelling_TwoPointCorrelation_wedges.cpp");
  
  vector<int> dataset_order_original = m_data_model->dataset_order;

  vector<int> new_dataset_order;
  vector<double> new_xx;
  
  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nWedges; n++)
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }

  m_data_model->dataset_order = new_dataset_order;
  m_likelihood->write_model(output_dir, output_file, parameter, new_xx);
  m_data_model->dataset_order = dataset_order_original;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run maximize_posterior first!", "write_model_at_bestfit", "Modelling_TwoPointCorrelation_wedges.cpp");

  vector<int> dataset_order_original = m_data_model->dataset_order;

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nWedges; n++)
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }

  m_data_model->dataset_order = new_dataset_order;
  m_posterior->write_model_at_bestfit(output_dir, output_file, new_xx);

  m_data_model->dataset_order = dataset_order_original;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start, const int thin)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run sample_posterior first!", "write_model_from_chains", "Modelling_TwoPointCorrelation_wedges.cpp");

  vector<int> dataset_order_original = m_data_model->dataset_order;

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nWedges; n++)
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }

  m_data_model->dataset_order = new_dataset_order;
  m_posterior->write_model_from_chain(output_dir, output_file, new_xx, {}, start, thin);

  m_data_model->dataset_order = dataset_order_original;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_fullShape_DeWiggled (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution SigmaNL_perpendicular_prior, const statistics::PriorDistribution SigmaNL_parallel_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaS_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_dewiggled";

  // compute the fiducial dark matter power spectrum terms used to construct the model
  if (compute_PkDM) set_fiducial_PkDM();
  
  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;
  
  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;
  
  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 7;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "SigmaNL_perpendicular";
  parameterName[3] = "SigmaNL_parallel";
  parameterName[4] = "f*sigma8";
  parameterName[5] = "b*sigma8";
  parameterName[6] = "Sigma_S";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, SigmaNL_perpendicular_prior, SigmaNL_parallel_prior, fsigma8_prior, bsigma8_prior, SigmaS_prior};

  //set the priors
  m_set_prior(priors);
  
  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}

// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_fullShape_ModeCoupling (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaV_prior, const statistics::PriorDistribution AMC_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_modecoupling";

  // compute the fiducial dark matter power spectrum terms used to construct the model
  if (compute_PkDM) set_fiducial_PkDM();

  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;
  
  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 6;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "f*sigma8";
  parameterName[3] = "b*sigma8";
  parameterName[4] = "sigma_v";
  parameterName[5] = "AMC";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, fsigma8_prior, bsigma8_prior, SigmaV_prior, AMC_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution Bperp_prior, const statistics::PriorDistribution Bpar_prior, const statistics::PriorDistribution Aperp0_prior, const statistics::PriorDistribution Apar0_prior, const statistics::PriorDistribution Aperp1_prior, const statistics::PriorDistribution Apar1_prior, const statistics::PriorDistribution Aperp2_prior, const statistics::PriorDistribution Apar2_prior, const bool compute_XiTemplate, const bool isRealSpace)
{
  // compute the fiducial dark matter two-point correlation function

  if (m_data_model->nWedges>2)
    ErrorCBL("BAO modelling can be done only with two wedges!", "set_model_BAO", "Modelling_TwoPointCorrelation_wedges");

  if (isRealSpace) {
    double lgf = m_data_model->linear_growth_rate_z;
    m_data_model->linear_growth_rate_z = 0.;

    set_fiducial_xiDM();

    m_data_model->linear_growth_rate_z = lgf;
  }
  else if (compute_XiTemplate)
    set_fiducial_xiDM();

  m_data_model->dataset_order = m_wedges_order;

  // set the model parameters
  const int nparameters = 10;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "Bperp";
  parameterName[3] = "Bpar";
  parameterName[4] = "Aperp0";
  parameterName[5] = "Apar0";
  parameterName[6] = "Aperp1";
  parameterName[7] = "Apar1";
  parameterName[8] = "Aperp2";
  parameterName[9] = "Apar2";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, Bperp_prior, Bpar_prior, Aperp0_prior, Apar0_prior, Aperp1_prior, Apar1_prior, Aperp2_prior, Apar2_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges_BAO, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_dispersion (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "dispersion_Gauss";
  else m_data_model->Pk_mu_model = "dispersion_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;

  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;

  m_data_model->nmultipoles = 3;
  
  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "alpha_perpendicular";
  parameterName[4] = "alpha_parallel";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}

// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_Scoccimarro (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "Scoccimarro_Gauss";
  else m_data_model->Pk_mu_model = "Scoccimarro_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();
  
  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;
  
  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;

  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "alpha_perpendicular";
  parameterName[4] = "alpha_parallel";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_Scoccimarro_fitPezzotta (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution kd_prior, const statistics::PriorDistribution kt_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "Scoccimarro_Pezzotta_Gauss";
  else m_data_model->Pk_mu_model = "Scoccimarro_Pezzotta_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();
  
  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;

  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;

  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 7;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "kd";
  parameterName[4] = "kt";
  parameterName[5] = "alpha_perpendicular";
  parameterName[6] = "alpha_parallel";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, kd_prior, kt_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_Scoccimarro_fitBel (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution kd_prior, const statistics::PriorDistribution bb_prior, const statistics::PriorDistribution a1_prior, const statistics::PriorDistribution a2_prior, const statistics::PriorDistribution a3_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "Scoccimarro_Bel_Gauss";
  else m_data_model->Pk_mu_model = "Scoccimarro_Bel_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;
  
  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;
  
  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 10;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "kd";
  parameterName[4] = "bb";
  parameterName[5] = "a1";
  parameterName[6] = "a2";
  parameterName[7] = "a3";
  parameterName[8] = "alpha_perpendicular";
  parameterName[9] = "alpha_parallel";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, kd_prior, bb_prior, a1_prior, a2_prior, a3_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_TNS (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "TNS_Gauss";
  else m_data_model->Pk_mu_model = "TNS_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;
  
  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;
  
  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "alpha_perpendicular";
  parameterName[4] = "alpha_parallel";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_eTNS (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution b1sigma8_prior, const statistics::PriorDistribution b2sigma8_prior, const statistics::PriorDistribution sigmav_prior, statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution N_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "eTNS_Gauss";
  else m_data_model->Pk_mu_model = "eTNS_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  // number of wedges to be used
  m_data_model->nWedges = m_nWedges;

  // integral limits used to measure the wedges
  m_data_model->mu_integral_limits = m_mu_integral_limits;
  
  // scales to be used for each wedges
  m_data_model->dataset_order = m_wedges_order;
  
  m_data_model->nmultipoles = 3;

  // set the model parameters
  const int nparameters = 7;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b1*sigma8";
  parameterName[2] = "b2*sigma8";
  parameterName[3] = "sigmav";
  parameterName[4] = "alpha_perpendicular";
  parameterName[5] = "alpha_parallel";
  parameterName[6] = "N";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, b1sigma8_prior, b2sigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior, N_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


