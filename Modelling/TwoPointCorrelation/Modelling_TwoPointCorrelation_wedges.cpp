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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Data1D.h"
#include "ModelFunction_TwoPointCorrelation_wedges.h"
#include "Modelling_TwoPointCorrelation_wedges.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nwedges(2), m_nwedges_fit(2), m_deltamu(0.5)
{ 
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;

  for (int j=0; j<m_nwedges; j++)
    for (int i=0; i<size; i++)
      m_wedges_order.push_back(j);
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const shared_ptr<data::Data> twop_dataset, const int nwedges)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nwedges(nwedges), m_nwedges_fit(m_nwedges), m_deltamu(1./m_nwedges)
{
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;

  for (int j=0; j<m_nwedges; j++)
    for (int i=0; i<size; i++)
      m_wedges_order.push_back(j);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const double xmin, const double xmax, const int nwedges)
{
  vector<vector<double>> fr(m_nwedges, vector<double>(2, -1.));

  int mp = (nwedges > 0 && nwedges < m_nwedges) ? nwedges : m_nwedges;

  for (int i=0; i<mp; i++) {
    fr[i][0] = xmin;
    fr[i][1] = xmax;
  }

  set_fit_range(fr);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const vector<vector<double>> fit_range)
{
  if ((int)fit_range.size() != m_nwedges)
    ErrorCBL("Error in set_fit_range of :Modelling_TwoPointCorrelation_wedges.cpp, wrong number of wedges provided!");

  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;
  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;
  vector<int> use_w(m_nwedges, 0);

  for (int j=0; j<m_nwedges; j++){
    for (int i=0; i<size; i++){
      if (m_data->xx(i+j*size) < fit_range[j][1] && m_data->xx(i+j*size) > fit_range[j][0]){
	m_wedges_order.push_back(j);
	xx.push_back(m_data->xx(i+j*size));
	use_w[j]=1;
	mask[i+j*size] = true;
      }
    }
  }

  vector<double> data, error;
  vector<vector<double>> covariance;
  m_data->cut(mask, data, error, covariance);

  m_data_fit = make_shared<cbl::data::Data1D>(cbl::data::Data1D(xx, data, covariance));
  m_fit_range = true; 

  m_nwedges_fit=0;
  for (size_t i =0; i<use_w.size(); i++)
    m_nwedges_fit += use_w[i];
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_PkDM ()
{
  m_data_model.nmultipoles = 3;
  m_data_model.nwedges = 2;

  if (m_data_model.sigmaNL==0) {    
    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0);

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk, "Spline"));
  }
  else {
    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0), PkNW(m_data_model.step,0);
    for (size_t i=0; i<kk.size(); i++){
      Pk[i] =  m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] =  m_data_model.cosmology->Pk(kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
    }

    m_data_model.kk = kk;
    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, PkNW, "Spline"));
  }
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_xiDM ()
{

  m_data_model.nmultipoles = 3;
  m_data_model.nwedges = 2;

  cout << endl; coutCBL << "Setting up the fiducial dark matter two-point correlation function model" << endl;

  const vector<double> rad = linear_bin_vector(m_data_model.step, m_data_model.r_min, m_data_model.r_max);

  m_data_model.rr = rad;

  if (m_data_model.sigmaNL==0) {    

    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0);

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk, "Spline"));

  }

  else {

    vector<double> Pk(m_data_model.step,0), PkNW(m_data_model.step,0);
    m_data_model.kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++) {
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
    }

    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model.kk, PkNW, "Spline"));

  }

  vector<vector<double>> xil = Xi_l(rad, m_data_model.nmultipoles, 1., 1., m_data_model.sigmaNL_perp, m_data_model.sigmaNL_par, m_data_model.bias, m_data_model.linear_growth_rate_z, 0., 0., m_data_model.kk, m_data_model.func_Pk, m_data_model.func_Pk_NW, m_data_model.prec);

  m_data_model.func_multipoles.erase(m_data_model.func_multipoles.begin(), m_data_model.func_multipoles.end());

  for(int i=0; i< m_data_model.nmultipoles; i++)
    m_data_model.func_multipoles.push_back(make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(rad, xil[i], "Spline")));

}



// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_fullShape (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, statistics::PriorDistribution fsigma8_prior, statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaS_prior, const bool compute_PkDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model.nwedges = m_nwedges_fit;
  m_data_model.dataset_order = m_wedges_order;
  m_data_model.nmultipoles = 2;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "f*sigma8";
  parameterName[3] = "b*sigma8";
  parameterName[4] = "Sigma_S";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, fsigma8_prior, bsigma8_prior, SigmaS_prior};

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, nparameters, parameterType, parameterName, inputs));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution Bperp_prior, const statistics::PriorDistribution Bpar_prior, const statistics::PriorDistribution Aperp0_prior, const statistics::PriorDistribution Apar0_prior, const statistics::PriorDistribution Aperp1_prior, const statistics::PriorDistribution Apar1_prior, const statistics::PriorDistribution Aperp2_prior, const statistics::PriorDistribution Apar2_prior, const bool compute_XiDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_XiDM) set_fiducial_xiDM();

  m_data_model.nwedges = m_nwedges_fit;
  m_data_model.dataset_order = m_wedges_order;

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

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges_BAO, nparameters, parameterType, parameterName, inputs));
}
