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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_multipoles.cpp
 *
 *  @brief Methods of the class
 *  Modelling_TwoPointCorrelation_multipoles
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_multipoles, used to model the
 *  multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Data1D.h"
#include "Modelling_TwoPointCorrelation_multipoles.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nmultipoles(3), m_nmultipoles_fit(3)
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;

  for (int j=0; j<m_nmultipoles; j++)
    for (int i=0; i<size; i++)
      m_multipoles_order.push_back(j);

}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const shared_ptr<data::Data> twop_dataset, const int nmultipoles)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nmultipoles(nmultipoles), m_nmultipoles_fit(nmultipoles)
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;

  for (int j=0; j<m_nmultipoles; j++)
    for (int i=0; i<size; i++)
      m_multipoles_order.push_back(j);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const double xmin, const double xmax, const int nmultipoles)
{
  vector<vector<double>> fr(m_nmultipoles, vector<double>(2, -1.));

  int mp = (nmultipoles > 0 && nmultipoles < m_nmultipoles) ? nmultipoles : m_nmultipoles;

  for (int i=0; i<mp; i++) {
    fr[i][0] = xmin;
    fr[i][1] = xmax;
  }

  set_fit_range(fr);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const vector<vector<double>> fit_range)
{
  if ((int)fit_range.size() != m_nmultipoles)
    ErrorCBL("Error in set_fit_range of :Modelling_TwoPointCorrelation_multipoles.cpp, wrong number of multipoles provided!");

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;
  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;
  vector<int> use_pole(m_nmultipoles, 0);

  for (int j=0; j<m_nmultipoles; j++){
    for (int i=0; i<size; i++){
      if (m_data->xx(i+j*size) < fit_range[j][1] && m_data->xx(i+j*size) > fit_range[j][0]){
	m_multipoles_order.push_back(j);
	xx.push_back(m_data->xx(i+j*size));
	use_pole[j]=1;
	mask[i+j*size] = true;
      }
    }
  }

  vector<double> data, error;
  vector<vector<double>> covariance;
  m_data->cut(mask, data, error, covariance);

  m_data_fit = make_shared<cbl::data::Data1D>(cbl::data::Data1D(xx, data, covariance));
  m_fit_range = true; 

  m_nmultipoles_fit = 0;
  for (size_t i =0; i<use_pole.size(); i++)
    m_nmultipoles_fit += use_pole[i];

  if (m_ModelIsSet) {
    m_data_model.dataset_order = m_multipoles_order;
    auto inputs = make_shared<STR_data_model>(m_data_model);
    m_model->set_inputs(inputs);
  }
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_PkDM ()
{
  m_data_model.nmultipoles = m_nmultipoles;

  if (m_data_model.sigmaNL==0) {    

    vector<double> Pk(m_data_model.step,0);
    m_data_model.kk= logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));
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
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial two-point correlation function model" << endl;

  m_data_model.nmultipoles = 3;

  const vector<double> rad = linear_bin_vector(m_data_model.step, m_data_model.r_min, m_data_model.r_max);

  if (m_data_model.sigmaNL==0) {    

    vector<double> Pk(m_data_model.step,0);
    m_data_model.kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));

  }

  else {

    vector<double> Pk(m_data_model.step, 0), PkNW(m_data_model.step, 0);
    m_data_model.kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++){
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
    }

    m_data_model.func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model.kk, PkNW, "Spline"));
  }


  vector<vector<double>> xil = Xi_l(rad, m_data_model.nmultipoles, 1., 1., m_data_model.sigmaNL_perp, m_data_model.sigmaNL_par, m_data_model.bias, m_data_model.linear_growth_rate_z, 0., 0., m_data_model.kk, m_data_model.func_Pk, m_data_model.func_Pk_NW, m_data_model.prec);

  m_data_model.func_multipoles.erase(m_data_model.func_multipoles.begin(), m_data_model.func_multipoles.end());
  for (int i=0; i< m_data_model.nmultipoles; i++)
    m_data_model.func_multipoles.push_back(make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(rad, xil[i], "Spline")));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, statistics::PriorDistribution fsigma8_prior, statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaS_prior, const bool compute_PkDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model.nmultipoles = m_nmultipoles_fit;
  m_data_model.dataset_order = m_multipoles_order;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType = {statistics::ParameterType::_Base_};

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
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, inputs));
  m_ModelIsSet = true;

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_sigma8_bias (const statistics::PriorDistribution sigma8_prior, const statistics::PriorDistribution bias_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_PkDM();

  m_data_model.nmultipoles = m_nmultipoles_fit;
  m_data_model.dataset_order = m_multipoles_order;

  // set the model parameters
  const int nparameters = 2;

  vector<statistics::ParameterType> parameterType = {statistics::ParameterType::_Base_};

  vector<string> parameterName(nparameters);
  parameterName[0] = "sigma8";
  parameterName[1] = "bias";

  vector<statistics::PriorDistribution> priors = {sigma8_prior, bias_prior};
  
  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_sigma8_bias, nparameters, parameterType, parameterName, inputs));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution B0_prior, const statistics::PriorDistribution B2_prior, const statistics::PriorDistribution A00_prior, const statistics::PriorDistribution A20_prior, const statistics::PriorDistribution A01_prior, const statistics::PriorDistribution A21_prior, const statistics::PriorDistribution A02_prior, const statistics::PriorDistribution A22_prior, const bool compute_XiDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_XiDM) set_fiducial_xiDM();

  m_data_model.nmultipoles = m_nmultipoles_fit;
  m_data_model.dataset_order = m_multipoles_order;

  // set the model parameters
  const int nparameters = 10;

  vector<statistics::ParameterType> parameterType = {statistics::ParameterType::_Base_};

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "B0";
  parameterName[3] = "B2";
  parameterName[4] = "A00";
  parameterName[5] = "A20";
  parameterName[6] = "A01";
  parameterName[7] = "A21";
  parameterName[8] = "A02";
  parameterName[9] = "A22";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, B0_prior, B2_prior, A00_prior, A20_prior, A01_prior, A21_prior, A02_prior, A22_prior};

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_BAO, nparameters, parameterType, parameterName, inputs));
  m_ModelIsSet = true;
}
