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
 *  @file Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_multipoles.cpp
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

using namespace cosmobl;

// ============================================================================================


cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nmultipoles(3), m_nmultipoles_fit(3)
{
  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;

  for (int j=0; j<m_nmultipoles; j++)
    for (int i=0; i<size; i++)
      m_multipoles_order.push_back(j);

}


// ============================================================================================


cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const shared_ptr<data::Data> twop_dataset, const int nmultipoles)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nmultipoles(nmultipoles), m_nmultipoles_fit(nmultipoles)
{
  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;

  for (int j=0; j<m_nmultipoles; j++)
    for (int i=0; i<size; i++)
      m_multipoles_order.push_back(j);
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const double xmin, const double xmax, const int nmultipoles)
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


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const vector<vector<double>> fit_range)
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

  m_data_fit = make_shared<cosmobl::data::Data1D>(cosmobl::data::Data1D(xx, data, covariance));
  m_fit_range = true; 

  m_nmultipoles_fit=0;
  for (size_t i =0; i<use_pole.size(); i++)
    m_nmultipoles_fit += use_pole[i];
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_PkDM ()
{
  m_data_model.nmultipoles = m_nmultipoles;

  if (m_data_model.sigmaNL==0) {    

    vector<double> Pk(m_data_model.step,0);
    m_data_model.kk= logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));
  }
  
  else {
    vector<double> Pk(m_data_model.step,0), PkNW(m_data_model.step,0);
    m_data_model.kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++) {
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
    }

    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(m_data_model.kk, PkNW, "Spline"));
  }
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial two-point correlation function model" << endl;

  m_data_model.nmultipoles = 3;

  const vector<double> rad = linear_bin_vector(m_data_model.step, m_data_model.r_min, m_data_model.r_max);

  if (m_data_model.sigmaNL==0) {    

    vector<double> Pk(m_data_model.step,0);
    m_data_model.kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));

  }

  else {

    vector<double> Pk(m_data_model.step, 0), PkNW(m_data_model.step, 0);
    m_data_model.kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model.step; i++){
      Pk[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] =  m_data_model.cosmology->Pk(m_data_model.kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
    }

    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(m_data_model.kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(m_data_model.kk, PkNW, "Spline"));
  }


  vector<vector<double>> xil = Xi_l(rad, m_data_model.nmultipoles, 1., 1., m_data_model.sigmaNL_perp, m_data_model.sigmaNL_par, m_data_model.bias, m_data_model.linear_growth_rate_z, 0., m_data_model.kk, m_data_model.func_Pk, m_data_model.func_Pk_NW, m_data_model.prec);

  m_data_model.func_multipoles.erase(m_data_model.func_multipoles.begin(), m_data_model.func_multipoles.end());
  for (int i=0; i< m_data_model.nmultipoles; i++)
    m_data_model.func_multipoles.push_back(make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(rad, xil[i], "Spline")));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape (const statistics::Prior alpha_perpendicular_prior, const statistics::Prior alpha_parallel_prior, statistics::Prior fsigma8_prior, statistics::Prior bsigma8_prior, const statistics::Prior SigmaS_prior, const bool compute_PkDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model.nmultipoles = m_nmultipoles_fit;
  m_data_model.dataset_order = m_multipoles_order;

  m_alpha_perpendicular = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_perpendicular_prior, "alpha_perpendicular"));
  m_alpha_parallel = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_parallel_prior, "alpha_parallel"));
  m_fsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(fsigma8_prior, "f*sigma8"));
  m_bsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bsigma8_prior, "b*sigma8"));
  m_SigmaS = make_shared<statistics::BaseParameter>(statistics::BaseParameter(SigmaS_prior, "Sigma_S"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_alpha_perpendicular, m_alpha_parallel, m_fsigma8, m_bsigma8, m_SigmaS};

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_BAO (const statistics::Prior alpha_perpendicular_prior, const statistics::Prior alpha_parallel_prior, const statistics::Prior B0_prior, const statistics::Prior B2_prior, const statistics::Prior A00_prior, const statistics::Prior A20_prior, const statistics::Prior A01_prior, const statistics::Prior A21_prior, const statistics::Prior A02_prior, const statistics::Prior A22_prior, const bool compute_XiDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_XiDM) set_fiducial_xiDM();

  m_data_model.nmultipoles = m_nmultipoles_fit;
  m_data_model.dataset_order = m_multipoles_order;

  m_alpha_perpendicular = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_perpendicular_prior, "alpha_perpendicular"));
  m_alpha_parallel = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_parallel_prior, "alpha_parallel"));
  m_B0 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(B0_prior, "B0"));
  m_B2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(B2_prior, "B2"));

  m_A00 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A00_prior, "A00"));
  m_A20 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A20_prior, "A20"));
  
  m_A01 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A01_prior, "A01"));
  m_A21 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A21_prior, "A21"));
  
  m_A02 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A02_prior, "A02"));
  m_A22 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A22_prior, "A22"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_alpha_perpendicular, m_alpha_parallel, m_B0, m_B2, m_A00, m_A20, m_A01, m_A21, m_A02, m_A22};

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_BAO, inputs));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_sigma8_bias (const statistics::Prior sigma8_prior, const statistics::Prior bias_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_PkDM();

  m_data_model.nmultipoles = m_nmultipoles_fit;
  m_data_model.dataset_order = m_multipoles_order;

  // model parameters
  m_sigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigma8_prior, "sigma8"));
  m_bias = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias_prior, "bias"));

  // set model parameters
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_sigma8, m_bias};
  set_parameters(ll_parameters);
  
  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_sigma8_bias, inputs));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model (const string dir, const string file, const vector<double> xx, const vector<double> parameter, const int start, const int thin)
{

  auto stored_inputs = m_model->fixed_parameters();

  vector<double> rad = (xx.size()==0) ? logarithmic_bin_vector(100, 0.1, 100.) : xx;

  // compute the model with best-fit parameters
  
  if (parameter.size()==0) {
    
    string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}
    string file_out = dir+file;
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
    
    vector<double> median_model, low_model, up_model;
    compute_model_from_chains(rad, start, thin, median_model, low_model, up_model); // check!!!
  
    for (size_t i=0; i<rad.size(); i++){
      fout << rad[i] << "  ";
      for (int j=0; j<m_nmultipoles_fit; j++)
	fout << median_model[i+j*rad.size()] << "  " << low_model[i+j*rad.size()] << "  " << up_model[i+j*rad.size()] << " ";
      fout << endl;
    }
    
    fout.clear(); fout.close();
    coutCBL << "I wrote the file: " << dir+file << endl;
  }

  // compute the model with input parameters
  
  else {

    vector<double> new_rad;

    for (int i=0; i<m_nmultipoles_fit; i++)
      for (size_t j=0; j<rad.size(); j++)
	new_rad.push_back(rad[j]);

    vector<int> multipoles_order;

    for (int j=0; j<m_nmultipoles_fit; j++)
      for (size_t i=0; i<rad.size(); i++)
	multipoles_order.push_back(j);

    m_data_model.dataset_order = multipoles_order;


    // input data used to construct the model
    auto inputs = make_shared<STR_data_model>(m_data_model);

    m_model->set_fixed_parameters(inputs);

    vector<double> par = parameter;
    vector<double> model = m_model->operator()(new_rad, par);

    string file_out = dir+file;
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  
    for (size_t i=0; i<rad.size(); i++){
      fout << rad[i] << "  ";
      for (int j=0; j<m_nmultipoles_fit; j++)
	fout << model[i+j*rad.size()] << " ";
      fout << endl;
    }

    fout.clear(); fout.close();
    coutCBL << "I wrote the file: " << dir+file << endl;

    m_data_model.dataset_order = m_multipoles_order;
    m_model->set_fixed_parameters(stored_inputs);
  }

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::compute_model_from_chains (const vector<double> xx, const int start, const int thin, vector<double> &median_model, vector<double> &low_model, vector<double> &up_model)
{
  auto stored_inputs = m_model->fixed_parameters();

  vector<double> rad = (xx.size()==0) ? logarithmic_bin_vector(100, 0.1, 100.) : xx;

  vector<double> new_rad;

  for (int i=0; i<m_nmultipoles_fit; i++)
    for (size_t j=0; j<rad.size(); j++)
      new_rad.push_back(rad[j]);
  
  vector<int> multipoles_order;

  for (int j=0; j<m_nmultipoles_fit; j++)
    for (size_t i=0; i<rad.size(); i++)
      multipoles_order.push_back(j);

  m_data_model.dataset_order = multipoles_order;
  
  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  m_model->set_fixed_parameters(inputs);

  // compute the model with best-fit parameters

  Modelling_TwoPointCorrelation1D_monopole::compute_model_from_chains(new_rad, start, thin, median_model, low_model, up_model); // check!!!

  m_data_model.dataset_order = m_multipoles_order;
  m_model->set_fixed_parameters(stored_inputs);
}

