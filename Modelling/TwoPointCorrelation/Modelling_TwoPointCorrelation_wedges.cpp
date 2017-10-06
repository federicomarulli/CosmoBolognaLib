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

using namespace cosmobl;


// ============================================================================================


cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nwedges(2), m_nwedges_fit(2), m_deltamu(0.5)
{ 
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;

  for (int j=0; j<m_nwedges; j++)
    for (int i=0; i<size; i++)
      m_wedges_order.push_back(j);
}


// ============================================================================================


cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const shared_ptr<data::Data> twop_dataset, const int nwedges)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nwedges(nwedges), m_nwedges_fit(m_nwedges), m_deltamu(1./m_nwedges)
{
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;

  for (int j=0; j<m_nwedges; j++)
    for (int i=0; i<size; i++)
      m_wedges_order.push_back(j);
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const double xmin, const double xmax, const int nwedges)
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


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const vector<vector<double>> fit_range)
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

  m_data_fit = make_shared<cosmobl::data::Data1D>(cosmobl::data::Data1D(xx, data, covariance));
  m_fit_range = true; 

  m_nwedges_fit=0;
  for (size_t i =0; i<use_w.size(); i++)
    m_nwedges_fit += use_w[i];
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_PkDM ()
{
  m_data_model.nmultipoles = 3;
  m_data_model.nwedges = 2;

  if (m_data_model.sigmaNL==0) {    
    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0);

    for (size_t i=0; i<(size_t)m_data_model.step; i++) 
      Pk[i] =  m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(kk, Pk, "Spline"));
  }
  else {
    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0), PkNW(m_data_model.step,0);
    for (size_t i=0; i<kk.size(); i++){
      Pk[i] =  m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] =  m_data_model.cosmology->Pk(kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
    }

    m_data_model.kk = kk;
    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(kk, PkNW, "Spline"));
  }
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_xiDM ()
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

    m_data_model.func_Pk = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(kk, Pk, "Spline"));

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

  vector<vector<double>> xil = Xi_l(rad, m_data_model.nmultipoles, 1., 1., m_data_model.sigmaNL_perp, m_data_model.sigmaNL_par, m_data_model.bias, m_data_model.linear_growth_rate_z, 0., m_data_model.kk, m_data_model.func_Pk, m_data_model.func_Pk_NW, m_data_model.prec);

  m_data_model.func_multipoles.erase(m_data_model.func_multipoles.begin(), m_data_model.func_multipoles.end());

  for(int i=0; i< m_data_model.nmultipoles; i++)
    m_data_model.func_multipoles.push_back(make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(rad, xil[i], "Spline")));

}



// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_fullShape (const statistics::Prior alpha_perpendicular_prior, const statistics::Prior alpha_parallel_prior, statistics::Prior fsigma8_prior, statistics::Prior bsigma8_prior, const statistics::Prior SigmaS_prior, const bool compute_PkDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model.nwedges = m_nwedges_fit;
  m_data_model.dataset_order = m_wedges_order;
  m_data_model.nmultipoles = 2;

  m_alpha_perpendicular = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_perpendicular_prior, "alpha_perpendicular"));
  m_alpha_parallel = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_parallel_prior, "alpha_parallel"));
  m_fsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(fsigma8_prior, "f*sigma8"));
  m_bsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bsigma8_prior, "b*sigma8"));
  m_SigmaS = make_shared<statistics::BaseParameter>(statistics::BaseParameter(SigmaS_prior, "b*sigma8"));
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_alpha_perpendicular, m_alpha_parallel, m_fsigma8, m_bsigma8, m_SigmaS};

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_BAO (const statistics::Prior alpha_perpendicular_prior, const statistics::Prior alpha_parallel_prior, const statistics::Prior Bperp_prior, const statistics::Prior Bpar_prior, const statistics::Prior Aperp0_prior, const statistics::Prior Apar0_prior, const statistics::Prior Aperp1_prior, const statistics::Prior Apar1_prior, const statistics::Prior Aperp2_prior, const statistics::Prior Apar2_prior, const bool compute_XiDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_XiDM) set_fiducial_xiDM();

  m_data_model.nwedges = m_nwedges_fit;
  m_data_model.dataset_order = m_wedges_order;

  m_alpha_perpendicular = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_perpendicular_prior, "alpha_perpendicular"));
  m_alpha_parallel = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_parallel_prior, "alpha_parallel"));
  m_Bperp = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Bperp_prior, "Bperp"));
  m_Bpar = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Bpar_prior, "Bpar"));
  m_Aperp0 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Aperp0_prior, "Aperp0"));
  m_Apar0 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Apar0_prior, "Apar0"));
  m_Aperp1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Aperp1_prior, "Aperp1"));
  m_Apar1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Apar1_prior, "Apar1"));
  m_Apar2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Apar2_prior, "Apar2"));
  m_Aperp2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Aperp2_prior, "Aperp2"));
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_alpha_perpendicular, m_alpha_parallel, m_Bperp, m_Bpar, m_Aperp0, m_Apar0, m_Aperp1, m_Apar1, m_Aperp2, m_Apar2};

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiWedges_BAO, inputs));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model (const string dir, const string file, const vector<double> xx, const vector<double> parameter, const int start, const int thin)
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
      for(int j=0; j<m_nwedges_fit; j++)
	fout << median_model[i+j*rad.size()] << "  " << low_model[i+j*rad.size()] << "  " << up_model[i+j*rad.size()] << " ";
      fout << endl;
    }
    
    fout.clear(); fout.close();
    coutCBL << "I wrote the file: " << dir+file << endl;
  }

  
  // compute the model with input parameters
  
  else{

    vector<double> new_rad;

    for(int i=0; i<m_nwedges_fit; i++)
      for(size_t j=0; j<rad.size(); j++)
	new_rad.push_back(rad[j]);

    vector<int> wedges_order;

    for (int j=0; j<m_nwedges_fit; j++)
      for (size_t i=0; i<rad.size(); i++)
	wedges_order.push_back(j);

    m_data_model.dataset_order = wedges_order;

    // input data used to construct the model
    auto inputs = make_shared<STR_data_model>(m_data_model);

    m_model->set_fixed_parameters(inputs);

    vector<double> par = parameter;
    vector<double> model = m_model->operator()(new_rad, par);

    string file_out = dir+file;
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  
    for (size_t i=0; i<rad.size(); i++){
      fout << rad[i] << "  ";
      for(int j=0; j<m_nwedges_fit; j++)
	fout << model[i+j*rad.size()] << " ";
      fout << endl;
    }

    fout.clear(); fout.close();
    coutCBL << "I wrote the file: " << dir+file << endl;

    m_data_model.dataset_order = m_wedges_order;
    m_model->set_fixed_parameters(stored_inputs);
  }

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::compute_model_from_chains (const vector<double> xx, const int start, const int thin, vector<double> &median_model, vector<double> &low_model, vector<double> &up_model)
{
  auto stored_inputs = m_model->fixed_parameters();

  vector<double> rad = (xx.size()==0) ? logarithmic_bin_vector(100, 0.1, 100.) : xx;

  vector<double> new_rad;

  for(int i=0; i<m_nwedges_fit; i++)
    for(size_t j=0; j<rad.size(); j++)
      new_rad.push_back(rad[j]);
  
  vector<int> wedges_order;

  for (int j=0; j<m_nwedges_fit; j++)
    for (size_t i=0; i<rad.size(); i++)
      wedges_order.push_back(j);

  m_data_model.dataset_order = wedges_order;
  
  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  m_model->set_fixed_parameters(inputs);

  // compute the model with best-fit parameters

  Modelling_TwoPointCorrelation1D_monopole::compute_model_from_chains(new_rad, start, thin, median_model, low_model, up_model); // check!!!

  m_data_model.dataset_order = m_wedges_order;
  m_model->set_fixed_parameters(stored_inputs);
}
