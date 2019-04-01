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


cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nwedges(2), m_nwedges_fit(2), m_deltamu(0.5)
{ 
  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;

  for (int j=0; j<m_nwedges; j++)
    for (int i=0; i<size; i++)
      m_wedges_order.push_back(j);
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::Modelling_TwoPointCorrelation_wedges (const std::shared_ptr<data::Data> twop_dataset, const int nwedges)
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

  int mp = (nwedges>0 && nwedges<m_nwedges) ? nwedges : m_nwedges;

  for (int i=0; i<mp; i++) {
    fr[i][0] = xmin;
    fr[i][1] = xmax;
  }

  set_fit_range(fr);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range (const std::vector<std::vector<double>> fit_range)
{
  if ((int)fit_range.size()!=m_nwedges)
    ErrorCBL("Error in cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fit_range() of Modelling_TwoPointCorrelation_wedges.cpp: wrong number of wedges provided, "+conv(fit_range.size(), cbl::par::fINT)+" instead of "+conv(m_nwedges, cbl::par::fINT)+"!");

  m_wedges_order.erase(m_wedges_order.begin(), m_wedges_order.end());

  int size = m_data->ndata()/m_nwedges;
  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;
  vector<int> use_w(m_nwedges, 0);

  for (int j=0; j<m_nwedges; j++) {
    for (int i=0; i<size; i++) {
      if (fit_range[j][0]<m_data->xx(i+j*size) && m_data->xx(i+j*size)<fit_range[j][1]) {
	m_wedges_order.push_back(j);
	xx.push_back(m_data->xx(i+j*size));
	use_w[j] = 1;
	mask[i+j*size] = true;
      }
    }
  }

  vector<double> data, error;
  vector<vector<double>> covariance;
  m_data->cut(mask, data, error, covariance);

  m_data_fit = make_shared<cbl::data::Data1D>(cbl::data::Data1D(xx, data, covariance));
  m_fit_range = true; 

  m_nwedges_fit = 0;
  for (size_t i=0; i<use_w.size(); i++)
    m_nwedges_fit += use_w[i];

  m_nwedges_fit = 2;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_PkDM ()
{
  m_data_model->nmultipoles = 3;
  m_data_model->nwedges = 2;

  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));
  vector<double> Pk(m_data_model->step, 0);

  for (size_t i=0; i<(size_t)m_data_model->step; i++) 
    Pk[i] = m_data_model->cosmology->Pk(m_data_model->kk[i], m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);

  m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));

  if (m_data_model->Pk_mu_model=="dispersion_dewiggled") {    
    vector<double> PkNW(m_data_model->step, 0);
    for (size_t i=0; i<(size_t)m_data_model->step; i++) 
      PkNW[i] = m_data_model->cosmology->Pk(m_data_model->kk[i], "EisensteinHu", false, m_data_model->redshift, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    
    m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));
  }
  
  else if (m_data_model->Pk_mu_model=="dispersion_modecoupling") {
    vector<double> kk_1loop, Pk_1loop;
    for (size_t i=0; i<(size_t)m_data_model->step; i++) 
      if (m_data_model->kk[i] < par::pi) {
	kk_1loop.push_back(m_data_model->kk[i]);
	Pk_1loop.push_back(m_data_model->cosmology->Pk_1loop(m_data_model->kk[i], m_data_model->func_Pk, 0, m_data_model->k_min, 5., m_data_model->prec)); 
      }
    
    m_data_model->func_Pk1loop = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk_1loop, Pk_1loop, "Spline"));
  }

  else ErrorCBL("Error in cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_PkDM() of Modelling_TwoPointCorrelation_wedges.cpp: the chosen model ("+m_data_model->Pk_mu_model+") is not currently implemented!");
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_fiducial_xiDM ()
{
  m_data_model->nmultipoles = 3;
  m_data_model->nwedges = 2;

  cout << endl; coutCBL << "Setting up the fiducial dark matter two-point correlation function model" << endl;

  const vector<double> rad = linear_bin_vector(m_data_model->step, m_data_model->r_min, m_data_model->r_max);

  m_data_model->rr = rad;

  if (m_data_model->sigmaNL==0) {    

    vector<double> kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.)), Pk(m_data_model->step,0);

    for (size_t i=0; i<(size_t)m_data_model->step; i++) 
      Pk[i] =  m_data_model->cosmology->Pk(kk[i], m_data_model->method_Pk, m_data_model->NL, m_data_model->redshift, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);

    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk, "Spline"));

  }

  else {

    vector<double> Pk(m_data_model->step,0), PkNW(m_data_model->step,0);
    m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));

    for (size_t i=0; i<(size_t)m_data_model->step; i++) {
      Pk[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
      PkNW[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], "EisensteinHu", false, m_data_model->redshift, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    }

    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

  }

  vector<vector<double>> xil = Xi_l(rad, m_data_model->nmultipoles, 0, {1., 1., m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->bias, m_data_model->linear_growth_rate_z, 0.}, {m_data_model->func_Pk, m_data_model->func_Pk_NW}, m_data_model->prec);

  m_data_model->func_multipoles.erase(m_data_model->func_multipoles.begin(), m_data_model->func_multipoles.end());

  for(int i=0; i< m_data_model->nmultipoles; i++)
    m_data_model->func_multipoles.push_back(make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(rad, xil[i], "Spline")));

}



// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_fullShape_DeWiggled (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution SigmaNL_perpendicular_prior, const statistics::PriorDistribution SigmaNL_parallel_prior, statistics::PriorDistribution fsigma8_prior, statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaS_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_dewiggled";

  // compute the fiducial dark matter power spectrum terms used to construct the model
  if (compute_PkDM) set_fiducial_PkDM();
  
  // number of wedges to be used
  m_data_model->nwedges = m_nwedges_fit;
  
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
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_fullShape_ModeCoupling (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, statistics::PriorDistribution fsigma8_prior, statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaV_prior, const statistics::PriorDistribution AMC_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_modecoupling";

  // compute the fiducial dark matter power spectrum terms used to construct the model
  if (compute_PkDM) set_fiducial_PkDM();
  
  // number of wedges to be used
  m_data_model->nwedges = m_nwedges_fit;

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
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution Bperp_prior, const statistics::PriorDistribution Bpar_prior, const statistics::PriorDistribution Aperp0_prior, const statistics::PriorDistribution Apar0_prior, const statistics::PriorDistribution Aperp1_prior, const statistics::PriorDistribution Apar1_prior, const statistics::PriorDistribution Aperp2_prior, const statistics::PriorDistribution Apar2_prior, const bool compute_XiDM)
{
  // compute the fiducial dark matter two-point correlation function
  if (compute_XiDM) set_fiducial_xiDM();

  m_data_model->nwedges = m_nwedges_fit;
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
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters)
{
  if (m_likelihood==NULL) ErrorCBL("Error in cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model() of Modelling_TwoPointCorrelation_wedges.cpp: this function requires the likelihood to be defined (with the function set_likelihood)!");
  
  vector<int> dataset_order_original = m_data_model->dataset_order;

  vector<int> new_dataset_order;
  vector<double> new_xx;
  
  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nwedges; n++)
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }

  m_data_model->dataset_order = new_dataset_order;
  m_likelihood->write_model(output_dir, output_file, parameters, new_xx);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_wedges::write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx)
{
  if(m_posterior==NULL)
    ErrorCBL("Error in write_model_at_bestfit of Modelling_TwoPointCorrelation_wedges.cpp. No posterior found! Run maximize_posterior() first");

  vector<int> dataset_order_original = m_data_model->dataset_order;

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nwedges; n++)
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
  if(m_posterior==NULL)
    ErrorCBL("Error in write_model_from_chains of Modelling_TwoPointCorrelation_wedges.cpp. No posterior found! Run sample_posterior() first");

  vector<int> dataset_order_original = m_data_model->dataset_order;

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nwedges; n++)
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }

  m_data_model->dataset_order = new_dataset_order;
  m_posterior->write_model_from_chain(output_dir, output_file, new_xx, {}, start, thin);

  m_data_model->dataset_order = dataset_order_original;
}
