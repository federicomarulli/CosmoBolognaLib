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


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nmultipoles(3)
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;

  for (int j=0; j<m_nmultipoles; j++)
    for (int i=0; i<size; i++)
      m_multipoles_order.push_back(j);

  m_use_pole.resize(3, true);
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<data::Data> twop_dataset, const int nmultipoles)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nmultipoles(nmultipoles)
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());
  m_use_pole.resize(3, false);

  for (int j=0; j<m_nmultipoles; j++) {
    m_use_pole[j] = true;
    for (int i=0; i<m_data->ndata()/m_nmultipoles; i++)
      m_multipoles_order.push_back(j);
  }

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const double xmin, const double xmax, const int nmultipoles)
{
  vector<vector<double>> fr(m_nmultipoles, vector<double>(2, -1.));

  int mp = (0<nmultipoles && nmultipoles<m_nmultipoles) ? nmultipoles : m_nmultipoles;

  for (int i=0; i<mp; i++) {
    fr[i][0] = xmin;
    fr[i][1] = xmax;
  }

  set_fit_range(fr);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const std::vector<std::vector<double>> fit_range)
{
  if ((int)fit_range.size()!=m_nmultipoles)
    ErrorCBL("the dimension input matrix must be equal to the number of multipoles to be fitted, i.e."+conv(m_nmultipoles, par::fINT)+"!", "set_fit_range", "Modelling_TwoPointCorrelation_multipoles.cpp");

  m_use_pole = {false, false, false};

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  const int size = m_data->ndata()/m_nmultipoles;
  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;

  for (int j=0; j<m_nmultipoles; j++) {
    for (int i=0; i<size; i++) {
      if (fit_range[j][0]<m_data->xx(i+j*size) && m_data->xx(i+j*size)<fit_range[j][1]) {
	m_multipoles_order.push_back(j);
	xx.push_back(m_data->xx(i+j*size));
	m_use_pole[j] = true;
	mask[i+j*size] = true;
      }
    }
  }

  m_data_fit = m_data->cut(mask);
  
  m_fit_range = true; 

  if (m_ModelIsSet) 
    m_data_model->dataset_order = m_multipoles_order;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_PkDM ()
{
  m_data_model->nmultipoles = m_nmultipoles;

  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));

  vector<double> Pk(m_data_model->step, 0);

  for (size_t i=0; i<(size_t)m_data_model->step; i++)
    Pk[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);

  if (m_data_model->Pk_mu_model=="dispersion_dewiggled") {
    vector<double> PkNW(m_data_model->step,0);
    for (size_t i=0; i<(size_t)m_data_model->step; i++)
      PkNW[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_NW);
  }

  else if (m_data_model->Pk_mu_model=="dispersion_modecoupling") {
    vector<double> kk_1loop, Pk_1loop;
    for (size_t i=0; i<(size_t)m_data_model->step; i++) {
      if (m_data_model->kk[i] < par::pi) {
	kk_1loop.push_back(m_data_model->kk[i]);
	Pk_1loop.push_back(m_data_model->cosmology->Pk_1loop(m_data_model->kk[i], m_data_model->func_Pk, 0,  m_data_model->k_min, 5., m_data_model->prec));
      }
    }
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk1loop = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk_1loop, Pk_1loop, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk1loop);
  }

  else if (m_data_model->Pk_mu_model=="DispersionGauss"  || m_data_model->Pk_mu_model=="DispersionLorentz") {
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
  }

  else if (m_data_model->Pk_mu_model=="ScoccimarroPezzottaGauss"  || m_data_model->Pk_mu_model=="ScoccimarroPezzottaLorentz" || m_data_model->Pk_mu_model=="ScoccimarroBelGauss"  || m_data_model->Pk_mu_model=="ScoccimarroBelLorentz") {
    vector<double> Pknonlin(m_data_model->step,0);
    for (size_t i=0; i<(size_t)m_data_model->step; i++)
      Pknonlin[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], m_data_model->method_Pk, true, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_nonlin = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pknonlin, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_nonlin);
  }

  else if (m_data_model->Pk_mu_model=="ScoccimarroGauss"  || m_data_model->Pk_mu_model=="ScoccimarroLorentz") {
    vector<vector<double>> Pk_terms = m_data_model->cosmology->Pk_TNS_dd_dt_tt(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->output_dir, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[2], "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
  }

  else if (m_data_model->Pk_mu_model=="TaruyaGauss"  || m_data_model->Pk_mu_model=="TaruyaLorentz") {
    vector<vector<double>> Pk_terms = m_data_model->cosmology->Pk_TNS_dd_dt_tt(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->output_dir, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
    vector<vector<double>> Pk_AB    = m_data_model->cosmology->Pk_TNS_AB_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, m_data_model->output_dir, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

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

  else ErrorCBL("the chosen model ("+m_data_model->Pk_mu_model+") is not currently implemented!", "set_fiducial_PkDM", "Modelling_TwoPointCorrelation_multipoles.cpp");
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial two-point correlation function model" << endl;

  m_data_model->nmultipoles = 3;
  
  const vector<double> rad = linear_bin_vector(m_data_model->step, m_data_model->r_min, m_data_model->r_max);

  m_data_model->rr = rad;

  vector<double> Pk(m_data_model->step, 0), PkNW(m_data_model->step, 0);
  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));

  for (size_t i=0; i<(size_t)m_data_model->step; i++) {
    Pk[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    PkNW[i] =  m_data_model->cosmology->Pk(m_data_model->kk[i], "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
  }

  m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
  m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

  std::vector<double> template_parameters = {m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->linear_growth_rate_z, m_data_model->bias, 0.};
  cbl::Print(template_parameters);

  vector<vector<double>> xil = Xi_l(rad, m_data_model->nmultipoles, "dispersion_dewiggled", {m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->linear_growth_rate_z, m_data_model->bias, 0.}, {m_data_model->func_Pk, m_data_model->func_Pk_NW}, m_data_model->prec, 1, 1);

  m_data_model->func_multipoles.erase(m_data_model->func_multipoles.begin(), m_data_model->func_multipoles.end());

  for (int i=0; i< m_data_model->nmultipoles; i++)
    m_data_model->func_multipoles.push_back(make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(rad, xil[i], "Spline")));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_DeWiggled (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution SigmaNL_perpendicular_prior, const statistics::PriorDistribution SigmaNL_parallel_prior, statistics::PriorDistribution fsigma8_prior, statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaS_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_dewiggled";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

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
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Dispersion (statistics::PriorDistribution fsigma8_prior, statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigma12_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "DispersionGauss";
  else m_data_model->Pk_mu_model = "DispersionLorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 3;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigma12";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigma12_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Scoccimarro_fitPezzotta (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigma12_prior, const statistics::PriorDistribution kd_prior, const statistics::PriorDistribution kt_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "ScoccimarroPezzottaGauss";
  else m_data_model->Pk_mu_model = "ScoccimarroPezzottaLorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigma12";
  parameterName[3] = "kd";
  parameterName[4] = "kt";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigma12_prior, kd_prior, kt_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Scoccimarro_fitBel (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigma12_prior, const statistics::PriorDistribution kd_prior, const statistics::PriorDistribution bb_prior, const statistics::PriorDistribution a1_prior, const statistics::PriorDistribution a2_prior, const statistics::PriorDistribution a3_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "ScoccimarroBelGauss";
  else m_data_model->Pk_mu_model = "ScoccimarroBelLorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 8;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigma12";
  parameterName[3] = "kd";
  parameterName[4] = "bb";
  parameterName[5] = "a1";
  parameterName[6] = "a2";
  parameterName[7] = "a3";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigma12_prior, kd_prior, bb_prior, a1_prior, a2_prior, a3_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_TNS (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigma12_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "TaruyaGauss";
  else m_data_model->Pk_mu_model = "TaruyaLorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 3;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigma12";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigma12_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Scoccimarro (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigma12_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "ScoccimarroGauss";
  else m_data_model->Pk_mu_model = "ScoccimarroLorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 3;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigma12";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigma12_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_ModeCoupling (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaV_prior, const statistics::PriorDistribution AMC_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_modecoupling";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

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
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_sigma8_bias (const statistics::PriorDistribution sigma8_prior, const statistics::PriorDistribution bias_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 2;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "sigma8";
  parameterName[1] = "bias";

  vector<statistics::PriorDistribution> priors = {sigma8_prior, bias_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_sigma8_bias, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution B0_prior, const statistics::PriorDistribution B2_prior, const statistics::PriorDistribution A00_prior, const statistics::PriorDistribution A20_prior, const statistics::PriorDistribution A01_prior, const statistics::PriorDistribution A21_prior, const statistics::PriorDistribution A02_prior, const statistics::PriorDistribution A22_prior, const bool compute_XiTemplate, const bool isRealSpace)
{
    
  // compute the fiducial dark matter two-point correlation function

  if (m_nmultipoles>2)
    ErrorCBL("BAO modelling can be done only with two multipoles!", "set_model_BAO", "Modelling_TwoPointCorrelation_multipoles");

  if (isRealSpace) {
    double lgf = m_data_model->linear_growth_rate_z;
    m_data_model->linear_growth_rate_z = 0.;

    set_fiducial_xiDM();

    m_data_model->linear_growth_rate_z = lgf;
  }
  else if (compute_XiTemplate)
    set_fiducial_xiDM();
  
  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;

  // set the model parameters
  const int nparameters = 10;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

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

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_BAO, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model (const std::string output_dir, const std::string output_file, const int nmultipoles, const std::vector<double> xx, const std::vector<double> parameters)
{
  if (m_likelihood==NULL) ErrorCBL("this function requires the likelihood to be defined (with the function set_likelihood)!", "write_model", "Modelling_TwoPointCorrelation_multipoles.cpp");
  
  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles = nmultipoles;

  vector<bool> new_use_pole(3, false);
  vector<int> new_dataset_order;
  vector<double> new_xx;
  
  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      new_use_pole[n] = true;
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }
    }

  m_data_model->dataset_order = new_dataset_order;
  m_data_model->use_pole = new_use_pole;

  m_likelihood->write_model(output_dir, output_file, parameters, new_xx);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model_at_bestfit (const std::string output_dir, const std::string output_file, const int nmultipoles, const std::vector<double> xx)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run maximize_posterior() first!", "write_model_at_best_fit", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles=nmultipoles;

  vector<bool> new_use_pole(3, false);
  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      new_use_pole[n] = true;
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }
    }

  m_data_model->dataset_order = new_dataset_order;
  m_data_model->use_pole = new_use_pole;
  m_posterior->write_model_at_bestfit(output_dir, output_file, new_xx);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model_from_chains (const std::string output_dir, const std::string output_file, const int nmultipoles, const std::vector<double> xx, const int start, const int thin)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run sample_posterior() first!", "write_model_from_chains", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles=nmultipoles;

  vector<bool> new_use_pole(3, false);
  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      new_use_pole[n] = true;
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }
    }

  m_data_model->dataset_order = new_dataset_order;
  m_data_model->use_pole = new_use_pole;
  m_posterior->write_model_from_chain(output_dir, output_file, new_xx, {}, start, thin);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}
