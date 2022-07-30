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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D_monopole.cpp
 *
 *  @brief Methods of the class
 *  Modelling_TwoPointCorrelation1D_monopole
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation1D_monopole, used to model the
 *  monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation1D_monopole.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial dark matter two-point correlation function model..." << endl;

  const vector<double> rad = linear_bin_vector(m_data_model->step, m_data_model->r_min, m_data_model->r_max);
  vector<double> xi0;

  if (m_data_model->sigmaNL==0) {    

    vector<double> kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.)), Pk(m_data_model->step,0);

    Pk = m_data_model->cosmology->Pk_matter(kk, m_data_model->method_Pk, m_data_model->NL, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);

    m_data_model->kk = kk;
    m_data_model->func_Pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, Pk, "Spline"));
    xi0 = wrapper::fftlog::transform_FFTlog(rad, 1, kk, Pk, 0);
  }

  else {

    vector<double> kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.)), Pk, PkNW, PkDW(m_data_model->step,0);
    Pk = m_data_model->cosmology->Pk_matter(kk, m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    PkNW = m_data_model->cosmology->Pk_matter(kk, "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    
    for (size_t i=0; i<kk.size(); i++)
      PkDW[i] = PkNW[i]*(1+(Pk[i]/PkNW[i]-1)*exp(-0.5*pow(kk[i]*m_data_model->sigmaNL, 2)));

    m_data_model->kk = kk;
    m_data_model->func_Pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, Pk, "Spline"));
    m_data_model->func_Pk_NW = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, PkNW, "Spline"));

    xi0 = wrapper::fftlog::transform_FFTlog(rad, 1, kk, PkDW);
  }

  m_data_model->func_xi = make_shared<glob::FuncGrid>(glob::FuncGrid(rad, xi0, "Spline"));
  
  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_sigma_data_model ()
{
  // create the grid file if it doesn't exist yet
  const string file_grid = m_data_model->cosmology->create_grid_sigmaM(m_data_model->method_Pk, 0., true, m_data_model->output_root, "Spline", m_data_model->k_max);

  
  // read the grid file
  
  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 

  double MMass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>MMass>>Sigma>>Dln_Sigma) {
    mass.push_back(MMass);
    sigma.push_back(Sigma);
  }


  // create the function to interpolate sigma(M) and dlg(sigma(M)

  m_data_model->func_sigma = make_shared<glob::FuncGrid>(glob::FuncGrid(mass, sigma, "Linear", BinType::_logarithmic_));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_PkDM ()
{
  coutCBL << "Setting up the fiducial matter power spectrum model" << endl;

  const vector<double> kk = logarithmic_bin_vector(m_data_HOD->step, max(m_data_HOD->k_min, 1.e-4), min(m_data_HOD->k_max, 500.));
  vector<double> PkDM(kk.size());
  
  PkDM = m_data_HOD->cosmology->Pk_matter(kk, m_data_HOD->method_Pk, m_data_HOD->NL, m_data_HOD->redshift, m_data_HOD->store_output, m_data_HOD->output_root, m_data_HOD->norm, m_data_HOD->k_min, m_data_HOD->k_max, m_data_HOD->prec, m_data_HOD->input_file);

  m_data_HOD->func_Pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, PkDM, "Spline"));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_sigma ()
{
  // create the grid file if it doesn't exist yet
  const string file_grid = m_data_HOD->cosmology->create_grid_sigmaM(m_data_HOD->method_Pk, 0., true, m_data_HOD->output_root, m_data_HOD->interpType, m_data_HOD->k_max, m_data_HOD->input_file, m_data_HOD->is_parameter_file);

  
  // read the grid file
  
  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 

  double MMass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dln_sigma;

  while (fin >>MMass>>Sigma>>Dln_Sigma) {
    mass.push_back(MMass);
    sigma.push_back(Sigma);
    dln_sigma.push_back(Dln_Sigma);
  }


  // create the function to interpolate sigma(M) and dlg(sigma(M)

  m_data_HOD->func_sigma = make_shared<glob::FuncGrid>(glob::FuncGrid(mass, sigma, "Spline"));
  
  m_data_HOD->func_dlnsigma = make_shared<glob::FuncGrid>(glob::FuncGrid(mass, dln_sigma, "Spline"));
  
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_bias_eff_grid (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<double> min_par, const std::vector<double> max_par, const std::vector<int> nbins_par, const std::string dir, const std::string file_grid_bias)
{
  if (m_data_model->cluster_mass_proxy->ndata()==0) ErrorCBL("m_data_model->cluster_mass_proxy->ndata() is not defined!", "set_bias_eff_grid", "Modelling_TwoPointCorrelation1D_monopole.cpp");
  
  const int npar = cosmo_param.size(); 
 
  const string file = dir+file_grid_bias;
  ifstream fin(file.c_str());
  
  if (!fin) {

    if (npar==1) {

      fin.clear(); fin.close();
      vector<double> mass_grid = logarithmic_bin_vector(m_data_model->cluster_mass_proxy->ndata()/10, Min(m_data_model->cluster_mass_proxy->data()), Max(m_data_model->cluster_mass_proxy->data()));
      vector<double> parameter, bias_eff;
     
      m_data_model->cosmology->generate_bias_eff_grid_one_cosmopar(parameter, bias_eff, dir, file_grid_bias, cosmo_param[0], min_par[0], max_par[0], nbins_par[0], m_data_model->cluster_mass_proxy->data(), mass_grid, m_data_model->cluster_mass_proxy->xx(), m_data_model->model_bias, m_data_model->method_Pk, m_data_model->meanType, true, m_data_model->output_root, m_data_model->Delta, 1., "Spline", m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, "NULL", false, m_data_model->cosmology_mass, m_data_model->redshift_source);

      m_data_model->cosmopar_bias_interp_1D = bind(interpolated, placeholders::_1, parameter, bias_eff, "Spline");

    }
    
    else if (npar==2) {
      
      vector<double> mass_grid = logarithmic_bin_vector(m_data_model->cluster_mass_proxy->ndata()/10, Min(m_data_model->cluster_mass_proxy->data()), Max(m_data_model->cluster_mass_proxy->data()));

      vector<double> parameter1, parameter2;
      vector<vector<double>> bias_eff;
      m_data_model->cosmology->generate_bias_eff_grid_two_cosmopars(parameter1, parameter2, bias_eff, dir, file_grid_bias, cosmo_param[0], min_par[0], max_par[0], nbins_par[0], cosmo_param[1], min_par[1], max_par[1], nbins_par[1], m_data_model->cluster_mass_proxy->data(), mass_grid, m_data_model->cluster_mass_proxy->xx(), m_data_model->model_bias, m_data_model->method_Pk, m_data_model->meanType, true, m_data_model->output_root, m_data_model->Delta, 1., "Spline", m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, par::defaultString, false, m_data_model->cosmology_mass, m_data_model->redshift_source);
      
      m_data_model->cosmopar_bias_interp_2D = bind(interpolated_2D, placeholders::_1, placeholders::_2, parameter1, parameter2, bias_eff, "Cubic");

    }
    
    else 
      ErrorCBL("this function works with 1 or 2 cosmological parameters!", "set_bias_eff_grid", "Modelling_TwoPointCorrelation1D_monopole.cpp");
  }
  
  else {
    
    if (npar==1) {
      vector<double> parameter, bias_eff;

      string line;
      while (getline(fin, line)) {
	stringstream SS(line); double _p, _b;
	SS >> _p >> _b;
	parameter.push_back(_p);
	bias_eff.push_back(_b);
      }
      fin.clear(); fin.close();

      m_data_model->cosmopar_bias_interp_1D = bind(interpolated, placeholders::_1, parameter, bias_eff, "Spline");
    }
    
    else if (npar==2) {
      fin.clear(); fin.close();
      vector<double> parameter1, parameter2;
      vector<vector<double>> bias_eff;
      
      read_matrix(file, parameter1, parameter2, bias_eff);

      m_data_model->cosmopar_bias_interp_2D = bind(interpolated_2D, placeholders::_1, placeholders::_2, parameter1, parameter2, bias_eff, "Cubic");
    }
    
    else 
      ErrorCBL("this function works with 1 or 2 cosmological parameters!", "set_bias_eff_grid", "Modelling_TwoPointCorrelation1D_monopole.cpp");
  }
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_bias_eff_grid (const std::string file_selection_function, const std::vector<int> column, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<double> min_par, const std::vector<double> max_par, const std::vector<int> nbins_par, const std::string dir, const std::string file_grid_bias)
{
  const int npar = cosmo_param.size(); 
  
  if (npar==1) {
    vector<double> parameter, bias_eff;
    m_data_model->cosmology->generate_bias_eff_grid_one_cosmopar(parameter, bias_eff, dir, file_grid_bias, cosmo_param[0], min_par[0], max_par[0], nbins_par[0], m_data_model->redshift, m_data_model->Mass_min, m_data_model->Mass_max, m_data_model->model_bias, m_data_model->model_MF, m_data_model->method_Pk, file_selection_function, column, 1., true, m_data_model->output_root, m_data_model->Delta, 1., "Spline", m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
      
    m_data_model->cosmopar_bias_interp_1D = bind(interpolated, placeholders::_1, parameter, bias_eff, "Spline");
  }
  
  else if (npar==2) {
    ErrorCBL("", "set_bias_eff_grid", "Modelling_TwoPointCorrelation1D_monopole.cpp", glob::ExitCode::_workInProgress_);
  }
  
  else 
    ErrorCBL("this function works with 1 or 2 cosmological parameters!", "set_bias_eff_grid", "Modelling_TwoPointCorrelation1D_monopole.cpp");
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_BAO_sigmaNL (const statistics::PriorDistribution sigmaNL_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution B_prior, const statistics::PriorDistribution A0_prior, const statistics::PriorDistribution A1_prior, const statistics::PriorDistribution A2_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  // set the model parameters
  const int nparameters = 6;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "sigmaNL";
  parameterName[1] = "alpha";
  parameterName[2] = "B";
  parameterName[3] = "A0";
  parameterName[4] = "A1";
  parameterName[5] = "A2";

  vector<statistics::PriorDistribution> priors = {sigmaNL_prior, alpha_prior, B_prior, A0_prior, A1_prior, A2_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_BAO_sigmaNL, nparameters, parameterType, parameterName, m_data_model));

  coutCBL << "Done!" << endl << endl;

}

// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const std::vector<statistics::PriorDistribution> polynomial_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_data_model->poly_order = polynomial_prior.size(); 

  // set the model parameters
  const int nparameters = 3+polynomial_prior.size();

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters);
  priors[0] = alpha_prior;
  priors[1] = fsigma8_prior;
  priors[2] = bsigma8_prior;

  parameterName[0] = "alpha";
  parameterName[1] = "f*sigma8";
  parameterName[2] = "b*sigma8";

  for (size_t i=0; i<polynomial_prior.size(); i++) {
    parameterName[i+3] = "A"+conv(i, par::fINT);
    priors[i+3] = polynomial_prior[i];
  }

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear, nparameters, parameterType, parameterName, m_data_model));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_LinearPoint (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const std::vector<statistics::PriorDistribution> polynomial_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_data_model->poly_order = polynomial_prior.size(); 

  // set the model parameters
  const int nparameters = 6+polynomial_prior.size();

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[0] = statistics::ParameterType::_Derived_;
  parameterType[1] = statistics::ParameterType::_Derived_;
  parameterType[2] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters-3);
  priors[0] = alpha_prior;
  priors[1] = fsigma8_prior;
  priors[2] = bsigma8_prior;
  
  parameterName[0] = "peak";
  parameterName[1] = "dip";
  parameterName[2] = "linear_point";
  parameterName[3] = "alpha";
  parameterName[4] = "f*sigma8";
  parameterName[5] = "b*sigma8";

  for (size_t i=0; i<polynomial_prior.size(); i++) {
    parameterName[i+6] = "A"+conv(i, par::fINT);
    priors[i+3] = polynomial_prior[i];
  }

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_LinearPoint, nparameters, parameterType, parameterName, m_data_model));

}


// =========================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_polynomial_LinearPoint (const std::vector<statistics::PriorDistribution> polynomial_prior)
{
  // set the model parameters
  const int nparameters = polynomial_prior.size()+3;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[0] = statistics::ParameterType::_Derived_;
  parameterType[1] = statistics::ParameterType::_Derived_;
  parameterType[2] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters-3);
  
  parameterName[0] = "peak";
  parameterName[1] = "dip";
  parameterName[2] = "linear_point";

  for (size_t i=0; i<polynomial_prior.size(); i++) {
    parameterName[i+3] = "a"+conv(i, par::fINT);
    priors[i] = polynomial_prior[i];
  }

  // set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_polynomial_LinearPoint, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_Kaiser (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior)
{
  set_model_linear(statistics::PriorDistribution(glob::DistributionType::_Constant_, 1), fsigma8_prior, bsigma8_prior, {});
}


// ============================================================================================
	

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_bias_cosmology (const statistics::PriorDistribution bias_prior, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior)
{
  // set the model parameters
  const int nparameters = cosmo_param.size()+1;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters-1);
  
  parameterName[0] = "bias";
  priors[0] = bias_prior;

  for (size_t i=0; i<cosmo_param.size(); i++) {
    parameterName[i+1] = cosmology::CosmologicalParameter_name(cosmo_param[i]);
    priors[i+1] = cosmo_param_prior[i];
  }

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_bias_cosmology, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================
	

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_sigma8_clusters (const statistics::PriorDistribution sigma8_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  // set the model parameters
  const int nparameters = 3;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[1] = statistics::ParameterType::_Derived_;
  parameterType[2] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);

  parameterName[0] = "sigma8";
  parameterName[1] = "bias";
  parameterName[2] = "sigma8(z)";

  vector<statistics::PriorDistribution> priors = {sigma8_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_sigma8_clusters, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================
	

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters_grid (const cbl::cosmology::CosmologicalParameter cosmo_param, const statistics::PriorDistribution cosmo_param_prior, const std::string dir, const std::string file_grid_bias, const double min_par, const double max_par, const int nbins_par, const std::string file_selection_function, const std::vector<int> column)
{
  // set the free cosmological parameter
  m_data_model->Cpar = {cosmo_param};
  
  // compute the bias on a grid, using a selection function if provided
  if (file_selection_function!=par::defaultString)
    set_bias_eff_grid(file_selection_function, column, {cosmo_param}, {min_par}, {max_par}, {nbins_par}, dir, file_grid_bias);
  else
    set_bias_eff_grid({cosmo_param}, {min_par}, {max_par}, {nbins_par}, dir, file_grid_bias);

  // set the model parameters
  const int nparameters = 2;

  vector<statistics::ParameterType> parameterType(nparameters);
  parameterType[0] = statistics::ParameterType::_Base_;
  parameterType[1] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);

  parameterName[0] = cosmology::CosmologicalParameter_name(cosmo_param);
  parameterName[1] = "bias";

  vector<statistics::PriorDistribution> priors = {cosmo_param_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_one_cosmo_par_clusters, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================
	

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters_grid (const cbl::cosmology::CosmologicalParameter cosmo_param1, const statistics::PriorDistribution cosmo_param_prior1, const cbl::cosmology::CosmologicalParameter cosmo_param2, const statistics::PriorDistribution cosmo_param_prior2, const std::string dir, const std::string file_grid_bias, const double min_par1, const double max_par1, const int nbins_par1, const double min_par2, const double max_par2, const int nbins_par2, const std::string file_selection_function, const std::vector<int> column)
{
  // set the two free cosmological parameters
  m_data_model->Cpar = {cosmo_param1, cosmo_param2};
  
  // compute the bias on a grid, using a selection function if provided
  (void)column;
  if (file_selection_function!=par::defaultString)
    ErrorCBL("", "set_model_linear_cosmology_clusters_grid", "Modelling_TwoPointCorrelation1D_monopole.cpp", glob::ExitCode::_workInProgress_);
  else 
    set_bias_eff_grid({cosmo_param1, cosmo_param2}, {min_par1, min_par2}, {max_par1, max_par2}, {nbins_par1, nbins_par2}, dir, file_grid_bias);

  // set the model parameters
  const int nparameters = 4;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[2] = statistics::ParameterType::_Derived_;
  parameterType[3] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);

  parameterName[0] = cosmology::CosmologicalParameter_name(cosmo_param1);
  parameterName[1] = cosmology::CosmologicalParameter_name(cosmo_param2);
  parameterName[2] = "bias";
  parameterName[3] = "alpha";

  vector<statistics::PriorDistribution> priors = {cosmo_param_prior1, cosmo_param_prior2};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_two_cosmo_pars_clusters, nparameters, parameterType, parameterName, m_data_model));

}


// ============================================================================================
	

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior)
{
  set_fiducial_xiDM();

  m_data_model->Cpar = cosmo_param;

  // set the model parameters
  const int nparameters = (cosmo_param.size());

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);

  for (size_t i=0; i<cosmo_param.size(); i++) 
    parameterName[i] = cosmology::CosmologicalParameter_name(cosmo_param[i]);

  vector<statistics::PriorDistribution> priors = cosmo_param_prior;

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_cosmology_clusters, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_sigma8_bias (const statistics::PriorDistribution sigma8_prior, const statistics::PriorDistribution bias_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

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
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_sigma8_bias, nparameters, parameterType, parameterName, m_data_model));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_BAO (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution bs8_prior, const statistics::PriorDistribution A0_prior, const statistics::PriorDistribution A1_prior, const statistics::PriorDistribution A2_prior)
{
  vector<statistics::PriorDistribution>  polynomial_prior = {A0_prior, A1_prior, A2_prior};
  
  vector<double> polynomial_value(polynomial_prior.size(), par::defaultDouble);

  set_model_linear(alpha_prior, statistics::PriorDistribution(glob::DistributionType::_Constant_, 0), bs8_prior, polynomial_prior);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_BAO_LinearPoint (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution BB_prior, const statistics::PriorDistribution A0_prior, const statistics::PriorDistribution A1_prior, const statistics::PriorDistribution A2_prior)
{
  vector<statistics::PriorDistribution> polynomial_prior = {A0_prior, A1_prior, A2_prior};
  
  vector<double> polynomial_value(polynomial_prior.size(), par::defaultDouble);

  set_model_linear_LinearPoint(alpha_prior, statistics::PriorDistribution(glob::DistributionType::_Constant_, 0), BB_prior, polynomial_prior);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_scaling_relation_sigmaz (const statistics::PriorDistribution M0_prior, const statistics::PriorDistribution slope_prior, const statistics::PriorDistribution scatter_prior, const statistics::PriorDistribution sigmaz_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();
  set_fiducial_sigma_data_model();

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[nparameters-1] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);

  parameterName[0] = "M0";
  parameterName[1] = "slope";
  parameterName[2] = "scatter"; 
  parameterName[3] = "sigmaz";
  parameterName[4] = "bias";

  vector<statistics::PriorDistribution> priors = {M0_prior, slope_prior, scatter_prior, sigmaz_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_damped_scaling_relation_sigmaz, nparameters, parameterType, parameterName, m_data_model));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_scaling_relation_sigmaz_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution sigmaz_prior, const std::string z_evo)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  // Set the scaling relation Modelling object
  if ((m_data_model->scaling_relation)->data()->xx().size() == 0)
    ErrorCBL("The mass-observable relation is not set! Use the correct set_data_model().", "set_model_scaling_relation_sigmaz_cosmology", "Modelling_TwoPointCorrelation1D_monopole.cpp");

  (m_data_model->scaling_relation)->set_model_MassObservableRelation_cosmology(z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  (m_data_model->scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {}); // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)
    
  // Set the parameter names and priors
  m_data_model->Cpar = cosmo_param;

  const size_t nParams_scaling_relation = (m_data_model->scaling_relation)->likelihood()->parameters()->nparameters();
  const size_t nParams_base = 1 + nParams_scaling_relation; // sigmaz + scaling relation parameters (which include the cosmological parameters)
  const size_t nParams_derived = 1; // the effective bias is derived from the scaling relation
  
  const size_t nParams = nParams_base + nParams_derived;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[nParams-1] = statistics::ParameterType::_Derived_;
  
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams_base);

  // Cosmological and scaling relation parameters
  for (size_t i=0; i<nParams_base-1; i++) {
    Par_string[i] = (m_data_model->scaling_relation)->likelihood()->parameters()->name(i);
    param_prior[i] = *(m_data_model->scaling_relation)->get_prior(i);
  }

  // sigmaz
  Par_string[nParams_base-1] = "sigmaz";
  param_prior[nParams_base-1] = sigmaz_prior;

  // bias
  Par_string[nParams-1] = "bias";

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_damped_scaling_relation_sigmaz_cosmology, nParams, Par_type, Par_string, m_data_model));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_bias_sigmaz (const statistics::PriorDistribution bias_prior, const statistics::PriorDistribution sigmaz_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();
  set_fiducial_sigma_data_model();

  // set the model parameters
  const int nparameters = 2;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);

  parameterName[0] = "bias";
  parameterName[1] = "sigmaz";

  vector<statistics::PriorDistribution> priors = {bias_prior, sigmaz_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_damped_bias_sigmaz, nparameters, parameterType, parameterName, m_data_model));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_HOD (const statistics::PriorDistribution Mmin_prior, const statistics::PriorDistribution sigmalgM_prior, const statistics::PriorDistribution M0_prior, const statistics::PriorDistribution M1_prior, const statistics::PriorDistribution alpha_prior)
{
  // compute the fiducial dark matter power spectrum
  set_fiducial_PkDM();

  // compute the fiducial mass variance and its logarithmic derivative
  set_fiducial_sigma();

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);

  parameterName[0] = "Mmin";
  parameterName[1] = "sigmalgM";
  parameterName[2] = "M0";
  parameterName[3] = "M1";
  parameterName[4] = "alpha";

  vector<statistics::PriorDistribution> priors = {Mmin_prior, sigmalgM_prior, M0_prior, M1_prior, alpha_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi_HOD, nparameters, parameterType, parameterName, m_data_HOD));
}


// ============================================================================================
	

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_cluster_selection_function (const statistics::PriorDistribution alpha_prior, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior)
{
  // set the model parameters
  const int nparameters = (cosmo_param.size()+2);

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[nparameters-1] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);
  vector<statistics::PriorDistribution> priors(nparameters-1);

  for (size_t i=0; i<cosmo_param.size(); i++) {
    parameterName[i] = cosmology::CosmologicalParameter_name(cosmo_param[i]);
    priors[i] = cosmo_param_prior[i];
  }
  
  priors[nparameters-2] = alpha_prior;

  parameterName[nparameters-2] = "alpha";
  parameterName[nparameters-1] = "bias";

  // add the cosmological parameters
  m_data_model->Cpar = cosmo_param;

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_cosmology_clusters_selection_function, nparameters, parameterType, parameterName, m_data_model));
}
