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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation1D_monopole.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial dark matter two-point correlation function model..." << endl;

  const vector<double> rad = linear_bin_vector(m_data_model.step, m_data_model.r_min, m_data_model.r_max);
  vector<double> xi0;

  if (m_data_model.sigmaNL==0) {    

    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0);

    for (size_t i=0; i<kk.size(); i++)
      Pk[i] = m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);

    m_data_model.kk = kk;
    m_data_model.func_Pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, Pk, "Spline"));
    xi0 = fftlog::transform_FFTlog(rad, 1, kk, Pk, 0);
  }

  else {

    vector<double> kk = logarithmic_bin_vector(m_data_model.step, max(m_data_model.k_min, 1.e-4), min(m_data_model.k_max, 500.)), Pk(m_data_model.step,0), PkNW(m_data_model.step,0), PkDW(m_data_model.step,0);
    for (size_t i=0; i<kk.size(); i++) {
      Pk[i] = m_data_model.cosmology->Pk(kk[i], m_data_model.method_Pk, false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkNW[i] = m_data_model.cosmology->Pk(kk[i], "EisensteinHu", false, m_data_model.redshift, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec, m_data_model.file_par);
      PkDW[i] = PkNW[i]*(1+(Pk[i]/PkNW[i]-1)*exp(-0.5*pow(kk[i]*m_data_model.sigmaNL, 2)));
    }

    m_data_model.kk = kk;
    m_data_model.func_Pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, Pk, "Spline"));
    m_data_model.func_Pk_NW = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, PkNW, "Spline"));

    xi0 = fftlog::transform_FFTlog(rad, 1, kk, PkDW);
  }

  m_data_model.func_xi = make_shared<glob::FuncGrid>(glob::FuncGrid(rad, xi0, "Spline"));
  
  coutCBL << "Done!" << endl << endl;
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_sigma_data_model ()
{
  // create the grid file if it doesn't exist yet
  const string file_grid = m_data_model.cosmology->create_grid_sigmaM(m_data_model.method_Pk, 0., m_data_model.output_root, "Spline", m_data_model.k_max);

  
  // read the grid file
  
  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 

  double MMass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>MMass>>Sigma>>Dln_Sigma) {
    mass.push_back(MMass);
    sigma.push_back(Sigma);
  }


  // create the function to interpolate sigma(M) and dlg(sigma(M)

  m_data_model.func_sigma = make_shared<glob::FuncGrid>(glob::FuncGrid(mass, sigma, "Linear", binType::_logarithmic_));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_PkDM ()
{
  coutCBL << "Setting up the fiducial matter power spectrum model" << endl;

  const vector<double> kk = logarithmic_bin_vector(m_data_HOD.step, max(m_data_HOD.k_min, 1.e-4), min(m_data_HOD.k_max, 500.));
  vector<double> PkDM(kk.size());
  
  for (size_t i=0; i<kk.size(); i++)
    PkDM[i] = m_data_HOD.cosmology->Pk(kk[i], m_data_HOD.method_Pk, m_data_HOD.NL, m_data_HOD.redshift, m_data_HOD.output_root, m_data_HOD.norm, m_data_HOD.k_min, m_data_HOD.k_max, m_data_HOD.prec, m_data_HOD.input_file);

  m_data_HOD.func_Pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, PkDM, "Spline"));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_fiducial_sigma ()
{
  // create the grid file if it doesn't exist yet
  
  const string file_grid = m_data_HOD.cosmology->create_grid_sigmaM(m_data_HOD.method_Pk, 0., m_data_HOD.output_root, m_data_HOD.interpType, m_data_HOD.k_max, m_data_HOD.input_file, m_data_HOD.is_parameter_file);

  
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

  m_data_HOD.func_sigma = make_shared<glob::FuncGrid>(glob::FuncGrid(mass, sigma, "Spline"));
  
  m_data_HOD.func_dlnsigma = make_shared<glob::FuncGrid>(glob::FuncGrid(mass, dln_sigma, "Spline"));
  
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_bias_eff_grid (const vector<cosmobl::cosmology::CosmoPar> cosmo_param, const vector<double> min_par, const vector<double> max_par, const vector<int> nbins_par, const string dir, const string file_grid_bias)
{
  if (m_data_model.cluster_mass_proxy->ndata()==0) ErrorCBL("Error in cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_bias_eff_grid of Modelling_TwoPointCorrelation1D_monopole.cpp: m_data_model.cluster_mass_proxy->ndata() is not defined!");
  
  const int npar = cosmo_param.size(); 

  const string file = dir+file_grid_bias;
  ifstream fin(file.c_str());
  
  if (!fin) {

    if (npar==1) {

      fin.clear(); fin.close();
      vector<double> mass_grid = logarithmic_bin_vector(m_data_model.cluster_mass_proxy->ndata()/10, Min(m_data_model.cluster_mass_proxy->data()), Max(m_data_model.cluster_mass_proxy->data()));
      vector<double> parameter, bias_eff;
      m_data_model.cosmology->generate_bias_eff_grid_one_cosmopar(parameter, bias_eff, dir, file_grid_bias, cosmo_param[0], min_par[0], max_par[0], nbins_par[0], m_data_model.cluster_mass_proxy->data(), mass_grid, m_data_model.cluster_mass_proxy->xx(), m_data_model.model_bias, m_data_model.method_Pk, m_data_model.output_root, m_data_model.Delta, 1., "Spline", m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec);

      m_data_model.cosmopar_bias_interp_1D = bind(interpolated , placeholders::_1, parameter, bias_eff, "Spline");

    }
    
    else if (npar==2) {
      
      vector<double> mass_grid = logarithmic_bin_vector(m_data_model.cluster_mass_proxy->ndata()/10, Min(m_data_model.cluster_mass_proxy->data()), Max(m_data_model.cluster_mass_proxy->data()));

      vector<double> parameter1, parameter2;
      vector<vector<double>> bias_eff;
      m_data_model.cosmology->generate_bias_eff_grid_two_cosmopars(parameter1, parameter2, bias_eff, dir, file_grid_bias, cosmo_param[0], min_par[0], max_par[0], nbins_par[0], cosmo_param[1], min_par[1], max_par[1], nbins_par[1], m_data_model.cluster_mass_proxy->data(), mass_grid, m_data_model.cluster_mass_proxy->xx(), m_data_model.model_bias, m_data_model.method_Pk, m_data_model.output_root, m_data_model.Delta, 1., "Spline", m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec);

      m_data_model.cosmopar_bias_interp_2D = bind(interpolated_2D , placeholders::_1, placeholders::_2, parameter1, parameter2, bias_eff, "Cubic");
    }
    
    else 
      ErrorCBL("Error in set_bias_eff_grid om ModellingTwoPointCorrelation1D_monopole.cpp, this function works with 1 or 2 cosmological parameters.");
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

      m_data_model.cosmopar_bias_interp_1D = bind(interpolated, placeholders::_1, parameter, bias_eff, "Spline");
    }
    
    else if (npar==2) {
      fin.clear(); fin.close();
      vector<double> parameter1, parameter2;
      vector<vector<double>> bias_eff;
      
      read_matrix(file, parameter1, parameter2, bias_eff);

      m_data_model.cosmopar_bias_interp_2D = bind(interpolated_2D , placeholders::_1, placeholders::_2, parameter1, parameter2, bias_eff, "Cubic");
    }
    
    else 
      ErrorCBL("Error in set_bias_eff_grid om ModellingTwoPointCorrelation1D_monopole.cpp, this function works with 1 or 2 cosmological parameters.");
  }
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_bias_eff_grid (const string file_selection_function, const vector<int> column, const vector<cosmobl::cosmology::CosmoPar> cosmo_param, const vector<double> min_par, const vector<double> max_par, const vector<int> nbins_par, const string dir, const string file_grid_bias)
{
  const int npar = cosmo_param.size(); 

  if (npar==1) {
    vector<double> parameter, bias_eff;
    m_data_model.cosmology->generate_bias_eff_grid_one_cosmopar(parameter, bias_eff, dir, file_grid_bias, cosmo_param[0], min_par[0], max_par[0], nbins_par[0], m_data_model.redshift,  m_data_model.Mass_min, m_data_model.Mass_max, m_data_model.model_bias, m_data_model.model_MF, m_data_model.method_Pk, file_selection_function, column, 1., m_data_model.output_root, m_data_model.Delta, 1., "Spline", m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec);
      
    m_data_model.cosmopar_bias_interp_1D = bind(interpolated , placeholders::_1, parameter, bias_eff, "Spline");
  }
  
  else if (npar==2) {
    ErrorCBL("Work in progress in cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_bias_eff_grid of Modelling_TwoPointCorrelation1D_monopole.cpp!", glob::ExitCode::_workInProgress_);
  }
  
  else 
    ErrorCBL("Error in set_bias_eff_grid om ModellingTwoPointCorrelation1D_monopole.cpp, this function works with 1 or 2 cosmological parameters.");
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_BAO_sigmaNL (const statistics::Prior sigmaNL_prior, const statistics::Prior alpha_prior, const statistics::Prior B_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_sigmaNL = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigmaNL_prior, "sigmaNL"));

  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));

  m_bias = make_shared<statistics::BaseParameter>(statistics::BaseParameter(B_prior, "B"));
  
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_sigmaNL, m_alpha, m_bias};
  
  m_polynomial.erase(m_polynomial.begin(), m_polynomial.end());

  m_data_model.poly_order = 3; 

  m_polynomial.resize(m_data_model.poly_order);

  m_polynomial[0] = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A0_prior, "A0"));
  m_polynomial[1] = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A1_prior, "A1"));
  m_polynomial[2] = make_shared<statistics::BaseParameter>(statistics::BaseParameter(A2_prior, "A2"));

  for (size_t i=0; i<m_polynomial.size(); i++) 
    ll_parameters.push_back(m_polynomial[i]);
  

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_BAO_sigmaNL, inputs));

}

// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear (const statistics::Prior alpha_prior, const statistics::Prior fsigma8_prior, const statistics::Prior bsigma8_prior, const vector<statistics::Prior> polynomial_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));

  m_fsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(fsigma8_prior, "f*sigma8"));

  m_bsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bsigma8_prior, "b*sigma8"));
  
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_alpha, m_fsigma8, m_bsigma8};
  
  m_polynomial.erase(m_polynomial.begin(), m_polynomial.end());

  if (polynomial_prior.size()>0) {

    m_data_model.poly_order = polynomial_prior.size(); 

    m_polynomial.resize(polynomial_prior.size());

    for (size_t i=0; i<m_polynomial.size(); i++) {
      m_polynomial[i] = make_shared<statistics::BaseParameter>(statistics::BaseParameter(polynomial_prior[i], "A"+conv(i, par::fINT)));

      ll_parameters.push_back(m_polynomial[i]);
    }
  }
  else
    m_data_model.poly_order = 0; 

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_LinearPoint (const statistics::Prior alpha_prior, const statistics::Prior fsigma8_prior, const statistics::Prior bsigma8_prior, const vector<statistics::Prior> polynomial_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_peak = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("peak"));
  
  m_dip = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("dip"));

  m_linear_point = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("linear_point"));

  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));

  m_fsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(fsigma8_prior, "f*sigma8"));

  m_bsigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bsigma8_prior, "b*sigma8"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_peak, m_dip, m_linear_point, m_alpha, m_fsigma8, m_bsigma8};

  m_polynomial.erase(m_polynomial.begin(), m_polynomial.end());

  if (polynomial_prior.size()>0) {

    m_data_model.poly_order=polynomial_prior.size(); 

    m_polynomial.resize(polynomial_prior.size());

    for (size_t i=0; i<m_polynomial.size(); i++) {
      m_polynomial[i] = make_shared<statistics::BaseParameter>(statistics::BaseParameter(polynomial_prior[i], "A"+conv(i, par::fINT)));

      ll_parameters.push_back(m_polynomial[i]);
    }
  }
  else
    m_data_model.poly_order = 0; 
 
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_LinearPoint, inputs));

}


// =========================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_polynomial_LinearPoint (const vector<statistics::Prior> polynomial_prior)
{

  m_peak = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("peak"));
  
  m_dip = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("dip"));

  m_linear_point = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("linear_point"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_peak, m_dip, m_linear_point};

  m_polynomial.erase(m_polynomial.begin(), m_polynomial.end());

  if (polynomial_prior.size()>0) {

    m_data_model.poly_order=polynomial_prior.size(); 

    m_polynomial.resize(polynomial_prior.size());

    for (size_t i =0; i<m_polynomial.size(); i++) {
      m_polynomial[i] = make_shared<statistics::BaseParameter>(statistics::BaseParameter(polynomial_prior[i], "A"+conv(i, par::fINT)));

      ll_parameters.push_back(m_polynomial[i]);
    }
  }

  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_polynomial_LinearPoint, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_Kaiser (const statistics::Prior fsigma8_prior, const statistics::Prior bsigma8_prior)
{
  set_model_linear(statistics::Prior(glob::DistributionType::_ConstantDistribution_, 1), fsigma8_prior, bsigma8_prior, {});
}


// ============================================================================================
	

void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_bias_cosmology (const statistics::Prior bias_prior, const vector<cosmobl::cosmology::CosmoPar> cosmo_param, const vector<statistics::Prior> cosmo_param_prior)
{
  auto bias = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias_prior, "bias"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {bias};

  m_data_model.Cpar = cosmo_param;
  
  for (size_t i=0; i<m_data_model.Cpar.size(); i++) {
    auto cosmopar = make_shared<statistics::BaseParameter>(statistics::BaseParameter(cosmo_param_prior[i], cosmology::CosmoPar_name(m_data_model.Cpar[i])));
    ll_parameters.push_back(cosmopar);
  }
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_bias_cosmology, inputs));
}


// ============================================================================================
	

void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_sigma8_clusters (const statistics::Prior sigma8_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_sigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigma8_prior, "sigma8"));

  m_bias = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("bias"));

  auto sigma8_z = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("sigma8(z)"));
  
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_sigma8, m_bias, sigma8_z};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_sigma8_clusters, inputs));
}


// ============================================================================================
	

void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters_grid (const cosmobl::cosmology::CosmoPar cosmo_param, const statistics::Prior cosmo_param_prior, const string dir, const string file_grid_bias, const double min_par, const double max_par, const int nbins_par, const string file_selection_function, const vector<int> column)
{
  // set the free cosmological parameter
  m_data_model.Cpar = {cosmo_param};

  // compute the bias on a grid, using a selection function if provided
  if (file_selection_function!=par::defaultString)
    set_bias_eff_grid(file_selection_function, column, {cosmo_param}, {min_par}, {max_par}, {nbins_par}, dir, file_grid_bias);
  else
    set_bias_eff_grid({cosmo_param}, {min_par}, {max_par}, {nbins_par}, dir, file_grid_bias);

  // set the base and derived model parameters
  auto cosmo_par = make_shared<statistics::BaseParameter>(statistics::BaseParameter(cosmo_param_prior, cosmology::CosmoPar_name(cosmo_param)));
  m_bias = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("bias"));

  // add the model parameters
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {cosmo_par, m_bias};
  set_parameters(ll_parameters);

  // set the input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_one_cosmo_par_clusters, inputs));
}


// ============================================================================================
	

void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters_grid (const cosmobl::cosmology::CosmoPar cosmo_param1, const statistics::Prior cosmo_param_prior1, const cosmobl::cosmology::CosmoPar cosmo_param2, const statistics::Prior cosmo_param_prior2, const string dir, const string file_grid_bias, const double min_par1, const double max_par1, const int nbins_par1, const double min_par2, const double max_par2, const int nbins_par2, const string file_selection_function, const vector<int> column)
{
  // set the two free cosmological parameters
  m_data_model.Cpar = {cosmo_param1, cosmo_param2};

  // compute the bias on a grid, using a selection function if provided
  (void)column;
  if (file_selection_function!=par::defaultString)
    ErrorCBL("Work in progress in cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters_grid of Modelling_TwoPointCorrelation1D_monopole.cpp", glob::ExitCode::_workInProgress_);
  else 
    set_bias_eff_grid({cosmo_param1, cosmo_param2}, {min_par1, min_par2}, {max_par1, max_par2}, {nbins_par1, nbins_par2}, dir, file_grid_bias);
  
  // set the base and derived model parameters
  auto cosmo_par1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(cosmo_param_prior1, cosmology::CosmoPar_name(cosmo_param1)));
  auto cosmo_par2 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(cosmo_param_prior2, cosmology::CosmoPar_name(cosmo_param2)));
  m_bias = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("bias"));

  // add the model parameters
  vector<shared_ptr<statistics::Parameter>> ll_parameters = {cosmo_par1, cosmo_par2, m_bias};
  set_parameters(ll_parameters);

  // set the input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_two_cosmo_pars_clusters, inputs));
}


// ============================================================================================
	

void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_clusters (const vector<cosmobl::cosmology::CosmoPar> cosmo_param, const vector<statistics::Prior> cosmo_param_prior)
{
  set_fiducial_xiDM();

  m_data_model.Cpar = cosmo_param;

  vector<shared_ptr<statistics::Parameter>> ll_parameters;

  for (size_t i=0; i<m_data_model.Cpar.size(); i++) {
    auto cosmopar = make_shared<statistics::BaseParameter>(statistics::BaseParameter(cosmo_param_prior[i], cosmology::CosmoPar_name(m_data_model.Cpar[i])));
    ll_parameters.push_back(cosmopar);
  }
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_cosmology_clusters, inputs));
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_sigma8_bias (const statistics::Prior sigma8_prior, const statistics::Prior bias_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();

  m_sigma8 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigma8_prior, "sigma8"));

  m_bias = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias_prior, "bias"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_sigma8, m_bias};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_sigma8_bias, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_BAO (const statistics::Prior alpha_prior, const statistics::Prior BB_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior)
{
  vector<statistics::Prior>  polynomial_prior = {A0_prior, A1_prior, A2_prior};
  
  vector<double> polynomial_value(polynomial_prior.size(), par::defaultDouble);

  set_model_linear(alpha_prior, statistics::Prior(glob::DistributionType::_ConstantDistribution_, 0), BB_prior, polynomial_prior);
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_BAO_LinearPoint (const statistics::Prior alpha_prior, const statistics::Prior BB_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior)
{
  vector<statistics::Prior> polynomial_prior = {A0_prior, A1_prior, A2_prior};
  
  vector<double> polynomial_value(polynomial_prior.size(), par::defaultDouble);

  set_model_linear_LinearPoint(alpha_prior, statistics::Prior(glob::DistributionType::_ConstantDistribution_, 0), BB_prior, polynomial_prior);
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_scaling_relation_sigmaz (const statistics::Prior M0_prior, const statistics::Prior slope_prior, const statistics::Prior scatter_prior, const statistics::Prior sigmaz_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();
  set_fiducial_sigma_data_model();

  auto M0 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(M0_prior, "M0"));

  auto slope = make_shared<statistics::BaseParameter>(statistics::BaseParameter(slope_prior, "slope"));
  
  auto scatter = make_shared<statistics::BaseParameter>(statistics::BaseParameter(scatter_prior, "scatter"));

  auto sigmaz = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigmaz_prior, "sigmaz"));

  auto bias = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("bias"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {M0, slope, scatter, sigmaz, bias};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_damped_scaling_relation_sigmaz, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_bias_sigmaz (const statistics::Prior bias_prior, const statistics::Prior sigmaz_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_xiDM();
  set_fiducial_sigma_data_model();

  m_bias = make_shared<statistics::BaseParameter>(statistics::BaseParameter(bias_prior, "bias"));

  m_sigmaz = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigmaz_prior, "sigmaz"));

  vector<shared_ptr<statistics::Parameter>> ll_parameters = {m_bias, m_sigmaz};
  
  set_parameters(ll_parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_damped_bias_sigmaz, inputs));

}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_HOD (const statistics::Prior Mmin_prior, const statistics::Prior sigmalgM_prior, const statistics::Prior M0_prior, const statistics::Prior M1_prior, const statistics::Prior alpha_prior)
{
  // compute the fiducial dark matter power spectrum
  set_fiducial_PkDM();
  
  // compute the fiducial mass variance and its logarithmic derivative
  set_fiducial_sigma();
  
    
  // set the parameter of the HOD model

  m_Mmin = make_shared<statistics::BaseParameter>(statistics::BaseParameter(Mmin_prior, "Mmin"));

  m_sigmalgM = make_shared<statistics::BaseParameter>(statistics::BaseParameter(sigmalgM_prior, "sigmalgM"));

  m_M0 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(M0_prior, "M0"));
  
  m_M1 = make_shared<statistics::BaseParameter>(statistics::BaseParameter(M1_prior, "M1"));

  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));

  set_parameters({m_Mmin, m_sigmalgM, m_M0, m_M1, m_alpha});

  
  // input data used to construct the model
  auto inputs = make_shared<STR_data_HOD>(m_data_HOD);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi_HOD, inputs));
}


// ============================================================================================
	

void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation1D_monopole::set_model_linear_cosmology_cluster_selection_function (const statistics::Prior alpha_prior, const vector<cosmobl::cosmology::CosmoPar> cosmo_param, const vector<statistics::Prior> cosmo_param_prior)
{
  // vector of model parameters
  vector<shared_ptr<statistics::Parameter>> parameters;

  // add the cosmological parameters
  m_data_model.Cpar = cosmo_param;
  for (size_t i=0; i<m_data_model.Cpar.size(); i++) {
    auto cosmopar = make_shared<statistics::BaseParameter>(statistics::BaseParameter(cosmo_param_prior[i], cosmology::CosmoPar_name(m_data_model.Cpar[i])));
    parameters.push_back(cosmopar);
  }

  // add the alpha parameter (of the mass scaling relation)
  m_alpha = make_shared<statistics::BaseParameter>(statistics::BaseParameter(alpha_prior, "alpha"));
  parameters.push_back(m_alpha);

  // add the bias parameter
  m_bias = make_shared<statistics::DerivedParameter>(statistics::DerivedParameter("bias"));
  parameters.push_back(m_bias);

  // set the model parameters
  set_parameters(parameters);

  // input data used to construct the model
  auto inputs = make_shared<STR_data_model>(m_data_model);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xi0_linear_cosmology_clusters_selection_function, inputs));
}
