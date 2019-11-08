/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Cosmology/Lib/Bias.cpp
 *
 *  @brief Methods of the class Cosmology used to model the bias
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the bias
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Cosmology::bias_halo (const double Mass, const double redshift, const std::string author, const std::string method_SS, const bool store_output_CAMB, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  const double SSS = sigma2M(Mass, method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file); 
  const double Sigma = sqrt(SSS); 
  
  double bias = m_bias_halo_generator(Sigma, redshift, author, Delta); 
  
  if (m_fNL!=0) 
    bias += bias_correction(kk, Mass, method_SS, store_output_CAMB, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*SSS*pow(bias-1, 2);

  return bias; 
}


// =====================================================================================


double cbl::cosmology::Cosmology::bias_halo (const double Mass, const double Sigma, const double redshift, const std::string model_bias, const bool store_output_CAMB, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double bias = m_bias_halo_generator(Sigma, redshift, model_bias, Delta); 

  if (m_fNL!=0)
    bias += bias_correction(kk, Mass, method_SS, store_output_CAMB, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*Sigma*pow(bias-1, 2); // check!!!

  return bias; 
}


// =====================================================================================


double cbl::cosmology::Cosmology::m_bias_halo_generator (const double Sigma, const double redshift, const std::string author, const double Delta) const
{
  const double deltacz = deltac(redshift);
  const double sigmaz = Sigma*DD(redshift)/DD(0.);
  
  double bias = -1000.;

  if (author=="ST99") {
    double aa = 0.707;
    double pp = 0.3;
    double ni = pow(deltacz/sigmaz, 2); 
    bias = 1.+(aa*ni-1.)/deltacz+(2.*pp/deltacz)/(1.+pow(aa*ni,pp));
  }

  else if (author=="SMT01") {
    double aa = 0.707;
    double bb = 0.5;
    double cc = 0.6;
    double ni = deltacz/sigmaz; 
    bias = 1.+1./(sqrt(aa)*deltacz)*(sqrt(aa)*aa*pow(ni,2.)+sqrt(aa)*bb*pow(aa*pow(ni,2.),1.-cc)-pow(aa*pow(ni,2.),cc)/(pow(aa*pow(ni,2.),cc)+bb*(1.-cc)*(1.-cc*0.5)));
  }
  
  else if (author=="SMT01_WL04") {
    double aa = 0.707;
    double bb = 0.5;
    double cc = 0.6;
    double ni = deltacz/sigmaz; 
    double niI = sqrt(aa)*ni;
    bias = 1.+1./deltacz*(pow(niI,2.)+bb*pow(niI,2.*(1.-cc))-pow(niI,2.*cc)/sqrt(aa)/(pow(niI,2.*cc)+bb*(1.-cc)*(1.-cc*0.5)));
  }
  
  else if (author=="Tinker") { // Tinker et al. (2010)
    double yy = log10(Delta);
    double AA = 1.+0.24*yy*exp(-pow(4./yy,4));
    double aa = 0.44*yy-0.88;
    double BB = 0.183;
    double bb = 1.5;
    double CC = 0.019+0.107*yy+0.19*exp(-pow(4./yy,4));
    double ccc = 2.4;
    double ni = 1.686/sigmaz;
    bias = 1.-AA*pow(ni,aa)/(pow(ni,aa)+pow(1.686,aa))+BB*pow(ni,bb)+CC*pow(ni,ccc);
  }
  
  else
    ErrorCBL("author = " + author + "!", "m_bias_halo_generator", "Bias.cpp");
  
  return bias;
}


// =====================================================================================


double cbl::cosmology::Cosmology::bias_eff (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output_CAMB, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  // ---------- create/read the grid file with sigma(M) and its derivative ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file);
  
  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  
  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Mass_min<Mass && Mass<Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }
  
  if (mass.size()==0)
    ErrorCBL("mass.size()=0!", "bias_eff", "Bias.cpp");
  

  // ---------- compute the effective bias ---------- 
  
  double Bias_eff = 0., Norm = 0.;

  for (size_t mm=0; mm<mass.size()-1; mm++) {

    const double MF = mass_function(mass[mm], sigma[mm], dlnsigma[mm], redshift, model_MF, store_output_CAMB, output_root, Delta, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    
    Bias_eff += bias_halo(mass[mm], sigma[mm], redshift, model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*MF*(mass[mm+1]-mass[mm]);

    Norm += MF*(mass[mm+1]-mass[mm]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


double cbl::cosmology::Cosmology::bias_eff (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output_CAMB, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  // ---------- create/read the grid file with sigma(M) and its derivative ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Min(MM)<Mass && Mass<Max(MM)) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
    }
  }
  
  if (mass.size()==0) 
    ErrorCBL("mass.size()=0, Min(MM) = " + conv(Min(MM),par::fDP3) + ", Max(MM) = " + conv(Max(MM),par::fDP3) + ", file_grid = " + file_grid, "bias_eff", "Bias.cpp");
  

  // ---------- compute the effective bias ---------- 
  
  double Bias_eff = 0., Norm = 0.;
  double mf, sig, err = -1;
  
  for (size_t k=0; k<MM.size()-1; k++) {
    mf = MF[k];
    sig = interpolated(MM[k], mass, sigma, "Linear");
    
    if (err/sig>0.1)
      ErrorCBL("err/sig = " + conv(err/sig, par::fDP3) + "!", "bias_eff", "Bias.cpp");

    Bias_eff += bias_halo(MM[k], sig, redshift, model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*mf*(MM[k+1]-MM[k]);
    Norm += mf*(MM[k+1]-MM[k]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_eff_mass_grid (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output_CAMB, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  // ---------- create/read the grid file with sigma(M) and its derivative ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Min(MM)<Mass && Mass<Max(MM)) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
    }
  }
  
  if (mass.size()==0) 
    ErrorCBL("mass.size()=0, Min(MM) = " + conv(Min(MM),par::fDP3) + ", Max(MM) = " + conv(Max(MM),par::fDP3) + ", file_grid = " + file_grid, "bias_eff_mass_grid", "Bias.cpp");
  

  // ---------- compute the effective bias ---------- 

  if (meanType!="mean_bias" && meanType!="pair_mean_bias")
    ErrorCBL("the chosen meanType is not allowed!", "bias_eff_mass_grid", "Bias.cpp");
  
  if (meanType=="mean_bias") {
    vector<double> bias(MM.size());
    
    for (size_t k=0; k<MM.size(); k++) {
      const double zz = (redshift.size()>1) ? redshift[k] : redshift[0];
      bias[k] = bias_halo(MM[k], interpolated(MM[k], mass, sigma, "Linear"), zz, model_bias, store_output_CAMB, output_root, interpType, Delta_vir(Delta_crit, zz), kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    }
    
    return {Average(bias), cbl::Sigma(bias)/sqrt(MM.size())};
  }

  else {
    vector<double> bias2(MM.size());
    
    for (size_t k=0; k<MM.size(); ++k) {
      const double z1 = (redshift.size()>1) ? redshift[k] : redshift[0];
      for (size_t l=k+1; l<MM.size(); ++l) {
	const double z2 = (redshift.size()>1) ? redshift[l] : redshift[0];
	bias2[k] = bias_halo(MM[k], interpolated(MM[k], mass, sigma, "Linear"), z1, model_bias, store_output_CAMB, output_root, interpType, Delta_vir(Delta_crit, z1), kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*bias_halo(MM[l], interpolated(MM[l], mass, sigma, "Linear"), z2, model_bias, store_output_CAMB, output_root, interpType, Delta_vir(Delta_crit, z2), kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      }
    }
    
    return {sqrt(Average(bias2)), sqrt(cbl::Sigma(bias2)/sqrt(MM.size()))};
  }
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_eff_mass (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output_CAMB, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  if (meanType!="mean_bias" && meanType!="pair_mean_bias")
    ErrorCBL("the chosen meanType is not allowed!", "bias_eff_mass", "Bias.cpp");
  
  if (meanType=="mean_bias") {

    vector<double> bias(MM.size());

#pragma omp parallel num_threads(omp_get_max_threads())
    {
      
#pragma omp for schedule(static, 2)
      for (size_t k=0; k<MM.size(); k++) {
	const double sigma = sqrt(sigma2M(MM[k], method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file, true));
	bias[k] = bias_halo(MM[k], sigma, (redshift.size()>1) ? redshift[k] : redshift[0], model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      }

    }
    
    return {Average(bias), cbl::Sigma(bias)/sqrt(MM.size())};
  }

  else {
    vector<double> bias2(MM.size());

#pragma omp parallel num_threads(omp_get_max_threads())
    {
      
#pragma omp for schedule(static, 2)
      for (size_t k=0; k<MM.size(); ++k)  {
	const double z1 = (redshift.size()>1) ? redshift[k] : redshift[0];
	const double sigma1 = sqrt(sigma2M(MM[k], method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file, true));
	for (size_t l=k+1; l<MM.size(); ++l) {
	  const double z2 = (redshift.size()>1) ? redshift[l] : redshift[0];
	  const double sigma2 = sqrt(sigma2M(MM[l], method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file, true));
	  bias2[k] = bias_halo(MM[k], sigma1, z1, model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*bias_halo(MM[l], sigma2, z2, model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
	}
      }

    }
    
    return {sqrt(Average(bias2)), sqrt(cbl::Sigma(bias2)/sqrt(MM.size()))};
  }

}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_eff_mass (const std::vector<double> mass,  const std::vector<double> mass_grid,  const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output_CAMB, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  if (meanType!="mean_bias" && meanType!="pair_mean_bias")
    ErrorCBL("the chosen meanType is not allowed!", "bias_eff_mass", "Bias.cpp");
  
  vector<double> Sigma;
  
  for (size_t k=0; k<mass_grid.size(); k++)
    Sigma.emplace_back(sqrt(sigma2M(mass_grid[k], method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file, true)));
  glob::FuncGrid sigma_interp(mass_grid, Sigma, "Spline");

  if (meanType=="mean_bias") {
    vector<double> bias(mass.size());
    for (size_t k=0; k<mass.size(); k++) 
      bias[k] = bias_halo(mass[k], sigma_interp(mass[k]), (redshift.size()>1) ? redshift[k] : redshift[0], model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    return {Average(bias), cbl::Sigma(bias)/sqrt(mass.size())};
  }

  else {
    vector<double> bias2(mass.size());
    for (size_t k=0; k<mass.size(); k++) {
      const double z1 = (redshift.size()>1) ? redshift[k] : redshift[0];
      for (size_t l=k+1; l<mass.size(); ++l) {
	const double z2 = (redshift.size()>1) ? redshift[l] : redshift[0];
	bias2[k] = bias_halo(mass[k], sigma_interp(mass[k]), z1, model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*bias_halo(mass[l], sigma_interp(mass[l]), z2, model_bias, store_output_CAMB, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      }
    }
    
    return {sqrt(Average(bias2)), sqrt(cbl::Sigma(bias2)/sqrt(mass.size()))};
  }
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DlnSigma, const glob::FuncGrid interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha, const bool store_output_CAMB, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{ 
  vector<double> Bias_eff(redshift.size(), 0.);

  for (size_t i=0; i<redshift.size(); ++i) {

    auto integrand_num = [&] (const double mass)
      {
	const double sigma = interp_sigma(mass);
	const double dlnsigma = interp_DlnSigma(mass);
	
	const double SF = interp_SF(mass/alpha);
	
	const double DD = Delta_vir(Delta_crit, redshift[i]);

	const double BH = bias_halo(mass, sigma, redshift[i], model_bias, store_output_CAMB, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
	
	const double MF = mass_function(mass, sigma, dlnsigma, redshift[i], model_MF, store_output_CAMB, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
	
	return SF*BH*MF;
      };

    auto integrand_denom = [&] (const double mass)
      {
	const double sigma = interp_sigma(mass);
	const double dlnsigma = interp_DlnSigma(mass);
	
	const double SF = interp_SF(mass/alpha); 
	
	const double DD = Delta_vir(Delta_crit, redshift[i]);
	
	const double MF = mass_function(mass, sigma, dlnsigma, redshift[i], model_MF, store_output_CAMB, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);

	return SF*MF;
      };

    Bias_eff[i] = wrapper::gsl::GSL_integrate_qag(integrand_num, Mass_min, Mass_max)/wrapper::gsl::GSL_integrate_qag(integrand_denom, Mass_min, Mass_max);
  }


  return Bias_eff;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DlnSigma, const glob::FuncGrid2D interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha, const bool store_output_CAMB, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{ 
  vector<double> Bias_eff(redshift.size(), 0.);

  for (size_t i=0; i<redshift.size(); ++i) {

    auto integrand_num = [&] (const double mass)
      {
	const double sigma = interp_sigma(mass);
	const double dlnsigma = interp_DlnSigma(mass);
	
	const double SF = interp_SF(mass/alpha, redshift[i]);
	
	const double DD = Delta_vir(Delta_crit, redshift[i]);

	const double BH = bias_halo(mass, sigma, redshift[i], model_bias, store_output_CAMB, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
	
	const double MF = mass_function(mass, sigma, dlnsigma, redshift[i], model_MF, store_output_CAMB, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
	
	return SF*BH*MF;
      };

    auto integrand_denom = [&] (const double mass)
      {
	const double sigma = interp_sigma(mass);
	const double dlnsigma = interp_DlnSigma(mass);
	
	const double SF = interp_SF(mass/alpha, redshift[i]); 
	
	const double DD = Delta_vir(Delta_crit, redshift[i]);
	
	const double MF = mass_function(mass, sigma, dlnsigma, redshift[i], model_MF, store_output_CAMB, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);

	return SF*MF;
      };

    Bias_eff[i] = wrapper::gsl::GSL_integrate_qag(integrand_num, Mass_min, Mass_max)/wrapper::gsl::GSL_integrate_qag(integrand_denom, Mass_min, Mass_max);
  }


  return Bias_eff;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_eff_selection_function (const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column, const double alpha, const bool store_output_CAMB, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{  
  // ---------- create/read the grid file with sigmaM, dlnsigmaM ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., store_output_CAMB, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;

  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Mass_min<Mass && Mass<Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }
  
  if (mass.size()==0)
    ErrorCBL("mass.size()=0, Mass_min = " + conv(Mass_min,par::fDP3) + ", Mass_max = " + conv(Mass_max,par::fDP3) + ", file_grid = " + file_grid, "bias_eff_selection_function", "Bias.cpp");
  
  const glob::FuncGrid interp_sigma(mass, sigma, "Spline");
  const glob::FuncGrid interp_DlnSigma(mass, dlnsigma, "Spline");
  
  
  // ---------- read the selection function ----------
  
  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, column);

  const glob::FuncGrid2D interp_SF(mass_SF, redshift_SF, selection_function, "Linear");

 
  // ---------- compute the effective bias the given redshifts ---------- 
  
  return bias_eff_selection_function(interp_sigma, interp_DlnSigma, interp_SF, Mass_min, Mass_max, redshift, model_bias, model_MF, method_SS, alpha, store_output_CAMB, output_root, Delta_crit, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);

}


// =====================================================================================


void cbl::cosmology::Cosmology::generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const cbl::cosmology::CosmologicalParameter cosmoPar, const double min_par, const double max_par, const int nbin_par, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output_CAMB, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const cbl::cosmology::Cosmology cosmology_mass, const std::vector<double> redshift_source)
{
  const double defaultValue = value(cosmoPar);
  
  const string file = dir_output+file_bias_eff_grid;
 
  ifstream fin(file.c_str());
  
  if (!fin) {
    
    vector<double> pp = linear_bin_vector(nbin_par, min_par, max_par);
   
    ofstream fout(file.c_str()); checkIO(fout, file);

    for (int i=0; i<nbin_par; i++) {
      set_parameter(cosmoPar, pp[i]);
      
      (void)cosmology_mass;
      (void)redshift_source;
      /*
      // convert the masses in the new cosmology, if they were computed in a different cosmology
      vector<double> _mass(mass.size());
      for (size_t mm=0; mm<mass.size(); mm++) 
	_mass[mm] = converted_mass(mass[mm], cosmology_mass, redshift[mm], (redshift_source.size()==redshift.size()) ? redshift_source[mm] : 0.);
      */
      
      fout << pp[i] << "  " << bias_eff_mass(/*_*/mass, mass_grid, redshift, model_bias, method_SS, meanType, store_output_CAMB, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0] << endl;
      
    }
    fout.clear(); fout.close();
  }
 
  fin.clear(); fin.close();

  fin.open(file.c_str()); checkIO(fin, file);
  
  parameter.erase(parameter.begin(), parameter.end());
  bias_eff.erase(bias_eff.begin(), bias_eff.end());
  
  string line;
  while (getline(fin, line)) {
    stringstream SS(line); double _p, _b;
    SS >> _p >> _b;
    parameter.push_back(_p);
    bias_eff.push_back(_b);
  }
  fin.clear(); fin.close();

  if (parameter.size()<2) ErrorCBL("parameter.size()<2; check the grid file: "+file+"!", "generate_bias_eff_grid_one_cosmopar", "Bias.cpp");
  
  set_parameter(cosmoPar, defaultValue);
}


// =====================================================================================


void cbl::cosmology::Cosmology::generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const cbl::cosmology::CosmologicalParameter cosmoPar, const double min_par, const double max_par, const int nbin_par, const double redshift, const double Mass_min, const double Mass_max, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column, const double alpha, const bool store_output_CAMB, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{
  const double defaultValue = value(cosmoPar);

  const string file = dir_output+file_bias_eff_grid;

  ifstream fin(file.c_str());
  
  if (!fin) {
    const vector<double> pp = linear_bin_vector(nbin_par, min_par, max_par);

    ofstream fout(file.c_str());

    if ((int)pp.size()<nbin_par) ErrorCBL(conv(pp.size(), par::fINT)+" < nbin_par!", "generate_bias_eff_grid_one_cosmopar", "Bias.cpp");

    for (int i=0; i<nbin_par; i++) {
      set_parameter(cosmoPar, pp[i]);

      fout << pp[i] << "  " << bias_eff_selection_function(Mass_min, Mass_max, {redshift}, model_bias, model_MF, method_SS, selection_function_file, column, alpha, store_output_CAMB, output_root, Delta_crit, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0] << endl;

    }
    
    fout.clear(); fout.close();
  }

  fin.clear(); fin.close();
  fin.open(file.c_str());

  parameter.erase(parameter.begin(), parameter.end());
  bias_eff.erase(bias_eff.begin(), bias_eff.end());
  
  string line;
  while (getline(fin, line)) {
    stringstream SS(line); double _p, _b;
    SS >> _p >> _b;

    parameter.push_back(_p);
    bias_eff.push_back(_b);
  }
  fin.clear(); fin.close();
  
  if (parameter.size()<2) ErrorCBL("parameter.size()<2; check the grid file: "+file+"!", "generate_bias_eff_grid_one_cosmopar", "Bias.cpp");
  
  set_parameter(cosmoPar, defaultValue);
}


// =====================================================================================


void cbl::cosmology::Cosmology::generate_bias_eff_grid_two_cosmopars (vector<double> &parameter1, vector<double> &parameter2, vector<vector<double>> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const cbl::cosmology::CosmologicalParameter cosmoPar1, const double min_par1, const double max_par1, const int nbin_par1, const cbl::cosmology::CosmologicalParameter cosmoPar2, const double min_par2, const double max_par2, const int nbin_par2, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output_CAMB, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const cbl::cosmology::Cosmology cosmology_mass, const std::vector<double> redshift_source)
{ 
  double defaultValue1 = value(cosmoPar1);
  double defaultValue2 = value(cosmoPar2);
  
  string file = dir_output+file_bias_eff_grid;

  ifstream fin(file.c_str());

  if (!fin) {
    vector<double> pp1 = linear_bin_vector(nbin_par1, min_par1, max_par1);
    vector<double> pp2 = linear_bin_vector(nbin_par2, min_par2, max_par2);

    ofstream fout(file.c_str()); checkIO(fout, file);

    for (int i=0; i<nbin_par1; i++) {
      set_parameter(cosmoPar1, pp1[i]);
      for (int j=0; j<nbin_par2; j++) {
	set_parameter(cosmoPar2, pp2[j]);

	// convert the masses in the new cosmology, if they were computed in a different cosmology
	vector<double> _mass(mass.size());
	for (size_t mm=0; mm<mass.size(); mm++) 
	  _mass[mm] = converted_mass(mass[mm], cosmology_mass, redshift[mm], (redshift_source.size()==redshift.size()) ? redshift_source[mm] : 0.);
	
	const double bias = bias_eff_mass(/*_*/mass, mass_grid, redshift, model_bias, method_SS, meanType, store_output_CAMB, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0];
	
	fout << pp1[i] << "  " << pp2[j] <<  "  " << bias << endl;
	coutCBL << "parameter1 = " << pp1[i] << ",  parameter2 = " << pp2[j] <<  ", bias = " << bias << endl;
      }
      fout << endl;
    }
    fout.clear(); fout.close();
  }

  read_matrix(file, parameter1, parameter2, bias_eff);
  
  set_parameter(cosmoPar1, defaultValue1);
  set_parameter(cosmoPar2, defaultValue2);
}


