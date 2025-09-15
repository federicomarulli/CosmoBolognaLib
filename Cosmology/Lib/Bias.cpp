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
 *  @brief Methods of the class Bias, used to model the halo bias
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  Bias used to model the halo bias function
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Bias.h"
#include "MassFunction.h"

using namespace std;
using namespace cbl;


// =====================================================================================


double cbl::cosmology::Bias::bias_halo (const double Mass, const double redshift, const std::string author, const std::string method_SS, const bool store_output, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  Sigma SG(m_cosmology);
  
  const double SSS = SG.sigma2M({Mass}, method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0]; 
  const double Sigma = sqrt(SSS); 

  double bias;
  if (author=="Castro") {
    const double RHO = m_cosmology->rho_m(0., true);
    double Radius = pow(3*Mass*0.25/cbl::par::pi/RHO, 1./3);
    double DlnSigma_DlnR = SG.dnsigma2R(1, {Radius}, method_SS, 0)[0]*(Radius / (2. * SSS));
    bias = m_bias_halo_generator(Sigma, redshift, author, Delta, DlnSigma_DlnR); 
  }
  else
    bias = m_bias_halo_generator(Sigma, redshift, author, Delta); 
  
  if (m_cosmology->fNL()!=0) {
    NG ng(m_cosmology);
    bias += ng.bias_correction(kk, Mass, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*SSS*pow(bias-1, 2);
  }
  
  return bias; 
}


// =====================================================================================


double cbl::cosmology::Bias::bias_halo (const double Mass, const double Sigma, const double redshift, const std::string model_bias, const bool store_output, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double bias;
  
  if (model_bias=="Castro") {
    cbl::cosmology::Sigma SG(m_cosmology);
    const double RHO = m_cosmology->rho_m(0., true);
    double Radius = pow(3*Mass*0.25/cbl::par::pi/RHO, 1./3);
    double DlnSigma_DlnR = SG.dnsigma2R(1, {Radius}, method_SS, 0)[0]*(Radius / (2. * Sigma*Sigma));
    bias = m_bias_halo_generator(Sigma, redshift, model_bias, Delta, DlnSigma_DlnR); 
  }
  else
    bias = m_bias_halo_generator(Sigma, redshift, model_bias, Delta); 

  if (m_cosmology->fNL()!=0) {
    NG ng(m_cosmology);
    bias += ng.bias_correction(kk, Mass, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*Sigma*pow(bias-1, 2); // check!!!
  }
  
  return bias;
}


// =====================================================================================


double cbl::cosmology::Bias::bias_halo (const double Mass, const double Sigma, const double redshift, const double DN, const std::string model_bias, const bool store_output, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double bias; 
  if (model_bias=="Castro") {
    cbl::cosmology::Sigma SG(m_cosmology);
    const double RHO = m_cosmology->rho_m(0., true);
    double Radius = pow(3*Mass*0.25/cbl::par::pi/RHO, 1./3);
    double DlnSigma_DlnR = SG.dnsigma2R(1, {Radius}, method_SS, 0)[0]*(Radius / (2. * Sigma*Sigma));
    bias = m_bias_halo_generator(Sigma, redshift, DN, model_bias, Delta, DlnSigma_DlnR); 
  }
  else
    bias = m_bias_halo_generator(Sigma, redshift, DN, model_bias, Delta); 

  if (m_cosmology->fNL()!=0) {
    NG ng(m_cosmology);
    bias += ng.bias_correction(kk, Mass, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*Sigma*pow(bias-1, 2); // check!!!
  }
  
  return bias; 
}


// =====================================================================================


double cbl::cosmology::Bias::m_bias_halo_generator (const double Sigma, const double redshift, const std::string author, const double Delta, const double dlnsigma_dlnR) const
{
  const double D_N = m_cosmology->DN(redshift);
  
  return m_bias_halo_generator(Sigma, redshift, D_N, author, Delta, dlnsigma_dlnR);
}


// =====================================================================================


double cbl::cosmology::Bias::m_bias_halo_generator (const double Sigma, const double redshift, const double D_N, const std::string author, const double Delta, const double dlnsigma_dlnR) const
{
  const double deltacz = m_cosmology->deltac(redshift);
  const double sigmaz = Sigma*D_N;
  
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

  else if (author=="Castro") { // Castro et al. (2023) 
    
    WarningMsgCBL("Check! The multiplicity function, used to compute the bias, differs for more than 1% from the original one at z>0.5 and M>10^14 solar masses. The difference is due to the integrand of the function sigma2M", "m_MF_generator", "MassFunction.cpp");
    
    // values for ROCKSTAR
    const double a1 = 0.7962217;
    const double a2 = 0.1448873;
    const double az = -0.06577418;
    const double p1 = -0.5612250;
    const double p2 = -0.4742516;    
    const double q1 = 0.3688495;
    const double q2 = -0.2803717;
    const double qz = 0.02514857;

    const double ar = a1 + a2*pow((dlnsigma_dlnR+0.6125) , 2);
    const double qr = q1 + q2*(dlnsigma_dlnR+0.5);
    const double pp = p1 + p2*(dlnsigma_dlnR+0.5);

    const double aa = ar*pow(m_cosmology->Omega_matter(redshift), az);
    const double qq = qr*pow(m_cosmology->Omega_matter(redshift), qz);

    const double AA = pow(pow(2, -0.5-pp+qq*0.5)/sqrt(par::pi)*(pow(2, pp)*std::tgamma(qq/2)+std::tgamma(-pp+qq*0.5)), -1);

    const double nu = deltacz / sigmaz;

    auto ln_nu_fnu = [AA, aa, pp, qq](double nu) -> double {
      const double BB = sqrt(2*aa*nu*nu / par::pi) * exp(-nu * nu * 0.5 * aa) * (1 + 1. / pow((nu * nu * aa), pp)) * pow(nu * sqrt(aa), qq-1.);
      return log(AA * BB);
    };
    
    double dln_nu_fnu_dnu = cbl::wrapper::gsl::GSL_derivative(ln_nu_fnu, nu, 1.e-1, 1.e-2);  
	     
    double bias_pbs =   1. - 1. / deltacz * nu* dln_nu_fnu_dnu;

    double S8 =   m_cosmology->sigma8(redshift)*sqrt(m_cosmology->Omega_matter(redshift)/0.3);
    double A0 = 1.150;
    double f0 = 1 + 0.0929 * m_cosmology->Omega_matter(redshift);
    double f1 = 1 + 0.256 * dlnsigma_dlnR + 0.173 * pow(dlnsigma_dlnR, 2);
    double f2 = 1-0.0372 * S8; 
    
    bias = bias_pbs*A0*f0*f1*f2; 
  }
 
  
  else
    ErrorCBL("author = " + author + "!", "m_bias_halo_generator", "Bias.cpp");
  
  return bias;
}


// =====================================================================================


double cbl::cosmology::Bias::bias_eff (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  // ---------- create/read the grid file with sigma(M) and its derivative ---------- 

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
  
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

    const double MF = m_mass_function.mass_function(mass[mm], sigma[mm], dlnsigma[mm], redshift, model_MF, store_output, output_root, Delta, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    
    Bias_eff += bias_halo(mass[mm], sigma[mm], redshift, model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*MF*(mass[mm+1]-mass[mm]);

    Norm += MF*(mass[mm+1]-mass[mm]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


double cbl::cosmology::Bias::bias_eff (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  // ---------- create/read the grid file with sigma(M) and its derivative ---------- 

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
  
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

    Bias_eff += bias_halo(MM[k], sig, redshift, model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*mf*(MM[k+1]-MM[k]);
    Norm += mf*(MM[k+1]-MM[k]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_eff_mass_grid (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  // ---------- create/read the grid file with sigma(M) and its derivative ---------- 

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

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
      bias[k] = bias_halo(MM[k], interpolated(MM[k], mass, sigma, "Linear"), zz, model_bias, store_output, output_root, interpType, Delta_crit/m_cosmology->Omega_matter(zz), kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    }
    
    return {Average(bias), cbl::Sigma(bias)/sqrt(MM.size())};
  }

  else {
    vector<double> bias2(MM.size());
    
    for (size_t k=0; k<MM.size(); ++k) {
      const double z1 = (redshift.size()>1) ? redshift[k] : redshift[0];
      for (size_t l=k+1; l<MM.size(); ++l) {
	const double z2 = (redshift.size()>1) ? redshift[l] : redshift[0];
	bias2[k] = bias_halo(MM[k], interpolated(MM[k], mass, sigma, "Linear"), z1, model_bias, store_output, output_root, interpType, Delta_crit/m_cosmology->Omega_matter(z1), kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*bias_halo(MM[l], interpolated(MM[l], mass, sigma, "Linear"), z2, model_bias, store_output, output_root, interpType, Delta_crit/m_cosmology->Omega_matter(z2), kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      }
    }
    
    return {sqrt(Average(bias2)), sqrt(cbl::Sigma(bias2)/sqrt(MM.size()))};
  }
}


// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_eff_mass (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  if (meanType!="mean_bias" && meanType!="pair_mean_bias")
    ErrorCBL("the chosen meanType is not allowed!", "bias_eff_mass", "Bias.cpp");

  Sigma SG(m_cosmology);
  
  vector<double> sigma2M = SG.sigma2M(MM, method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file, true);
  
  if (meanType=="mean_bias") {

    vector<double> bias(MM.size());

#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for schedule(static, 2)
      for (size_t k=0; k<MM.size(); k++) {
        const double sigma = sqrt(sigma2M[k]);
        bias[k] = bias_halo(MM[k], sigma, (redshift.size() > 1) ? redshift[k] : redshift[0], model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      }
    }

    return {Average(bias), cbl::Sigma(bias)/sqrt(MM.size())};
  }

  else {
    
    vector<double> bias2(MM.size());

#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for schedule(static, 2)
      for (size_t k=0; k<MM.size(); ++k) {
        const double z1 = (redshift.size() > 1) ? redshift[k] : redshift[0];
        const double sigma1 = sqrt(sigma2M[k]);
        for (size_t l=k+1; l<MM.size(); ++l) {
          const double z2 = (redshift.size() > 1) ? redshift[l] : redshift[0];
          const double sigma2 = sqrt(sigma2M[l]);
          bias2[k] = bias_halo(MM[k], sigma1, z1, model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file) * bias_halo(MM[l], sigma2, z2, model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
        }
      }
    }
    
    return {sqrt(Average(bias2)), sqrt(cbl::Sigma(bias2)/sqrt(MM.size()))};
  }
}

// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_eff_mass (const std::vector<double> mass,  const std::vector<double> mass_grid,  const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  if (meanType != "mean_bias" && meanType != "pair_mean_bias")
    ErrorCBL("the chosen meanType is not allowed!", "bias_eff_mass", "Bias.cpp");

  Sigma SG(m_cosmology);
 
  vector<double> Sigma2M = SG.sigma2M(mass_grid, method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file, true);
  vector<double> Sigma(mass_grid.size());

  for (size_t k = 0; k < mass_grid.size(); k++)
    Sigma[k] = sqrt(Sigma2M[k]);

  glob::FuncGrid sigma_interp(mass_grid, Sigma, "Spline");

  if (meanType=="mean_bias") {
    vector<double> bias(mass.size());
    for (size_t k = 0; k < mass.size(); k++)
      bias[k] = bias_halo(mass[k], sigma_interp(mass[k]), (redshift.size() > 1) ? redshift[k] : redshift[0], model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    return {Average(bias), cbl::Sigma(bias) / sqrt(mass.size())};
  }

  else {
    vector<double> bias2(mass.size());
    for (size_t k=0; k<mass.size(); k++) {
      const double z1 = (redshift.size() > 1) ? redshift[k] : redshift[0];
      for (size_t l=k+1; l<mass.size(); ++l) {
        const double z2 = (redshift.size() > 1) ? redshift[l] : redshift[0];
        bias2[k] = bias_halo(mass[k], sigma_interp(mass[k]), z1, model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file) * bias_halo(mass[l], sigma_interp(mass[l]), z2, model_bias, store_output, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      }
    }

    return {sqrt(Average(bias2)), sqrt(cbl::Sigma(bias2) / sqrt(mass.size()))};
  }
}


// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DlnSigma, const glob::FuncGrid interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha, const bool store_output, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{ 
  vector<double> Bias_eff(redshift.size(), 0.);
  
  vector<double> mass_vec = logarithmic_bin_vector(50, Mass_min, Mass_max);
  
  for (size_t i=0; i<redshift.size(); ++i) {

    const double DD = Delta_crit/m_cosmology->Omega_matter(redshift[i]);
    vector<double> sigma = interp_sigma.eval_func(mass_vec);
    vector<double> dlnsigma = interp_DlnSigma.eval_func(mass_vec);
    vector<double> mass_func = m_mass_function.mass_function(mass_vec, sigma, dlnsigma, redshift[i], model_MF, store_output, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    vector<double> bias_func = bias_halo(mass_vec, sigma, redshift[i], model_bias, store_output, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);

    cbl::glob::FuncGrid interp_MF(mass_vec, mass_func, "Spline");
    cbl::glob::FuncGrid interp_BH(mass_vec, bias_func, "Spline");

    auto integrand_num = [&] (const double lg_mass) {
      const double mass = exp(lg_mass);
      const double SF = interp_SF(mass/alpha);
      const double BH = interp_BH(mass);
      const double MF = interp_MF(mass);
	
      return SF*BH*MF*mass;
    };

    auto integrand_denom = [&] (const double lg_mass) {
      const double mass = exp(lg_mass);
      const double SF = interp_SF(mass/alpha); 
      const double MF = interp_MF(mass);

      return SF*MF*mass;
    };

    Bias_eff[i] = wrapper::gsl::GSL_integrate_qag(integrand_num, log(Mass_min), log(Mass_max))/wrapper::gsl::GSL_integrate_qag(integrand_denom, log(Mass_min), log(Mass_max));
  }


  return Bias_eff;
}


// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DlnSigma, const glob::FuncGrid2D interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha, const bool store_output, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{ 
  vector<double> Bias_eff(redshift.size(), 0.);

  for (size_t i=0; i<redshift.size(); ++i) {

    const double DD = Delta_crit/m_cosmology->Omega_matter(redshift[i]);

    auto integrand_num = [&] (const double mass) {
      const double sigma = interp_sigma(mass);
      const double dlnsigma = interp_DlnSigma(mass);
      const double SF = interp_SF(mass/alpha, redshift[i]);
      const double BH = bias_halo(mass, sigma, redshift[i], model_bias, store_output, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      const double MF = m_mass_function.mass_function(mass, sigma, dlnsigma, redshift[i], model_MF, store_output, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
	
      return SF*BH*MF;
    };

    auto integrand_denom = [&] (const double mass) {
      const double sigma = interp_sigma(mass);
      const double dlnsigma = interp_DlnSigma(mass);
      const double SF = interp_SF(mass/alpha, redshift[i]); 
      const double MF = m_mass_function.mass_function(mass, sigma, dlnsigma, redshift[i], model_MF, store_output, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);

      return SF*MF;
    };

    Bias_eff[i] = wrapper::gsl::GSL_integrate_qag(integrand_num, Mass_min, Mass_max)/wrapper::gsl::GSL_integrate_qag(integrand_denom, Mass_min, Mass_max);
  }

  return Bias_eff;
}


// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_eff_selection_function (const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column, const double alpha, const bool store_output, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{  
  // ---------- create/read the grid file with sigmaM, dlnsigmaM ---------- 

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

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
  
  return bias_eff_selection_function(interp_sigma, interp_DlnSigma, interp_SF, Mass_min, Mass_max, redshift, model_bias, model_MF, method_SS, alpha, store_output, output_root, Delta_crit, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);

}


// =====================================================================================


void cbl::cosmology::Bias::generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const std::string cosmoPar, const double min_par, const double max_par, const int nbin_par, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const std::shared_ptr<cosmology::Cosmology> cosmology_mass, const std::vector<double> redshift_source)
{
  const double defaultValue = m_cosmology->value(cosmoPar);
  const string file = dir_output+file_bias_eff_grid;
 
  ifstream fin(file.c_str());
  
  if (!fin) {
    
    vector<double> pp = linear_bin_vector(nbin_par, min_par, max_par);
   
    ofstream fout(file.c_str()); checkIO(fout, file);

    for (int i=0; i<nbin_par; i++) {
      m_cosmology->set_parameter(cosmoPar, pp[i]);
      
      (void)cosmology_mass;
      (void)redshift_source;
      /*
      // convert the masses in the new cosmology, if they were computed in a different cosmology
      vector<double> _mass(mass.size());
      for (size_t mm=0; mm<mass.size(); mm++) 
      _mass[mm] = m_mass_function.converted_mass(mass[mm], cosmology_mass, redshift[mm], (redshift_source.size()==redshift.size()) ? redshift_source[mm] : 0.);
      */
      
      fout << pp[i] << "  " << bias_eff_mass(/*_*/mass, mass_grid, redshift, model_bias, method_SS, meanType, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0] << endl;
      
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
  
  m_cosmology->set_parameter(cosmoPar, defaultValue);
}


// =====================================================================================


void cbl::cosmology::Bias::generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const std::string cosmoPar, const double min_par, const double max_par, const int nbin_par, const double redshift, const double Mass_min, const double Mass_max, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column, const double alpha, const bool store_output, const std::string output_root, const double Delta_crit, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{
  const double defaultValue = m_cosmology->value(cosmoPar);
  const string file = dir_output+file_bias_eff_grid;

  ifstream fin(file.c_str());
  
  if (!fin) {
    const vector<double> pp = linear_bin_vector(nbin_par, min_par, max_par);

    ofstream fout(file.c_str());

    if ((int)pp.size()<nbin_par) ErrorCBL(conv(pp.size(), par::fINT)+" < nbin_par!", "generate_bias_eff_grid_one_cosmopar", "Bias.cpp");

    for (int i=0; i<nbin_par; i++) {
      m_cosmology->set_parameter(cosmoPar, pp[i]);

      fout << pp[i] << "  " << bias_eff_selection_function(Mass_min, Mass_max, {redshift}, model_bias, model_MF, method_SS, selection_function_file, column, alpha, store_output, output_root, Delta_crit, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0] << endl;

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
  
  m_cosmology->set_parameter(cosmoPar, defaultValue);
}


// =====================================================================================


void cbl::cosmology::Bias::generate_bias_eff_grid_two_cosmopars (vector<double> &parameter1, vector<double> &parameter2, vector<vector<double>> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const std::string cosmoPar1, const double min_par1, const double max_par1, const int nbin_par1, const std::string cosmoPar2, const double min_par2, const double max_par2, const int nbin_par2, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const std::shared_ptr<cosmology::Cosmology> cosmology_mass, const std::vector<double> redshift_source)
{ 
  double defaultValue1 = m_cosmology->value(cosmoPar1);
  double defaultValue2 = m_cosmology->value(cosmoPar2);
  
  string file = dir_output+file_bias_eff_grid;

  ifstream fin(file.c_str());

  if (!fin) {
    vector<double> pp1 = linear_bin_vector(nbin_par1, min_par1, max_par1);
    vector<double> pp2 = linear_bin_vector(nbin_par2, min_par2, max_par2);

    ofstream fout(file.c_str()); checkIO(fout, file);

    for (int i=0; i<nbin_par1; i++) {
      m_cosmology->set_parameter(cosmoPar1, pp1[i]);
      for (int j=0; j<nbin_par2; j++) {
	m_cosmology->set_parameter(cosmoPar2, pp2[j]);

	// convert the masses in the new cosmology, if they were computed in a different cosmology
	vector<double> _mass(mass.size());
	for (size_t mm=0; mm<mass.size(); mm++) 
	  _mass[mm] = m_mass_function.converted_mass(mass[mm], cosmology_mass, redshift[mm], (redshift_source.size()==redshift.size()) ? redshift_source[mm] : 0.);
	
	const double bias = bias_eff_mass(mass, mass_grid, redshift, model_bias, method_SS, meanType, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0];
	
	fout << pp1[i] << "  " << pp2[j] <<  "  " << bias << endl;
	coutCBL << "parameter1 = " << pp1[i] << ",  parameter2 = " << pp2[j] <<  ", bias = " << bias << endl;
      }
      fout << endl;
    }
    fout.clear(); fout.close();
  }

  read_matrix(file, parameter1, parameter2, bias_eff);
  
  m_cosmology->set_parameter(cosmoPar1, defaultValue1);
  m_cosmology->set_parameter(cosmoPar2, defaultValue2);
}


// =====================================================================================


vector<double> cbl::cosmology::Bias::bias_halo (const std::vector<double> Mass, const std::vector<double> Sigma, const double redshift, const std::string model_bias, const bool store_output, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double D_N = m_cosmology->DN(redshift);
  vector<double> bias(Mass.size()), Radius(Mass.size()), DlnSigma_DlnR(Mass.size());

  if (model_bias=="Castro") {
    cbl::cosmology::Sigma SG(m_cosmology);
    const double RHO = m_cosmology->rho_m(0., true);
    for (size_t i=0; i<Mass.size(); i++) 
      Radius[i] = pow(3*Mass[i]*0.25/cbl::par::pi/RHO, 1./3);
    
    vector<double> DlnSigma_DlnR = SG.dnsigma2R(1, Radius, method_SS, 0); 

    for (size_t i=0; i<Mass.size(); i++) {
      DlnSigma_DlnR[i]*=(Radius[i] / (2.*Sigma[i]*Sigma[i]));
      bias[i] = m_bias_halo_generator(Sigma[i], redshift, model_bias, Delta, DlnSigma_DlnR[i]); 
    }
  }
  else
    for (size_t i=0; i<Mass.size(); i++)
      bias[i] = m_bias_halo_generator(Sigma[i], redshift, D_N, model_bias, Delta); 

  if (m_cosmology->fNL()!=0) {
    NG ng(m_cosmology);
    for (size_t i=0; i<Mass.size(); i++)
      bias[i] += ng.bias_correction(kk, Mass[i], method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*Sigma[i]*pow(bias[i]-1, 2); // check!!!
  }
  
  return bias; 
}


