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

using namespace cosmobl;


// =====================================================================================


double cosmobl::cosmology::Cosmology::bias_halo (const double Mass, const double redshift, const string author, const string method_SS, const string output_root, const string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file) 
{
  const double SSS = sigma2M(Mass, method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file); 
  const double Sigma = sqrt(SSS); 
  
  double bias = m_bias_halo_generator(Sigma, redshift, author, Delta); 
  
  if (m_fNL!=0) 
    bias += bias_correction(kk, Mass, method_SS, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file) * SSS * (bias-1) * (bias-1);

  return bias; 
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::bias_halo (const double Mass, const double Sigma, const double redshift, const string author_bias, const string output_root, const string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const string method_SS, const string input_file, const bool is_parameter_file) 
{
  double bias = m_bias_halo_generator(Sigma, redshift, author_bias, Delta); 
  
  if (m_fNL!=0)
    bias += bias_correction(kk, Mass, method_SS, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file) * sigma2M(Mass, method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file) * (bias-1) * (bias-1); 

  return bias; 
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::m_bias_halo_generator (const double Sigma, const double redshift, const string author, const double Delta) const
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
  
  if (author=="SMT01") {
    double aa = 0.707;
    double bb = 0.5;
    double cc = 0.6;
    double ni = deltacz/sigmaz; 
    bias = 1.+1./(sqrt(aa)*deltacz)*(sqrt(aa)*aa*pow(ni,2.)+sqrt(aa)*bb*pow(aa*pow(ni,2.),1.-cc)-pow(aa*pow(ni,2.),cc)/(pow(aa*pow(ni,2.),cc)+bb*(1.-cc)*(1.-cc*0.5)));
  }
  
  if (author=="SMT01_WL04") {
    double aa = 0.707;
    double bb = 0.5;
    double cc = 0.6;
    double ni = deltacz/sigmaz; 
    double niI = sqrt(aa)*ni;
    bias = 1.+1./deltacz*(pow(niI,2.)+bb*pow(niI,2.*(1.-cc))-pow(niI,2.*cc)/sqrt(aa)/(pow(niI,2.*cc)+bb*(1.-cc)*(1.-cc*0.5)));
  }
  
  if (author=="Tinker") { // Tinker et al. (2010)
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
  
  if (bias<-100) {
    string Err = "Error in cosmobl::cosmology::Cosmology::bias_halo of Bias.cpp: author = " + author + "!";
    ErrorCBL(Err);
  }
  
  return bias;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::bias_eff (const double Mass_min, const double Mass_max, const double redshift, const string author_bias, const string author_MF, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file) 
{
  // ---------- read the grid file ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file);
  
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


  // ---------- compute the effective bias ---------- 
  
  double Bias_eff = 0., Norm = 0.;

  for (unsigned int mm=0; mm<mass.size()-1; mm++) {

    double MF = mass_function(mass[mm], sigma[mm], dlnsigma[mm], redshift, author_MF, output_root, Delta, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    
    Bias_eff += bias_halo(mass[mm], sigma[mm], redshift, author_bias, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*MF*(mass[mm+1]-mass[mm]);

    Norm += MF*(mass[mm+1]-mass[mm]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::bias_eff (const vector<double> MM, const vector<double> MF, const double redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file) 
{
  // ---------- read the grid file ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Min(MM)<Mass && Mass<Max(MM)) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
    }
  }
  if (mass.size()==0) { 
    string Err = "Error in cosmobl::cosmology::Cosmology::bias_eff of Bias.cpp: mass.size()=0, Min(MM) = " + conv(Min(MM),par::fDP3) + ", Max(MM) = " + conv(Max(MM),par::fDP3) + ", file_grid = " + file_grid;
    ErrorCBL(Err);
  }
  

  // ---------- compute the effective bias ---------- 
  
  double Bias_eff = 0., Norm = 0.;
  double mf, sig, err = -1;
  
  for (size_t k=0; k<MM.size()-1; k++) {
    mf = MF[k];
    sig = interpolated(MM[k], mass, sigma, "Linear");
    
    if (err/sig>0.1) { 
      string Err = "Error in cosmobl::cosmology::Cosmology::bias_eff of Bias.cpp: err/sig = " + conv(err/sig, par::fDP3) + "!";
      ErrorCBL(Err);
    } 

    Bias_eff += bias_halo(MM[k], sig, redshift, author_bias, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file)*mf*(MM[k+1]-MM[k]);
    Norm += mf*(MM[k+1]-MM[k]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::bias_eff_mass_grid (const vector<double> MM, const vector<double> redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file) 
{
  // ---------- read the grid file ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str()); checkIO(fin, file_grid); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Min(MM)<Mass && Mass<Max(MM)) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
    }
  }
  
  if (mass.size()==0) { 
    string Err = "Error in cosmobl::cosmology::Cosmology::bias_eff of Bias.cpp: mass.size()=0, Min(MM) = " + conv(Min(MM),par::fDP3) + ", Max(MM) = " + conv(Max(MM),par::fDP3) + ", file_grid = " + file_grid;
    ErrorCBL(Err);
  }
  

  // ---------- compute the effective bias ---------- 
  
  vector<double> bias(MM.size());

  for (size_t k=0; k<MM.size(); k++)
    bias[k] = bias_halo(MM[k], interpolated(MM[k], mass, sigma, "Linear"), (redshift.size()>1) ? redshift[k] : redshift[0], author_bias, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    
  return {Average(bias), cosmobl::Sigma(bias)/sqrt(MM.size())};
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::bias_eff_mass (const vector<double> MM, const vector<double> redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file) 
{
  vector<double> bias(MM.size());

  for (size_t k=0; k<MM.size(); k++) {

    const double sigma = sqrt(sigma2M(MM[k], method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file, true));
			      
    bias[k] = bias_halo(MM[k], sigma, (redshift.size()>1) ? redshift[k] : redshift[0], author_bias, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    
  }
			
  return {Average(bias), cosmobl::Sigma(bias)/sqrt(MM.size())};
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::bias_eff_mass (const vector<double> mass,  const vector<double> mass_grid,  const vector<double> redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file) 
{
  vector<double> bias(mass.size()), Sigma;
    
  for (size_t k=0; k<mass_grid.size(); k++) 
    Sigma.emplace_back(sqrt(sigma2M(mass_grid[k], method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file, true)));

  glob::FuncGrid sigma_interp(mass_grid, Sigma, "Spline");

  for (size_t k=0; k<mass.size(); k++) 
    bias[k] = bias_halo(mass[k], sigma_interp(mass[k]), (redshift.size()>1) ? redshift[k] : redshift[0], author_bias, output_root, interpType, Delta, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
			
  return {Average(bias), cosmobl::Sigma(bias)/sqrt(mass.size())};
}



// =====================================================================================


void cosmobl::cosmology::Cosmology::generate_bias_eff_grid_one_cosmopar (vector<double> &parameter, vector<double> &bias_eff, const string dir_output, const string file_bias_eff_grid, const cosmobl::cosmology::CosmoPar cosmoPar, const double min_par, const double max_par, const int nbin_par, const vector<double> mass, const vector<double> mass_grid, const vector<double> redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file)
{
  double defaultValue = value(cosmoPar);

  string file = dir_output+file_bias_eff_grid;

  ifstream fin(file.c_str());

  if(!fin)
  {
    vector<double> pp = linear_bin_vector(nbin_par, min_par, max_par);
    fin.clear(); fin.close();
    ofstream fout(file.c_str());

    for (int i=0; i<nbin_par; i++){
      set_parameter (cosmoPar, pp[i]);
      double bias_eff = bias_eff_mass (mass, mass_grid, redshift, author_bias, method_SS, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0];

      fout << pp[i] << " " <<bias_eff << endl;
    }
    fout.clear(); fout.close();
  }
  fin.clear(); fin.close();

  fin.open(file.c_str());

  parameter.erase(parameter.begin(), parameter.end());
  bias_eff.erase(bias_eff.begin(), bias_eff.end());
  
  string line;
  while(getline(fin, line)){
    stringstream SS(line); double _p, _b;
    SS >> _p >> _b;
    parameter.push_back(_p);
    bias_eff.push_back(_b);
  }
  fin.clear(); fin.close();

  set_parameter(cosmoPar, defaultValue);
}


// =====================================================================================


void cosmobl::cosmology::Cosmology::generate_bias_eff_grid_two_cosmopars (vector<double> &parameter1, vector<double> &parameter2, vector<vector<double>> &bias_eff, const string dir_output, const string file_bias_eff_grid, const cosmobl::cosmology::CosmoPar cosmoPar1, const double min_par1, const double max_par1, const int nbin_par1, const cosmobl::cosmology::CosmoPar cosmoPar2, const double min_par2, const double max_par2, const int nbin_par2, const vector<double> mass, const vector<double> mass_grid, const vector<double> redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file)
{
  double defaultValue1 = value(cosmoPar1);
  double defaultValue2 = value(cosmoPar2);

  string file = dir_output+file_bias_eff_grid;

  ifstream fin(file.c_str());

  if(!fin)
  {
    vector<double> pp1 = linear_bin_vector(nbin_par1, min_par1, max_par1);
    vector<double> pp2 = linear_bin_vector(nbin_par2, min_par2, max_par2);
    fin.clear(); fin.close();
    ofstream fout(file.c_str());

    for (int i=0; i<nbin_par1; i++){
      set_parameter (cosmoPar1, pp1[i]);
      for (int j=0; j<nbin_par2; j++){
	set_parameter (cosmoPar2, pp2[j]);
	double bias_eff = bias_eff_mass (mass, mass_grid, redshift, author_bias, method_SS, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)[0];
	fout << pp1[i] << " " << pp2[j] <<  "  " << bias_eff << endl;
      }
      fout << endl;
    }
    fout.clear(); fout.close();
  }
  fin.clear(); fin.close();

  read_matrix(file, parameter1, parameter2, bias_eff);
  
  set_parameter(cosmoPar1, defaultValue1);
  set_parameter(cosmoPar2, defaultValue2);
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::bias_eff_selection_function (const vector<double> redshift, const double Mass_min, const double Mass_max, const string author_bias, const string author_MF, const string method_SS, const string selection_function_file, const vector<int> cols, const string output_root, const double Delta, const bool isDelta_vir, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file)
{
  // ---------- read the grid file ---------- 
  
  const string file_grid = create_grid_sigmaM(method_SS, 0., output_root, interpType, k_max, input_file, is_parameter_file);

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
  if (mass.size()==0) { 
    string Err = "Error in cosmobl::cosmology::Cosmology::bias_eff of Bias.cpp: mass.size()=0, Mass_min = " + conv(Mass_min,par::fDP3) + ", Mass_max = " + conv(Mass_max,par::fDP3) + ", file_grid = " + file_grid;
    ErrorCBL(Err);
  }

  glob::FuncGrid interp_sigma(mass, sigma, "Spline");
  glob::FuncGrid interp_DlnSigma(mass, dlnsigma, "Spline");
  
  // ---------- read the selection function ---------- 

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);

  glob::FuncGrid2D interp_SF( mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the effective bias ---------- 
  
  const int redshift_bins = redshift.size();
  vector<double> Bias_eff(redshift_bins, 0.);

  for(int i=0; i<redshift_bins; i++){
    auto integrand_num = [&] (const double mass) {
      double sigma = interp_sigma(mass);
      double dlnsigma = interp_DlnSigma(mass);
      double SF = interp_SF(mass, redshift[i]);

      double DD = (isDelta_vir) ? Delta_vir(Delta, redshift[i]) : Delta;
      double MF = mass_function(mass, sigma, dlnsigma, redshift[i], author_MF, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      double BH = bias_halo(mass, sigma, redshift[i], author_bias, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      return BH*MF*SF;
    };

    auto integrand_denom = [&] (const double mass) {
      double DD = (isDelta_vir) ? Delta_vir(Delta, redshift[i]) : Delta;
      double sigma = interp_sigma(mass);
      double dlnsigma = interp_DlnSigma(mass);
      double SF = interp_SF(mass, redshift[i]);
      double MF = mass_function(mass, sigma, dlnsigma, redshift[i], author_MF, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
      return MF*SF;
    };

    Bias_eff[i] = gsl::GSL_integrate_qag(integrand_num, Mass_min, Mass_max)/gsl::GSL_integrate_qag(integrand_denom, Mass_min, Mass_max);
  }


  return Bias_eff;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::bias_eff_selection_function (const double z_min, const double z_max, const double Mass_min, const double Mass_max, const string author_bias, const string author_MF, const string method_SS, const string selection_function_file, const vector<int> cols, const string output_root, const double Delta, const bool isDelta_vir, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file)
{
  vector<double> _kk = logarithmic_bin_vector(200, ((k_min==0) ? 1.e-4 : k_min), k_max);

  vector<double> _Pk = Pk(_kk, method_SS, false, 0., output_root, -1, k_min, k_max, prec);

  glob::FuncGrid interp_Pk(_kk, _Pk, "Spline");

  // compute sigma, dlnsigma

  vector<double> mass = logarithmic_bin_vector(50, Mass_min, Mass_max);
  double rho = rho_m(0.);

  vector<double> _sigma, _dnsigma;

  for(size_t i=0; i<mass.size(); i++){

    const double RR = cosmobl::Radius(mass[i], rho); 

    // Compute sigma
    auto func_sigma = [&] (double _k){
      return pow(TopHat_WF(_k*RR)*_k, 2)*interp_Pk(_k); 
    };
    _sigma.emplace_back(sqrt(1./(2.*pow(par::pi, 2)) *gsl::GSL_integrate_qag (func_sigma, k_min, k_max, prec)));

    // Compute dlnsigma
    const double dRdM = pow(3./(4.*cosmobl::par::pi*rho), 1./3.)*pow(mass[i], -2./3.)/3.;

    auto func_dnsigma = [&] (const double _k)
    {
      double filter = 2*cosmobl::TopHat_WF(_k*RR)*cosmobl::TopHat_WF_D1(_k*RR)*_k*dRdM;
      return filter*pow(_k, 2)*interp_Pk(_k); 
    };

    _dnsigma.emplace_back(1./(2.*pow(par::pi, 2)) *gsl::GSL_integrate_qag (func_dnsigma, k_min, k_max, prec));
  }

  glob::FuncGrid interp_sigma(mass, _sigma, "Spline");
  glob::FuncGrid interp_dnsigma(mass, _dnsigma, "Spline");


  // ---------- read the selection function ---------- 

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);

  glob::FuncGrid2D interp_SF( mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the effective bias ---------- 
        
  auto integrand_num = [&] (const vector<double> parameters) {
    double redshift = parameters[0];
    double mass = parameters[1];

    double sigma = interp_sigma(mass);
    double dlnsigma = interp_dnsigma(mass)*(mass/(2.*sigma*sigma));
    
    double SF = interp_SF(mass, redshift);

    double DD = (isDelta_vir) ? Delta_vir(Delta, redshift) : Delta;
    double MF = mass_function(mass, sigma, dlnsigma, redshift, author_MF, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    double BH = bias_halo(mass, sigma, redshift, author_bias, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    return BH*MF*SF;
  };

  auto integrand_denom = [&] (const vector<double> parameters) {
    double redshift = parameters[0];
    double mass = parameters[1];

    double sigma = interp_sigma(mass);
    double dlnsigma = interp_dnsigma(mass)*(mass/(2.*sigma*sigma));
    double SF = interp_SF(mass, redshift);

    double DD = (isDelta_vir) ? Delta_vir(Delta, redshift) : Delta;
    double MF = mass_function(mass, sigma, dlnsigma, redshift, author_MF, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    return MF*SF;
  };

  cosmobl::cuba::CUBAwrapper IntegrandNum(integrand_num, 2);
  cosmobl::cuba::CUBAwrapper IntegrandDenom(integrand_denom, 2);

  vector<vector<double>> limits = {{z_min, z_max}, {Mass_min, Mass_max}};

  return IntegrandNum.IntegrateVegas(limits)/IntegrandDenom.IntegrateVegas(limits);
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::bias_eff_selection_function_2 (const double redshift, const double Mass_min, const double Mass_max, const string author_bias, const string author_MF, const string method_SS, const string selection_function_file, const vector<int> cols, const string output_root, const double Delta, const bool isDelta_vir, const double kk, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file)
{
  vector<double> _kk = logarithmic_bin_vector(200, ((k_min==0) ? 1.e-4 : k_min), k_max);

  vector<double> _Pk = Pk(_kk, method_SS, false, 0., output_root, -1, k_min, k_max, prec);

  glob::FuncGrid interp_Pk(_kk, _Pk, "Spline");

  // compute sigma, dlnsigma

  vector<double> mass = logarithmic_bin_vector(50, Mass_min, Mass_max);
  double rho = rho_m(0.);

  vector<double> _sigma, _dnsigma;

  for(size_t i=0; i<mass.size(); i++){

    const double RR = cosmobl::Radius(mass[i], rho); 

    // Compute sigma
    auto func_sigma = [&] (double _k){
      return pow(TopHat_WF(_k*RR)*_k, 2)*interp_Pk(_k); 
    };
    _sigma.emplace_back(sqrt(1./(2.*pow(par::pi, 2)) *gsl::GSL_integrate_qag (func_sigma, k_min, k_max, prec)));

    // Compute dlnsigma
    const double dRdM = pow(3./(4.*cosmobl::par::pi*rho), 1./3.)*pow(mass[i], -2./3.)/3.;

    auto func_dnsigma = [&] (const double _k)
    {
      double filter = 2*cosmobl::TopHat_WF(_k*RR)*cosmobl::TopHat_WF_D1(_k*RR)*_k*dRdM;
      return filter*pow(_k, 2)*interp_Pk(_k); 
    };

    _dnsigma.emplace_back(1./(2.*pow(par::pi, 2)) *gsl::GSL_integrate_qag (func_dnsigma, k_min, k_max, prec));
  }

  glob::FuncGrid interp_sigma(mass, _sigma, "Spline");
  glob::FuncGrid interp_dnsigma(mass, _dnsigma, "Spline");

  // ---------- read the selection function ---------- 

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);

  glob::FuncGrid2D interp_SF( mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the effective bias ---------- 

  double DD = (isDelta_vir) ? Delta_vir(Delta, redshift) : Delta;

  auto funcN = [&] (const double mass){
    double sigma = interp_sigma(mass);
    double dlnsigma = interp_dnsigma(mass)*(mass/(2.*sigma*sigma));

    double SF = interp_SF(mass, redshift);

    double MF = mass_function(mass, sigma, dlnsigma, redshift, author_MF, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    double BH = bias_halo(mass, sigma, redshift, author_bias, output_root, interpType, DD, kk, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    return BH*MF*SF;
  };

  auto funcD = [&] (const double mass){
    double sigma = interp_sigma(mass);
    double dlnsigma = interp_dnsigma(mass)*(mass/(2.*sigma*sigma));

    double SF = interp_SF(mass, redshift);

    double MF = mass_function(mass, sigma, dlnsigma, redshift, author_MF, output_root, DD, interpType, norm, k_min, k_max, prec, method_SS, input_file, is_parameter_file);
    return MF*SF;
  };



  return gsl::GSL_integrate_qag(funcN, Mass_min, Mass_max)/gsl::GSL_integrate_qag(funcD, Mass_min, Mass_max);
}
