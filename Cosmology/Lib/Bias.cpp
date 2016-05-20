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


double cosmobl::Cosmology::bias_halo (const double Mass, const double redshift, const string author, const string method_SS, const string output_root, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  double SSS = SSM_norm(Mass, method_SS, redshift, output_root, k_max, file_par); 
  double Sigma = sqrt(SSS); 
  
  double bias = bias_halo_generator (Sigma, redshift, author, Delta); 

  if (m_fNL!=0) {
    double fact = bias_correction(kk, Mass, method_SS, output_root, norm, k_min, k_max, GSL, prec, file_par) * SSS * (bias-1) * (bias-1);
    bias += fact; 
  }

  return bias; 
}


// =====================================================================================


double cosmobl::Cosmology::bias_halo (const double Mass, const double Sigma, const double redshift, const string author_bias, const string output_root, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string method_SS, const string file_par) 
{
  double bias = bias_halo_generator(Sigma, redshift, author_bias, Delta); 
  
  if (m_fNL!=0) {
    double z0 = 0.;
    double fact = bias_correction(kk, Mass, method_SS, output_root, norm, k_min, k_max, GSL, prec, file_par) * SSM_norm(Mass, method_SS, z0, output_root, k_max, file_par) * (bias-1) * (bias-1);
    bias += fact; 
  }

  return bias; 
}


// =====================================================================================


double cosmobl::Cosmology::bias_halo_generator (const double Sigma, const double redshift, const string author, const double Delta) const
{
  double Z0 = 0.;
  double deltacz = deltac(redshift);
  double sigmaz = Sigma*DD(redshift)/DD(Z0);
  
  double bias = -1000.;

  if (author=="ST99") {
    double aa = 0.707;
    double pp = 0.3;
    double ni = pow(deltacz/sigmaz,2); 
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
    double ni = deltacz/sigmaz; 
    bias = 1.-AA*pow(ni,aa)/(pow(ni,aa)+pow(deltacz,aa))+BB*pow(ni,bb)+CC*pow(ni,ccc);
  }
  
  if (bias<-100) {
    string Err = "Error in cosmobl::Cosmology::bias_halo of Bias.cpp: author = " + author + "!";
    ErrorMsg(Err);
  }
  
  return bias;
}


// =====================================================================================


double cosmobl::Cosmology::bias_eff (const double Mass_min, const double Mass_max, const double redshift, const string author_bias, const string author_MF, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int Num, const double stepsize, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  // ---------- read the grid file ---------- 
  
  double zero = 0.;
  string file_grid = create_grid_sigmaM(method_SS, zero, output_root, interpType, Num, stepsize, k_max, file_par);
  //cout <<"grid file: "<<file_grid<<endl;
  ifstream fin (file_grid.c_str()); checkIO (file_grid,1); 
  
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

    double MF = mass_function(mass[mm], sigma[mm], dlnsigma[mm], redshift, author_MF, output_root, Delta, interpType, Num, stepsize, norm, k_min, k_max, GSL, prec, method_SS, file_par);
    
    Bias_eff += bias_halo(mass[mm], sigma[mm], redshift, author_bias, output_root, Delta, kk, norm, k_min, k_max, GSL, prec, method_SS, file_par)*MF*(mass[mm+1]-mass[mm]);

    Norm += MF*(mass[mm+1]-mass[mm]);
  }

  return Bias_eff/Norm;
}


// =====================================================================================


double cosmobl::Cosmology::bias_eff (const vector<double> MM, const vector<double> MF, const double redshift, const string author_bias, const string method_SS, const string output_root, const double Delta, const double kk, const string interpType, const int Num, const double stepsize, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  // ---------- read the grid file ---------- 
  
  double zero = 0.;
  string file_grid = create_grid_sigmaM(method_SS, zero, output_root, interpType, Num, stepsize, k_max, file_par);
  ifstream fin(file_grid.c_str()); checkIO(file_grid, 1); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma;

  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Min(MM)<Mass && Mass<Max(MM)) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
    }
  }
  if (mass.size()==0) { 
    string Err = "Error in cosmobl::Cosmology::bias_eff of Bias.cpp: mass.size()=0, Min(MM) = " + conv(Min(MM),par::fDP3) + ", Max(MM) = " + conv(Max(MM),par::fDP3) + ", file_grid = " + file_grid;
    ErrorMsg(Err);
  }
  

  // ---------- compute the effective bias ---------- 
  
  double Bias_eff = 0., Norm = 0.;
  double mf, sig, err = -1;
  
  for (size_t k=0; k<MM.size()-1; k++) {
    mf = MF[k];
    sig = interpolated(MM[k], mass, sigma, "Linear");
    
    if (err/sig>0.1) { 
      string Err = "Error in cosmobl::Cosmology::bias_eff of Bias.cpp: err/sig = " + conv(err/sig, par::fDP3) + "!";
      ErrorMsg(Err);
    } 

    Bias_eff += bias_halo(MM[k], sig, redshift, author_bias, output_root, Delta, kk, norm, k_min, k_max, GSL, prec, method_SS, file_par)*mf*(MM[k+1]-MM[k]);
    Norm += mf*(MM[k+1]-MM[k]);
  }

  return Bias_eff/Norm;
}
