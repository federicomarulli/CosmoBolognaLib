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
 *  @file Cosmology/Lib/RSD.cpp
 *
 *  @brief Methods of the class Cosmology used to model redshift-space
 *  distortions
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the dynamic redshift-space distortions of
 *  two-point statistics
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Cosmology::linear_growth_rate (const double redshift, const double kk) const
{
  if (m_wa!=0) WarningMsgCBL("the current implementation is not correct with m_wa="+conv(m_wa, par::fDP3), "linear_growth_rate", "RSD.cpp");
  
  if (redshift>10) WarningMsgCBL("the approximation to compute the linear growth rate is not very accurate at z>10 (see Kiakotou, Elgaroy & Lahav 2008)!", "linear_growth_rate", "RSD.cpp");
  

  // Wang & Steinhardt 1998

  double alpha0 = 3./(5.-m_w0/(1.-m_w0));
  double alpha1 = 3./125.*(1.-m_w0)*(1.-3.*m_w0/2.)/pow(1.-6.*m_w0/5.,3);
  double alpha = alpha0+alpha1*(1-OmegaM(redshift));

  
  // Kiakotou, Elgarøy & Lahav 2008
  
  double fnu = m_Omega_neutrinos/m_Omega_matter;
  double ff = -1.;

  if (fnu>0) {
    if (kk<0) ErrorCBL("kk<0!", "linear_growth_rate", "RSD.cpp");

    double kkk[] = {0.001, 0.01, 0.05, 0.07, 0.1, 0.5};
    double aa[] = {0., 0.132, 0.613, 0.733, 0.786, 0.813};
    double bb[] = {0., 1.62, 5.59, 6., 5.09, 0.803};
    double cc[] = {0., 7.13, 21.13, 21.45, 15.5, -0.844};
    
    vector<double> lgKK;
//lgKK.push_back(log10(0.005)); // to improve the interpolation
    for (int i=0; i<6; i++) lgKK.push_back(log10(kkk[i]));

    //vector<double> KK;
    //KK.push_back(0.005); // to improve the interpolation
    //for (int i=0; i<6; i++) KK.push_back(kkk[i]);

    vector<double> MU;

    for (int i=0; i<6; i++) MU.push_back(1.-aa[i]*m_Omega_DE*fnu+bb[i]*fnu*fnu-cc[i]*fnu*fnu*fnu); 
    
    double mu;

    if (kk<0.001) mu = interpolated(log10(kk), lgKK, MU, "Poly");
    else if (0.001<=kk && kk<=0.5) mu = interpolated(log10(kk), lgKK, MU, "Poly");

    else mu = pow(1-fnu, alpha0); 
    
    ff = mu*pow(OmegaM(redshift), alpha);
  }

  else 
    ff = pow(OmegaM(redshift), alpha);
  
  ff += (alpha-4./7.)*m_Omega_k; // Gong et al. 2009 

  return ff;
}


// =====================================================================================


double cbl::cosmology::Cosmology::fsigma8 (const double redshift, const string method_Pk, const bool store_output, const string output_root, const double kk, const bool NL, const double k_min, const double k_max, const double prec, const string file_par) const
{
  return linear_growth_rate(redshift, kk)*sigma8_Pk(method_Pk, redshift, store_output, output_root, NL, k_min, k_max, prec, file_par);
}


// =====================================================================================


double cbl::cosmology::Cosmology::beta (const double redshift, const double bias, const double kk) const
{
  return linear_growth_rate(redshift, kk)/bias;
}


// =====================================================================================


double cbl::cosmology::Cosmology::error_beta (const double redshift, const double bias, const double err_bias, const double kk) const
{
  return linear_growth_rate(redshift, kk)/pow(bias,2)*err_bias;
}


// =====================================================================================


double cbl::cosmology::Cosmology::beta (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  return linear_growth_rate(redshift)/bias_eff(Mass_min, Mass_max, redshift, model_bias, model_MF, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
}


// =====================================================================================


double cbl::cosmology::Cosmology::error_beta (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double err_bias, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{
  return linear_growth_rate(redshift)/pow(bias_eff(Mass_min, Mass_max, redshift, model_bias, model_MF, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file),2)*err_bias;
}


// =====================================================================================


double cbl::cosmology::Cosmology::beta (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  return linear_growth_rate(redshift)/bias_eff(MM, MF, redshift, model_bias, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
}


// =====================================================================================


double cbl::cosmology::Cosmology::error_beta (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const double err_bias, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  return linear_growth_rate(redshift)/pow(bias_eff(MM, MF, redshift, model_bias, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file),2)*err_bias;
}


// =====================================================================================


double cbl::cosmology::Cosmology::error_beta_measured (const double Volume, const double density, const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{ // from Eq. 20 of Bianchi et al. 2012
  
  double bias = bias_eff(Mass_min, Mass_max, redshift, model_bias, model_MF, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);

  return relative_error_beta (bias, Volume, density);
}


// =====================================================================================


double cbl::cosmology::Cosmology::quadrupole (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  double Beta = beta(Mass_min, Mass_max, redshift, model_bias, model_MF, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  return (4./3.*Beta+4./7.*Beta*Beta)/(1.+2./3*Beta+1./5.*Beta*Beta);
}


// =====================================================================================


double cbl::cosmology::Cosmology::quadrupole (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  double Beta = beta(MM, MF, redshift, model_bias, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  return (4./3.*Beta+4./7.*Beta*Beta)/(1.+2./3*Beta+1./5.*Beta*Beta);
}


// =====================================================================================


double cbl::cosmology::Cosmology::Pk_DeltaDelta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1)
{
  double Pkdd = 0;
  if(author == "Pezzotta" || author == "Bel")
  Pkdd = Pk_DM(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);
  else WarningMsgCBL("the current implementation is not correct with author = " + author, "Pk_DeltaDelta_fitting_function", "RSD.cpp");
  return Pkdd;
}


// =====================================================================================


double cbl::cosmology::Cosmology::Pk_DeltaTheta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const bool NL, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1)
{
  double sigma8_z = sigma8_Pk(method_Pk, redshift, store_output, output_root, NL, k_min, k_max, prec, file_par);
  double kd = 1./(-0.017 + 1.496*pow(sigma8_z, 2.));
  double Pkdt = 0;
  if (author == "Pezzotta"){
    Pkdt = sqrt(Pk_DM(kk, method_Pk, true, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*Pk_DM(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1))*exp(-kk/kd);
  }
  else if (author == "Bel"){
    double b = 0.091 + 0.702*sigma8_z*sigma8_z;
    Pkdt = sqrt(Pk_DM(kk, method_Pk, true, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*Pk_DM(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1))*exp(-kk/kd-b*pow(kk,6.0));
  }
  else WarningMsgCBL("the current implementation is not correct with author = " + author, "Pk_DeltaTheta_fitting_function", "RSD.cpp");
    return Pkdt;
}


// =====================================================================================


double cbl::cosmology::Cosmology::Pk_ThetaTheta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const bool NL, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1)
{
  double sigma8_z = sigma8_Pk(method_Pk, redshift, store_output, output_root, NL, k_min, k_max, prec, file_par);
  double Pktt = 0;
  if (author == "Pezzotta"){
    double kt = 1./(-0.048 + 1.917*sigma8_z*sigma8_z);
    Pktt = Pk_DM(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*exp(-kk/kt);
  }

  else if (author == "Bel"){
    double a1 = -0.817 + 3.198*sigma8_z;
    double a2 = 0.877 - 4.191*sigma8_z;
    double a3 = -1.199 + 4.629*sigma8_z;
    Pktt = Pk_DM(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*exp(-kk*(a1 + a2*kk + a3*kk*kk));
  }
  
  else WarningMsgCBL("the current implementation is not correct with author = " + author, "Pk_ThetaTheta_fitting_function", "RSD.cpp");

  return Pktt;
}
