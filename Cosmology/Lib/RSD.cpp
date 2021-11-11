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
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Cosmology::fsigma8 (const double redshift, const string method_Pk, const bool store_output, const string output_root, const bool NL, const double k_min, const double k_max, const double prec, const string file_par) const
{
  return linear_growth_rate(redshift)*sigma8_Pk(method_Pk, redshift, store_output, output_root, NL, k_min, k_max, prec, file_par);
}


// =====================================================================================


double cbl::cosmology::Cosmology::beta (const double redshift, const double bias) const
{
  return linear_growth_rate(redshift)/bias;
}


// =====================================================================================


double cbl::cosmology::Cosmology::error_beta (const double redshift, const double bias, const double err_bias) const
{
  return linear_growth_rate(redshift)/pow(bias, 2)*err_bias;
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
  return (4./3.*Beta+4./7.*Beta*Beta)/(1.+2./3.*Beta+1./5.*Beta*Beta);
}


// =====================================================================================


double cbl::cosmology::Cosmology::quadrupole (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const double kk, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) 
{
  double Beta = beta(MM, MF, redshift, model_bias, method_SS, store_output, output_root, Delta, kk, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  return (4./3.*Beta+4./7.*Beta*Beta)/(1.+2./3.*Beta+1./5.*Beta*Beta);
}


// =====================================================================================


double cbl::cosmology::Cosmology::Pk_DeltaDelta_fitting_function (const double kk, const std::string method_Pk, const double redshift, const std::string author, const bool store_output, const std::string output_root, const int norm, double k_min, double k_max, const double prec, const std::string file_par, const bool unit1)
{
  double Pkdd = 0;
  if(author == "Pezzotta" || author == "Bel")
  Pkdd = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);
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
    Pkdt = sqrt(Pk_matter(kk, method_Pk, true, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1))*exp(-kk/kd);
  }
  else if (author == "Bel"){
    double b = 0.091 + 0.702*sigma8_z*sigma8_z;
    Pkdt = sqrt(Pk_matter(kk, method_Pk, true, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1))*exp(-kk/kd-b*pow(kk,6.0));
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
    Pktt = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*exp(-kk/kt);
  }

  else if (author == "Bel"){
    double a1 = -0.817 + 3.198*sigma8_z;
    double a2 = 0.877 - 4.191*sigma8_z;
    double a3 = -1.199 + 4.629*sigma8_z;
    Pktt = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1)*exp(-kk*(a1 + a2*kk + a3*kk*kk));
  }
  
  else WarningMsgCBL("the current implementation is not correct with author = " + author, "Pk_ThetaTheta_fitting_function", "RSD.cpp");

  return Pktt;
}


// =====================================================================================


double cbl::cosmology::Cosmology::sigma_v (const double redshift, const std::string method_Pk, const bool store_output, const std::string output_root, const int norm, const double k_min, const double k_max, const int bin_k, const double prec, const std::string file_par, const bool unit1)
{
  const vector<double> kk = logarithmic_bin_vector(bin_k, k_min, k_max);
  const vector<double> Pk = Pk_matter(kk, method_Pk, false, redshift, store_output, output_root, norm, k_min, k_max, prec, file_par, unit1);

  auto interp_Pk = glob::FuncGrid(kk, Pk, "Spline");
  
  auto integrand = [&] (const double log_q)
		   {
		     const double qq = exp(log_q);
		     return qq*interp_Pk(qq);
		   };

  return sqrt(1./(6.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_cquad(integrand, log(k_min), log(k_max), 1.e-5, 0., 10000));
}
