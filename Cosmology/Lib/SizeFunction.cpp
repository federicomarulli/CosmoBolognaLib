/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Tommaso Ronconi      *
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
 *  @file Cosmology/Lib/SizeFunction.cpp
 *
 *  @brief Methods of the class SizeFunction used to model the void
 *  size function
 *
 *  This file contains the implementation of the methods of the class
 *  SizeFunction used to model the size function of cosmic voids
 *
 *  @authors Federico Marulli, Tommaso Ronconi
 *
 *  @authors federico.marulli3@unibo.it, tommaso.ronconi@studio.unibo.it
 */

#include "SizeFunction.h"

using namespace std;

using namespace cbl;

// =====================================================================================

double cbl::cosmology::SizeFunction::deltav_L (const double deltav_NL, const double b_eff, double slope, double offset) const
{
  if (b_eff == 1.) {
    slope = 1.;
    offset = 0.;
  }

  return 1.594 * (1. - pow(1 + deltav_NL / (slope * b_eff + offset), (-1. / 1.594)));
}

// =====================================================================================

double cbl::cosmology::SizeFunction::deltav_NL (const double deltav_L) const
{
  return pow((1. - deltav_L / 1.594), -1.594) - 1.;
}

// =====================================================================================

double cbl::cosmology::SizeFunction::r_rL (const double deltav) const
{
  return pow(pow((1. - deltav / 1.594), -1.594), -1. / 3.);
}

// =====================================================================================

double cbl::cosmology::SizeFunction::f_nu (const double SS, const double del_v, const double del_c) const
{
  double radnu = fabs(del_v) / SS;
  double nu = pow(radnu, 2.);
  double DDD = fabs(del_v) / (del_c + fabs(del_v));
  double xx = DDD / radnu;

  if (xx <= 0.276)
    return sqrt(2. / par::pi) * radnu * exp(-0.5 * nu);

  else {
    double ff = 0.;
    for (int j=1; j<=4; j++)
      ff += 2. * exp(-pow((j * par::pi * xx), 2.) / 2.) * j * par::pi * pow(xx, 2.) * sin(j * par::pi * DDD);
    return ff;
  }
}


// =====================================================================================


vector<double> cbl::cosmology::SizeFunction::AP_corr (const std::shared_ptr<cosmology::Cosmology> cosm_true, const std::vector<double> redshift)
{
  std::vector<double> APcorr_vect(redshift.size(), 0.);

  for (size_t ii=0; ii<redshift.size(); ii++) {
    double q_par = m_cosmology->Hubble(redshift[ii]) / cosm_true->Hubble(redshift[ii]);
    double q_perp = cosm_true->D_M(redshift[ii]) / m_cosmology->D_M(redshift[ii]);
    APcorr_vect[ii] = pow(q_par * q_perp * q_perp, (-1. / 3.));
  }
  return APcorr_vect;
}

// =====================================================================================

double cbl::cosmology::SizeFunction::size_function (const double RV, const double redshift, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{
  double del_v = deltav_L(deltav_NL, b_eff, slope, offset);
  double RL;

  if ((model == "Vdn") || (model == "SvdW"))
    RL = RV / r_rL(del_v);

  else if (model == "linear")
    RL = RV;

  else {
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "size_function", "SizeFunction.cpp");
    return 0;
  }

  double fact = m_cosmology->DN(redshift, 0.);
  
  PkXi PX(m_cosmology);
  Sigma SG(m_cosmology);
  
  function<double(double)> Pk0 = PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
  
  double sigmaR = sqrt(SG.sigma2R({RL}, Pk0, k_max)[0]);
  double sigmaRz = sigmaR * fact;
  double SSSR = sigmaRz * sigmaRz;
  double Dln_SigmaR = SG.dnsigma2R(1, {RL}, Pk0, k_max)[0] * (RL / (2. * SSSR)) * fact * fact;

  if (model == "Vdn")
    return f_nu(sigmaRz, del_v, del_c) / volume_sphere(RV) * fabs(Dln_SigmaR);
  else
    return f_nu(sigmaRz, del_v, del_c) / volume_sphere(RL) * fabs(Dln_SigmaR);
}

// =====================================================================================

std::vector<double> cbl::cosmology::SizeFunction::size_function (const std::vector<double> RV, const double redshift, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{
  double del_v = deltav_L(deltav_NL, b_eff, slope, offset);

  std::vector<double> RL(RV.size());

  if ((model == "Vdn") || (model == "SvdW"))
    for (size_t i=0; i<RV.size(); i++)
      RL[i] = RV[i] / r_rL(del_v);

  else if (model == "linear")
    for (size_t i=0; i<RV.size(); i++)
      RL[i] = RV[i];

  else
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "size_function", "SizeFunction.cpp");

  
  double fact = m_cosmology->DN(redshift);

  PkXi PX(m_cosmology);
  Sigma SG(m_cosmology);
  
  function<double(double)> Pk0 = PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  std::vector<double> sigma2R = SG.sigma2R(RL, Pk0, k_max);
  std::vector<double> Dln_SigmaR = SG.dnsigma2R(1, RL, Pk0, k_max);
  std::vector<double> sigmaR(RV.size());
  std::vector<double> sigmaRz(RV.size());
  std::vector<double> SSSR(RV.size());

  for (size_t i=0; i<RV.size(); i++) {
    sigmaR[i] = sqrt(sigma2R[i]);
    sigmaRz[i] = sigmaR[i] * fact;
    SSSR[i] = sigmaRz[i] * sigmaRz[i];
    Dln_SigmaR[i] *= (RL[i] / (2. * SSSR[i])) * fact * fact;
  }

  std::vector<double> result(RV.size());

  if (model == "Vdn")
    for (size_t i=0; i<RV.size(); i++)
      result[i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RV[i]) * fabs(Dln_SigmaR[i]);
  else
    for (size_t i=0; i<RV.size(); i++)
      result[i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RL[i]) * fabs(Dln_SigmaR[i]);

  return result;
}

// =====================================================================================

std::vector<std::vector<double>> cbl::cosmology::SizeFunction::Nvoids (const double min_r, const double max_r, const int num_bins, const double mean_z, const double Volume, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{
  double del_v = deltav_L(deltav_NL, b_eff, slope, offset);
  double NLcorr;

  if ((model == "Vdn") || (model == "SvdW"))
    NLcorr = r_rL(del_v);
  else if (model == "linear")
    NLcorr = 1.;
  else
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "Nvoids", "SizeFunction.cpp");

  double fact = m_cosmology->DN(mean_z);

  std::vector<double> r_bins = cbl::logarithmic_bin_vector(num_bins + 1, min_r, max_r);

  std::vector<double> RV(num_bins);
  std::vector<double> RL(num_bins);

  for (int i = 0; i < num_bins; i++) {
    RV[i] = pow(10., log10(r_bins[i]) + 0.5 * (log10(r_bins[i + 1] / r_bins[i])));
    RL[i] = RV[i] / NLcorr;
  }

  PkXi PX(m_cosmology);
  Sigma SG(m_cosmology);
  
  function<double(double)> Pk0 = PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  std::vector<double> sigma2R = SG.sigma2R(RL, Pk0, k_max);
  std::vector<double> Dln_SigmaR = SG.dnsigma2R(1, RL, Pk0, k_max);
  std::vector<double> sigmaR(RV.size());
  std::vector<double> sigmaRz(RV.size());
  std::vector<double> SSSR(RV.size());

  for (int i=0; i<num_bins; i++) {
    sigmaR[i] = sqrt(sigma2R[i]);
    sigmaRz[i] = sigmaR[i] * fact;
    SSSR[i] = sigmaRz[i] * sigmaRz[i];
    Dln_SigmaR[i] *= (RL[i] / (2. * SSSR[i])) * fact * fact;
  }

  std::vector<std::vector<double>> result(2, std::vector<double>(num_bins));

  if (model == "Vdn")
    for (int i=0; i<num_bins; i++) {
      result[0][i] = RV[i];
      result[1][i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RV[i]) * fabs(Dln_SigmaR[i]) * (r_bins[i + 1] - r_bins[i]) / RV[i] * Volume;
    }

  else
    for (int i=0; i<num_bins; i++) {
      result[0][i] = RV[i];
      result[1][i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RL[i]) * fabs(Dln_SigmaR[i]) * (r_bins[i + 1] - r_bins[i]) / RV[i] * Volume;
    }
  
  return result;
}

// =====================================================================================

std::vector<std::vector<double>> cbl::cosmology::SizeFunction::Nvoids (const double min_r, const double max_r, const int num_bins, const double min_z, const double max_z, const double mean_z, const double Area, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{
  double del_v = deltav_L(deltav_NL, b_eff, slope, offset);
  double NLcorr;

  if ((model == "Vdn") || (model == "SvdW"))
    NLcorr = r_rL(del_v);
  else if (model == "linear")
    NLcorr = 1.;
  else
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "Nvoids", "SizeFunction.cpp");

  double fact = m_cosmology->DN(mean_z);

  std::vector<double> r_bins = cbl::logarithmic_bin_vector(num_bins + 1, min_r, max_r);

  std::vector<double> RV(num_bins);
  std::vector<double> RL(num_bins);

  for (int i=0; i<num_bins; i++) {
    RV[i] = pow(10., log10(r_bins[i]) + 0.5 * (log10(r_bins[i + 1] / r_bins[i])));
    RL[i] = RV[i] / NLcorr;
  }

  PkXi PX(m_cosmology);
  Sigma SG(m_cosmology);

  function<double(double)> Pk0 = PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  std::vector<double> sigma2R = SG.sigma2R(RL, Pk0, k_max);
  std::vector<double> Dln_SigmaR = SG.dnsigma2R(1, RL, Pk0, k_max);
  std::vector<double> sigmaR(RV.size());
  std::vector<double> sigmaRz(RV.size());
  std::vector<double> SSSR(RV.size());

  for (int i=0; i<num_bins; i++) {
    sigmaR[i] = sqrt(sigma2R[i]);
    sigmaRz[i] = sigmaR[i] * fact;
    SSSR[i] = sigmaRz[i] * sigmaRz[i];
    Dln_SigmaR[i] *= (RL[i] / (2. * SSSR[i])) * fact * fact;
  }

  std::vector<std::vector<double>> result(2, std::vector<double>(num_bins));
  double volume = m_cosmology->Volume(min_z, max_z, Area);

  if (model == "Vdn")
    for (int i=0; i<num_bins; i++) {
      result[0][i] = RV[i];
      result[1][i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RV[i]) * fabs(Dln_SigmaR[i]) * (r_bins[i + 1] - r_bins[i]) / RV[i] * volume;
    }

  else
    for (int i=0; i<num_bins; i++) {
      result[0][i] = RV[i];
      result[1][i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RL[i]) * fabs(Dln_SigmaR[i]) * (r_bins[i + 1] - r_bins[i]) / RV[i] * volume;
    }
  
  return result;
}

// =====================================================================================

std::vector<std::vector<double>> cbl::cosmology::SizeFunction::Nvoids (const double min_r, const double max_r, const int num_bins, const double min_z, const double max_z, const double mean_z, const double Area, const std::vector<std::vector<double>> Sigma2R, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c) const
{
  double del_v = deltav_L(deltav_NL, b_eff, slope, offset);
  double NLcorr;

  if ((model == "Vdn") || (model == "SvdW"))
    NLcorr = r_rL(del_v);
  else if (model == "linear")
    NLcorr = 1.;
  else
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "Nvoids", "SizeFunction.cpp");

  double fact = m_cosmology->DN(mean_z);

  std::vector<double> r_bins = cbl::logarithmic_bin_vector(num_bins + 1, min_r, max_r);

  std::vector<double> RV(num_bins);
  std::vector<double> RL(num_bins);

  for (int i=0; i <num_bins; i++) {
    RV[i] = pow(10., log10(r_bins[i]) + 0.5 * (log10(r_bins[i + 1] / r_bins[i])));
    RL[i] = RV[i] / NLcorr;
  }

  std::vector<double> sigmaR(RV.size());
  std::vector<double> sigmaRz(RV.size());
  std::vector<double> SSSR(RV.size());
  std::vector<double> Dln_SigmaR(RV.size());

  function<double(double)> sigma2R = glob::FuncGrid(Sigma2R[0], Sigma2R[1], "Spline", cbl::BinType::_linear_);

  for (int i=0; i<num_bins; i++) {
    sigmaR[i] = sqrt(sigma2R(RL[i]));
    sigmaRz[i] = sigmaR[i] * fact;
    SSSR[i] = sigmaRz[i] * sigmaRz[i];
    Dln_SigmaR[i] = wrapper::gsl::GSL_derivative(sigma2R, RL[i], 1.e-3) * (RL[i] / (2. * SSSR[i])) * fact * fact;
  }

  std::vector<std::vector<double>> result(2, std::vector<double>(num_bins));
  double volume = m_cosmology->Volume(min_z, max_z, Area);

  if (model == "Vdn")
    for (int i=0; i <num_bins; i++) {
      result[0][i] = RV[i];
      result[1][i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RV[i]) * fabs(Dln_SigmaR[i]) * (r_bins[i + 1] - r_bins[i]) / RV[i] * volume;
    }

  else
    for (int i=0; i<num_bins; i++) {
      result[0][i] = RV[i];
      result[1][i] = f_nu(sigmaRz[i], del_v, del_c) / volume_sphere(RL[i]) * fabs(Dln_SigmaR[i]) * (r_bins[i + 1] - r_bins[i]) / RV[i] * volume;
    }
  
  return result;
}

// =====================================================================================

double cbl::cosmology::SizeFunction::size_function (const double RV, const double redshift, const std::string model_mf, const double del_v, const std::string model_sf, const std::string method_Pk, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file)
{
  double RL;

  if ((model_sf == "Vdn") || (model_sf == "SvdW"))
    RL = RV / r_rL(del_v);

  else if (model_sf == "linear")
    RL = RV;

  else {
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "size_function", "SizeFunction.cpp");
    return 0;
  }

  double RHO = m_cosmology->rho_m(redshift, true);
  double MM = Mass(RL, RHO);

  MassFunction MF(m_cosmology);

  if (model_sf == "Vdn")
    return 3. * MM * MF.mass_function(MM, redshift, model_mf, method_Pk, store_output, output_root, Delta, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file, false, del_v) * deltav_NL(del_v);

  else if ((model_sf == "SvdW") || (model_sf == "linear"))
    return 3. * MM * MF.mass_function(MM, redshift, model_mf, method_Pk, store_output, output_root, Delta, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file, false, del_v);

  else {
    ErrorCBL("the model name is not allowed; the allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)", "size_function", "SizeFunction.cpp");    
    return 0;
  }
}
