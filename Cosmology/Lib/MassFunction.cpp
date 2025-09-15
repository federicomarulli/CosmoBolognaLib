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
 *  @file Cosmology/Lib/MassFunction.cpp
 *
 *  @brief Methods of the class MassFunction used to model the halo
 *  mass function
 *
 *  This file contains the implementation of the methods of the class
 *  MassFunction used to model the mass function of dark matter haloes
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "MassFunction.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::MassFunction::mass_function (const double Mass, const double redshift, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const bool default_delta, const double delta_t) const 
{
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  double MASS = Mass * fact;
  double RADIUS;

  const bool unit1 = true;

  Sigma SG(m_cosmology);
  
  double SSS = SG.sigma2M({MASS}, method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file, unit1)[0];
  double Sigma = sqrt(SSS);
  double Dln_Sigma = SG.dnsigma2M(1, {MASS}, method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file, unit1)[0]*(MASS/(2.*SSS));
  double DlnSigma_DlnR;

  const double RHO = m_cosmology->rho_m(0., true);
  double MF;
  if (model_MF=="Castro") {
      RADIUS = pow(3*MASS*0.25/cbl::par::pi/RHO, 1./3);
      DlnSigma_DlnR = SG.dnsigma2R(1, {RADIUS},  method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file, unit1)[0]*(RADIUS / (2. * SSS));
      MF = RHO/(MASS*MASS)*DlnSigma_DlnR/(-3) * m_MF_generator(Sigma, DlnSigma_DlnR, redshift, model_MF, Delta, default_delta, delta_t) * pow(fact, 4.);
  }
  else
    MF = RHO/(MASS*MASS)*m_MF_generator(Sigma, Dln_Sigma, redshift, model_MF, Delta, default_delta, delta_t)*pow(fact, 4.);
  
  if (m_cosmology->fNL() != 0 && default_delta) {
    NG ng(m_cosmology);
    MF *= ng.MF_correction(MASS, redshift, model_MF, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }
  
  else if (m_cosmology->fNL() != 0 && !default_delta)
    ErrorCBL("Non Gaussianity still not available for user-defined density contrast threshold!", "mass_function", "MassFunction.cpp");

  return MF;
}


// =====================================================================================


std::vector<double> cbl::cosmology::MassFunction::mass_function (const std::vector<double> Mass, const double redshift, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const bool default_delta, const double delta_t) const
{
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  vector<double> MASS(Mass.size()), RADII(Mass.size());
  for (size_t ii = 0; ii < Mass.size(); ii++)
    MASS[ii] = Mass[ii] * fact;

  const bool unit1 = true;

  PkXi PX(m_cosmology);
  Sigma SG(m_cosmology);
  
  function<double(double)> Pk0 = PX.Pk_interp_for_sigma(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
  
  vector<double> SSS = SG.sigma2M(MASS, Pk0, k_max, unit1);
  SSS = SG.sigma2M(MASS, method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file); 
  
  vector<double> Dln_Sigma = SG.dnsigma2M(1, MASS, Pk0, k_max, unit1);
  vector<double> DlnSigma_DlnR;
  const double RHO = m_cosmology->rho_m(0., true);
  if (model_MF=="Castro") {
    for (size_t ii = 0; ii < Mass.size(); ii++)
      RADII[ii] = pow(3*MASS[ii]*0.25/cbl::par::pi/RHO, 1./3);
    DlnSigma_DlnR = SG.dnsigma2R(1, RADII, Pk0, k_max, unit1);
  }

  vector<double> MF(Mass.size());
  double D_N = m_cosmology->DN(redshift);
  double Sigma;
  for (size_t ii = 0; ii < Mass.size(); ii++) {
    Sigma = sqrt(SSS[ii]);
    if (model_MF=="Castro") {
      DlnSigma_DlnR[ii] *=(RADII[ii] / (2. * SSS[ii]));
      MF[ii] = RHO / (MASS[ii] * MASS[ii]) * DlnSigma_DlnR[ii]/(-3) * m_MF_generator(Sigma, DlnSigma_DlnR[ii], redshift, D_N, model_MF, Delta, default_delta, delta_t) * pow(fact, 4.);
    }
    else{
      Dln_Sigma[ii] *= (MASS[ii] / (2. * SSS[ii]));
      MF[ii] = RHO / (MASS[ii] * MASS[ii]) * m_MF_generator(Sigma, Dln_Sigma[ii], redshift, D_N, model_MF, Delta, default_delta, delta_t) * pow(fact, 4.);
    }
  }

  if (m_cosmology->fNL() != 0 && default_delta) {
    NG ng(m_cosmology);
    for (size_t ii = 0; ii < Mass.size(); ii++)
      MF[ii] *= ng.MF_correction(MASS[ii], redshift, model_MF, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }
  else if (m_cosmology->fNL() != 0 && !default_delta)
    ErrorCBL("Non Gaussianity still not available for user-defined density contrast threshold!", "mass_function", "MassFunction.cpp");

  return MF;
}


// ================================================================================


vector<double> cbl::cosmology::MassFunction::mass_function_fR (const vector<double> Mass, const double redshift, const std::string model_MF, const double f_R0, const bool store_output, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const bool default_delta, const double delta_t) const
{
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  vector<double> MASS(Mass.size());
  for (size_t ii = 0; ii < Mass.size(); ii++)
    MASS[ii] = Mass[ii] * fact;

  const bool unit1 = true;

  cbl::Path path;
  string pathToMGCAMB = path.DirCosmo() + "/External/MGCAMB/";

  if (system(("cp " + pathToMGCAMB + "params_MG_backup.ini " + pathToMGCAMB + "params_MG.ini").c_str())) {}
  int MG_flag = (f_R0 == 0.) ? 0 : 3;
  string cosmoTag = (f_R0 == 0.) ? "LCDM" : "fR" + to_string(f_R0);

  string sed;
  sed = "sed -i 's/MG_flag = 0/MG_flag = " + cbl::conv(MG_flag, cbl::par::fINT) + "/g' " + pathToMGCAMB + "params_MG.ini";
  if (system(sed.c_str())) {}
  sed = "sed -i 's/F_R0 = 0.001d0/F_R0 = " + cbl::conv(f_R0, cbl::par::fDP8) + "d0/g' " + pathToMGCAMB + "params_MG.ini";
  if (system(sed.c_str())) {}

  Sigma SG(m_cosmology);
  
  vector<double> SSS = SG.sigma2M(MASS, "MGCAMB", redshift, store_output, cosmoTag, interpType, k_max, input_file, is_parameter_file, unit1);
  vector<double> Dln_Sigma = SG.dnsigma2M(1, MASS, "MGCAMB", redshift, store_output, cosmoTag, interpType, k_max, input_file, is_parameter_file, unit1);
  vector<double> MF(Mass.size(), pow(fact, 4.));

  const double RHO = m_cosmology->rho_m(0., true);
  double Sigma;

  for (size_t ii = 0; ii < Mass.size(); ii++) {
    Sigma = sqrt(SSS[ii]);
    Dln_Sigma[ii] *= (MASS[ii] / (2. * SSS[ii]));
    MF[ii] *= RHO / (MASS[ii] * MASS[ii]) * m_MF_generator(Sigma, Dln_Sigma[ii], redshift, 1., model_MF, Delta, default_delta, delta_t);
  }

  if (m_cosmology->fNL() != 0 && default_delta) {
    NG ng(m_cosmology);
    for (size_t ii = 0; ii < Mass.size(); ii++)
      MF[ii] *= ng.MF_correction(MASS[ii], redshift, model_MF, store_output, cosmoTag, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }
  else if (m_cosmology->fNL() != 0 && !default_delta)
    ErrorCBL("Non Gaussianity still not available for user-defined density contrast threshold!", "mass_function", "MassFunction.cpp");

  return MF;
}


// =====================================================================================


double cbl::cosmology::MassFunction::m_mass_function (const double Mass, std::shared_ptr<void> mass_function_params) const
{
  std::shared_ptr<glob::STR_MF> pp = static_pointer_cast<glob::STR_MF>(mass_function_params);

  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  double MASS = Mass * fact;
  const double RHO = m_cosmology->rho_m(0., true);

  const bool unit1 = true;

  Sigma SG(m_cosmology);
 
  double SSS = SG.sigma2M({MASS}, pp->method_SS, 0., pp->store_output, pp->output_root, pp->interpType, pp->k_max, pp->input_file, pp->is_parameter_file, unit1)[0];
  double Sigma = sqrt(SSS);
  double Dln_Sigma = SG.dnsigma2M(1, {MASS}, pp->method_SS, 0., pp->store_output, pp->output_root, pp->interpType, pp->k_max, pp->input_file, pp->is_parameter_file, unit1)[0] * (MASS / (2. * SSS));
  double DlnSigma_DlnR, RADIUS, MF;

  if (pp->model_MF=="Castro") {
      RADIUS = pow(3*MASS*0.25/cbl::par::pi/RHO, 1./3);
      DlnSigma_DlnR = SG.dnsigma2R(1, {RADIUS},  pp->method_SS, 0., pp->store_output, pp->output_root, pp->interpType, pp->k_max, pp->input_file, pp->is_parameter_file, unit1)[0]*(RADIUS / (2. * SSS));
      MF = RHO/(MASS*MASS)*DlnSigma_DlnR/(-3) * m_MF_generator(Sigma, DlnSigma_DlnR, pp->redshift, pp->model_MF, pp->Delta, pp->default_delta, pp->delta_t) * pow(fact, 4.);
    
  }
  else
    MF = RHO / (MASS * MASS) * m_MF_generator(Sigma, Dln_Sigma, pp->redshift, pp->model_MF, pp->Delta, pp->default_delta, pp->delta_t) * pow(fact, 4.);

  if (m_cosmology->fNL() != 0 && pp->default_delta) {
    NG ng(m_cosmology);
    MF *= ng.MF_correction(MASS, pp->redshift, pp->model_MF, pp->store_output, pp->output_root, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->input_file, pp->is_parameter_file);
  }
  else if (m_cosmology->fNL() != 0 && !pp->default_delta)
    ErrorCBL("Non Gaussianity still not available for user-defined density contrast threshold!", "m_mass_function", "MassFunction.cpp");

  return MF;
}


// =====================================================================================


double cbl::cosmology::MassFunction::mass_function_fast (const double Mass, const double redshift, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file) const
{
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  double MASS = Mass * fact;

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double MMass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dln_sigma;
  while (fin >> MMass >> Sigma >> Dln_Sigma) {
    mass.push_back(MMass);
    sigma.push_back(Sigma);
    dln_sigma.push_back(Dln_Sigma);
  }

  // ---------- compute the MF ----------

  double sig = interpolated(MASS, mass, sigma, "Steffen");
  double dlsig = interpolated(MASS, mass, dln_sigma, "Steffen");
  const double RHO = m_cosmology->rho_m(0., true);

  double MF = RHO / (MASS * MASS) * m_MF_generator(sig, dlsig, redshift, model_MF, Delta) * pow(fact, 4.);
  if (std::isnan(MF))
    ErrorCBL("MF = " + conv(MF, par::fDP3) + "!", "", "MassFunction.cpp");

  if (m_cosmology->fNL() != 0) {
    NG ng(m_cosmology);
    MF *= ng.MF_correction(MASS, redshift, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }

  return MF;
}


// =====================================================================================


double cbl::cosmology::MassFunction::mass_function (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const std::string model_MF, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) const
{
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  double MASS = Mass * fact;
  const double RHO = m_cosmology->rho_m(0., true);

  double MF;
  if (model_MF=="Castro"){
    cbl::cosmology::Sigma SG(m_cosmology);
    double RADIUS = pow(3*MASS*0.25/cbl::par::pi/RHO, 1./3);
    double DlnSigma_DlnR = SG.dnsigma2R(1, {RADIUS},  method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0]*(RADIUS / (2. * Sigma*Sigma));
    MF = RHO/(MASS*MASS)*DlnSigma_DlnR/(-3) * m_MF_generator(Sigma, DlnSigma_DlnR, redshift, model_MF, Delta) * pow(fact, 4.);
  } else
    MF = RHO / (MASS * MASS) * m_MF_generator(Sigma, Dln_Sigma, redshift, model_MF, Delta) * pow(fact, 4.);

  if (m_cosmology->fNL() != 0) {
    NG ng(m_cosmology);
    MF *= ng.MF_correction(MASS, redshift, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }
  
  return MF;
}


// =====================================================================================


double cbl::cosmology::MassFunction::mass_function (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const double D_N, const std::string model_MF, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) const
{ 
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  double MASS = Mass * fact;
  const double RHO = m_cosmology->rho_m(0., true);

  double MF;
  if (model_MF=="Castro"){
    cbl::cosmology::Sigma SG(m_cosmology);
    double RADIUS = pow(3*MASS*0.25/cbl::par::pi/RHO, 1./3);
    double DlnSigma_DlnR = SG.dnsigma2R(1, {RADIUS},  method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0]*(RADIUS / (2. * Sigma*Sigma));
    MF = RHO/(MASS*MASS)*DlnSigma_DlnR/(-3) * m_MF_generator(Sigma, DlnSigma_DlnR, redshift, D_N, model_MF, Delta) * pow(fact, 4.);
  } 
  else
    MF = RHO / (MASS * MASS) * m_MF_generator(Sigma, Dln_Sigma, redshift, D_N, model_MF, Delta) * pow(fact, 4.);

  if (m_cosmology->fNL() != 0) {
    NG ng(m_cosmology);
    MF *= ng.MF_correction(MASS, redshift, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }

  return MF;
}


// =====================================================================================


double cbl::cosmology::MassFunction::m_MF_generator (const double Sigma, const double Dln_Sigma, const double redshift, const std::string model_MF, const double Delta, const bool default_delta, const double delta_t) const
{
  double D_N = m_cosmology->DN(redshift);

  return m_MF_generator(Sigma, Dln_Sigma, redshift, D_N, model_MF, Delta, default_delta, delta_t);
}


// =====================================================================================


double cbl::cosmology::MassFunction::m_MF_generator (const double Sigma, const double Dln_Sigma, const double redshift, const double D_N, const std::string model_MF, const double Delta, const bool default_delta, const double delta_t) const
{
  const double deltacz = (default_delta) ? m_cosmology->deltac(redshift) : fabs(delta_t * D_N);
  const double sigmaz = Sigma * D_N;

  if (model_MF == "PS") // Press & Schechter
    return sqrt(2. / par::pi) * deltacz / sigmaz * fabs(Dln_Sigma) * exp(-(deltacz * deltacz) * 0.5 / (sigmaz * sigmaz));

  else if (model_MF == "ST") { // Sheth & Tormen
    const double aa = 0.707;
    const double qq = 0.3;
    const double AA = 0.3222;
    return AA * sqrt(2. * aa / par::pi) * deltacz / sigmaz * (1 + pow(sigmaz / (sqrt(aa) * deltacz), 2. * qq)) * fabs(Dln_Sigma) * exp(-(aa * deltacz * deltacz) * 0.5 / (sigmaz * sigmaz));
  }

  else if (model_MF == "Jenkins") // Jenkins et al. (2001)
    return 0.315 * exp(-pow(fabs(log(1. / sigmaz) + 0.61), 3.8)) * fabs(Dln_Sigma);

  else if (model_MF == "Warren") { // Warren et al. (2006)
    const double AA = 0.7234;
    const double aa = 1.625;
    const double bb = 0.2538;
    const double cc = 1.1982;
    return AA * (pow(sigmaz, -aa) + bb) * exp(-cc / (sigmaz * sigmaz)) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Reed") { // Reed et al. (2007)
    const double cc = 1.08;
    const double aa = 0.764 / cc;
    const double pp = 0.3;
    const double AA = 0.3222;
    const double Fatt = AA * sqrt(2. * aa / par::pi);
    const double n_eff = -6. * Dln_Sigma - 3.;
    const double G1 = exp(-pow(log(1. / sigmaz) - 0.4, 2) / 0.72);
    const double G2 = exp(-pow(log(1. / sigmaz) - 0.75, 2) / 0.08);
    return Fatt * deltacz / sigmaz * (1 + pow(sigmaz * sigmaz / (aa * deltacz * deltacz), pp) + 0.6 * G1 + 0.4 * G2) * fabs(Dln_Sigma) * exp(-(cc * aa * deltacz * deltacz) * 0.5 / (sigmaz * sigmaz) - (0.03 / pow(n_eff + 3., 2) * pow(deltacz / sigmaz, 0.6)));
  }

  else if (model_MF == "Pan") { // Pan (2007)
    const double alpha = 0.435; // check!!!
    return 4. * alpha / sqrt(2. * par::pi) * deltacz / pow(sigmaz, 2. * alpha) * fabs(Dln_Sigma) * exp(-(deltacz * deltacz) * 0.5 / pow(sigmaz, 4. * alpha));
  }

  else if (model_MF == "ShenH") { // halo MF, Shen et al. (2006) // check!!!
    const double alpha = -0.55;
    const double beta = -0.56;
    const double ni = pow(deltacz / sigmaz, 2);
    const double ni_fni = sqrt(ni / (2. * par::pi)) * exp(-ni * pow(1. - beta * pow(ni, alpha), 2) * 0.5) * (1. - beta / pow(ni, -alpha) * (1 + alpha + (alpha * (alpha + 1.)) * 0.5));
    return ni_fni * fabs(Dln_Sigma) * 2;
  }

  else if (model_MF == "ShenF") { // filament MF, Shen et al. (2006) // check!!!
    const double alpha = -0.28;
    const double beta = -0.012;
    const double ni = pow(deltacz / sigmaz, 2);
    const double ni_fni = sqrt(ni / (2. * par::pi)) * exp(-ni * pow(1. - beta * pow(ni, alpha), 2) * 0.5) * (1. - beta / pow(ni, -alpha) * (1 + alpha + (alpha * (alpha + 1.)) * 0.5));
    return ni_fni * fabs(Dln_Sigma) * 2;
  }

  else if (model_MF == "ShenS") { // sheet MF, Shen et al. (2006) // check!!!
    const double alpha = -0.61;
    const double beta = 0.45;
    const double ni = pow(deltacz / sigmaz, 2);
    const double ni_fni = sqrt(ni / (2. * par::pi)) * exp(-ni * pow(1. - beta * pow(ni, alpha), 2) * 0.5) * (1. - beta / pow(ni, -alpha) * (1 + alpha + (alpha * (alpha + 1.)) * 0.5));
    return fabs(ni_fni * fabs(Dln_Sigma) * 2); // check!!!
  }

  else if (model_MF == "Tinker") { // Tinker et al. (2008)

    //if (redshift>2) WarningMsgCBL("the Tinker mass function has been tested for z<~2!", "m_MF_generator", "MassFunction.cpp");

    double A0, a0, b0, c0;

    if (Delta == 200) {
      A0 = 0.186;
      a0 = 1.47;
      b0 = 2.57;
      c0 = 1.19;
    }
    else if (Delta == 300) {
      A0 = 0.200;
      a0 = 1.52;
      b0 = 2.25;
      c0 = 1.27;
    }
    else if (Delta == 400) {
      A0 = 0.212;
      a0 = 1.56;
      b0 = 2.05;
      c0 = 1.34;
    }
    else if (Delta == 600)  {
      A0 = 0.218;
      a0 = 1.61;
      b0 = 1.87;
      c0 = 1.45;
    }
    else if (Delta == 800) {
      A0 = 0.248;
      a0 = 1.87;
      b0 = 1.59;
      c0 = 1.58;
    }
    else if (Delta == 1200) {
      A0 = 0.255;
      a0 = 2.13;
      b0 = 1.51;
      c0 = 1.80;
    }
    else if (Delta == 1600) {
      A0 = 0.260;
      a0 = 2.30;
      b0 = 1.46;
      c0 = 1.97;
    }
    else if (Delta == 2400) {
      A0 = 0.260;
      a0 = 2.53;
      b0 = 1.44;
      c0 = 2.24;
    }
    else if (Delta == 3200) {
      A0 = 0.260;
      a0 = 2.66;
      b0 = 1.41;
      c0 = 2.44;
    }
    else {
      A0 = (Delta < 1600) ? 0.1 * log10(Delta) - 0.05 : 0.26;
      a0 = 1.43 + pow(max(log10(Delta) - 2.3, 0.), 1.5); // check!!!
      b0 = 1. + pow(log10(Delta) - 1.6, -1.5);
      c0 = 1.2 + pow(max(log10(Delta) - 2.35, 0.), 1.6); // check!!!
    }

    const double alpha = pow(10., -pow(max(0.75 / log10(Delta / 75.), 0.), 1.2)); // check!!!
    const double AA = A0 * pow(1. + redshift, -0.14);
    const double aa = a0 * pow(1. + redshift, -0.06);
    const double bb = b0 * pow(1. + redshift, -alpha);
    const double cc = c0;

    return AA * (pow(sigmaz / bb, -aa) + 1.) * exp(-cc / (sigmaz * sigmaz)) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Tinker08_Interp") { // Tinker et al. (2008)

    //if (redshift>2) WarningMsgCBL("the Tinker mass function has been tested for z<~2!", "m_MF_generator", "MassFunction.cpp");

    vector<double> _Delta = {200, 300, 400, 600, 800, 1200, 1600, 2400, 3200};
    vector<double> _A0 = {1.858659e-01, 1.995973e-01, 2.115659e-01, 2.184113e-01, 2.480968e-01, 2.546053e-01, 2.600000e-01, 2.600000e-01, 2.600000e-01};
    vector<double> _a0 = {1.466904, 1.521782, 1.559186, 1.614585, 1.869936, 2.128056, 2.301275, 2.529241, 2.661983};
    vector<double> _b0 = {2.571104, 2.254217, 2.048674, 1.869559, 1.588649, 1.507134, 1.464374, 1.436827, 1.405210};
    vector<double> _c0 = {1.193958, 1.270316, 1.335191, 1.446266, 1.581345, 1.795050, 1.965613, 2.237466, 2.439729};

    cbl::glob::FuncGrid interp_A0(_Delta, _A0, "Spline");
    cbl::glob::FuncGrid interp_a0(_Delta, _a0, "Spline");
    cbl::glob::FuncGrid interp_b0(_Delta, _b0, "Spline");
    cbl::glob::FuncGrid interp_c0(_Delta, _c0, "Spline");

    const double alpha = pow(10., -pow(max(0.75 / log10(Delta / 75.), 0.), 1.2)); // check!!!
    const double AA = interp_A0(Delta) * pow(1. + redshift, -0.14);
    const double aa = interp_a0(Delta) * pow(1. + redshift, -0.06);
    const double bb = interp_b0(Delta) * pow(1. + redshift, -alpha);
    const double cc = interp_c0(Delta);

    return AA * (pow(sigmaz / bb, -aa) + 1.) * exp(-cc / (sigmaz * sigmaz)) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Crocce") { // Crocce et al. (2010)
    const double AA = 0.58 * pow(1. + redshift, -0.13);
    const double aa = 1.37 * pow(1. + redshift, -0.15);
    const double bb = 0.3 * pow(1. + redshift, -0.084);
    const double cc = 1.036 * pow(1. + redshift, -0.024);
    return AA * (pow(sigmaz, -aa) + bb) * exp(-cc * pow(sigmaz, -2)) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Angulo_FOF") // FoF MF, Angulo et al. (2012)
    return 0.201 * (pow(2.08 / sigmaz, 1.7) + 1) * exp(-1.172 / (sigmaz * sigmaz)) * fabs(Dln_Sigma);

  else if (model_MF == "Angulo_Sub") // SUBFIND MF, Angulo et al. (2012)
    return 0.265 * (pow(1.675 / sigmaz, 1.9) + 1) * exp(-1.4 / (sigmaz * sigmaz)) * fabs(Dln_Sigma);

  else if (model_MF == "Bhattacharya") { // Bhattacharya et al. (2011) // check!!!
    const double aa = 0.788 / pow(1 + redshift, 0.01);
    const double qq = 1.795;
    const double AA = 0.333 / pow(1 + redshift, 0.11);
    const double pp = 0.807;
    const double nu = deltacz / sigmaz;
    //if (redshift>=2) aa = 0.799/pow(1+redshift,0.024);
    return fabs(Dln_Sigma) * AA * sqrt(2 / par::pi) * exp(-nu * nu * 0.5 * aa) * (1 + pow(1 / (nu * nu * aa), pp)) * pow(nu * sqrt(aa), qq);
  }

  else if (model_MF == "Courtin") { // Courtin et al. (2010)
    const double aa = 0.695;
    const double qq = 0.1;
    const double AA = 0.348;
    return AA * sqrt(2. * aa / par::pi) * deltacz / sigmaz * (1 + pow(sigmaz / (sqrt(aa) * deltacz), 2. * qq)) * fabs(Dln_Sigma) * exp(-(aa * deltacz * deltacz) * 0.5 / (sigmaz * sigmaz));
  }

  else if (model_MF == "Manera") { // Manera et al. (2010)

    double aa, qq;
    if (fabs(redshift - 0.) < 1.e-30) {
      aa = 0.709;
      qq = 0.248;
    }
    else if (fabs(redshift - 0.5) < 1.e-30) {
      aa = 0.724;
      qq = 0.241;
    }
    else
      ErrorCBL("the Manera et al. mass function has been tested only for z=0 and 0.5!", "m_MF_generator", "MassFunction.cpp");

    return 0.3222 * sqrt(2. * aa / par::pi) * deltacz / sigmaz * (1 + pow(sigmaz / (sqrt(aa) * deltacz), 2. * qq)) * fabs(Dln_Sigma) * exp(-(aa * deltacz * deltacz) * 0.5 / (sigmaz * sigmaz));
  }

  else if (model_MF == "Watson_FOF") { // FoF MF, Watson et al. (2012)
    const double a = 0.282;
    const double b = 1.406;
    const double c = 2.163;
    const double d = 1.210;
    return a * (pow((b / sigmaz), c) + 1) * exp(-d / (sigmaz * sigmaz)) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Watson_SOH") { // Spherical Overdensity halo MF, Watson et al. (2012)
    const double Omega = m_cosmology->Omega_matter(redshift);
    const double C = exp(0.023 * (Delta / 178 - 1));
    const double d = -0.456 * Omega - 0.139;
    const double p = 0.072;
    const double q = 2.13;
    const double Gamma = C * pow((Delta / 178), d) * exp(p * (1 - Delta / 178) / pow(sigmaz, q));
    double A = 0.194;
    double alpha = 2.267;
    double beta = 1.805;
    double gamma = 1.287;
    if (redshift > 0 && redshift < 6) {
      A = Omega * (1.097 * pow((1 + redshift), -3.216) + 0.074);
      alpha = Omega * (3.136 * pow((1 + redshift), -3.056) + 2.349);
      beta = Omega * (5.907 * pow((1 + redshift), -3.599) + 2.344);
      gamma = 1.318;
    }
    if (redshift >= 6) {
      A = 0.563;
      alpha = 0.874;
      beta = 3.810;
      gamma = 1.453;
    }
    return fabs(Dln_Sigma) * Gamma * A * (pow((beta / sigmaz), alpha) + 1.) * exp(-gamma / (sigmaz * sigmaz));
  }

  else if (model_MF == "Peacock") { // Peacock et al. (2007)
    const double aa = 1.529;
    const double bb = 0.704;
    const double cc = 0.412;
    const double nu = deltacz / sigmaz;
    return nu * exp(-cc * nu * nu) / pow(1 + aa * pow(nu, bb), 2) * (bb * aa * pow(nu, bb - 1) + 2 * cc * nu * (1 + aa * pow(nu, bb))) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Despali_Z0") {  // Despali et al. 2016
    const double aa = 0.794; // 0.794 ± 0.005, z=0 fit
    const double pp = 0.247; // 0.247 ± 0.009, z=0 fit
    const double AA = 0.333; // 0.333 ± 0.001, z=0 fit
    const double nu = pow(deltacz / sigmaz, 2);
    const double nup = aa * nu;
    return 2. * AA * (1. + pow(nup, -pp)) * sqrt(nup / (2. * par::pi)) * exp(-0.5 * nup) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Despali_AllZ") { // Despali et al. 2016
    const double aa = 0.7663; // 0.7663 ± 0.0013, all z fit
    const double pp = 0.2579; // 0.2579 ± 0.0026, all z fit
    const double AA = 0.3298; // 0.3298 ± 0.0003, all z fit
    const double nu = pow(deltacz / sigmaz, 2);
    const double nup = aa * nu;
    return 2. * AA * (1. + pow(nup, -pp)) * sqrt(nup / (2. * par::pi)) * exp(-0.5 * nup) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Despali_AllZAllCosmo") { // Despali et al. 2016
    const double aa = 0.7689; // 0.7689 ± 0.0011, all z fit
    const double pp = 0.2536; // 0.2536 ± 0.0026, all z fit
    const double AA = 0.3295; // 0.3295 ± 0.0003, all z fit
    const double nu = pow(deltacz / sigmaz, 2);
    const double nup = aa * nu;
    return 2. * AA * (1. + pow(nup, -pp)) * sqrt(nup / (2. * par::pi)) * exp(-0.5 * nup) * fabs(Dln_Sigma);
  }

  else if (model_MF == "Despali_HighM") {  // Despali et al. 2016
    const double aa = 0.8199; // 0.8199 ± 0.0010, M>3e13
    const double pp = 0.0;    // 0, M>3e13
    const double AA = 0.3141; // 0.3141 ± 0.0006, M>3e13
    const double nu = pow(deltacz / sigmaz, 2);
    const double nup = aa * nu;
    return 2. * AA * (1. + pow(nup, -pp)) * sqrt(nup / (2. * par::pi)) * exp(-0.5 * nup) * fabs(Dln_Sigma);
  }

  else if (model_MF=="Castro") { // Castro et al. (2023) 

    cbl::WarningMsgCBL("Check! The implemented MF differs for more than 1% from the original one at z>0.5 and M>10^14 solar masses. The difference is due to the integrand of the function sigma2M", "m_MF_generator", "MassFunction.cpp");  
    
    //values for ROCKSTAR
    const double a1 = 0.7962217;
    const double a2 = 0.1448873;
    const double az = -0.06577418;
    const double p1 = -0.5612250;
    const double p2 = -0.4742516;
    const double q1 = 0.3688495;
    const double q2 = -0.2803717;
    const double qz = 0.02514857;

    const double dlnsigma_dlnR = Dln_Sigma;
    //here Dln_Sigma is the derivative with respect to R
    const double ar = a1+a2*pow((dlnsigma_dlnR+0.6125), 2);
    const double qr = q1+q2*(dlnsigma_dlnR+0.5);
    const double pp = p1+p2*(dlnsigma_dlnR+0.5);

    const double aa = ar*pow(m_cosmology->Omega_matter(redshift), az);
    const double qq = qr*pow(m_cosmology->Omega_matter(redshift), qz);

    const double nu = deltacz/sigmaz;
    const double AA = 1./(pow(2, -0.5-pp+qq*0.5)/sqrt(par::pi)*(pow(2, pp)*std::tgamma(qq*0.5)+std::tgamma(-pp+qq*0.5)));
    const double BB = sqrt(2*aa*nu*nu/par::pi)*exp(-nu*nu*0.5*aa)*(1+1./pow((nu*nu*aa), pp))*pow(nu*sqrt(aa), qq-1.);
    return AA*BB;
  }

  else
    return ErrorCBL("model_MF not allowed!", "m_MF_generator", "MassFunction.cpp");
}


// =====================================================================================


double cbl::cosmology::MassFunction::n_haloes (const double Mass_min, const double Mass_max, const double z_min, const double z_max, const bool angle_rad, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
 
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma) {
    if (Mass_min < Mass && Mass < Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }

  // ---------- compute the number of haloes ----------

  int step_z = 100;
  double delta_z = (z_max - z_min) / step_z;
  double redshift = z_min;
  double N_haloes = 0.;

  for (int Red = 0; Red < step_z; Red++) {
    double Int = 0.;
    for (unsigned int k = 0; k < mass.size() - 1; k++)
      Int += mass_function(mass[k], sigma[k], dlnsigma[k], redshift, model_MF, store_output, output_root, Delta) * m_cosmology->dV_dZdOmega(redshift, angle_rad) * (mass[k + 1] - mass[k]);
    
    N_haloes += Int * delta_z;
    
    redshift += delta_z;
  }

  return N_haloes;
}


// =====================================================================================


double cbl::cosmology::MassFunction::n_haloes (const double Mass_min, const double Mass_max, const double Volume, const double redshift, const std::string model_MF, const std::string method_SS, const int nbin_mass, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string input_file, const bool is_parameter_file, const bool default_delta, const double delta_t) const
{
  glob::STR_MF mass_function_par;

  mass_function_par.redshift = redshift;
  mass_function_par.model_MF = model_MF;
  mass_function_par.method_SS = method_SS;
  mass_function_par.store_output = store_output;
  mass_function_par.output_root = output_root;
  mass_function_par.Delta = Delta;
  mass_function_par.interpType = interpType;
  mass_function_par.norm = norm;
  mass_function_par.k_min = k_min;
  mass_function_par.k_max = k_max;
  mass_function_par.prec = prec;
  mass_function_par.input_file = input_file;
  mass_function_par.is_parameter_file = is_parameter_file;
  mass_function_par.default_delta = default_delta;
  mass_function_par.delta_t = delta_t;

  auto pp = make_shared<glob::STR_MF>(mass_function_par);

  if (nbin_mass == 0) {
    function<double(double)> func = bind(&MassFunction::m_mass_function, this, std::placeholders::_1, pp);
    return Volume * wrapper::gsl::GSL_integrate_qag(func, Mass_min, Mass_max, prec);
  }

  else if (nbin_mass > 0) {
    vector<double> mass = logarithmic_bin_vector(nbin_mass, Mass_min, Mass_max);
    vector<double> dndm(nbin_mass);

    for (int i = 0; i < nbin_mass; i++)
      dndm[i] = m_mass_function(mass[i], pp);

    glob::FuncGrid mass_function_interpolator(mass, dndm, interpType);

    return Volume * mass_function_interpolator.integrate_qag(Mass_min, Mass_max, prec);
  }

  else
    ErrorCBL("nbinmass must be >=0", "n_haloes", "MassFunction.cpp");

  return 0;
}


// =====================================================================================


double cbl::cosmology::MassFunction::MhaloMin (const int n_halo, const double Area, const bool angle_rad, const double z_min, const double z_max, const double Mmax, const double lgM1_guess, const double lgM2_guess, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{ 
  function<double(double)> func_MhaloMin = [&] (const double lgMmin) 
  {
    double Mmin = pow(10., lgMmin);

    double n_halo_expected = n_haloes(Mmin, Mmax, z_min, z_max, angle_rad, model_MF, method_SS, store_output, output_root, Delta, interpType, k_max, input_file, is_parameter_file)*Area;
   
    return n_halo-n_halo_expected;
  };

  return pow(10., wrapper::gsl::GSL_root_brent(func_MhaloMin, 0., lgM1_guess, lgM2_guess, 0.0001));

}


// =====================================================================================


double cbl::cosmology::MassFunction::unevolved_mass_function (const double mass_accr) const
{
  const double aa = -0.8;
  const double mm = mass_accr / aa; // mass_accr = m/M
  const double yn = log(0.3);       // log(0.21) check!!!

  return aa * log(-mm) + log(exp(2. * par::pi * pow(mm, 3.))) + yn;
}


// =====================================================================================


double cbl::cosmology::MassFunction::n_haloes_selection_function (const double Mass_min, const double Mass_max, const double z_min, const double z_max, const bool angle_rad, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> cols, const bool store_output, const std::string output_root, const double Delta, const bool isDelta_critical, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma) {
    if (Mass_min < Mass && Mass < Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }
  fin.clear();
  fin.close();

  glob::FuncGrid interp_sigma(mass, sigma, "Spline");
  glob::FuncGrid interp_dlnsigma(mass, dlnsigma, "Spline");

  // ---------- read the selection function ----------

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);
  function<double(const double, const double)> interp_SF = bind(interpolated_2D, std::placeholders::_1, std::placeholders::_2, mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the number of haloes ----------

  vector<vector<double>> limits = {{Mass_min, Mass_max}, {z_min, z_max}};

  int step_z = 100;
  double delta_z = (z_max - z_min) / step_z;
  double redshift = z_min;
  double N_haloes = 0.;
  //double norm_SF = 0.;

  for (int Red = 0; Red < step_z; Red++) {
    double Int = 0.;
    for (unsigned int k = 0; k < mass.size() - 1; k++) {
      double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
      double SF = interp_SF(mass[k], redshift);
      Int += SF * mass_function(mass[k], sigma[k], dlnsigma[k], redshift, model_MF, store_output, output_root, DD) * m_cosmology->dV_dZdOmega(redshift, angle_rad) * (mass[k + 1] - mass[k]);
    }

    N_haloes += Int * delta_z;

    redshift += delta_z;
  }

  return N_haloes;
}


// =====================================================================================


std::vector<double> cbl::cosmology::MassFunction::mass_function (const std::vector<double> mass, const double z_min, const double z_max, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const bool isDelta_critical, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const 
{
  vector<double> MF(mass.size());

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> _mass, _sigma, _dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma)
    {
      _mass.push_back(Mass);
      _sigma.push_back(Sigma);
      _dlnsigma.push_back(Dln_Sigma);
    }
  fin.clear();
  fin.close();

  glob::FuncGrid interp_sigma(_mass, _sigma, "Spline");
  glob::FuncGrid interp_dlnsigma(_mass, _dlnsigma, "Spline");

  // ---------- compute the number of haloes ----------

  //int step_z = 100;
  //double delta_z = (z_max-z_min)/step_z;

  double VV = m_cosmology->Volume(z_min, z_max, 1.);

  for (unsigned int k = 0; k < mass.size(); k++)
    {
      double sigma = interp_sigma(mass[k]);
      double dlnsigma = interp_dlnsigma(mass[k]);

      auto func = [&](const double redshift) {
	double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
	return mass_function(mass[k], sigma, dlnsigma, redshift, model_MF, store_output, output_root, DD) * m_cosmology->dV_dZdOmega(redshift, false);
      };
      MF[k] = wrapper::gsl::GSL_integrate_qag(func, z_min, z_max) / VV;
    }

  return MF;
}


// =====================================================================================


std::vector<double> cbl::cosmology::MassFunction::mass_function_selection_function_vector (const std::vector<double> mass, const double z_min, const double z_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> cols, const bool store_output, const std::string output_root, const double Delta, const bool isDelta_critical, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{
  vector<double> MF(mass.size());

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> _mass, _sigma, _dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma)
    {
      _mass.push_back(Mass);
      _sigma.push_back(Sigma);
      _dlnsigma.push_back(Dln_Sigma);
    }
  fin.clear();
  fin.close();

  glob::FuncGrid interp_sigma(_mass, _sigma, "Spline");
  glob::FuncGrid interp_dlnsigma(_mass, _dlnsigma, "Spline");

  // ---------- read the selection function ----------

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);
  function<double(const double, const double)> interp_SF = bind(interpolated_2D, std::placeholders::_1, std::placeholders::_2, mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the number of haloes ----------

  //int step_z = 100;
  //double delta_z = (z_max-z_min)/step_z;

  double VV = m_cosmology->Volume(z_min, z_max, 1.);

  for (unsigned int k = 0; k < mass.size(); k++)
    {
      double sigma = interp_sigma(mass[k]);
      double dlnsigma = interp_dlnsigma(mass[k]);

      auto func = [&](const double redshift) {
	double SF = interp_SF(mass[k], redshift);
	double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
	return SF * mass_function(mass[k], sigma, dlnsigma, redshift, model_MF, store_output, output_root, DD) * m_cosmology->dV_dZdOmega(redshift, false);
      };

      MF[k] = wrapper::gsl::GSL_integrate_qag(func, z_min, z_max) / VV; ///wrapper::gsl::GSL_integrate_qag(funcD, z_min, z_max);
    }

  return MF;
}


// =====================================================================================


std::vector<double> cbl::cosmology::MassFunction::redshift_distribution_haloes (const double z_min, const double z_max, const int step_z, const double Area_degrees, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const bool store_output, const std::string output_root, const double Delta, const bool isDelta_critical, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma) {
    if (Mass_min < Mass && Mass < Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }
  fin.clear();
  fin.close();

  glob::FuncGrid interp_sigma(mass, sigma, "Spline");
  glob::FuncGrid interp_dlnsigma(mass, dlnsigma, "Spline");

  // ---------- compute the number of haloes ----------

  vector<double> redshift_distribution(step_z, 0);
  double delta_z = (z_max - z_min) / step_z;
  double redshift = z_min;

  for (int Red = 0; Red < step_z; Red++) {
    vector<vector<double>> integration_limits(2);
    integration_limits[0] = {redshift, redshift + delta_z};
    integration_limits[1] = {Mass_min, Mass_max};

    auto integrand = [&](const vector<double> parameters) {
      double redshift = parameters[0];
      double mass = parameters[1];

      double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
      double sigma = interp_sigma(mass);
      double dlnsigma = interp_dlnsigma(mass);
      return Area_degrees * m_cosmology->dV_dZdOmega(redshift, false) * mass_function(mass, sigma, dlnsigma, redshift, model_MF, store_output, output_root, DD);
    };

    cbl::wrapper::cuba::CUBAwrapper CubaIntegrand(integrand, 2);

    redshift_distribution[Red] = CubaIntegrand.IntegrateVegas(integration_limits);
    redshift += delta_z;
  }

  return redshift_distribution;
}


// =====================================================================================


std::vector<double> cbl::cosmology::MassFunction::redshift_distribution_haloes_selection_function (const std::vector<double> redshift, const double Area_degrees, const double Mass_min, const double Mass_max, const string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> cols, const bool store_output, const std::string output_root, const double Delta, const bool isDelta_critical, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma) {
    if (Mass_min < Mass && Mass < Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }
  fin.clear();
  fin.close();

  glob::FuncGrid interp_sigma(mass, sigma, "Spline");
  glob::FuncGrid interp_dlnsigma(mass, dlnsigma, "Spline");

  // ---------- read the selection function ----------

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);
  //function<double(const double, const double)> interp_SF = bind(interpolated_2D, std::placeholders::_1, std::placeholders::_2, mass_SF, redshift_SF, selection_function, "Linear");

  glob::FuncGrid2D interp_SF(mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the number of haloes ----------

  vector<double> redshift_distribution(redshift.size(), 0.);

  for (size_t zz = 0; zz < redshift.size() - 1; ++zz) {

    vector<vector<double>> integration_limits(2);
    integration_limits[0] = {redshift[zz], redshift[zz + 1]};
    integration_limits[1] = {Mass_min, Mass_max};

    auto integrand = [&](const vector<double> parameters) {
      const double redshift = parameters[0];
      const double mass = parameters[1];

      const double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
      const double sigma = interp_sigma(mass);
      const double dlnsigma = interp_dlnsigma(mass);
      const double SF = interp_SF(mass, redshift);

      return SF * Area_degrees * m_cosmology->dV_dZdOmega(redshift, false) * mass_function(mass, sigma, dlnsigma, redshift, model_MF, store_output, output_root, DD);
    };

    wrapper::cuba::CUBAwrapper CubaIntegrand(integrand, 2);
    const double distr = CubaIntegrand.IntegrateVegas(integration_limits);
    redshift_distribution[zz] = (distr != distr) ? 0. : distr;
  }

  return redshift_distribution;
}


// =====================================================================================


double cbl::cosmology::MassFunction::mean_redshift_haloes_selection_function (const double z_min, const double z_max, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> cols, const bool store_output, const std::string output_root, const double Delta, const bool isDelta_critical, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file) const
{

  // ---------- read the grid file ----------

  Sigma SG(m_cosmology);
  
  const string file_grid = SG.create_grid_sigmaM(method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  ifstream fin(file_grid.c_str());
  checkIO(fin, file_grid);

  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  while (fin >> Mass >> Sigma >> Dln_Sigma)
    {
      if (Mass_min < Mass && Mass < Mass_max)
	{
	  mass.push_back(Mass);
	  sigma.push_back(Sigma);
	  dlnsigma.push_back(Dln_Sigma);
	}
    }
  fin.clear();
  fin.close();

  glob::FuncGrid interp_sigma(mass, sigma, "Spline");
  glob::FuncGrid interp_dlnsigma(mass, dlnsigma, "Spline");

  // ---------- read the selection function ----------

  vector<double> mass_SF, redshift_SF;
  vector<vector<double>> selection_function;

  read_matrix(selection_function_file, mass_SF, redshift_SF, selection_function, cols);

  glob::FuncGrid2D interp_SF(mass_SF, redshift_SF, selection_function, "Linear");

  // ---------- compute the mean redshift ----------

  auto integrand_num = [&](const vector<double> parameters) {
    double redshift = parameters[0];
    double mass = parameters[1];
    double sigma = interp_sigma(mass);
    double dlnsigma = interp_dlnsigma(mass);
    double SF = interp_SF(mass, redshift);

    double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
    double MF = mass_function(mass, sigma, dlnsigma, redshift, model_MF, store_output, output_root, DD);
    return redshift * MF * SF * m_cosmology->dV_dZdOmega(redshift, false);
  };

  auto integrand_denom = [&](const vector<double> parameters) {
    double redshift = parameters[0];
    double mass = parameters[1];
    double sigma = interp_sigma(mass);
    double dlnsigma = interp_dlnsigma(mass);
    double SF = interp_SF(mass, redshift);

    double DD = (isDelta_critical) ? Delta / m_cosmology->Omega_matter(redshift) : Delta;
    double MF = mass_function(mass, sigma, dlnsigma, redshift, model_MF, store_output, output_root, DD);
    return MF * SF * m_cosmology->dV_dZdOmega(redshift, false);
  };

  cbl::wrapper::cuba::CUBAwrapper IntegrandNum(integrand_num, 2);
  cbl::wrapper::cuba::CUBAwrapper IntegrandDenom(integrand_denom, 2);

  vector<vector<double>> limits = {{z_min, z_max}, {Mass_min, Mass_max}};

  return IntegrandNum.IntegrateVegas(limits) / IntegrandDenom.IntegrateVegas(limits);
}


// =====================================================================================


double cbl::cosmology::MassFunction::converted_mass (const double mass, const std::shared_ptr<cbl::cosmology::Cosmology> assumed_cosmology, const double redshift, const double redshift_source) const
{
  const double Dds2 = assumed_cosmology->D_A(redshift, redshift_source);
  const double Dds1 = m_cosmology->D_A(redshift, redshift_source);

  const double Ds2 = assumed_cosmology->D_A(redshift);
  const double Ds1 = m_cosmology->D_A(redshift);

  const double H2 = assumed_cosmology->Hubble(redshift);
  const double H1 = m_cosmology->Hubble(redshift);

  return mass * pow(Dds2 / Ds2, 1.5) * H2 * pow(Dds1 / Ds1, -1.5) / H1;
}


// =====================================================================================


vector<double> cbl::cosmology::MassFunction::mass_function (const std::vector<double> Mass, const std::vector<double> Sigma, const std::vector<double> Dln_Sigma, const double redshift, const std::string model_MF, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) const 
{
  double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  double D_N = m_cosmology->DN(redshift);
  const double RHO = m_cosmology->rho_m(0., true);

  vector<double> MASS(Mass.size(), 0), MF(Mass.size(), 0), RADIUS(Mass.size(), 0), DlnSigma_DlnR(Mass.size(), 0);
  
  for (size_t i = 0; i < Mass.size(); i++){ 
    MASS[i] = Mass[i] * fact;
    if (model_MF=="Castro")
      RADIUS[i] = pow(3*MASS[i]*0.25/cbl::par::pi/RHO, 1./3);
  }
  
  if (model_MF=="Castro") {
    cbl::cosmology::Sigma SG(m_cosmology);
    DlnSigma_DlnR = SG.dnsigma2R(1, RADIUS,  method_SS, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file); 
    for (size_t i = 0; i < Mass.size(); i++) {
      DlnSigma_DlnR[i]*= (RADIUS[i] / (2. * Sigma[i]*Sigma[i]));
      MF[i] = RHO/(MASS[i]*MASS[i])*DlnSigma_DlnR[i]/(-3) * m_MF_generator(Sigma[i], DlnSigma_DlnR[i], redshift, model_MF, Delta) * pow(fact, 4.);
    }
  }
  else
    for (size_t i = 0; i < Mass.size(); i++)
      MF[i] = RHO / (MASS[i] * MASS[i]) * m_MF_generator(Sigma[i], Dln_Sigma[i], redshift, D_N, model_MF, Delta) * pow(fact, 4.);
  
  
  if (m_cosmology->fNL() != 0) {
  NG ng(m_cosmology);
  for (size_t i = 0; i < Mass.size(); i++)
    MF[i] *= ng.MF_correction(MASS[i], redshift, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);
  }
  
  return MF;
}
