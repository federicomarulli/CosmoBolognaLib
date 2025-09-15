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
 *  @file Cosmology/Lib/Sigma.cpp
 *
 *  @brief Methods of the class Cosmology used to model the amplitude
 *  of the matter power spectrum
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the rms fluctuations in the matter mass
 *  density
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Sigma.h"

using namespace std;
using namespace cbl;


// =====================================================================================


vector<double> cbl::cosmology::Sigma::m_func_sigma (const function<double(double)> pkInterp, const function<double(double,int)> filter, const vector<double> radius, const double k_max, const bool unit1) const
{
  double fact = (m_cosmology->unit() || unit1) ? 1. : m_cosmology->little_h();

  vector<double> sigma(radius.size(), pow(fact, -3.)/(2.*pow(par::pi, 2)));
  
  for (size_t index=0; index<radius.size(); index++) {
    
    auto ff = [&] (const double kk)
      {
	return pkInterp(kk/fact)*kk*kk*filter(kk,index);
      };
    
    sigma[index] *= wrapper::gsl::GSL_integrate_qag(ff, 1.e-4, k_max, 1.e-3);
  }

  return sigma;
}


// =====================================================================================


vector<double> cbl::cosmology::Sigma::m_sigma2R (const vector<double> radius, const function<double(double)> pkInterp, const double k_max, const double fact, const bool unit1) const
{
  function<double(double,int)> filter = [&] (const double k, const int index)
  {  
    return pow(TopHat_WF(k*radius[index]), 2)*fact;
  };

  return cosmology::Sigma::m_func_sigma(pkInterp, filter, radius, k_max, unit1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::sigma2R (const std::vector<double> radius, const std::function<double(double)> pk0Interp, const double k_max, const bool unit1) const
{
  for (auto &&rr : radius)
    if (rr<0.)
      ErrorCBL("the radius must be >0!", "dnsigma2R", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  if (m_cosmology->sigma8()>0.)
    // sigma_8 = sigma(8Mpc/h)
    fact = pow(m_cosmology->sigma8(), 2)/m_sigma2R({8.}, pk0Interp, k_max, 1., true)[0]; // normalization factor

  return m_sigma2R(radius, pk0Interp, k_max, fact, unit1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::sigma2R (const std::vector<double> radius, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file, const bool unit1) const
{
  for (auto &&rr : radius)
    if (rr<0.)
      ErrorCBL("the radius must be >0!", "sigma2R", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  PkXi PX(m_cosmology);
  
  auto Pk = PX.Pk_interp_for_sigma(method_Pk, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  // if the power spectrum is read from file or m_cosmology->sigma8()<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2R=sigma2R_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file==par::defaultString || is_parameter_file) {
    if (m_cosmology->sigma8()>0.) {
      auto Pk_0 = (redshift == 0.) ? Pk : PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
      // sigma_8 = sigma(8Mpc/h)
      fact = pow(m_cosmology->sigma8(), 2)/m_sigma2R({8.}, Pk_0, k_max, 1., true)[0]; // normalization factor
    }
  }

  return m_sigma2R(radius, Pk, k_max, fact, unit1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::dnsigma2R (const int nd, const std::vector<double> radius, std::function<double(double)> pk0Interp, const double k_max, const bool unit1) const
{
  for (auto &&rr : radius)
    if (rr < 0)
      ErrorCBL("the radius must be >0!", "dnsigma2R", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  if (m_cosmology->sigma8()>0.)
    // sigma_8 = sigma(8Mpc/h)
    fact = pow(m_cosmology->sigma8(), 2)/m_sigma2R({8.}, pk0Interp, k_max, 1., true)[0]; // normalization factor

  if (nd == 1) {
    
    function<double(double,int)> filter = [&](const double k, const int index) {
      return cbl::TopHat_WF(k*radius[index])*cbl::TopHat_WF_D1(k*radius[index])*k*2.*fact;
    };

    return cbl::cosmology::Sigma::m_func_sigma(pk0Interp, filter, radius, k_max, unit1);
  }

  else {
    ErrorCBL("", "dnsigma2R", "Sigma.cpp", glob::ExitCode::_workInProgress_);
    return vector<double>();
  }
}


// =====================================================================================

std::vector<double> cbl::cosmology::Sigma::dnsigma2R (const int nd, const std::vector<double> radius, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file, const bool unit1) const
{
  for (auto &&rr : radius)
    if (rr < 0)
      ErrorCBL("the radius must be >0!", "dnsigma2R", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  PkXi PX(m_cosmology);
  
  function<double(double)> Pk = PX.Pk_interp_for_sigma(method_Pk, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  // if the power spectrum is read from file or m_cosmology->sigma8()<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2R=sigma2R_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file==par::defaultString || is_parameter_file) {
    if (m_cosmology->sigma8()>0.) {
      auto Pk_0 = (redshift == 0.) ? Pk : PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
      // sigma_8 = sigma(8Mpc/h)
      fact = pow(m_cosmology->sigma8(), 2)/m_sigma2R({8.}, Pk_0, k_max, 1., true)[0]; // normalization factor
    }
  }

  if (nd==1) {
    
    function<double(double,int)> filter = [&](const double k, const int index) {
      return cbl::TopHat_WF(k*radius[index])*cbl::TopHat_WF_D1(k*radius[index])*k*2.*fact;
    };

    return cbl::cosmology::Sigma::m_func_sigma(Pk, filter, radius, k_max, unit1);
  }

  else {
    ErrorCBL("", "dnsigma2R", "Sigma.cpp", glob::ExitCode::_workInProgress_);
    return vector<double>();
  }
}


// =====================================================================================


vector<double> cbl::cosmology::Sigma::m_sigma2M (const vector<double> mass, const function<double(double)> pkInterp, const double k_max, const double fact, const bool unit1) const
{
  vector<double> radius(mass.size());
  
  for (size_t ii=0; ii<mass.size(); ii++)
    radius[ii] = cbl::Radius(mass[ii], m_cosmology->rho_m(0., unit1));

  function<double(double,int)> filter = [&](const double k, const int index) {
    return pow(TopHat_WF(k*radius[index]), 2)*fact;
  };

  return cosmology::Sigma::m_func_sigma(pkInterp, filter, radius, k_max, unit1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::sigma2M (const std::vector<double> mass, const std::function<double(double)> pk0Interp, const double k_max, const bool unit1) const
{
  for (auto &&mm : mass)
    if (mm<0.)
      ErrorCBL("the mass must be >0!", "sigma2M", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  if (m_cosmology->sigma8()>0.)
    // (sigma8 = sigma(8Mpc/h))
    fact = pow(m_cosmology->sigma8(), 2)/m_sigma2M({Mass(8., m_cosmology->rho_m(0., true))}, pk0Interp, k_max, 1., true)[0];

  return m_sigma2M(mass, pk0Interp, k_max, fact, unit1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::sigma2M (const std::vector<double> mass, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file, const bool unit1) const
{
  for (auto &&mm : mass)
    if (mm<0.)
      ErrorCBL("the mass must be >0!", "sigma2M", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  PkXi PX(m_cosmology);
  
  const function<double(double)> Pk = PX.Pk_interp_for_sigma(method_Pk, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  // if the power spectrum is read from file or m_cosmology->sigma8()<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2M=sigma2M_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file == par::defaultString || is_parameter_file) {
    if (m_cosmology->sigma8() > 0) {
      function<double(double)> Pk_0 = (redshift == 0.) ? Pk : PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
      // (sigma8 = sigma(8Mpc/h))
      fact = pow(m_cosmology->sigma8(), 2)/m_sigma2M({Mass(8., m_cosmology->rho_m(0., true))}, Pk_0, k_max, 1., true)[0];
    }
  }

  return m_sigma2M(mass, Pk, k_max, fact, unit1);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::dnsigma2M (const int nd, const std::vector<double> mass, const std::function<double(double)> pk0Interp, const double k_max, const bool unit1) const
{
  for (auto &&mm : mass)
    if (mm<0.)
      ErrorCBL("the mass must be >0!", "dnsigma2M", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  if (m_cosmology->sigma8()>0.)
    // (sigma8 = sigma(8Mpc/h))
    fact = pow(m_cosmology->sigma8(), 2)/m_sigma2M({Mass(8., m_cosmology->rho_m(0., true))}, pk0Interp, k_max, 1., true)[0];

  if (nd==1) {
    
    vector<double> radius(mass.size());
    vector<double> dRdM(mass.size(), pow(3./(4.*cbl::par::pi*m_cosmology->rho_m(0., unit1)), 1./3.)/3.);
    
    for (size_t ii=0; ii<mass.size(); ii++) {
      radius[ii] = cbl::Radius(mass[ii], m_cosmology->rho_m(0., unit1));
      dRdM[ii] *= pow(mass[ii], -2./3.);
    }

    function<double(double,int)> filter = [&](const double k, const int index) {
      return cbl::TopHat_WF(k*radius[index])*cbl::TopHat_WF_D1(k*radius[index])*k*dRdM[index]*2.*fact;
    };

    return cosmology::Sigma::m_func_sigma(pk0Interp, filter, radius, k_max, unit1);
  }

  else {
    ErrorCBL("", "dnsigma2M", "Sigma.cpp", glob::ExitCode::_workInProgress_);
    return vector<double>();
  }
}


// =====================================================================================


std::vector<double> cbl::cosmology::Sigma::dnsigma2M (const int nd, const std::vector<double> mass, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file, const bool unit1) const
{
  for (auto &&mm : mass)
    if (mm < 0.)
      ErrorCBL("the mass must be >0!", "dnsigma2M", "Sigma.cpp");

  // the normalisation factor
  double fact = 1.;

  PkXi PX(m_cosmology);
  
  auto Pk = PX.Pk_interp_for_sigma(method_Pk, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  // if the power spectrum is read from file or m_cosmology->sigma8()<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2M=sigma2M_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file == par::defaultString || is_parameter_file) {
    if (m_cosmology->sigma8()>0.) {
      auto Pk_0 = (redshift == 0.) ? Pk : PX.Pk_interp_for_sigma(method_Pk, 0., store_output, output_root, interpType, k_max, input_file, is_parameter_file);
      // (sigma8 = sigma(8Mpc/h))
      fact = pow(m_cosmology->sigma8(), 2)/m_sigma2M({Mass(8., m_cosmology->rho_m(0., true))}, Pk_0, k_max, 1., true)[0];
    }
  }

  if (nd==1) {
    
    const double rho = (input_file!=par::defaultString && !is_parameter_file) ? m_cosmology->rho_m(0.) : m_cosmology->rho_m(0., unit1);

    vector<double> radius(mass.size());
    vector<double> dRdM(mass.size(), pow(3./(4.*cbl::par::pi*rho), 1./3.)/3.);
      
    for (size_t ii=0; ii<mass.size(); ii++) {
      radius[ii] = cbl::Radius(mass[ii], m_cosmology->rho_m(0., unit1));
      dRdM[ii] *= pow(mass[ii], -2./3.);
    }

    function<double(double,int)> filter = [&](const double k, const int index) {
      return cbl::TopHat_WF(k*radius[index])*cbl::TopHat_WF_D1(k*radius[index])*k*dRdM[index]*2.*fact;
    };

    return cosmology::Sigma::m_func_sigma(Pk, filter, radius, k_max, unit1);
  }

  else {
    ErrorCBL("", "dnsigma2M", "Sigma.cpp", glob::ExitCode::_workInProgress_);
    return vector<double>();
  }
}


// =====================================================================================


string cbl::cosmology::Sigma::create_grid_sigmaM (const string method_SS, const double redshift, const bool store_output, const string output_root, const string interpType, const double k_max, const string input_file, const bool is_parameter_file) const
{
  string norm = (m_cosmology->sigma8()>0.) ? "_sigma8" + conv(m_cosmology->sigma8(), par::fDP3) : "_scalar_amp" + conv(m_cosmology->scalar_amp(), par::ee3);

  cbl::Path path;
  string dir_grid = path.DirCosmo() + "/Cosmology/Tables/grid_SigmaM/unit" + conv(m_cosmology->unit(), par::fINT) + "/";
  string MK = "mkdir -p " + dir_grid;
  if (system(MK.c_str())) {};

  string file_grid = dir_grid + "grid_" + method_SS + norm + "_h" + conv(m_cosmology->little_h(), par::fDP6) + "_OmB" + conv(m_cosmology->Omega_baryon(), par::fDP6) + "_OmCDM" + conv(m_cosmology->Omega_CDM(), par::fDP6) + "_OmL" + conv(m_cosmology->Omega_DE(), par::fDP6) + "_OmN" + conv(m_cosmology->Omega_neutrinos(), par::fDP6) + "_Z" + conv(redshift, par::fDP6) + "_scalar_amp" + conv(m_cosmology->scalar_amp(), par::ee3) + "_scalar_pivot" + conv(m_cosmology->scalar_pivot(), par::fDP6) + "_n" + conv(m_cosmology->n_spec(), par::fDP6) + "_w0" + conv(m_cosmology->w0(), par::fDP6) + "_wa" + conv(m_cosmology->wa(), par::fDP6) + ".dat";

  ifstream fin(file_grid.c_str());

  if (!fin) {
    coutCBL << endl << "I'm creating the grid file with sigma(M): " << file_grid.c_str() << "..." << endl;

    ofstream fout(file_grid.c_str());
    checkIO(fout, file_grid);

    vector<double> MM = logarithmic_bin_vector(1000, 1.e6, 3.e16);

    vector<double> Sigma2M = sigma2M(MM, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file, true);
    vector<double> Dln_Sigma = dnsigma2M(1, MM, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file, true);
    double Sigma;

    int dp = cout.precision();
    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout.precision(1);

    for (size_t k=0; k<MM.size(); k++) {
      coutCBL << "\r............." << double(k)/double(MM.size())*100. << "% completed \r";
      cout.flush();

      Sigma = sqrt(Sigma2M[k]);
      Dln_Sigma[k] *= MM[k]/(2.*Sigma2M[k]);
      fout << MM[k] << "   " << Sigma << "   " << Dln_Sigma[k] << endl;
    }

    cout.unsetf(ios::fixed);
    cout.unsetf(ios::showpoint);
    cout.precision(dp);
    fout.clear();
    fout.close();
    coutCBL << endl  << "I wrote the file: " << file_grid << endl;
  }

  fin.clear();
  fin.close();

  return file_grid;
}


// =====================================================================================


double cbl::cosmology::Sigma::Sn_PT (const int nn, const double RR, const string method_SS, const bool store_output, const string output_root, const string interpType, const double k_max, const string input_file, const bool is_parameter_file) const
{
  if (3>nn || nn>5) ErrorCBL("nn = " + conv(nn, par::fINT), "Sn_PT", "PkXi.cpp");

  double redshift = 0.; // (the hierarchical moments predicted by the PT do not depend on the redshift)

  double gamma1 = 1., gamma2 = -1., gamma3 = -1., d2S = -1., d3S = -1., Sn = 1.;

  double RHO = m_cosmology->rho_m(0., true);
  double MASS = Mass(RR,RHO);
  double SSS = sigma2M({MASS}, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0];

  gamma1 = RR/SSS*dnsigma2R(1, {RR}, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0];

  if (nn>3) {
    d2S = dnsigma2R(2, {RR}, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0];
    gamma2 = gamma1+pow(RR, 2)/SSS*d2S;
  }

  if (nn>4) {
    d3S = dnsigma2R(3, {RR}, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0];
    gamma3 = gamma2+pow(RR, 2)/SSS*(2.*d2S+RR*d3S);
  }

  if (nn==3) Sn = 34./7.+gamma1;
  if (nn==4) Sn = 60712./1323.+62./3.*gamma1+7./3.*pow(gamma1, 2)+2./3.*gamma2;
  if (nn==5) Sn = 200575880./305613.+1847200./3969.*gamma1+6940./63.*pow(gamma1, 2)+235./27.*pow(gamma1, 3)+1490./63.*gamma2+50./9.*gamma1*gamma2+10./27.*gamma3;

  return Sn;
}


// =====================================================================================


double cbl::cosmology::Sigma::Sigman_PT (const int nn, const double RR, const string method_SS, const bool store_output, const string output_root, const string interpType, const double k_max, const string input_file, const bool is_parameter_file) const
{
  if (3>nn || nn>5) ErrorCBL("nn = " + conv(nn, par::fINT), "Sigman_PT", "PkXi.cpp");

  double redshift = 0.; // (the hierarchical moments predicted by the PT do not depend on the redshift)

  double RHO = m_cosmology->rho_m(redshift, true);
  double MASS = Mass(RR, RHO);
  double SSS = sigma2M({MASS}, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0];

  double gamma1 = RR/SSS*dnsigma2R(1, {RR}, method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file)[0];

  double Sn = -1.;

  if (nn==3) Sn = 36./7.+3./2.*(gamma1+1.);
  if (nn==4) Sn = 2540./49.+33.*(gamma1+1.)+21./4.*pow(gamma1+1.,2);
  if (nn==5) Sn = 793.+794.*(gamma1+1.)+265.*pow(gamma1+1.,2)+29.4*pow(gamma1+1.,3);

  return Sn;
}
