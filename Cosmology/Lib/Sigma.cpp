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

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Cosmology::m_func_sigma (const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, function<double(double)> filter, const bool unit1) const 
{
  function<double(double)> func;
 
  vector<double> kk, Pk;
  
  double fact = (m_unit || unit1) ? 1. : m_hh;
 
  
  // the power spectrum is read from file
  
  if (input_file!=par::defaultString && !is_parameter_file) {
    
    string line;
    double KK, PK;
    
    ifstream fin(input_file.c_str()); checkIO(fin, input_file);

    while (getline(fin, line)) {
      if (line.find("#")==string::npos) { // skip comments
	stringstream ss(line);
	vector<double> num;
	ss >> KK >> PK;
	if (KK<kmax) {
	  kk.emplace_back(KK);
	  Pk.emplace_back(PK);
	}
      }
    }

    fin.clear(); fin.close();

    func = glob::FuncGrid(kk, Pk, interpType, cbl::BinType::_logarithmic_);
  }
    

  // alternatively, the power spectrum is computed using either the internal cosmological parameters, or using a parameter file
    
  else if (method_Pk=="EisensteinHu" && input_file==par::defaultString) {
    
    auto ff = [&] (const double kk)
    {
      EisensteinHu eh;
      eh.TFmdm_set_cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massive_neutrinos, m_Omega_DE, m_hh, redshift, m_scalar_amp, m_scalar_pivot, m_n_spec);
      
      if (eh.Pk(kk)!=eh.Pk(kk)) ErrorCBL("eh.Pk=nan!", "m_func_sigma", "Sigma.cpp");

      return eh.Pk(kk*fact)*pow(fact, -3.);
    };
    
    func = ff;
  }

  else if (method_Pk=="CAMB" || method_Pk=="MGCAMB" || method_Pk=="CLASS") {
    
    vector<double> lgkk, lgPk; 
    Table_PkCodes(method_Pk, false, lgkk, lgPk, redshift, store_output, output_root, kmax, input_file);

    for (size_t i=0; i<lgkk.size(); i++) {
      const double KK = pow(10., lgkk[i])*fact;
      if (KK<kmax) {
	kk.emplace_back(KK);
	Pk.emplace_back(pow(10., lgPk[i])*pow(fact, -3.));
      }
    }

    func = glob::FuncGrid(kk, Pk, interpType, cbl::BinType::_linear_);
  }

  else if (method_Pk=="EisensteinHu" && input_file!=par::defaultString)
    ErrorCBL("in the EisensteiHu case, no input files can be read!", "m_func_sigma", "Sigma.cpp");
  
  else
    ErrorCBL("the chosen method_Pk is not available!", "m_func_sigma", "Sigma.cpp");

  auto ff = [&] (const double kk)
    {
      return func(kk)*kk*kk*filter(kk);
    };
  
  // compute the mass variance
  //return 1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag(ff, 0., 1., 1.e-4)+wrapper::gsl::GSL_integrate_qagiu(ff, 1., 1.e-5);
 
  return 1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag(ff, 1.e-4, kmax, 1.e-3);
}



// =====================================================================================


double cbl::cosmology::Cosmology::m_sigma2R_notNormalised (const double radius, const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, const bool unit1) const 
{
  auto filter = [&] (const double k)
  {
    return pow(TopHat_WF(k*radius),2);
  };

  return cosmology::Cosmology::m_func_sigma(method_Pk, redshift, store_output, output_root, interpType, kmax, input_file, is_parameter_file, filter, unit1);

}


// =====================================================================================


double cbl::cosmology::Cosmology::sigma2R (const double radius, const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, const bool unit1) const 
{
  if (radius<0) ErrorCBL("the radius must be >0!", "sigma2R", "Sigma.cpp");
  
  // the normalisation factor
  double fact = 1.;

  // if the power spectrum is read from file or m_sigma8<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2R=sigma2R_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file==par::defaultString || is_parameter_file) {
    if (m_sigma8>0) 
      // sigma_8 = sigma(8Mpc/h)
      fact = pow(m_sigma8, 2)/m_sigma2R_notNormalised(8., method_Pk, 0., store_output, output_root, interpType, kmax, input_file, is_parameter_file, true); // normalization factor
  }

  return m_sigma2R_notNormalised(radius, method_Pk, redshift, store_output, output_root, interpType, kmax, input_file, is_parameter_file, unit1)*fact;
}


// =====================================================================================


double cbl::cosmology::Cosmology::dnsigma2R (const int nd, const double radius, const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, const bool unit1) const 
{
  if (radius<0) ErrorCBL("the radius must be >0!", "dnsigma2R", "Sigma.cpp");
  
  // the normalisation factor
  double fact = 1.;

  // if the power spectrum is read from file or m_sigma8<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2R=sigma2R_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file==par::defaultString || is_parameter_file) {
    if (m_sigma8>0) 
      // sigma_8 = sigma(8Mpc/h)
      fact = pow(m_sigma8, 2)/m_sigma2R_notNormalised(8., method_Pk, 0., store_output, output_root, interpType, kmax, input_file, is_parameter_file, true); // normalization factor    
  }

  if (nd==1) {

    auto filter = [&] (const double k)
    {
      return 2.*cbl::TopHat_WF(k*radius)*cbl::TopHat_WF_D1(k*radius)*k;
    };

    return cbl::cosmology::Cosmology::m_func_sigma(method_Pk, redshift, store_output, output_root, interpType, kmax, input_file, is_parameter_file, filter, unit1)*fact;

  }

  else
    return ErrorCBL("", "dnsigma2R", "Sigma.cpp", glob::ExitCode::_workInProgress_);
}


// =====================================================================================


double cbl::cosmology::Cosmology::m_sigma2M_notNormalised (const double mass, const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, const bool unit1) const 
{
  if (mass<0) ErrorCBL("the mass must be >0!", "m_sigma2M_notNormalised", "Sigma.cpp");
  
  const double radius = (m_RhoZero>0) ? cbl::Radius(mass, m_RhoZero) : cbl::Radius(mass, rho_m(redshift, unit1)); 
  
  auto filter = [&] (const double k)
  {
    return pow(TopHat_WF(k*radius), 2);
  };
  
  return Cosmology::m_func_sigma(method_Pk, redshift, store_output, output_root, interpType, kmax, input_file, is_parameter_file, filter, unit1);
}


// =====================================================================================


double cbl::cosmology::Cosmology::sigma2M (const double mass, const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, const bool unit1) const 
{
  if (mass<0) ErrorCBL("the mass must be >0!", "sigma2M", "Sigma.cpp");
  
  // the normalisation factor
  double fact = 1.;
  
  // if the power spectrum is read from file or m_sigma8<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2M=sigma2M_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file==par::defaultString || is_parameter_file) {
    if (m_sigma8>0)
      // (sigma8 = sigma(8Mpc/h))
      fact = pow(m_sigma8, 2)/m_sigma2M_notNormalised(Mass(8., rho_m(0., true)), method_Pk, 0., store_output, output_root, interpType, kmax, input_file, is_parameter_file, true);
  }
  
  return m_sigma2M_notNormalised(mass, method_Pk, redshift, store_output, output_root, interpType, kmax, input_file, is_parameter_file, unit1)*fact;
}


// =====================================================================================


double cbl::cosmology::Cosmology::dnsigma2M (const int nd, const double mass, const string method_Pk, const double redshift, const bool store_output, const string output_root, const string interpType, const double kmax, const string input_file, const bool is_parameter_file, const bool unit1) const 
{
  if (mass<0) ErrorCBL("the mass must be >0!", "dnsigma2M", "Sigma.cpp");
  
  // the normalisation factor
  double fact = 1.;

  // if the power spectrum is read from file or m_sigma8<0, then the
  // mass variance does not need to be normalised (hence fact=1 and
  // sigma2M=sigma2M_unnormalised); otherwise the normalisation factor
  // is computed
  if (input_file==par::defaultString || is_parameter_file) {
    if (m_sigma8>0)
      // (sigma8 = sigma(8Mpc/h))
      fact = pow(m_sigma8, 2)/m_sigma2M_notNormalised(Mass(8., rho_m(0., true)), method_Pk, 0., store_output, output_root, interpType, kmax, input_file, is_parameter_file, true);
  }

  if (nd==1) {

    //return wrapper::gsl::GSL_derivative(bind(&Cosmology::sigma2M, this, placeholders::_1, method_Pk, redshift, output_root, interpType, kmax, input_file, is_parameter_file, unit1), mass, 1.e11);
    //return wrapper::gsl::GSL_derivative(bind(&Cosmology::sigma2M, this, placeholders::_1, method_Pk, redshift, output_root, interpType, kmax, input_file, is_parameter_file), mass, mass*1.e-3);

    const double rho = (m_RhoZero>0) ?  m_RhoZero : rho_m(redshift, unit1);

    const double radius = cbl::Radius(mass, rho); 

    const double dRdM = pow(3./(4.*cbl::par::pi*rho), 1./3.)*pow(mass, -2./3.)/3.;

    auto filter = [&] (const double k)
    {
      return 2.*cbl::TopHat_WF(k*radius)*cbl::TopHat_WF_D1(k*radius)*k*dRdM;
    };

    return cosmology::Cosmology::m_func_sigma(method_Pk, redshift, store_output, output_root, interpType, kmax, input_file, is_parameter_file, filter, unit1)*fact;

  }
  else
    return ErrorCBL("", "dnsigma2M", "Sigma.cpp", glob::ExitCode::_workInProgress_);
}


// =====================================================================================


std::string cbl::cosmology::Cosmology::create_grid_sigmaM (const string method_SS, const double redshift, const bool store_output,const string output_root, const string interpType, const double k_max, const string input_file, const bool is_parameter_file) const 
{
  string norm = (m_sigma8>0) ? "_sigma8"+conv(m_sigma8, par::fDP3) : "_scalar_amp"+conv(m_scalar_amp, par::ee3);

  cbl::Path path;
  string dir_grid = path.DirCosmo()+"/Cosmology/Tables/grid_SigmaM/unit"+conv(m_unit,par::fINT)+"/";
  string MK = "mkdir -p "+dir_grid; if (system (MK.c_str())) {};

  string file_grid = dir_grid+"grid_"+method_SS+norm+"_h"+conv(m_hh, par::fDP6)+"_OmB"+conv(m_Omega_baryon, par::fDP6)+"_OmCDM"+conv(m_Omega_CDM, par::fDP6)+"_OmL"+conv(m_Omega_DE, par::fDP6)+"_OmN"+conv(m_Omega_neutrinos, par::fDP6)+"_Z"+conv(redshift, par::fDP6)+"_scalar_amp"+conv(m_scalar_amp, par::ee3)+"_scalar_pivot"+conv(m_scalar_pivot, par::fDP6)+"_n"+conv(m_n_spec, par::fDP6)+"_w0"+conv(m_w0, par::fDP6)+"_wa"+conv(m_wa, par::fDP6)+".dat";

  ifstream fin(file_grid.c_str());

  if (!fin) {

    coutCBL << endl << "I'm creating the grid file with sigma(M): " << file_grid.c_str() << "..." << endl;
    
    ofstream fout(file_grid.c_str()); checkIO(fout, file_grid); 

    vector<double> MM = logarithmic_bin_vector(1000, 1.e6, 3.e16);
    
    double SSS, Sigma, Dln_Sigma;
    
    int dp = cout.precision();
    cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(1);
      
    for (size_t k=0; k<MM.size(); k++) {
      coutCBL << "\r............." << double(k)/double(MM.size())*100. << "% completed \r"; cout.flush(); 

      SSS = sigma2M(MM[k], method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file, true);
      Sigma = sqrt(SSS);
      Dln_Sigma = dnsigma2M(1, MM[k], method_SS, redshift, store_output, output_root, interpType, k_max, input_file, is_parameter_file, true)*(MM[k]/(2.*SSS));
      fout << MM[k] << "   " << Sigma << "   " << Dln_Sigma << endl;
    }
    
    cout.unsetf(ios::fixed); cout.unsetf(ios::showpoint); cout.precision(dp);
    fout.clear(); fout.close(); coutCBL << endl << "I wrote the file: " << file_grid << endl;
  }
  
  fin.clear(); fin.close();
  
  return file_grid;
}
