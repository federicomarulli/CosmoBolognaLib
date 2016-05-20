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
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::Cosmology::SSR (const double RR, const string method_Pk, const double redshift, const string output_root, const double kmax, const string file_par) const 
{
  double SS = -1;

  if (method_Pk=="EisensteinHu") {
    cosmobl::classfunc::func_SSR func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, RR, redshift); 
    Midpnt<cosmobl::classfunc::func_SSR> q1(func,0.,1.); 
    Midinf<cosmobl::classfunc::func_SSR> q2(func,1.,1.e99);
    double Int = qromo(q1)+qromo(q2);

    SS = 1./(2.*pow(par::pi,2))*Int;
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0;
    double RHO = Rho(m_Omega_matter,m_Omega_neutrinos); 
   
    Table_PkCodes (method_Pk, do_nonlinear, lgkk, lgPk, redshift, output_root, kmax, file_par);
        
    cosmobl::classfunc::func_SSR_Table func (m_hh, m_n_spec, RHO, m_unit, lgkk, lgPk, RR); 

    Midpnt<cosmobl::classfunc::func_SSR_Table> q1(func,0.,1.); 
    Midinf<cosmobl::classfunc::func_SSR_Table> q2(func,1.,1.e99);
    double Int = qromo(q1)+qromo(q2);

    SS = 1./(2.*pow(par::pi,2))*Int;
  }

  return SS;
}


// =====================================================================================


double cosmobl::Cosmology::SSR_norm (const double RR, const string method_Pk, const double redshift, const string output_root, const double kmax, const string file_par) const 
{
  double fact = 1.;

  if (m_sigma8>0) {
    double RRR = 8.; // sigma_8 = sigma(8Mpc/h)
    fact = (m_sigma8*m_sigma8)/SSR(RRR, method_Pk, redshift, output_root, kmax, file_par); // normalization factor
  }
  else if (method_Pk=="EisensteinHu") ErrorMsg("Error in cosmobl::Cosmology::SSR_norm of Sigma.cpp!");
  
  return SSR(RR, method_Pk, redshift, output_root, kmax, file_par)*fact;
}


// =====================================================================================


double cosmobl::Cosmology::dnSR (const int nd, const double RR, const string method_Pk, const double redshift, const string output_root, const string interpType, const int Num, const double stepsize, const double kmax, const string file_par) const 
{
  double dR = RR*1.e-7;
  double RRR = RR+dR;
  return (SSR_norm(RRR, method_Pk, redshift, output_root, kmax, file_par)-SSR_norm(RR, method_Pk, redshift, output_root, kmax, file_par))/dR;
  /*
    cosmobl::classfunc::func_SSRd SSRd(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, redshift, output_root, kmax, file_par);

    return Deriv(nd, RR, SSRd, interpType, Num, stepsize);
  */
}


// =====================================================================================


double cosmobl::Cosmology::SSM (const double MM, const string method_Pk, const double redshift, const string output_root, const double kmax, const string file_par) const 
{
  double SS = -1;

  if (method_Pk=="EisensteinHu") {
    cosmobl::classfunc::func_SSM func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, MM, redshift); 
    Midpnt<cosmobl::classfunc::func_SSM> q1(func,0.,1.); 
    Midinf<cosmobl::classfunc::func_SSM> q2(func,1.,1.e99);
    double Int = qromo(q1)+qromo(q2);

    SS = 1./(2.*pow(par::pi,2))*Int;
  }

  else if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0;
    double RHO = Rho(m_Omega_matter,m_Omega_neutrinos); 
    
    if (file_par!=par::defaultString) {
      string Warn = "Check that in the parameter file Omega_matter = " + conv(m_Omega_matter,par::fDP3) + " and Omega_neutrinos = " + conv(m_Omega_neutrinos,par::fDP3) + "!";
      WarningMsg(Warn);
    }

    Table_PkCodes (method_Pk, do_nonlinear, lgkk, lgPk, redshift, output_root, kmax, file_par);
    
    cosmobl::classfunc::func_SSM_Table func (m_hh, m_n_spec, RHO, m_unit, lgkk, lgPk, MM); 
    Midpnt<cosmobl::classfunc::func_SSM_Table> q1(func, 0., 1.); 
    Midinf<cosmobl::classfunc::func_SSM_Table> q2(func, 1., 1.e99);
    
    double Int = qromo(q1)+qromo(q2);

    SS = 1./(2.*pow(par::pi,2))*Int;
  }

  else ErrorMsg("Error in cosmobl::Cosmology::SSM of Sigma.cpp: method_Pk is wrong!");

  return SS;
}


// =====================================================================================


double cosmobl::Cosmology::SSM_norm (const double MM, const string method_Pk, const double redshift, const string output_root, const double kmax, const string file_par) const 
{
  double fact = 1.;

  if (m_sigma8>0) {
    double RR = 8.; // sigma_8 = sigma(8Mpc/h)
    double RHO = Rho(m_Omega_matter,m_Omega_neutrinos); 
    double Mss = Mass(RR, RHO);
    fact = (m_sigma8*m_sigma8)/SSM(Mss, method_Pk, redshift, output_root, kmax, file_par); // normalization factor
  }
  else if (method_Pk=="EisensteinHu") ErrorMsg("Error in cosmobl::Cosmology::SSM_norm of Sigma.cpp!");

  return SSM(MM, method_Pk, redshift, output_root, kmax, file_par)*fact;
}


// =====================================================================================


double cosmobl::Cosmology::dnSM (const int nd, const double MM, const string method_Pk, const double redshift, const string output_root, const string interpType, const int Num, const double stepsize, const double kmax, const string file_par) const 
{
 
  double dM = MM*1.e-7;
  double MMM = MM+dM;
  return (SSM_norm(MMM, method_Pk, redshift, output_root, kmax, file_par)-SSM_norm(MM, method_Pk, redshift, output_root, kmax, file_par))/dM;
  
  /*
  cosmobl::classfunc::func_SSMd SSMd(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, redshift, output_root, kmax, file_par);
  
  return Deriv(nd, MM, SSMd, interpType, Num, stepsize);
  */

  /*
  double RHO = Rho(m_Omega_matter, m_Omega_neutrinos); 
  double RR = Radius(MM,RHO);

  cosmobl::classfunc::func_SSRd SSRd(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, redshift, output_root, kmax, file_par);

  return Deriv(nd, RR, SSRd, interpType, Num, stepsize)/(4.*par::pi*RHO*RR*RR);
  */
}


// =====================================================================================


string cosmobl::Cosmology::create_grid_sigmaM (const string method_SS, const double redshift, const string output_root, const string interpType, const int Num, const double stepsize, const double kmax, const string file_par) const 
{ 
  string norm = (m_sigma8>0) ? "_sigma8"+conv(m_sigma8,par::fDP3) : "_scalar_amp"+conv(m_scalar_amp,par::ee3);
  string dir_cosmo=fullpath(par::DirCosmo);
  string dir_grid = dir_cosmo+"Cosmology/Tables/grid_SigmaM/unit"+conv(m_unit,par::fINT)+"/";
  string MK = "mkdir -p "+dir_grid; if (system (MK.c_str())) {};

  string file_grid = dir_grid+"grid_"+method_SS+norm+"_h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+".dat";

  ifstream fin(file_grid.c_str());

  if (!fin) {

    cout << endl << "I'm creating the grid file with sigma(M): " << file_grid.c_str() << "..." << endl;
    
    ofstream fout (file_grid.c_str()); checkIO (file_grid, 0); 

    vector<double> MM = logarithmic_bin_vector(1000, 1.e6, 3.e16);
    
    double SSS, Sigma, Dln_Sigma;
    
    for (size_t k=0; k<MM.size(); k++) {
      cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(1);
      cout << "\r..." << double(k)/double(MM.size())*100. << "% completed \r"; cout.flush(); 

      SSS = SSM_norm(MM[k], method_SS, redshift, output_root, kmax, file_par);
      Sigma = sqrt(SSS);
      Dln_Sigma = dnSM(1, MM[k], method_SS, redshift, output_root, interpType, Num, stepsize, kmax, file_par)*(MM[k]/(2.*SSS));
      fout << MM[k] << "   " << Sigma << "   " << Dln_Sigma << endl;
    }

    fout.clear(); fout.close(); cout << endl << "I wrote the file: " << file_grid << endl;
  }
  
  fin.clear(); fin.close();

  return file_grid;
}
