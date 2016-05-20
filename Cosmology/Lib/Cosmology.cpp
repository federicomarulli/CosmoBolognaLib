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
 *  @file Cosmology/Lib/Cosmology.cpp
 *
 *  @brief Generic methods of the class Cosmology  
 *
 *  This file contains the implementation of the \e generic methods of
 *  the class Cosmology
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


cosmobl::Cosmology::Cosmology (const double Omega_matter, const double Omega_baryon, const double Omega_neutrinos, const double massless_neutrinos, const int massive_neutrinos, const double Omega_DE, const double Omega_radiation, const double hh, const double scalar_amp, const double n_spec, const double w0, const double wa, const double fNL, const int type_NG, const string model, const bool unit)
  : m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_sigma8(-1.), m_scalar_amp(scalar_amp), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_model(model), m_unit(unit)              
{
  if (m_Omega_matter==0) ErrorMsg("Error in cosmobl::Cosmology::Cosmology of Cosmology.cpp: Omega_matter=0!");

  m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
  m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
  m_H0 = (m_unit) ? 100. : 100.*m_hh;
  m_t_H = 1./m_H0;
  m_D_H = par::cc*m_t_H;
  m_RhoZero = Rho(m_Omega_matter, m_Omega_neutrinos); 
  m_Pk0_EH = 1., m_Pk0_CAMB = 1., m_Pk0_MPTbreeze = 1., m_Pk0_CLASS = 1.;
} 


// =====================================================================================


double cosmobl::Cosmology::value (const CosmoPar parameter) const
{
  double param_value;
  
  switch (parameter) {

  case (_Omega_matter_LCDM_):
    param_value = m_Omega_matter;
    break;

  case (_Omega_matter_):
    return m_Omega_matter;
    break;

  case (_Omega_baryon_):        
    param_value = m_Omega_baryon;
    break;

  case (_Omega_neutrinos_):      
    param_value = m_Omega_neutrinos;
    break;

  case (_massless_neutrinos_):   
    param_value = m_massless_neutrinos;
    break;

  case (_massive_neutrinos_):    
    param_value = m_massive_neutrinos;
    break;

  case (_Omega_DE_):            
    param_value = m_Omega_DE;
    break;

  case (_Omega_radiation_):     
    param_value = m_Omega_radiation;
    break;

  case (_H0_):
    param_value = m_H0;
    break;

  case (_scalar_amp_):           
    param_value = m_scalar_amp;
    break;

  case (_n_spec_):               
    param_value = m_n_spec;
    break;

  case (_w0_):
    param_value = m_w0;
    break;

  case (_wa_):             
    param_value = m_wa;
    break;

  case (_fNL_):                  
    param_value = m_fNL;
    break;

  case (_sigma8_):
    param_value = m_sigma8;
    break;

  default:
    ErrorMsg("Error in cosmobl::Cosmology::value_CosmoPar of Cosmology.cpp: no such a variable in the list!");
  }
  
  return param_value;
}


// =====================================================================================


void cosmobl::Cosmology::set_parameter (const CosmoPar parameter, const double value)
{
  switch (parameter) {
    
  case (_Omega_matter_LCDM_):
    set_Omega(value);
    break;

  case (_Omega_matter_):
    set_OmegaM(value);
    break;
  
  case (_Omega_baryon_):        
    set_OmegaB(value);
    break;

  case (_Omega_neutrinos_):      
    set_OmegaNu (value, m_massless_neutrinos, m_massive_neutrinos);
    break;

  case (_massless_neutrinos_):   
    set_OmegaNu (m_Omega_neutrinos, value, m_massive_neutrinos);
    break;

  case (_massive_neutrinos_):    
    set_OmegaNu (m_Omega_neutrinos, m_massless_neutrinos, int(value));
    break;

  case (_Omega_DE_):            
    set_OmegaDE(value); 
    break;

  case (_Omega_radiation_):     
    set_Omega_radiation(value);
    break;

  case (_H0_):
    set_H0(value); 
    break;

  case (_scalar_amp_):           
    set_scalar_amp(value);
    break;

  case (_n_spec_):               
    set_n_spec(value);
    break;

  case (_w0_):
    set_w0(value);
    break;

  case (_wa_):             
    set_wa(value);
    break;

  case (_fNL_):                  
    set_fNL(value);
    break;

  case (_sigma8_):
    set_sigma8(value); 
    break;

  default:
    ErrorMsg("Error in cosmobl::Cosmology::set_CosmoPar of Cosmology.cpp: no such a variable in the list!");
  }
}


// =====================================================================================


void cosmobl::Cosmology::set_parameter (const vector<CosmoPar> parameter, const vector<double> value)
{
  for (size_t i=0; i<parameter.size(); i++)
    set_parameter(parameter[i], value[i]);
}


// =====================================================================================


double cosmobl::Cosmology::w_CPL (const double redshift) const 
{
  return m_w0+m_wa*redshift/(1.+redshift);
}


// =====================================================================================


double cosmobl::Cosmology::f_DE (const double redshift) const
{
  // CPL parameterisation (see e.g. Bassett & Hlozek 2010)
  
  return pow((1.+redshift),3.*(1.+m_w0+m_wa))*exp(-3.*m_wa*redshift/(1.+redshift)); 
  
  /*
  // direct calculation, useful for future generalizations:
  func_fDE Func(m_w0,m_wa);
  return exp(3.*qromb(Func,0.,redshift));
  */
}


// =====================================================================================


double cosmobl::Cosmology::EE (const double redshift) const
{
  return sqrt(m_Omega_matter*pow((1.+redshift),3)+m_Omega_DE*f_DE(redshift)+m_Omega_k*pow((1.+redshift),2)+m_Omega_radiation*pow(1.+redshift,4));
}


// =====================================================================================


double cosmobl::Cosmology::HH (const double redshift) const 
{
  return m_H0*EE(redshift);
}


// =====================================================================================


double cosmobl::Cosmology::OmegaM (const double redshift) const
{
  return m_Omega_matter/EE2(redshift)*1./(1.+redshift);
  //return pow(m_H0/HH(redshift),2)*m_Omega_matter*pow(1.+redshift,3);
}

double cosmobl::Cosmology::OmegaDE (const double redshift) const 
{
  if (m_wa!=0) ErrorMsg("Error in cosmobl::Cosmology::OmegaDE of Cosmology.cpp: w_a!=0 -> work in progress...");

  return m_Omega_DE/EE2(redshift)*pow(1./(1.+redshift),1.-3.*m_w0);
} 

double cosmobl::Cosmology::OmegaR (const double redshift) const 
{
  return m_Omega_radiation/EE2(redshift);
} 

double cosmobl::Cosmology::OmegaK (const double redshift) const 
{
  return m_Omega_k/EE2(redshift)*pow(1./(1.+redshift),2);
} 

double cosmobl::Cosmology::Omega (const double redshift) const 
{
  if (m_wa!=0) ErrorMsg("Error in cosmobl::Cosmology::Omega of Cosmology.cpp: w_a!=0 -> work in progress...");

  double aa = 1./(1.+redshift);
  return (m_Omega_radiation+m_Omega_matter*aa+m_Omega_DE*pow(aa,1.-3.*m_w0))/EE2(redshift);
} 


// =====================================================================================


double cosmobl::Cosmology::gg (const double redshift) const 
{
  return 2.5*OmegaM(redshift)*pow(pow(OmegaM(redshift),4./7.)-(1.-OmegaM(redshift))+(1.+OmegaM(redshift)*0.5)*(1.+((1.-OmegaM(redshift))/70.)),-1.);
}


// =====================================================================================


double cosmobl::Cosmology::DD (const double redshift) const 
{
  // by Carroll, Press, & Turner 1992
  //return 5.*OmegaM(redshift)/(2*(1+redshift))/(1./70.+209./140.*OmegaM(redshift)-pow(OmegaM(redshift),2)/140.+pow(OmegaM(redshift),4./7.)); 

  double zero = 0.;
  return 1./(1.+redshift)*gg(redshift)/gg(zero);
}


// =====================================================================================


double cosmobl::Cosmology::D_C (const double redshift) const 
{
  if (redshift<0) ErrorMsg ("Error in cosmobl::Cosmology::D_C of Cosmology.cpp: redshift have to be >=0!");
  
  double Dc;

  if (m_model=="LCDM") {
    cosmobl::classfunc::E_inv EE (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit);
    Dc = qromb(EE, 0., redshift);
  }
  
  else {
    string dir = par::DirCosmo+"Cosmology/Tables/dc_cDE/";
    string file_in;
    if (m_model=="LCDM_Baldi_wmap7") file_in = dir+"LCDM-wmap7-comovingdist.dat"; 
    else if (m_model=="EXP005_Baldi_wmap7") file_in = dir+"EXP005-wmap7-comovingdist.dat";
    else if (m_model=="EXP010e2_Baldi_wmap7") file_in = dir+"EXP010e2-wmap7-comovingdist.dat";
    else if (m_model=="LCDM_Baldi_CoDECS") file_in = dir+"LCDM_CoDECS-comovingdist.dat"; 
    else if (m_model=="EXP001_Baldi_CoDECS") file_in = dir+"EXP001_CoDECS-comovingdist.dat";
    else if (m_model=="EXP002_Baldi_CoDECS") file_in = dir+"EXP002_CoDECS-comovingdist.dat";
    else if (m_model=="EXP003_Baldi_CoDECS") file_in = dir+"EXP003_CoDECS-comovingdist.dat";
    else if (m_model=="EXP008e3_Baldi_CoDECS") file_in = dir+"EXP008e3_CoDECS-comovingdist.dat";
    else if (m_model=="EXP010e2_Baldi_CoDECS") file_in = dir+"EXP010e2_CoDECS-comovingdist.dat";
    else if (m_model=="SUGRA003_Baldi_CoDECS") file_in = dir+"SUGRA003_CoDECS-comovingdist.dat";
    else { string Err = "Error in cosmobl::Cosmology::D_C of Cosmology.cpp: model = " + m_model + "!"; ErrorMsg(Err); }
                     
    ifstream fin (file_in.c_str()); checkIO (file_in,1); 
    
    double Red, DC;
    vector<double> Redshift, dc;
    while (fin >>Red>>DC) {
      Redshift.push_back(Red);
      dc.push_back(DC);
    }
    fin.clear(); fin.close();

    double err = -1.;
    Dc = interpolated(redshift, Redshift, dc, "Rat");
    
    if (err/Dc>0.1) {
      string Err = "Error in cosmobl::Cosmology::D_C of Cosmology.cpp: " + conv(redshift,par::fDP3) + "   " + conv(Redshift.size(),par::fINT) + "   " + conv(dc.size(),par::fINT);
      ErrorMsg(Err);
    }
  }
  
  return m_D_H*Dc;
}


// =====================================================================================


void cosmobl::Cosmology::D_C_table (const string file_table, const double z_min, const double z_max, const int step, vector<double> &Redshift, vector<double> &dc) const
{
  string File_table = par::DirCosmo+"Cosmology/Tables/dc/"+file_table;
 
  ifstream fin;
  fin.open (File_table.c_str());
  if (!fin) {

    ofstream fout (File_table.c_str()); checkIO (File_table,0); 
    
    double delta_z = (z_max-z_min)/step;
    double z1 = z_min;
    double z2 = z_min+delta_z;

    for (int i=0; i<step; i++) {
      double zmean = (z1+z2)*0.5;
      fout <<zmean<<"   "<<D_C(zmean)<<endl;
      z1 = z2; z2 += delta_z;
    }
    
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<File_table<<endl;
  }
  fin.clear(); fin.close();
  
  fin.open(File_table.c_str());
  double Red, DC;
  while (fin >>Red>>DC) {
    Redshift.push_back(Red);
    dc.push_back(DC);
  }
  fin.clear(); fin.close();
}


// =====================================================================================


double cosmobl::Cosmology::D_M (const double redshift) const 
{
  if (m_Omega_k>1.e-10) 
    return m_D_H/sqrt(m_Omega_k)*sinh(sqrt(m_Omega_k)*D_C(redshift)/m_D_H);
					
  else if (fabs(m_Omega_k)<1.e-10) 
    return D_C(redshift);
  
  else 
    return m_D_H/sqrt(-m_Omega_k)*sin(sqrt(-m_Omega_k)*D_C(redshift)/m_D_H);
}


// =====================================================================================


double cosmobl::Cosmology::D_A (const double redshift) const 
{
  return D_M(redshift)/(1.+redshift);
}


// =====================================================================================


double cosmobl::Cosmology::D_L (const double redshift) const 
{
  return (1.+redshift)*D_M(redshift);
}


// =====================================================================================


double cosmobl::Cosmology::D_V (const double redshift) const 
{
  return pow(pow(D_M(redshift),2)*par::cc*redshift/HH(redshift),1./3.);
}


// =====================================================================================


double cosmobl::Cosmology::Distance (const double redshift, const string distance_type) const 
{
  if (distance_type=="DC")
    return D_C(redshift);

  else if (distance_type=="DL")
    return D_L(redshift);

  else if (distance_type=="DA")
    return D_A(redshift);

  else if (distance_type=="Dv")
    return D_V(redshift);

  else if (distance_type=="Dvrs")
    return D_V(redshift)/rs_CAMB();

  else if (distance_type=="rsDv")
    return rs_CAMB()/D_V(redshift);

  else {
    ErrorMsg("Error in Distance of Cosmology, Cosmology/Lib/Cosmology.cpp. No such a distance type");
    return -1;
  }
}


// =====================================================================================


double cosmobl::Cosmology::lookback_time (const double redshift) const 
{
  cosmobl::classfunc::E_inv2 func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit);

  double Mpc = par::mega*par::pc*1.e-3; // in Km;
  double Gyr = par::giga*par::yr; // in sec

  return 1./(m_hh*100.)*qromb(func,0.,redshift)*Mpc/Gyr;
}


// =====================================================================================


double cosmobl::Cosmology::cosmic_time (const double redshift) const
{
  cosmobl::classfunc::E_inv3 func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_model);
  
  double aa = 1./(1.+redshift);
  Midpnt<cosmobl::classfunc::E_inv3> q1(func,0.,aa); 

  double Mpc = par::mega*par::pc*1.e-3; // in Km;
  double Gyr = par::giga*par::yr; // in sec

  return 1./(m_hh*100.)*qromo(q1)*Mpc/Gyr;
}


// =====================================================================================


double cosmobl::Cosmology::EE2 (const double redshift) const // see e.g. de Araujo 2005
{
  double aa = 1./(1.+redshift);
  return m_Omega_radiation+m_Omega_matter*aa+m_Omega_k*aa*aa+m_Omega_DE*pow(aa,1.-3.*m_w0);

}


// =====================================================================================


double cosmobl::Cosmology::qq (const double redshift) const
{
  if (m_wa!=0) ErrorMsg("Error in cosmobl::Cosmology::qq of Cosmology.cpp: w_a!=0 -> work in progress...");

  double aa = 1./(1.+redshift);
  return (m_Omega_matter*aa+2.*m_Omega_radiation+(1.+3.*m_w0)*m_Omega_DE*pow(aa,1.-3.*m_w0))/(2.*EE2(redshift));
}


// =====================================================================================


double cosmobl::Cosmology::Hdot (const double redshift) const
{
  return -pow(HH(redshift),2)*(1.+qq(redshift));
}


// =====================================================================================


double cosmobl::Cosmology::z_acc () const
{
  if (m_wa!=0) ErrorMsg("Error in cosmobl::Cosmology::z_acc of Cosmology.cpp: w_a!=0 -> work in progress...");
  
  double z_acc = pow(-(1.+3.*m_w0)*m_Omega_DE/m_Omega_matter,-1./(3.*m_w0))-1.;
  if (std::isnan(z_acc)) ErrorMsg("Error in cosmobl::Cosmology::z_acc of Cosmology.cpp!");
  
  return z_acc;
}


// =====================================================================================


double cosmobl::Cosmology::z_eq () const
{
  if (m_wa!=0) ErrorMsg("Error in cosmobl::Cosmology::z_eq of Cosmology.cpp: w_a!=0 -> work in progress...");
  return pow(m_Omega_DE/m_Omega_matter,-1./(3.*m_w0))-1.;
}


// =====================================================================================


double cosmobl::Cosmology::Mag_Volume_limited (const double z_max, const double mag_lim) const
{
  return (mag_lim-5.*log10(D_L(z_max)))-25.;
}


// =====================================================================================


double cosmobl::Cosmology::Lum_bol (const double redshift, const double flux) const 
{
  return 4.*par::pi*flux*pow(D_L(redshift),2);
}


// =====================================================================================


double cosmobl::Cosmology::Redshift (const double d_c, const double z1_guess, const double z2_guess, const double prec) const
{
  double Prec = prec;
  if (Prec<1.e-5) {
    WarningMsg("Attention in cosmobl::Cosmology::Redshift of Cosmology.cpp: prec has been set to 1.e-5");
    Prec = 1.e-5;
  }

  double redshift = -1.;
  
  if (m_model=="LCDM") {
    cosmobl::classfunc::func_z func(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, d_c); 
    redshift = zbrent(func, z1_guess, z2_guess, prec);
  }
  
  else {
    WarningMsg("Attention in cosmobl::Cosmology::Redshift of Cosmology.cpp: the quantity prec is not used");

    string dir = par::DirCosmo+"Cosmology/Tables/dc_cDE/";
    string file_in;

    if (m_model=="LCDM_Baldi_wmap7") file_in = dir+"LCDM-wmap7-comovingdist.dat"; 
    else if (m_model=="EXP005_Baldi_wmap7") file_in = dir+"EXP005-wmap7-comovingdist.dat";
    else if (m_model=="EXP010e2_Baldi_wmap7") file_in = dir+"EXP010e2-wmap7-comovingdist.dat";
    else if (m_model=="LCDM_Baldi_CoDECS") file_in = dir+"LCDM_CoDECS-comovingdist.dat"; 
    else if (m_model=="EXP001_Baldi_CoDECS") file_in = dir+"EXP001_CoDECS-comovingdist.dat";
    else if (m_model=="EXP002_Baldi_CoDECS") file_in = dir+"EXP002_CoDECS-comovingdist.dat";
    else if (m_model=="EXP003_Baldi_CoDECS") file_in = dir+"EXP003_CoDECS-comovingdist.dat";
    else if (m_model=="EXP008e3_Baldi_CoDECS") file_in = dir+"EXP003_CoDECS-comovingdist.dat";
    else if (m_model=="EXP010e2_Baldi_CoDECS") file_in = dir+"EXP010e2_CoDECS-comovingdist.dat";
    else if (m_model=="SUGRA003_Baldi_CoDECS") file_in = dir+"SUGRA003_CoDECS-comovingdist.dat";
    else { string Err = "Error in cosmobl::Cosmology::Redshift of Cosmology.cpp: model = " + m_model + "!"; ErrorMsg(Err); }
    
    ifstream fin (file_in.c_str()); checkIO (file_in,1); 
    
    double Red, DC;
    vector<double> Redshift, dc;
    while (fin >>Red>>DC) {
      Redshift.push_back(Red);
      dc.push_back(m_D_H*DC);
    }
    fin.clear(); fin.close();
    
    double err = -1;
    redshift = interpolated(d_c, dc, Redshift, "Rat");
    if (err/redshift>0.1) ErrorMsg("Error in cosmobl::Cosmology::Redshift of Cosmology.cpp!");
  }
  
  return redshift;
}


// =====================================================================================


double cosmobl::Cosmology::Redshift_LCDM (const double d_c, const double z1_guess, const double z2_guess, const bool go_fast, const double prec) const 
{
  if (m_model!="LCDM") 
    ErrorMsg("Error in cosmobl::Cosmology::Redshift_LDCM of Cosmology.cpp: this"
	     " method works only for a LambdaCDM universe");

  double Prec = prec;
  if (Prec<1.e-5) {
    WarningMsg("Attention in cosmobl::Cosmology::Redshift_LDCM of Cosmology.cpp: prec"
	       " has been set to 1.e-5");
    Prec = 1.e-5;
  }

  double fact = 1./(par::cc/100.); // factor to convert [Mpc/h] into [c/H0]
  
  double D_c = d_c * fact; // [Mpc/h] -> [c/H0]


  // first interval definition: based on user definition and cosmology

  double z0 = max(max(z1_guess, D_c), 4./pow(sqrt(m_Omega_matter)*D_c-2.,2)-1.); // z is always > than these two numbers

  double z1 = (D_c<2.) ? min(z2_guess, 4./pow(D_c-2.,2)-1.) : z2_guess; // z is always < than this number, but valid only when d_c<2
  
  z1 = max(z1, z0);


  // choice of the method

  function<double(double)> DC = (go_fast) ? bind(&Cosmology::D_C_LCDM, this, std::placeholders::_1) : bind(&Cosmology::D_C, this, std::placeholders::_1);

  
  // interval check: needs z0<=z and z1>=z

  double d0 = DC(z0)*fact;
  
  if (d0>D_c) {
    z0 = D_c/d0*z0/EE(z0);
    d0 = DC(z0)*fact;
  }

  double d1 = DC(z1)*fact;
  if (d1<D_c) {
    z1 = EE(z1)*D_c/d1*z1;
    d1 = DC(z1)*fact;
  }
 

  // calculation

  double diffz = z1-z0;
  double diffd = d1-d0;
  double prec2 = 2.*Prec; // prec is the tolerance on the result, prec2 on the interval length
 
  while (diffz>prec2 && diffd>0.) {
    z1 = z0+(D_c-d0)*(diffz/diffd); // given the shape this estimate is always > z
    d1 = DC(z1)*fact;
    z0 = z0+EE(z0)*(D_c-d0); // given the shape this estimate is always < z
    d0 = DC(z0)*fact;
    diffz = z1-z0;
    diffd = d1-d0;
  }
  
  return 0.5*(z0+z1);
}


// =====================================================================================


double cosmobl::Cosmology::Redshift_time (const double time, const double z1_guess, const double z2_guess) const 
{
  if (m_model!="LCDM") ErrorMsg("Error in cosmobl::Cosmology::Redshift_time of Cosmology.cpp: model!=LCDM -> Work in progress...");

  double prec = 0.0001;

  cosmobl::classfunc::func_zt func(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, time); 

  return zbrent (func, z1_guess, z2_guess, prec);
}


// =====================================================================================
  

double cosmobl::Cosmology::Volume (const double z1, const double z2, const double Area) const 
{
  double Area_steradians = Area/pow(180./par::pi,2);

  return 4./3.*par::pi*fabs(pow(D_C(z1),3)-pow(D_C(z2),3))*Area_steradians/(4.*par::pi);
}


// =====================================================================================

/* Alfonso Veropalumbo */

// Total Comoving volume from z=0 to z, all sky (Hogg 2000, Eq 29)
double cosmobl::Cosmology::Volume (const double zz) const 
{
  double DDMM = D_M(zz);

  if (m_Omega_k>1.e-10) 
    return (4*par::pi*pow(m_D_H,3)/(2*m_Omega_k))*(DDMM/m_D_H*sqrt(1+m_Omega_k*pow(DDMM/m_D_H,2))-pow(fabs(m_Omega_k),-0.5)*asinh(sqrt(fabs(m_Omega_k))*DDMM/m_D_H));

  else if (fabs(m_Omega_k)<1.e-10) 
    return 4*par::pi*pow(DDMM,3)/3;

  else 
    return (4*par::pi*pow(m_D_H,3)/(2*m_Omega_k))*(DDMM/m_D_H*sqrt(1+m_Omega_k*pow(DDMM/m_D_H,2))-pow(fabs(m_Omega_k),-0.5)*asin(sqrt(fabs(m_Omega_k))*DDMM/m_D_H));
}


// =====================================================================================


double cosmobl::Cosmology::max_redshift (const double Volume, const double Area, const double z_min) const
{
  cosmobl::classfunc::func_V func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, z_min, Area, Volume);
  return zbrent (func, z_min, 10., 1.e-9);
}


// =====================================================================================


double cosmobl::Cosmology::dV_dZdOmega (const double redshift, const bool angle_rad) const 
{
  // angle_rad: 1 -> Omega in steradians; 0 -> Omega in square degrees
  double conv = (angle_rad) ? 1. : 3282.80635; 

  return m_D_H*pow((1.+redshift)*D_A(redshift),2)/EE(redshift)/conv;
}


// =====================================================================================
  

double cosmobl::Cosmology::deltac (const double redshift) const 
{ 
  return 1.686*(1.+0.012299*log10(OmegaM(redshift))); 
}


// =====================================================================================


double cosmobl::Cosmology::Rho (const double Omega_matter, const double Omega_neutrinos, const bool unit1) const  
{
  double fact = (m_unit && unit1) ? 1. : m_hh*m_hh;

  //return 2.778e11*Omega_matter*fact; 
  return 2.778e11*(Omega_matter-Omega_neutrinos)*fact; // check!!!!
}


// =====================================================================================


double cosmobl::Cosmology::DeltaR (const double Delta_crit, const double redshift) const 
{
  return Delta_crit/OmegaM(redshift);
}


// =====================================================================================


double cosmobl::Cosmology::D_C_LCDM (const double redshift) const 
{  
  if (m_model!="LCDM" || (1.-m_Omega_matter-m_Omega_DE)>1.e-30 || fabs(m_Omega_k)>1.e-30) 
    ErrorMsg("Error in cosmobl::Cosmology::D_C_LCDM of Cosmology.cpp: this method"
    	     " works only for a flat LambdaCDM universe; it does not work for"
	     " non-standard dark energy or non-flat models");

  if (redshift>10)
    WarningMsg("Attention: the method cosmobl::Cosmology::D_C_LCDM of Cosmology.cpp"
	       " works only at low enough redshifts where &Omega<SUB>r</SUB> is"
	       " negligible");

  double CC = 1./pow(3.,0.25)/(pow(m_Omega_matter,1./3.)*pow(1.-m_Omega_matter,1./6.));
  double f_m = (1.-sqrt(3.))*pow((1.-m_Omega_matter)/m_Omega_matter,1./3.);
  double f_p = (1.+sqrt(3.))*pow((1.-m_Omega_matter)/m_Omega_matter,1./3.);
  double phi0 = acos((1.+f_m)/(1.+f_p));
  double F_phi0 = elf_dz(phi0);

  double aa = 1./(1.+redshift);
  double phi1 = acos((1.+f_m*aa)/(1.+f_p*aa));
  
  return CC*(F_phi0-elf_dz(phi1))*2997.9199;
}  


// =====================================================================================


double cosmobl::Cosmology::elf_dz (const double phi) const 
{
  int jj = round(phi/M_PI);
  double phi0 = phi-jj*M_PI;

  // then, it has to be reduced again to [0,pi/2], taking the sign into account
  double ss = phi0/abs(phi0);
  phi0 = ss*phi0;

  // optimisation parameters
  double phiS = 1.249;
  double yS = 0.9;
  double phic = 1.5707963-phi0; // pi/2 - phi0

  if (phi0 < phiS) 
    return ss*asn_dz(sin(phi0))+jj*5.5361264;
  else {
    double cc = sin(phic);
    double xx = cc*cc;
    double d2 = 0.066987298 + 0.93301270 * xx;
    if (xx < yS*d2) 
      return ss*(2.7680632-asn_dz(cc/sqrt(d2)))+jj*5.5361264;
    else {
      double vv = 0.066987298 * (1.-xx);
      if (vv < xx*d2) return ss*acn_dz(cc);
      else return ss*acn_dz(sqrt(vv/d2)) + jj*5.5361264;
    }
  }
}


// =====================================================================================


double cosmobl::Cosmology::acn_dz (const double cc) const 
{
  double pp = 1.;
  double xx = cc*cc;
  for (int j=1; j<10; j++) {
    if (xx > 0.5) return pp*asn_dz(sqrt(1.-xx));
    double dd = sqrt(0.066987298+0.93301270*xx);
    xx = (sqrt(xx)+dd)/(1.+dd);
    pp *= 2.;
  }
  
  ErrorMsg("Error in cosmobl::Cosmology::acn_dz of Cosmology.cpp: too many half"
	   " argument transformations of cn");
  return 0;
}


// =====================================================================================


double cosmobl::Cosmology::asn_dz (const double ss) const 
{
  //double yA = 0.034856757;
  double yA = 0.153532;

  double yy = ss*ss;
  if (yy < yA) 
    return ss*serf_dz(yy);

  double pp = 1.;
  for (int j=1; j<10; j++) {
    yy /= ((1.+sqrt(1.-yy))*(1.+sqrt(1.-0.93301270*yy)));
    pp *= 2.;
    if (yy < yA)
      return pp*sqrt(yy)*serf_dz(yy);
  }
  
  ErrorMsg("Error in cosmobl::Cosmology::asn_dz: too many half argument"
	   " transformations of sn");
  return 0;
}


// =====================================================================================


double cosmobl::Cosmology::serf_dz (const double yy) const 
{
  return (1.+yy*(0.32216878+yy*(0.18693909+yy*(0.12921048+yy*(0.097305732+yy*(0.077131543+yy*0.063267775+yy*(0.053185339+yy*(0.045545557+yy*0.039573617))))))));
}
