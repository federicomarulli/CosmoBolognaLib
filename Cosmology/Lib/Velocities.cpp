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
 *  @file Cosmology/Lib/Velocities.cpp
 *
 *  @brief Methods of the class Cosmology used to model the peculiar
 *  velocities
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the cosmic peculiar velocities
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::Cosmology::square_bulk_flow (const double rr, const double k_int_min, const string method_Pk, const double redshift, const string output_root, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{
  double bulk = -1.;
  Pk_0(method_Pk, redshift, output_root, k_min, k_max, GSL, prec, file_par); 
  
  if (method_Pk=="EisensteinHu") {
    if (m_sigma8<0) ErrorMsg("Error in cosmobl::Cosmology::square_bulk_flow: sigma8 must be >0 using EisensteinHu!");
    cosmobl::classfunc::func_V2 func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, rr, redshift);
    Midpnt<cosmobl::classfunc::func_V2> q1(func,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_V2> q2(func,1.,1.e30);
    bulk = m_Pk0_EH*(qromo(q1)+qromo(q2));
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    Table_PkCodes (method_Pk, do_nonlinear, lgkk, lgPk, redshift, output_root, k_max, file_par);

    cosmobl::classfunc::func_V2_Table func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, lgkk, lgPk, rr, redshift);
   
    Midpnt<cosmobl::classfunc::func_V2_Table> q1(func,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_V2_Table> q2(func,1.,1.e30);

    double PP0 = (method_Pk=="CAMB") ? m_Pk0_CAMB : m_Pk0_CLASS;
    bulk = PP0*(qromo(q1)+qromo(q2));
  }
  
  return pow(HH(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*bulk;
}


// =====================================================================================


double cosmobl::Cosmology::square_bulk_flow_Table (const double rr, const double k_int_min, const vector<double> lgkk, const vector<double> lgPk, const double redshift) const 
{
  cosmobl::classfunc::func_V2_Table func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, lgkk, lgPk, rr, redshift);
  
  Midpnt<cosmobl::classfunc::func_V2_Table> q1(func,k_int_min,1.); 
  Midinf<cosmobl::classfunc::func_V2_Table> q2(func,1.,1.e30);

  return pow(HH(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*(qromo(q1)+qromo(q2));
}


// =====================================================================================


double cosmobl::Cosmology::square_velocity_dispersion (const double rr, const double k_int_min, const string method_Pk, const double redshift, const string output_root, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{
  double sigma2 = -1.;
  Pk_0(method_Pk, redshift, output_root, k_min, k_max, GSL, prec, file_par); 
  
  if (method_Pk=="EisensteinHu") {
    if (m_sigma8<0) ErrorMsg("Error in cosmobl::Cosmology::square_velocity_dispersion: sigma8 must be >0 using EisensteinHu!");
    cosmobl::classfunc::func_sigma2 func2 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, rr, redshift);
    Midpnt<cosmobl::classfunc::func_sigma2> q11(func2,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_sigma2> q22(func2,1.,100.);
    sigma2 = m_Pk0_EH*(qromo(q11)+qromo(q22));
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    Table_PkCodes (method_Pk, do_nonlinear, lgkk, lgPk, redshift, output_root, k_max, file_par);
 
    cosmobl::classfunc::func_sigma2_Table func2 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, lgkk, lgPk, rr, redshift);

    Midpnt<cosmobl::classfunc::func_sigma2_Table> q11(func2,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_sigma2_Table> q22(func2,1.,100.);

    double PP0 = (method_Pk=="CAMB") ? m_Pk0_CAMB : m_Pk0_CLASS;
    sigma2 = PP0*(qromo(q11)+qromo(q22));
  }

  return pow(HH(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*sigma2;
}


// =====================================================================================


double cosmobl::Cosmology::CMN (const double rr, const double k_int_min, const string method_Pk, const double redshift, const string output_root, const double k_max, const string file_par) const 
{
  double CMN = -1000.; 
  
  if (method_Pk=="EisensteinHu") {
    
    cosmobl::classfunc::func_V2 func1 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, rr, redshift);
    cosmobl::classfunc::func_sigma2 func2 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, method_Pk, rr, redshift);
    
    Midpnt<cosmobl::classfunc::func_V2> q1(func1,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_V2> q2(func1,1.,1.e30);
    double Int1 = qromo(q1)+qromo(q2);
    
    Midpnt<cosmobl::classfunc::func_sigma2> q11(func2,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_sigma2> q22(func2,1.,100.);
    double Int2 = qromo(q11)+qromo(q22);
    
    CMN = Int1/Int2;
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {

    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    Table_PkCodes (method_Pk, do_nonlinear, lgkk, lgPk, redshift, output_root, k_max, file_par);
 
    cosmobl::classfunc::func_V2_Table func1 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, lgkk, lgPk, rr, redshift);
    cosmobl::classfunc::func_sigma2_Table func2 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, lgkk, lgPk, rr, redshift);
    
    Midpnt<cosmobl::classfunc::func_V2_Table> q1(func1,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_V2_Table> q2(func1,1.,1.e30);
    double Int1 = qromo(q1)+qromo(q2);
    
    Midpnt<cosmobl::classfunc::func_sigma2_Table> q11(func2,k_int_min,1.); 
    Midinf<cosmobl::classfunc::func_sigma2_Table> q22(func2,1.,100.);
    double Int2 = qromo(q11)+qromo(q22);
    
    CMN = Int1/Int2;
  }
				       
  return sqrt(CMN);
}
