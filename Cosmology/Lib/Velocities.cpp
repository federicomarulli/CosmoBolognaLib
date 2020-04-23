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
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Cosmology::square_bulk_flow (const double rr, const double k_int_min, const string method_Pk, const double redshift, const bool store_output, const string output_root, const double k_min, const double k_max, const double prec, const string file_par)
{
  double bulk = -1.;
  Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par); 

  function<double(double)> ff;

  if (method_Pk=="EisensteinHu") {
    if (m_sigma8<0) ErrorCBL("sigma8 must be >0 using EisensteinHu!", "square_bulk_flow", "Velocities.cpp");
    cbl::classfunc::func_V2 func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, method_Pk, rr, redshift, store_output);

    ff = bind(&cbl::classfunc::func_V2::operator(), func, std::placeholders::_1);
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    cbl::classfunc::func_V2_Table func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, lgkk, lgPk, rr, redshift);

    ff = bind(&cbl::classfunc::func_V2_Table::operator(), func, std::placeholders::_1);

  }

  double Int1 = wrapper::gsl::GSL_integrate_qag(ff, k_int_min, 1., 1.e-3);
  double Int2 = wrapper::gsl::GSL_integrate_qag(ff, 1., 1.e30, 1.e-3);

  bulk = m_Pk0_EH*(Int1+Int2);
  return pow(HH(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*bulk;
}


// =====================================================================================


double cbl::cosmology::Cosmology::square_bulk_flow_Table (const double rr, const double k_int_min, const vector<double> lgkk, const vector<double> lgPk, const double redshift) const 
{
  cbl::classfunc::func_V2_Table func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, lgkk, lgPk, rr, redshift);

  function<double(double)> ff = bind(&cbl::classfunc::func_V2_Table::operator(), func, std::placeholders::_1);

  double Int1 = wrapper::gsl::GSL_integrate_qag(ff, k_int_min, 1., 1.e-3);
  double Int2 = wrapper::gsl::GSL_integrate_qag(ff, 1., 1.e30, 1.e-3);

  return pow(HH(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*(Int1+Int2);
}


// =====================================================================================


double cbl::cosmology::Cosmology::square_velocity_dispersion (const double rr, const double k_int_min, const string method_Pk, const double redshift, const bool store_output, const string output_root, const double k_min, const double k_max, const double prec, const string file_par)
{
  (void)k_int_min;
  
  double sigma2 = -1.;
  Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par); 
  function<double(double)> ff;

  if (method_Pk=="EisensteinHu") {
    if (m_sigma8<0) ErrorCBL("sigma8 must be >0 using EisensteinHu!", "square_velocity_dispersion", "Velocities.cpp");
    cbl::classfunc::func_sigma2 func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, method_Pk, rr, redshift, store_output);

    ff = bind(&cbl::classfunc::func_sigma2::operator(), func, std::placeholders::_1);

  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);
 
    cbl::classfunc::func_sigma2_Table func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, lgkk, lgPk, rr, redshift);

    ff = bind(&cbl::classfunc::func_sigma2_Table::operator(), func, std::placeholders::_1);
  }

  return pow(HH(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*sigma2;
}


// =====================================================================================


double cbl::cosmology::Cosmology::CMN (const double rr, const double k_int_min, const string method_Pk, const double redshift, const bool store_output, const string output_root, const double k_max, const string file_par) const 
{
  double CMN = -1000.; 

  function<double(double)> ff1;
  function<double(double)> ff2;

  if (method_Pk=="EisensteinHu") {

    cbl::classfunc::func_V2 func1 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, method_Pk, rr, redshift, store_output);
    cbl::classfunc::func_sigma2 func2 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, method_Pk, rr, redshift, store_output);

    ff1 = bind(&cbl::classfunc::func_V2::operator(), func1, std::placeholders::_1);
    ff2 = bind(&cbl::classfunc::func_sigma2::operator(), func2, std::placeholders::_1);
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {

    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    cbl::classfunc::func_V2_Table func1 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, lgkk, lgPk, rr, redshift);
    cbl::classfunc::func_sigma2_Table func2 (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, lgkk, lgPk, rr, redshift);

    ff1 = bind(&cbl::classfunc::func_V2_Table::operator(), func1, std::placeholders::_1);
    ff2 = bind(&cbl::classfunc::func_sigma2_Table::operator(), func2, std::placeholders::_1);

  }

  double i1 = wrapper::gsl::GSL_integrate_qag(ff1, k_int_min, 1.,1.e-3);
  double i2 = wrapper::gsl::GSL_integrate_qag(ff1, 1., 1.e30,1.e-3);
  double Int1 = i1+i2;

  i1 = wrapper::gsl::GSL_integrate_qag(ff2, k_int_min, 1.,1.e-3);
  i2 = wrapper::gsl::GSL_integrate_qag(ff2, 1., 1.e30,1.e-3);
  double Int2 = i1+i2;

  CMN = Int1/Int2;			       
  return sqrt(CMN);
}
