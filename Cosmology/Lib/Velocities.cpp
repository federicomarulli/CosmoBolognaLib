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
 *  @brief Methods of the class Velocities used to model the peculiar
 *  velocity statistics of cosmic structure
 *
 *  This file contains the implementation of the methods of the class
 *  Velocities used for calculations of the statistics of the cosmic
 *  velocity field
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Velocities.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Velocities::square_bulk_flow (const double rr, const double k_int_min, const string method_Pk, const double redshift, const bool store_output, const string output_root, const double k_min, const double k_max, const double prec, const string file_par)
{
  PkXi PX(m_cosmology);
  
  double bulk = -1.;
  PX.Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par); 

  function<double(double)> ff;

  if (method_Pk=="EisensteinHu") {
    if (m_cosmology->sigma8()<0) ErrorCBL("sigma8 must be >0 using EisensteinHu!", "square_bulk_flow", "Velocities.cpp");

    PkXi PX(m_cosmology);
    
    function<double(double)> func = [&] (const double kk)
    {
      return pow(m_cosmology->linear_growth_rate(redshift, kk),2)*PX.Pk_matter({kk}, method_Pk, false, redshift, store_output)[0]*pow(cbl::TopHat_WF(kk*rr), 2);   
    };

    ff = func;
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = 0; 

    PkXi PX(m_cosmology);
    
    PX.Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);
    
    function<double(double)> func = [&] (const double kk) 
    { 
      double fact = (m_cosmology->unit()) ? 1. : m_cosmology->little_h();
      double lgk = log10(kk/fact);

      double lgPkK = cbl::interpolated(lgk, lgkk, lgPk, "Linear");

      return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*pow(10., lgPkK)/pow(fact, m_cosmology->n_spec())*pow(cbl::TopHat_WF(kk*rr), 2);  
    };

    ff = func;
  }
  
  double Int1 = wrapper::gsl::GSL_integrate_qag(ff, k_int_min, 1., 1.e-3);
  double Int2 = wrapper::gsl::GSL_integrate_qag(ff, 1., 1.e30, 1.e-3);

  bulk = PX.Pk0_EH()*(Int1+Int2);

  return pow(m_cosmology->Hubble(redshift)/(1.+redshift), 2)/(2.*par::pi*par::pi)*bulk;
}


// =====================================================================================


double cbl::cosmology::Velocities::square_bulk_flow_Table (const double rr, const double k_int_min, const vector<double> lgkk, const vector<double> lgPk, const double redshift) const 
{
  const double fact = (m_cosmology->unit()) ? 1. : m_cosmology->little_h();

  function<double(double)> ff = [&] (const double kk)
  { 
    const double lgk = log10(kk/fact);
    const double lgPkK = cbl::interpolated(lgk, lgkk, lgPk, "Linear");

    return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*pow(10., lgPkK)/pow(fact, m_cosmology->n_spec())*pow(cbl::TopHat_WF(kk*rr), 2);  
  };

  double Int1 = wrapper::gsl::GSL_integrate_qag(ff, k_int_min, 1., 1.e-3);
  double Int2 = wrapper::gsl::GSL_integrate_qag(ff, 1., 1.e30, 1.e-3);

  return pow(m_cosmology->Hubble(redshift)/(1.+redshift),2)/(2.*par::pi*par::pi)*(Int1+Int2);
}


// =====================================================================================


double cbl::cosmology::Velocities::square_velocity_dispersion (const double rr, const double k_int_min, const string method_Pk, const double redshift, const bool store_output, const string output_root, const double k_min, const double k_max, const double prec, const string file_par)
{
  (void)k_int_min;
  
  double sigma2 = -1.;

  PkXi PX(m_cosmology);
  
  PX.Pk_0(method_Pk, redshift, store_output, output_root, k_min, k_max, prec, file_par); 
  function<double(double)> ff;
    
  if (method_Pk=="EisensteinHu") {
    if (m_cosmology->sigma8()<0) ErrorCBL("sigma8 must be >0 using EisensteinHu!", "square_velocity_dispersion", "Velocities.cpp");
    
    function<double(double)> func_sigma2 = [&] (const double kk) 
    {
      return pow(m_cosmology->linear_growth_rate(redshift, kk),2)*PX.Pk_matter({kk}, method_Pk, false, redshift, store_output)[0]*(1.-pow(cbl::TopHat_WF(kk*rr), 2));
    };

    ff = func_sigma2;
  }

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {
    vector<double> lgkk, lgPk;
    bool do_nonlinear = false; 
    
    PX.Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    const double fact = (m_cosmology->unit()) ? 1. : m_cosmology->little_h();
    
    function<double(double)> func_sigma2_Table = [&] (const double kk) 
    {
      const double lgk = log10(kk/fact);
      const double lgPkK = cbl::interpolated(lgk, lgkk, lgPk, "Linear");
      return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*pow(10., lgPkK)/pow(fact, m_cosmology->n_spec())*(1.-pow(cbl::TopHat_WF(kk*rr), 2));  
    };
    
    ff = func_sigma2_Table;
  }

  return pow(m_cosmology->Hubble(redshift)/(1.+redshift), 2)/(2.*par::pi*par::pi)*sigma2;
}


// =====================================================================================


double cbl::cosmology::Velocities::CMN (const double rr, const double k_int_min, const string method_Pk, const double redshift, const bool store_output, const string output_root, const double k_max, const string file_par) const 
{
  double CMN = -1000.; 

  PkXi PX(m_cosmology);

  function<double(double)> ff1 = [&] (const double kk) 
  { 
    return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*PX.Pk_matter({kk}, method_Pk, false, redshift, store_output)[0]*pow(cbl::TopHat_WF(kk*rr), 2);   
  };
  

  function<double(double)> ff2 = [&] (const double kk) 
  {
    return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*PX.Pk_matter({kk}, method_Pk, false, redshift, store_output)[0]*(1.-pow(cbl::TopHat_WF(kk*rr), 2));
  };
    

  if (method_Pk=="CAMB" || method_Pk=="CLASS") {

    vector<double> lgkk, lgPk;
    bool do_nonlinear = false; 

    PkXi PX(m_cosmology);
    
    PX.Table_PkCodes(method_Pk, do_nonlinear, lgkk, lgPk, redshift, store_output, output_root, k_max, file_par);

    const double fact = (m_cosmology->unit()) ? 1. : m_cosmology->little_h();
    
    function<double(double)> func_V2_Table = [&] (const double kk) 
    { 
      const double lgk = log10(kk/fact);
      const double lgPkK = cbl::interpolated(lgk, lgkk, lgPk, "Linear");
      
      return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*pow(10.,lgPkK)/pow(fact, m_cosmology->n_spec())*pow(cbl::TopHat_WF(kk*rr), 2);  
    };

    function<double(double)> func_sigma2_Table = [&] (const double kk) 
    { 
      const double lgk = log10(kk/fact);
      const double lgPkK = cbl::interpolated(lgk, lgkk, lgPk, "Linear");
      
      return pow(m_cosmology->linear_growth_rate(redshift, kk), 2)*pow(10., lgPkK)/pow(fact, m_cosmology->n_spec())*pow(cbl::TopHat_WF(kk*rr), 2);  
    };

    ff1 = func_V2_Table;
    ff2 = func_sigma2_Table;
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
