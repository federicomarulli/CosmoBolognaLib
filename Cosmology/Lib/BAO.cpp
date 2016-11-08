/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Cosmology/Lib/BAO.cpp
 *
 *  @brief Methods of the class Cosmology used to model the BAO
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the Baryon Acoustic Oscillations (BAO)
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


//redshift at wich occurs baryon photon decoupling, see Hu & Sugiyama (1996).
double cosmobl::cosmology::Cosmology::z_decoupling() const
{
  double ombh2 = m_Omega_baryon*m_hh*m_hh;
  double omdmh2 = m_Omega_CDM*m_hh*m_hh;
  double g1 = 0.0783*pow(ombh2,-0.238)/(1.+39.5*pow(ombh2,0.763));
  double g2 = 0.560/(1+21.1*pow(ombh2,1.81));
  double zdec = 1048*(1.+0.00124*pow(ombh2,-0.738))*(1.+g1*pow(ombh2+omdmh2,g2));
  return zdec;

}


// =====================================================================================


//redshift at wich occurs baryon photon decoupling, see Hu & Sugiyama (1996).
double cosmobl::cosmology::Cosmology::z_drag() const
{
  double wb = m_Omega_baryon*m_hh*m_hh;
  double wm = m_Omega_matter*m_hh*m_hh;

  double b1 = 0.313*pow(wm,-0.419)*(1+0.607*pow(wm,0.674));
  double b2 = 0.238*pow(wm,0.223);
  double zd = 1291*pow(wm,0.251)*(1+b1*pow(wb,b2))/(1+0.659*pow(wm,0.828));
  return zd;
}


// =====================================================================================

// Sound horizon at drag epoch

double cosmobl::cosmology::Cosmology::rs (const string method_Pk, const double T_CMB) const
{
  if (method_Pk=="EisensteinHu") 
    return rs_EH(T_CMB);
   
  else if (method_Pk=="CAMB")
    return rs_CAMB();

  else
    return ErrorCBL("Error in cosmobl::cosmology::Cosmology::rs of BAO.cpp: 'method_Pk' not allowed!");
}


// =====================================================================================

// Sound horizon at drag epoch (Eisentein & Hu 1998, Section 2.1)

double cosmobl::cosmology::Cosmology::rs_EH (const double T_CMB) const
{
  double Om0h2 = m_Omega_matter*pow(m_hh,2);
  double Ombh2 = m_Omega_baryon*pow(m_hh,2);
  double Tratio = T_CMB/2.7;

  double zeq = 2.5*pow(10.,4)*Om0h2*pow(Tratio,-4);    // Matter radiation Equivalence (explicit function?)
  double keq = 7.46*pow(10.,-2)*Om0h2*pow(Tratio,-2);  // Particle Horizon at zeq

  double b1 = 0.313 * pow(Om0h2,-0.419) * (1.+0.607*pow(Om0h2,0.674));
  double b2 = 0.238 * pow(Om0h2,0.223);

  double zdrag = 1291.*(pow(Om0h2,0.251)/(1.+0.659*pow(Om0h2,0.828)))*(1.+b1*pow(Ombh2,b2));

  double Rd = 31.5*Ombh2*pow(Tratio,-4)*pow(10.,3)/zdrag;
  double Req = 31.5*Ombh2*pow(Tratio,-4)*pow(10.,3)/zeq;

  return 2.*pow(3.*keq,-1)*pow(6./Req,0.5)*log((pow(1.+Rd,0.5)+pow(Rd+Req,0.5))/(1.+pow(Req,0.5)));
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::rs_CAMB () const
{
  double wcb= m_Omega_matter*pow(m_hh,2);
  double wb= m_Omega_baryon*pow(m_hh,2);
  double wnu = m_Omega_neutrinos*pow(m_hh,2);

  double rd =  55.154*exp(-72.3*pow(wnu+0.0006,2))/(pow(wcb,0.25351)*pow(wb,0.12807));
  return ((m_unit) ? rd*m_hh : rd);
}


// =====================================================================================

// Fiducial cosmology independent ratio rs/DV (rs,DV [Mpc])

double cosmobl::cosmology::Cosmology::ys (const double redshift, const string method_Pk, const double T_CMB) const
{
  return rs(method_Pk, T_CMB)/((m_unit) ? D_V(redshift)/m_hh : D_V(redshift));
}


// =====================================================================================

// Acoustic parameter (see Eisenstein 2005)

double cosmobl::cosmology::Cosmology::Az (const double redshift) const
{
  return ((m_unit) ? D_V(redshift)/m_hh : D_V(redshift))*1.e2*sqrt(m_Omega_matter*m_hh*m_hh)/(par::cc*redshift);
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::sound_speed(const double redshift, const double T_CMB) const
{
  double R = 31500.*m_Omega_baryon*m_hh*m_hh*pow(T_CMB/2.7,-4)/(1+redshift);
  double cs = 1./sqrt(3*(1+R));
  return par::cc*cs;
}


// =====================================================================================
// Sound horizon integrand

double cosmobl::cosmology::Cosmology::rs_integrand (const double a, const double T_CMB) const
{
  double redshift=1./a-1;

  double zeq = 2.5e4*m_Omega_matter*m_hh*m_hh*pow(T_CMB/2.7,-4);
  double a_eq = 1./(1+zeq);

  double factor;
  if(m_Omega_radiation ==0)
    factor = sqrt(m_Omega_matter*(a+a_eq)+m_Omega_k*a*a+m_Omega_DE*f_DE(redshift)*pow(a,4));
  else
    factor = a*a*EE(redshift);

  return sound_speed(redshift, T_CMB)/factor;
}


// =====================================================================================
// Sound horizon 

double cosmobl::cosmology::Cosmology::rs (const double redshift, const double T_CMB) const
{
  function<double(double)> integrand = bind(&Cosmology::rs_integrand, this, std::placeholders::_1, T_CMB);
  double a = 1./(1+redshift);
  return GSL_integrate_qag(integrand,0, a)/m_H0;
}
