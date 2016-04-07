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

// Sound horizon at drag epoch

double cosmobl::Cosmology::rs (const string method_Pk, const double T_CMB) const
{
  if (method_Pk=="EisensteinHu") 
    return rs_EH(T_CMB);
   
  else if (method_Pk=="CAMB")
    return rs_CAMB();

  else { ErrorMsg("Error in cosmobl::Cosmology::rs of BAO.cpp: 'method_Pk' not allowed!"); return 0; }
}


// =====================================================================================

// Sound horizon at drag epoch (Eisentein & Hu 1998, Section 2.1)

double cosmobl::Cosmology::rs_EH (const double T_CMB) const
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


double cosmobl::Cosmology::rs_CAMB () const
{
  double wcb= m_Omega_matter*pow(m_hh,2);
  double wb= m_Omega_baryon*pow(m_hh,2);
  double wnu = m_Omega_neutrinos*pow(m_hh,2);

  return 55.154*exp(-72.3*pow(wnu+0.0006,2))/(pow(wcb,0.25351)*pow(wb,0.12807));
}


// =====================================================================================

// Fiducial cosmology independent ratio rs/DV (rs,DV [Mpc])

double cosmobl::Cosmology::ys (const double redshift, const string method_Pk, const double T_CMB) const
{
  return rs(method_Pk, T_CMB)/((m_unit) ? D_V(redshift)/m_hh : D_V(redshift));
}


// =====================================================================================

// Acoustic parameter (see Eisenstein 2005)

double cosmobl::Cosmology::Az (const double redshift) const
{
  return ((m_unit) ? D_V(redshift)/m_hh : D_V(redshift))*1.e2*sqrt(m_Omega_matter*m_hh*m_hh)/(par::cc*redshift);
}
