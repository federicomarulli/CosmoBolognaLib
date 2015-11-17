/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Tommaso Ronconi      *
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
 *  @file Cosmology/Lib/SizeFunction.cpp
 *
 *  @brief Methods of the class Cosmology used to model the mass
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the mass function of dark matter haloes
 *
 *  @authors Federico Marulli, Tommaso Ronconi
 *
 *  @authors federico.marulli3@unbo.it, tommaso.ronconi@studio.unibo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::Cosmology::VolS (const double R) const
{
  return 4./3.*par::pi*pow(R,3);
}


// =====================================================================================


double cosmobl::Cosmology::deltav (const double rho_vm) const
{
  return 1.594*(1. - pow((rho_vm),(-1./1.594)));
}


// =====================================================================================


double cosmobl::Cosmology::r_rL (const double rho_vm) const
{
  return pow(rho_vm, -1./3.);
}


// =====================================================================================


double cosmobl::Cosmology::f_nu (const double SS, const double del_v, const double del_c) const
{	
  double radnu = fabs(del_v)/SS;
  double nu = pow(radnu, 2.);
  double DDD = fabs(del_v)/(del_c + fabs(del_v));
  double xx = DDD/radnu;
	
  if (xx <= 0.276) 
    return sqrt(2./par::pi)*radnu*exp(-0.5*nu);
  
  else {
    double ff = 0.;
    for (int j = 1; j <= 4; j++)
      ff += 2.*exp(-pow((j*par::pi*xx), 2.)/2.)*j*par::pi*pow(xx, 2.)*sin(j*par::pi*DDD);
    return ff;
  }
}


// =====================================================================================


double cosmobl::Cosmology::size_function (const double R, const double redshift, const double rho_vm, const double del_v, const double del_c, const string method_Pk, const string output_root, const string interpType, const int Num, const double stepsize, const double k_max, const string file_par) const
{
  double Z0 = 0.;
  double zero = 0.;
  double RR = R/r_rL(rho_vm);
  double sigmaR = sqrt(SSR_norm(RR, method_Pk, zero, output_root, k_max, file_par));
  double sigmaRz = sigmaR*DD(redshift)/DD(Z0);
  double SSSR = sigmaRz*sigmaRz;
	
  double Dln_SigmaR = dnSR(1, RR, method_Pk, redshift, output_root, interpType, Num, stepsize, k_max, file_par)*(RR/(2.*SSSR))*DD(redshift)/DD(Z0);

  return f_nu(sigmaRz, del_v, del_c)/VolS(R)*fabs(Dln_SigmaR);
}
