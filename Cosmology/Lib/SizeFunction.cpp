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

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::Cosmology::deltav_L (const double bias, const double rho_vm) const
{
  return 1.594*(1.-pow(1+(rho_vm-1.)/bias, (-1./1.594)));
}


// =====================================================================================


double cbl::cosmology::Cosmology::deltav_NL (const double deltav) const
{
  return pow((1.-deltav/1.594), -1.594) - 1.;
}


// =====================================================================================


double cbl::cosmology::Cosmology::r_rL (const double deltav) const
{
  return pow(pow((1.-deltav/1.594), -1.594), -1./3.);
}


// =====================================================================================


double cbl::cosmology::Cosmology::f_nu (const double SS, const double del_v, const double del_c) const
{	
  double radnu = fabs(del_v)/SS;
  double nu = pow(radnu, 2.);
  double DDD = fabs(del_v)/(del_c+fabs(del_v));
  double xx = DDD/radnu;
	
  if (xx<=0.276) 
    return sqrt(2./par::pi)*radnu*exp(-0.5*nu);
  
  else {
    double ff = 0.;
    for (int j = 1; j <= 4; j++)
      ff += 2.*exp(-pow((j*par::pi*xx), 2.)/2.)*j*par::pi*pow(xx, 2.)*sin(j*par::pi*DDD);
    return ff;
  }
}


// =====================================================================================


double cbl::cosmology::Cosmology::size_function (const double RV, const double redshift, const double del_v, const double del_c, const string model, const string method_Pk, const string output_root, const string interpType, const double k_max, const string input_file, const bool is_parameter_file) const
{
  double RL;
  if ((model == "Vdn") || (model == "SvdW"))
    RL = RV/r_rL(del_v);
  
  else if (model == "linear")
    RL = RV;
  
  else 
    { ErrorCBL("Error in cbl::cosmology::Cosmology::size_function of SizeFunction.cpp: model name not allowed! Allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)"); return 0; }
  
  double sigmaR = sqrt(sigma2R(RL, method_Pk, 0., output_root, interpType, k_max, input_file, is_parameter_file));
  double sigmaRz = sigmaR*DD(redshift)/DD(0.);
  double SSSR = sigmaRz*sigmaRz;
        
  double Dln_SigmaR = dnsigma2R(1, RL, method_Pk, 0., output_root, interpType, k_max, input_file, is_parameter_file)*(RL/(2.*SSSR))*pow(DD(redshift)/DD(0.), 2.);
  
  if (model == "Vdn")
    return f_nu(sigmaRz, del_v, del_c)/volume_sphere(RV)*fabs(Dln_SigmaR);
  
  else if ((model == "SvdW") || (model == "linear"))
    return f_nu(sigmaRz, del_v, del_c)/volume_sphere(RL)*fabs(Dln_SigmaR);
  
  else 
    { ErrorCBL("Error in cbl::cosmology::Cosmology::size_function of SizeFunction.cpp: model name not allowed! Allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)"); return 0; }
}


// =====================================================================================


double cbl::cosmology::Cosmology::size_function (const double RV, const double redshift, const string model_mf, const double del_v, const string model_sf, const string method_Pk, const string output_root, const double Delta, const string interpType, const int norm, const double k_min, const double k_max, const double prec, const string input_file, const bool is_parameter_file)
{
  double RL;
  
  if ((model_sf == "Vdn") || (model_sf == "SvdW"))
    RL = RV/r_rL(del_v);
  
  else if (model_sf == "linear")
    RL = RV;
  
  else 
    { ErrorCBL("Error in cbl::cosmology::Cosmology::size_function of SizeFunction.cpp: model name not allowed! Allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)"); return 0; }
  
  double RHO = rho_m(redshift, true);
  double MM = Mass(RL, RHO);
  
  if (model_sf == "Vdn")
    return 3.*MM*cosmology::Cosmology::mass_function(MM, redshift, model_mf, method_Pk, output_root, Delta, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file, false, del_v)*deltav_NL(del_v);
  
  else if ((model_sf == "SvdW") || (model_sf == "linear"))
    return 3.*MM*cosmology::Cosmology::mass_function(MM, redshift, model_mf, method_Pk, output_root, Delta, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file, false, del_v);
  
  else 
    { ErrorCBL("Error in cbl::cosmology::Cosmology::size_function of SizeFunction.cpp: model name not allowed! Allowed names are: SvdW (Sheth and van de Weygaert, 2004), linear/Vdn (Jennings, Li and Hu, 2013)"); return 0; }

}
