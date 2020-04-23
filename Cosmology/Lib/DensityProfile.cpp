/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  @file Cosmology/Lib/DensityProfile.cpp
 *
 *  @brief Methods of the class Cosmology used to model the halo
 *  density profile
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model halo density profile
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;
using namespace glob;


// =====================================================================================


double cbl::cosmology::Cosmology::concentration (const double Mass, const double redshift, const string author, const string profile, const string halo_def) const
{
  double AA, BB, CC;
  
  if (author=="Duffy") {

    if (redshift>2) ErrorCBL("the concentration-mass relation by Duffy et al. has been tested only at z<2", "concentration", "DensityProfile.cpp");
      
    if (profile=="NFW") {
      
      if (halo_def=="200") {
	AA = 5.71;
	BB = -0.084;
	CC = -0.47;
      }
      
      else if (halo_def=="vir") {
	AA = 7.85;
	BB = -0.081;
	CC = -0.71;
      }
      
      else if (halo_def=="mean") {
	AA = 10.14;
	BB = -0.081;
	CC = -1.01;
      }
      
      else return ErrorCBL("halo_def not allowed!", "concentration", "DensityProfile.cpp");

    }
    
    else if (profile=="Einasto") {

      if (halo_def=="200") {
	AA = 6.4;
	BB = -0.108;
	CC = -0.62;
      }
      
      else if (halo_def=="vir") {
	AA = 8.82;
	BB = -0.106;
	CC = -0.87;
      }
      
      else if (halo_def=="mean") {
	AA = 11.39;
	BB = -0.107;
	CC = -1.16;
      }
      
      else return ErrorCBL("halo_def not allowed!", "concentration", "DensityProfile.cpp");
      
    }

    else return ErrorCBL("profile not allowed!", "concentration", "DensityProfile.cpp");
    
  }
  
  else return ErrorCBL("author not allowed!", "concentration", "DensityProfile.cpp");

  
  const double Mpivot = 2.e12; // in Msun/h
  
  return AA*pow(Mass/Mpivot, BB)*pow(1.+redshift, CC);
}


// ============================================================================================


double cbl::cosmology::Cosmology::density_profile (const double rad, const double Mass, const double redshift, const string model_cM, const string profile, const string halo_def) const
{
  if (profile!="NFW")
    return ErrorCBL("profile not allowed!", "density_profile", "DensityProfile.cpp");
  
  const double conc = concentration(Mass, redshift, model_cM, profile, halo_def);

  const double rho_s = rho_crit(redshift)*Delta_c(redshift)/3.*pow(conc, 3)/(log(1.+conc)-conc/(1.+conc));

  const double r_s = r_vir(Mass, redshift)/conc; 
  
  return rho_s/((rad/r_s)*pow(1.+rad/r_s, 2));
}


// ============================================================================================


double cbl::cosmology::Cosmology::density_profile_FourierSpace (const double kk, const double Mass, const double redshift, const string model_cM, const string profile, const string halo_def) const
{
  const double conc = concentration(Mass, redshift, model_cM, profile, halo_def);

  const double rho_s = rho_crit(redshift)*Delta_c(redshift)/3.*pow(conc, 3)/(log(1.+conc)-conc/(1.+conc));

  const double r_s = r_vir(Mass, redshift)/conc;
  
  const double mu = kk*r_s;  
  
  return 4.*par::pi*rho_s*pow(r_s, 3)/Mass*(cos(mu)*(gsl_sf_Ci(mu+mu*conc)-gsl_sf_Ci(mu))+sin(mu)*(gsl_sf_Si(mu+mu*conc)-gsl_sf_Si(mu))-sin(mu*conc)/(mu+mu*conc));
}


// ============================================================================


double cbl::cosmology::Cosmology::concentration2 (const double Vmax, const double Rmax) const
{ 
  const int nn = 128;
  vector<double> xxi = linear_bin_vector(nn, 0.1, 50.);
  vector<double> yyi(nn);

  for (int i=0; i<nn; i++)
    // reset 200 for the spherical collapse model
    yyi[i] = 200./3.*pow(xxi[i],3.)/(log(1.+xxi[i])-xxi[i]/(1.+xxi[i]))-14.426*pow(Vmax/Rmax/m_H0,2.);
  
  return interpolated(0., yyi, xxi, "Poly");
}


// ============================================================================


double cbl::cosmology::Cosmology::Mass_Delta (const double Mass, const double Delta_in, const double Delta_out, const double conc, const bool is_input_conc, const double rRmin_guess, const double rRmax_guess) const
{
  auto func = [&] (const double xx)
    {
      const double c_in = (is_input_conc) ? conc : conc*xx;
      const double c_out = (!is_input_conc) ? conc : conc/xx;
      const double AA = log(1.+c_out)-c_out/(1.+c_out);
      const double BB = log(1.+c_in)-c_in/(1.+c_in);
      return fabs(Delta_in/Delta_out*AA*pow(xx, 3)-BB);
    };
  
  // R_[Delta_out]/R_[Delta_in]
  const double Rratio = wrapper::gsl::GSL_minimize_1D(func, 1., max(1.e-3, rRmin_guess), rRmax_guess);

  // M_[Delta_out]
  return Delta_out/Delta_in/pow(Rratio, 3)*Mass;
}
