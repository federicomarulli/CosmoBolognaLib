/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Catalogue/Cluster.cpp
 *
 *  @brief Methods of the class Cluster
 *
 *  This file contains the implementation of the methods of the class
 *  Cluster
 *
 *  @author Federico Marulli, Giorgio Lesci
 *
 *  @author federico.marulli3@unibo.it, giorgio.lesci2@unibo.it
 */

#include "Cluster.h"

using namespace std;

using namespace cbl;
using namespace glob;


// =====================================================================================


double cbl::catalogue::Cluster::m_concentration_Duffy(std::vector<double> cosmo_par, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names)
{
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *m_cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<cosmo_par_names.size(); ++i)
    cosmo.set_parameter(cosmo_par_names[i], cosmo_par[i]);

  const double Mpivot = 2.e12; // in Msun/h
  
  return m_A*pow(m_mass/Mpivot, m_B)*pow(1.+m_redshift, m_C);
}


// =====================================================================================


double cbl::catalogue::Cluster::m_NFW(const double conc, const double rad, std::vector<double> cosmo_par, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names)
{
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *m_cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<cosmo_par_names.size(); ++i)
    cosmo.set_parameter(cosmo_par_names[i], cosmo_par[i]);

  // compute rho_s and r_s
  const double rho_s = cosmo.rho_crit(m_redshift)*cosmo.Delta_c(m_redshift)/3.*pow(conc, 3)/(log(1.+conc)-conc/(1.+conc));
  const double r_s = cosmo.r_vir(m_mass, m_redshift)/conc; 
  
  return rho_s/((rad/r_s)*pow(1.+rad/r_s, 2));
}


// =====================================================================================


void cbl::catalogue::Cluster::set_profile (const cbl::cosmology::Cosmology cosmology, const double Mass, const double redshift, const std::string cM_author, const std::string profile_author, const std::string halo_def)
{
  m_isSet_cM_relation = true;
  m_isSet_profile = true;

  m_cosmology = make_shared<cosmology::Cosmology>(cosmology::Cosmology(move(cosmology)));
  m_mass = Mass;
  m_redshift = redshift;
  m_profile_author = profile_author;

  
  // Set the c-M relation
    
  if (cM_author=="Duffy") {
    if (redshift>2) ErrorCBL("the concentration-mass relation by Duffy et al. has been tested only at z<2", "set_profile", "Cluster.cpp");
      
    if (profile_author=="NFW") {      
      if (halo_def=="200") {
	m_A = 5.71;
	m_B = -0.084;
	m_C = -0.47;
      }      
      else if (halo_def=="vir") {
	m_A = 7.85;
	m_B = -0.081;
	m_C = -0.71;
      }      
      else if (halo_def=="mean") {
	m_A = 10.14;
	m_B = -0.081;
	m_C = -1.01;
      }      
      else ErrorCBL("halo_def not allowed!", "set_profile", "Cluster.cpp");
    }
    else if (profile_author=="Einasto") {
      if (halo_def=="200") {
	m_A = 6.4;
	m_B = -0.108;
	m_C = -0.62;
      }      
      else if (halo_def=="vir") {
	m_A = 8.82;
	m_B = -0.106;
	m_C = -0.87;
      }      
      else if (halo_def=="mean") {
	m_A = 11.39;
	m_B = -0.107;
	m_C = -1.16;
      }      
      else ErrorCBL("halo_def not allowed!", "set_profile", "Cluster.cpp");
    }
    else ErrorCBL("profile not allowed!", "set_profile", "Cluster.cpp");
    
    m_concentration_from_mass = &cbl::catalogue::Cluster::m_concentration_Duffy;
  }
  else ErrorCBL("concentration-mass relation author not allowed!", "set_profile", "Cluster.cpp");


  // Set the profile
  
  if (profile_author == "NFW")
    m_profile = &cbl::catalogue::Cluster::m_NFW;
  else
    ErrorCBL("density profile author not allowed!", "set_profile", "Cluster.cpp");
}


// =====================================================================================


void cbl::catalogue::Cluster::set_cosmology (const cbl::cosmology::Cosmology cosmology)
{
  m_cosmology = make_shared<cosmology::Cosmology>(cosmology::Cosmology(move(cosmology)));
}


// =====================================================================================


double cbl::catalogue::Cluster::concentration_from_mass(std::vector<double> cosmo_par, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names)
{
  if (m_isSet_cM_relation == false)
    return ErrorCBL("the concentration-mass relation is not set! Use set_profile().", "concentration", "Cluster.cpp");
  else
    return (this->*m_concentration_from_mass)(cosmo_par, cosmo_par_names);
}


// ============================================================================================


double cbl::catalogue::Cluster::density_profile (const double rad, std::vector<double> cosmo_par, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names)
{  
  double conc = concentration_from_mass(cosmo_par, cosmo_par_names);

  return (this->*m_profile)(conc, rad, cosmo_par, cosmo_par_names);
}


// ============================================================================================


double cbl::catalogue::Cluster::density_profile (const double rad, const double conc, std::vector<double> cosmo_par, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names)
{
  if (m_isSet_profile == false)
    return ErrorCBL("the density profile is not set! Use set_profile().", "density_profile", "Cluster.cpp");
  
  return (this->*m_profile)(conc, rad, cosmo_par, cosmo_par_names);
}


// ============================================================================================


double cbl::catalogue::Cluster::density_profile_FourierSpace (const double kk, std::vector<double> cosmo_par, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_par_names)
{
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *m_cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<cosmo_par_names.size(); ++i)
    cosmo.set_parameter(cosmo_par_names[i], cosmo_par[i]);
  
  double conc = concentration_from_mass(cosmo_par, cosmo_par_names);
  const double rho_s = cosmo.rho_crit(m_redshift)*cosmo.Delta_c(m_redshift)/3.*pow(conc, 3)/(log(1.+conc)-conc/(1.+conc));
  const double r_s = cosmo.r_vir(m_mass, m_redshift)/conc;  
  const double mu = kk*r_s;  
  
  return 4.*par::pi*rho_s*pow(r_s, 3)/m_mass*(cos(mu)*(gsl_sf_Ci(mu+mu*conc)-gsl_sf_Ci(mu))+sin(mu)*(gsl_sf_Si(mu+mu*conc)-gsl_sf_Si(mu))-sin(mu*conc)/(mu+mu*conc));
}


// ============================================================================


double cbl::catalogue::Cluster::concentration2 (const double Vmax, const double Rmax) const
{ 
  const int nn = 128;
  vector<double> xxi = linear_bin_vector(nn, 0.1, 50.);
  vector<double> yyi(nn);

  for (int i=0; i<nn; i++)
    // reset 200 for the spherical collapse model
    yyi[i] = 200./3.*pow(xxi[i],3.)/(log(1.+xxi[i])-xxi[i]/(1.+xxi[i]))-14.426*pow(Vmax/Rmax/m_cosmology->H0(),2.);
  
  return interpolated(0., yyi, xxi, "Poly");
}


// ============================================================================


double cbl::catalogue::Cluster::Mass_Delta (const double Mass, const double Delta_in, const double Delta_out, const double conc, const bool is_input_conc, const double rRmin_guess, const double rRmax_guess) const
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
