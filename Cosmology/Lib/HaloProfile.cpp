/********************************************************************
 *  Copyright (C) 2022 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Cosmology/Lib/HaloProfile.cpp
 *
 *  @brief Methods of the class HaloProfile
 *
 *  This file contains the implementation of the methods of the class
 *  HaloProfile
 *
 *  @author Giorgio Lesci, Federico Marulli
 *
 *  @author giorgio.lesci2@unibo.it, federico.marulli3@unibo.it
 */

#include "HaloProfile.h"

using namespace std;

using namespace cbl;
using namespace glob;


// =====================================================================================


cbl::cosmology::HaloProfile::HaloProfile (const cbl::cosmology::Cosmology cosmology, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact, const bool miscentering, const bool single_profile, const double sigma_off, const double f_off)
{
  m_set_profile(cosmology, redshift, conc, Mass, Delta, profile_author, halo_def, trunc_fact, miscentering, single_profile, sigma_off, f_off);
}


// =====================================================================================


cbl::cosmology::HaloProfile::HaloProfile (const cbl::cosmology::Cosmology cosmology, const double redshift, const std::string cM_author, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact, const bool miscentering, const bool single_profile, const double sigma_off, const double f_off)
{
  m_set_profile(cosmology, redshift, 0., Mass, Delta, profile_author, halo_def, trunc_fact, miscentering, single_profile, sigma_off, f_off);
  m_set_cM_relation(cM_author);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_concentration_Duffy ()
{  
  const double fact = (m_cosmology->unit()) ? 1 : m_cosmology->hh();
  const double Mpivot = 2.e12/fact; // in Msun/h
  
  return m_A*pow(m_mass/Mpivot, m_B)*pow(1.+m_redshift, m_C);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_rho_s_NFW ()
{
  // Compute rho_s
  const double conc = (this->*m_return_concentration)();
  const double rho_s =  pow(conc, 3) * m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift) / (3.*(log(1.+conc)-conc/(1.+conc)));
  
  return rho_s;
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_rho_s_NFW_trunc ()
{
  // compute rho_s
  const double conc = (this->*m_return_concentration)();  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double r_t = m_trunc_fact * r_s*conc;
  const double tau = r_t/r_s;
  const double fact = tau*tau / (2.*pow(tau*tau+1.,3)*(1.+conc)*(tau*tau+conc*conc)) * (conc*(tau*tau+1.) * (conc*(conc+1.) - tau*tau*(conc-1.)*(2.+3.*conc) - 2.*pow(tau,4)) + tau*(conc+1.)*(tau*tau+conc*conc) * (2.*(3.*tau*tau-1.) * std::atan(conc/tau) + tau*(tau*tau-3.)*log(tau*tau*(1.+conc)*(1.+conc)/(tau*tau+conc*conc)))); // Eq. 10 from Oguri & Hamana 2011
  
  const double rho_s =  pow(conc, 3) * m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift) / (3.*fact);
  
  return rho_s;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_rho_NFW (const std::vector<double> rad)
{
  // compute rho_s and r_s
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double rho_s = m_rho_s_NFW();

  // compute rho
  std::vector<double> rho(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    rho[i] = rho_s/((rad[i]/r_s)*pow(1.+rad[i]/r_s, 2));
  
  return rho;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_rho_NFW_trunc (const std::vector<double> rad)
{
  // compute rho_s and r_s
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double r_t = m_trunc_fact * r_s*conc;
  const double rho_s = m_rho_s_NFW_trunc();

  // compute rho
  std::vector<double> rho(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    rho[i] = rho_s/((rad[i]/r_s)*pow(1.+rad[i]/r_s, 2)) * pow(r_t*r_t/(rad[i]*rad[i]+r_t*r_t),2);
  
  return rho;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_NFW (const std::vector<double> rad)
{
  // Implementation of second part of eq. 4 from Golse et al. 2002 (https://ui.adsabs.harvard.edu/abs/2002A%26A...390..821G/abstract)
  std::function<double(double)> F = [] (double x) { 
				      if (x < 1.)
					return (1 - std::acosh(1./x)/std::sqrt(1.-x*x)) / (x*x-1.);
				      else if (x == 1.)
					return 1./3.;
				      else
					return (1 - std::acos(1./x)/std::sqrt(x*x-1.)) / (x*x-1.);
				    };

  // compute the density
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double rho_s = m_rho_s_NFW();

  std::vector<double> Sigma(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    Sigma[i] = 2. * rho_s * r_s * F(rad[i]/r_s) * 1.e-12;
  
  return Sigma;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_mean_NFW (const std::vector<double> rad)
{
  // Implementation of eq. 5 from Golse et al. 2002
  std::function<double(double)> G = [] (double x) { 
				      if (x < 1.)
					return log(x/2.) + std::acosh(1./x) / sqrt(1.-x*x);
				      else if (x == 1.)
					return 1+log(1./2.);
				      else
					return log(x/2.) + std::acos(1./x) / sqrt(x*x-1.);;
				    };

  // compute the density
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double rho_s = m_rho_s_NFW();

  std::vector<double> Sigma_mean(rad.size());
  for (size_t i=0; i<rad.size(); i++) {
    double x = rad[i]/r_s;
    Sigma_mean[i] = 4. * rho_s * r_s * G(x)/x/x * 1.e-12;
  }
  
  return Sigma_mean;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_NFW_trunc (const std::vector<double> rad)
{
  // Implementation of eq. A.5 from Baltz et al. 2009 (https://ui.adsabs.harvard.edu/abs/2009JCAP...01..015B/abstract)  
  const double conc = (this->*m_return_concentration)();
  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double r_t = m_trunc_fact * r_s*conc;
  const double tau = r_t/r_s;
  const double rho_s = m_rho_s_NFW_trunc();
  
  std::function<double(double)> F = [] (double x){
				      if (x < 1.)
					return std::acosh(1./x) / std::sqrt(1.-x*x);
				      else if (x == 1.)
					return 1.;
				      else
					return std::acos(1./x) / std::sqrt(x*x-1.);
				    };

  // Implementation of eq. A.28 from Baltz et al. 2009
  std::function<double(double)> G = [&F] (double x){
				      if (x < 1.)
					return (F(x)-1.) / (1.-x*x);
				      else if (x == 1.)
					return 1./3.;
				      else
					return (1.-F(x)) / (x*x-1.);
				    };

  // Implementation of second part of eq. A.6 from Baltz et al. 2009
  std::function<double(double)> L = [&tau] (double x){
				      return log(x/(sqrt(x*x+tau*tau)+tau));
				    };

  // compute the density
  std::vector<double> Sigma(rad.size());
  double tausq = tau*tau;
  double tau4th = tausq*tausq;
  
  for (size_t i=0; i<rad.size(); i++) {
    double x = rad[i]/r_s;
    double tausq_xsq = tausq+x*x;

    double prefact = rho_s*r_s*tau4th/(tausq+1.)/(tausq+1.)/(tausq+1.);
    double a = 2*(tausq+1.)*G(x);
    double b = 8*F(x);
    double c = (tau4th-1.)/tausq/tausq_xsq;
    double d = -cbl::par::pi*(4*tausq_xsq+tausq+1.)/pow(tausq_xsq,1.5);
    double e = (tausq*(tau4th-1.)+(tausq_xsq)*(3*tau4th-6*tausq-1))*L(x)/(tau*tau*tau)/pow(tausq_xsq,1.5);
    
    Sigma[i] = prefact*(a+b+c+d+e) * 1.e-12;
  }
  
  return Sigma;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_mean_NFW_trunc (const std::vector<double> rad)
{
  const double conc = (this->*m_return_concentration)();
  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta,*m_cosmology,m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. ) / conc;
  const double r_t = m_trunc_fact * r_s*conc;
  const double tau = r_t/r_s;
  const double rho_s = m_rho_s_NFW_trunc();
  
  std::function<double(double)> F = [] (double x) {
				      if (x < 1.)
					return std::acosh(1./x) / std::sqrt(1.-x*x);
				      else if (x == 1.)
					return 1.;
				      else
					return std::acos(1./x) / std::sqrt(x*x-1.);
				    };

  std::function<double(double)> L = [&tau] (double x) {
				      return log(x/(sqrt(x*x+tau*tau)+tau));
				    };

  // compute the density
  std::vector<double> Sigma_mean(rad.size());
  double tausq = tau*tau;
  double tau4th = tausq*tausq;
  
  for (size_t i=0; i<rad.size(); i++) {
    double x = rad[i]/r_s;
    double tausq_xsq = tausq+x*x;

    double prefact = 2. * cbl::par::pi * rho_s * pow(r_s,3);
    double term1 = tau4th / (tausq+1.) / (tausq+1.) / (tausq+1.);
    double term2 = 2. * (tausq + 1. + 4.*(x*x-1.)) * F(x);
    double term3 = (cbl::par::pi * (3.*tausq-1.) + 2. * tau * (tausq-3.) * log(tau)) / tau;
    double term4 = tausq * tau * sqrt(tausq_xsq);
    double term5 = - tausq * tau * cbl::par::pi * (4. * tausq_xsq - tausq - 1.);
    double term6 = - tausq * (tau4th - 1.) + tausq_xsq * (3. * tau4th - 6. * tausq - 1.);

    double M_proj = prefact * term1 * (term2 + term3 + (term5 + term6 * L(x)) / term4);
    
    Sigma_mean[i] = M_proj / (cbl::par::pi * rad[i]*rad[i]) * 1.e-12;
  }
  
  return Sigma_mean;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_mis (const std::vector<double> rad)
{
  // Define the integrands (where P(R_s) is a Rayleigh distribution)
  
  double dummy_R_s, dummy_R;
  
  std::function<double(double)> integrand = [&] (double theta) {
					      double new_rad2 = dummy_R*dummy_R + dummy_R_s*dummy_R_s + 2*dummy_R*dummy_R_s*std::cos(theta);					      
					      return (this->*m_Sigma_cen_ptr)({std::sqrt(new_rad2)})[0];
					    };


  double r=0, a_meantau=0, b_meantau=0, inv_mis_scale = 1./m_sigma_off/m_sigma_off;
  std::function<double(double)> tau_factor1 = [&a_meantau,&b_meantau,&inv_mis_scale] (double s) {
						return (exp(-0.5*(s*s+a_meantau)*inv_mis_scale))/sqrt(b_meantau-a_meantau-s*s);
					      };
  
  std::function<double(double)> tau_factor2 = [&a_meantau,&b_meantau,&inv_mis_scale] (double s){
						return (exp(-0.5*(b_meantau-s*s)*inv_mis_scale))/sqrt(b_meantau-a_meantau-s*s);
					      };
  
  std::function<double(double)> sigma_off_int = [&] (double tau) {
						  a_meantau = (tau-r)*(tau-r);
						  b_meantau = (tau+r)*(tau+r);
						  double c_meantau = (a_meantau+b_meantau)*0.5;
						  double upper_limit = std::min(sqrt(c_meantau-a_meantau),sqrt(25*m_sigma_off*m_sigma_off-a_meantau));
						  double sigma_off_int = r*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor1, 0., upper_limit);
						  if (c_meantau < 25*m_sigma_off*m_sigma_off)
						    {
						      double int_left = std::max(0.,sqrt(b_meantau-25*m_sigma_off*m_sigma_off));
						      double int_right = sqrt(b_meantau-c_meantau);						     
						      sigma_off_int += r*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor2, int_left, int_right);
						    }
						  return sigma_off_int * tau * (this->*m_Sigma_cen_ptr)({tau})[0];
						};

  // Compute Sigma_mis
  std::vector<double> Sigma_mis(rad.size());

  if (m_single_profile) {

    dummy_R_s = m_sigma_off;
    
    for (size_t i=0; i<rad.size(); i++) {
      dummy_R = rad[i];
      const double integral = cbl::wrapper::gsl::GSL_integrate_cquad(integrand, 0., cbl::par::pi);
      Sigma_mis[i] = integral / (cbl::par::pi);
    }
    
  }
  
  else {

    for (size_t i=0; i<rad.size(); i++) {
      
      r = rad[i];
      Sigma_mis[i] = cbl::wrapper::gsl::GSL_integrate_cquad(sigma_off_int, std::max(0., rad[i]-5*m_sigma_off), rad[i]+5.*m_sigma_off) /rad[i]/cbl::par::pi/m_sigma_off/m_sigma_off;
    
    }

  }
  
  return Sigma_mis;
  
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_including_miscentering (const std::vector<double> rad)
{
  if (m_sigma_off > m_sigma_off_threshold && m_f_off > m_f_off_threshold) {
    const std::vector<double> Sigma_cen = (this->*m_Sigma_cen_ptr)(rad);
    const std::vector<double> Sigma_mis = (this->*m_Sigma_mis_ptr)(rad);
  
    std::vector<double> Sigma_tot(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      Sigma_tot[i] = (1.-m_f_off)*Sigma_cen[i] + m_f_off*Sigma_mis[i];

    return Sigma_tot;
  }
  else
    return (this->*m_Sigma_cen_ptr)(rad);
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_DeltaSigma_cen (const std::vector<double> rad)
{
  const std::vector<double> Sigma = (this->*m_Sigma_cen_ptr)(rad);
  const std::vector<double> Sigma_mean = (this->*m_Sigma_mean_cen_ptr)(rad);  

  // Compute DeltaSigma
  std::vector<double> DeltaSigma(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    DeltaSigma[i] = Sigma_mean[i] - Sigma[i];
  
  return DeltaSigma;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_DeltaSigma_mis (const std::vector<double> rad)
{
  if (m_sigma_off > m_sigma_off_threshold) {
    
    // Interpolate Sigma(R)
    const std::vector<double> rad_forInterp = cbl::logarithmic_bin_vector((size_t)(100), m_min_Smis, m_max_Smis);
    const std::vector<double> Sigma = (this->*m_Sigma_mis_ptr)(rad_forInterp);
    
    cbl::glob::FuncGrid Sigma_interp (rad_forInterp, Sigma, "Spline");

    // Compute the mean surface density
    std::function<double(double)> mean_sigma_integrand = [&Sigma_interp] (double rad){
							   return Sigma_interp(rad) * rad;
							 };
  
    std::vector<double> Sigma_mean(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      Sigma_mean[i] = cbl::wrapper::gsl::GSL_integrate_cquad(mean_sigma_integrand, 0., rad[i]) * 2./rad[i]/rad[i];

    // Compute DeltaSigma
    std::vector<double> DeltaSigma(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      DeltaSigma[i] = Sigma_mean[i] - Sigma_interp(rad[i]);
  
    return DeltaSigma;

  }
  else
    return (this->*m_DeltaSigma_cen_ptr)(rad);
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_DeltaSigma_including_miscentering (const std::vector<double> rad)
{
  if (m_sigma_off > m_sigma_off_threshold && m_f_off > m_f_off_threshold) {
    const std::vector<double> DeltaSigma_cen = (this->*m_DeltaSigma_cen_ptr)(rad);
    const std::vector<double> DeltaSigma_mis = (this->*m_DeltaSigma_mis_ptr)(rad);
  
    std::vector<double> DeltaSigma_tot(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      DeltaSigma_tot[i] = (1.-m_f_off)*DeltaSigma_cen[i] + m_f_off*DeltaSigma_mis[i];

    return DeltaSigma_tot;
  }
  else
    return (this->*m_DeltaSigma_cen_ptr)(rad);
}


// =====================================================================================


void cbl::cosmology::HaloProfile::m_set_cM_relation (const std::string cM_author)
{
  m_isSet_cM_relation = true;

  if (cM_author=="Duffy") {
    if (m_redshift>2) ErrorCBL("the concentration-mass relation by Duffy et al. has been tested only at z<2", "set_cM_relation", "HaloProfile.cpp");

    if ( (m_halo_def == "critical" && m_Delta != 200.) || (m_halo_def == "mean" && 200.*m_cosmology->OmegaM(m_redshift) != m_Delta_func(m_Delta,*m_cosmology,m_redshift)) )
      ErrorCBL("the concentration-mass relation by Duffy et al. has been implemented only for critical/mean overdensity factors equal to 200.", "set_cM_relation", "HaloProfile.cpp");
      
    if (m_profile_author=="NFW" || m_profile_author=="NFW_trunc") {      
      if (m_halo_def=="critical") {
	m_A = 5.71;
	m_B = -0.084;
	m_C = -0.47;
      }      
      else if (m_halo_def=="vir") {
	m_A = 7.85;
	m_B = -0.081;
	m_C = -0.71;
      }      
      else if (m_halo_def=="mean") {
	m_A = 10.14;
	m_B = -0.081;
	m_C = -1.01;
      }      
      else ErrorCBL("halo_def not allowed!", "set_cM_relation", "HaloProfile.cpp");
    }
    else if (m_profile_author=="Einasto") {
      if (m_halo_def=="critical") {
	m_A = 6.4;
	m_B = -0.108;
	m_C = -0.62;
      }      
      else if (m_halo_def=="vir") {
	m_A = 8.82;
	m_B = -0.106;
	m_C = -0.87;
      }      
      else if (m_halo_def=="mean") {
	m_A = 11.39;
	m_B = -0.107;
	m_C = -1.16;
      }      
      else ErrorCBL("halo_def not allowed!", "set_cM_relation", "HaloProfile.cpp");
    }
    else ErrorCBL("profile not allowed!", "set_cM_relation", "HaloProfile.cpp");
    
    m_return_concentration = &cbl::cosmology::HaloProfile::m_concentration_Duffy;
  }
  else ErrorCBL("concentration-mass relation author not allowed!", "set_cM_relation", "HaloProfile.cpp");
}


// =====================================================================================


void cbl::cosmology::HaloProfile::m_set_profile (const cbl::cosmology::Cosmology cosmology, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact, const bool miscentering, const bool single_profile, const double sigma_off, const double f_off)
{
  if (miscentering && single_profile && f_off != 1.)
    ErrorCBL("if a single profile is considered in the miscentering model, then f_off must be set equal to 1.", "m_set_profile", "HaloProfile.cpp");

  m_sigma_off_threshold = 1.e-4;
  m_f_off_threshold = 1.e-4;
  m_min_Smis = 1.e-5;
  m_max_Smis = 300;

  m_cosmology = make_shared<cosmology::Cosmology>(cosmology::Cosmology(move(cosmology)));
  m_redshift = redshift;
  m_concentration = conc;
  m_mass = Mass;
  m_Delta = Delta;
  m_profile_author = profile_author;
  m_halo_def = halo_def;
  
  m_trunc_fact = trunc_fact;
  m_miscentering = miscentering;
  m_single_profile = single_profile;
  m_sigma_off = sigma_off;
  m_f_off = f_off;

  // Set the overdensity factor Delta
  if (halo_def == "critical")
    m_Delta_func = [] (const double delta, cbl::cosmology::Cosmology cosmo, const double redshift) {
		     (void)cosmo; (void)redshift;
		     return delta;};
  else if (halo_def == "mean")
    m_Delta_func = [] (const double delta, cbl::cosmology::Cosmology cosmo, const double redshift) {
		     return delta * cosmo.OmegaM(redshift);};
  else if (halo_def == "vir")
    m_Delta_func = [] (const double delta, cbl::cosmology::Cosmology cosmo, const double redshift) {
		     (void)delta;
		     return cosmo.Delta_c(redshift, "BryanNorman");};
  else
    ErrorCBL("wrong halo_def declaration!", "m_set_profile", "HaloProfile.cpp");

  // Set the function returning the concentration
  m_return_concentration = &cbl::cosmology::HaloProfile::m_return_set_concentration;

  // Set the profile  
  if (profile_author == "NFW") {
    
    m_rho_ptr = &cbl::cosmology::HaloProfile::m_rho_NFW;
    m_rho_s_ptr = &cbl::cosmology::HaloProfile::m_rho_s_NFW;
    
    if (miscentering) {
      m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW;
      m_Sigma_mean_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_mean_NFW;
      m_Sigma_mis_ptr = &cbl::cosmology::HaloProfile::m_Sigma_mis;
      m_Sigma_ptr = &cbl::cosmology::HaloProfile::m_Sigma_including_miscentering;
      m_DeltaSigma_cen_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_cen;
      m_DeltaSigma_mis_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_mis;
      m_DeltaSigma_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_including_miscentering;
    }
    else {
      m_Sigma_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW;
      m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW;
      m_Sigma_mean_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_mean_NFW;
      m_DeltaSigma_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_cen;
    }
    
  }
  
  else if (profile_author == "NFW_trunc") {
    
    m_rho_ptr = &cbl::cosmology::HaloProfile::m_rho_NFW_trunc;
    m_rho_s_ptr = &cbl::cosmology::HaloProfile::m_rho_s_NFW_trunc;

    if (miscentering) {
      m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW_trunc;
      m_Sigma_mean_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_mean_NFW_trunc;
      m_Sigma_mis_ptr = &cbl::cosmology::HaloProfile::m_Sigma_mis;
      m_Sigma_ptr = &cbl::cosmology::HaloProfile::m_Sigma_including_miscentering;
      m_DeltaSigma_cen_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_cen;
      m_DeltaSigma_mis_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_mis;
      m_DeltaSigma_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_including_miscentering;
    }
    else {
      m_Sigma_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW_trunc;
      m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW_trunc;
      m_Sigma_mean_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_mean_NFW_trunc;
      m_DeltaSigma_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_cen;
    }
    
  }
  
  else
    ErrorCBL("density profile author not allowed!", "m_set_profile", "HaloProfile.cpp");
}


// =====================================================================================


void cbl::cosmology::HaloProfile::set_cosmology (const cbl::cosmology::Cosmology cosmology)
{
  m_cosmology = make_shared<cosmology::Cosmology>(cosmology::Cosmology(move(cosmology)));
}


// ============================================================================================


double cbl::cosmology::HaloProfile::rho_s ()
{
  return (this->*m_rho_s_ptr)();
}

// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::density_profile_3D (const std::vector<double> rad)
{
  return (this->*m_rho_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::Sigma (const std::vector<double> rad)
{  
  return (this->*m_Sigma_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::Sigma_mis (const std::vector<double> rad)
{  
  if (m_miscentering == false)
    ErrorCBL("the miscentering is not set!", "Sigma_mis", "HaloProfile.cpp");

  if (m_sigma_off > m_sigma_off_threshold)
    return (this->*m_Sigma_mis_ptr)(rad);
  else
    return (this->*m_Sigma_cen_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::Sigma_cen (const std::vector<double> rad)
{  
  if (m_miscentering == false)
    ErrorCBL("the miscentering is not set!", "Sigma_cen", "HaloProfile.cpp");

  return (this->*m_Sigma_cen_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::Sigma_2h (const std::vector<double> rad, const double bias, const std::string method_Pk, const std::string interp_type, const bool NL)
{
  if (cbl::Max(rad) > 70. || cbl::Min(rad) < 1.e-2)
    ErrorCBL("model tested only in the range [1.e-2,70] Mpc.", "Sigma_2h", "HaloProfile.cpp");
  if (method_Pk != "EisensteinHu")
    ErrorCBL("model working only with EisensteinHu!", "Sigma_2h", "HaloProfile.cpp");
  
  double Rho_m = m_cosmology->rho_m(m_redshift);
  double Dl = m_cosmology->D_A(m_redshift);
  
  double theta = 0; // in radians

  // Define the minimum and maximum l
  const double l_min = 1.e-1;
  const double l_max = 1.e6;

  // Interpolate P(k)
  const std::vector<double> kl_forInterp = cbl::logarithmic_bin_vector((size_t)(500), l_min/((1+m_redshift)*Dl), l_max/((1+m_redshift)*Dl));
  std::vector<double> Pk = m_cosmology->Pk_matter(kl_forInterp, method_Pk, NL, m_redshift, false, "test", 1, 1.e-4, 100.);
  cbl::glob::FuncGrid Pk_interp (kl_forInterp, Pk, interp_type);

  // Integrand function of the integral defining the 2-halo density profile
  std::function<double(double)> Func = [&] (double l)    
				       {
					 l = pow(10,l);
					 
					 double X  = theta*l;
					 double kl = l/((1+m_redshift)*Dl);
					 double Pk_ = Pk_interp(kl);
					 double J0 = gsl_sf_bessel_J0(X);
					 return l*J0*Pk_ * l;
				       };
  
  std::vector<double> S_2h(rad.size());  
  for (size_t i=0; i<rad.size(); i++) {
    theta = rad[i]/Dl;
    const double integral = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(l_min), log10(l_max)) * log(10);
    S_2h[i] = 1.e-12 * bias*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3)) * integral;
  }

  return S_2h;
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::Sigma_2h (const std::vector<double> rad, const std::string bias_author, const std::string method_Pk, const std::string interp_type, const bool NL)
{  
  double delta_bkg = m_Delta_func(m_Delta,*m_cosmology,m_redshift)/m_cosmology->OmegaM(m_redshift);
  double bias = m_cosmology->bias_halo(m_mass, m_redshift, bias_author, method_Pk, false, "test", interp_type, delta_bkg, -1., -1, 0.001, 100);
  
  return Sigma_2h (rad, bias, method_Pk, interp_type, NL);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma (const std::vector<double> rad)
{  
  return (this->*m_DeltaSigma_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma_mis (const std::vector<double> rad)
{  
  if (m_miscentering == false)
    ErrorCBL("the miscentering is not set!", "DeltaSigma_mis", "HaloProfile.cpp");

  return (this->*m_DeltaSigma_mis_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma_cen (const std::vector<double> rad)
{  
  if (m_miscentering == false)
    ErrorCBL("the miscentering is not set!", "DeltaSigma_cen", "HaloProfile.cpp");

  return (this->*m_DeltaSigma_cen_ptr)(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma_2h (const std::vector<double> rad, const double bias, const std::string method_Pk, const std::string interp_type, const bool NL)
{
  if (cbl::Max(rad) > 70. || cbl::Min(rad) < 1.e-2)
    ErrorCBL("model tested only in the range [1.e-2,70] Mpc.", "DeltaSigma_2h", "HaloProfile.cpp");
  if (method_Pk != "EisensteinHu")
    ErrorCBL("model working only with EisensteinHu!", "DeltaSigma_2h", "HaloProfile.cpp");
  
  double Rho_m = m_cosmology->rho_m(m_redshift);
  double Dl = m_cosmology->D_A(m_redshift);
  
  double theta = 0; // in radians

  // Define the minimum and maximum l
  const double l_min = 1.e-1;
  const double l_max = 1.e6;

  // Interpolate P(k)  
  const std::vector<double> kl_forInterp = cbl::logarithmic_bin_vector((size_t)(500), l_min/((1+m_redshift)*Dl), l_max/((1+m_redshift)*Dl));
  std::vector<double> Pk = m_cosmology->Pk_matter(kl_forInterp, method_Pk, NL, m_redshift, false, "test", 1, 1.e-4, 100.);
  cbl::glob::FuncGrid Pk_interp (kl_forInterp, Pk, interp_type);
  
  // Integrand function of the integral defining the 2-halo density profile
  std::function<double(double)> Func = [&] (double l)    
				       {
					 l = pow(10,l);
					 
					 double X  = theta*l;
					 double kl = l/((1+m_redshift)*Dl);         
					 double Pk_ = Pk_interp(kl);
					 double J2 = gsl_sf_bessel_Jn(2,X);
					 return l*J2*Pk_ * l;
				       };
  
  std::vector<double> DS_2h(rad.size());
  for (size_t i=0; i<rad.size(); i++) {
    theta = rad[i]/Dl;
    const double integral = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(l_min), log10(l_max)) * log(10);    
    DS_2h[i] = 1.e-12 * bias*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3)) * integral;
  }

  return DS_2h;
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma_2h (const std::vector<double> rad, const std::string bias_author, const std::string method_Pk, const std::string interp_type, const bool NL)
{  
  double delta_bkg = m_Delta_func(m_Delta,*m_cosmology,m_redshift)/m_cosmology->OmegaM(m_redshift);
  double bias = m_cosmology->bias_halo(m_mass, m_redshift, bias_author, method_Pk, false, "test", interp_type, delta_bkg, -1., -1, 0.001, 100);
  
  return DeltaSigma_2h (rad, bias, method_Pk, interp_type, NL);
}


// ============================================================================================


double cbl::cosmology::HaloProfile::density_profile_FourierSpace (const double kk)
{  
  const double conc = concentration();
  const double rho_s = m_cosmology->rho_crit(m_redshift)*m_cosmology->Delta_c(m_redshift)/3.*pow(conc, 3)/(log(1.+conc)-conc/(1.+conc));
  const double r_s = m_cosmology->r_vir(m_mass, m_redshift)/conc;  
  const double mu = kk*r_s;  
  
  return 4.*par::pi*rho_s*pow(r_s, 3)/m_mass*(cos(mu)*(gsl_sf_Ci(mu+mu*conc)-gsl_sf_Ci(mu))+sin(mu)*(gsl_sf_Si(mu+mu*conc)-gsl_sf_Si(mu))-sin(mu*conc)/(mu+mu*conc));
}


// ============================================================================


double cbl::cosmology::HaloProfile::concentration2 (const double Vmax, const double Rmax) const
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


double cbl::cosmology::HaloProfile::Mass_Delta (const double Mass, const double Delta_in, const double Delta_out, const double conc, const bool is_input_conc, const double rRmin_guess, const double rRmax_guess) const
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
