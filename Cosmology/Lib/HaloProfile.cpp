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


cbl::cosmology::HaloProfile::HaloProfile (std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const bool two_halo, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact, const bool miscentering, const bool single_profile, const double sigma_off, const double f_off, const double b_e, const double s_e, const double AB_fact, const std::string bias, const std::string method_Pk, const double alpha_0, const double alpha_nu, const double beta, const double gamma_0)
  : m_cosmology(std::move(cosmology))
{
  m_set_profile(two_halo, redshift, conc, Mass, Delta, profile_author, halo_def, trunc_fact, miscentering, single_profile, sigma_off, f_off, b_e, s_e, AB_fact, bias, method_Pk, alpha_0, alpha_nu, beta, gamma_0);
}


//=====================================================================================


cbl::cosmology::HaloProfile::HaloProfile (std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const bool two_halo, const double redshift, const std::string cM_author, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact, const bool miscentering, const bool single_profile, const double sigma_off, const double f_off, const double b_e, const double s_e, const double AB_fact, const std::string bias, const std::string method_Pk, const double alpha_0, const double alpha_nu, const double beta, const double gamma_0)
: m_cosmology(std::move(cosmology))
{
  m_set_profile(two_halo, redshift, 0., Mass, Delta, profile_author, halo_def, trunc_fact, miscentering, single_profile, sigma_off, f_off, b_e, s_e, AB_fact, bias, method_Pk, alpha_0, alpha_nu, beta, gamma_0);
  m_set_cM_relation(cM_author);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_nu_vir ()
{
  double mass_vir = 0.;
  
  if (m_halo_def == "vir")
    mass_vir = m_mass;

  else if (m_halo_def == "critical" || m_halo_def == "mean") {
  
    // Define the mass integrand
    std::function<double(double)> Integrand_mass = [&] (double r) {
      return 4.*par::pi*r*r*m_rho_NFW({r})[0];
    };

    // Interpolate the enclosed mass
    std::vector<double> r_forInterp = cbl::logarithmic_bin_vector(50, 0.001, 100.);
    std::vector<double> mass_forInterp (r_forInterp.size(), 0.);
    
    for (size_t i=0; i<r_forInterp.size(); i++)
      mass_forInterp[i] = wrapper::gsl::GSL_integrate_qag(Integrand_mass, 0., r_forInterp[i]);
    
    cbl::glob::FuncGrid mass_interp (r_forInterp, mass_forInterp, "Spline");
    
    // Get the virial overdensity
    const double Delta_vir = m_cosmology->rho_crit(m_redshift)*m_cosmology->Delta_c(m_redshift, "BryanNorman");

    // Define the function for which we search the zero
    std::function<double(double)> fc = [&] (double r) {
      return mass_interp(r)/4.0/cbl::par::pi*3.0/pow(r,3.) - Delta_vir;
    };
    
    // Get the virial radius, then the virial mass
    const double R_vir = cbl::wrapper::gsl::GSL_root_brent(fc, Delta_vir, cbl::Min(r_forInterp), cbl::Max(r_forInterp));
    mass_vir = 4./3.*cbl::par::pi*Delta_vir*pow(R_vir, 3.);
    
  }
  
  else
    ErrorCBL("wrong halo_def declaration!", "m_nu_vir", "HaloProfile.cpp");

  // Compute the virial peak height
  const double lagrangian_R = pow( 3.*mass_vir/(4.*par::pi*m_cosmology->rho_m(0.)), 1./3. );

  
  cosmology::Sigma SG(m_cosmology);
  
  const double nu = 1.68647/sqrt(SG.sigma2R({lagrangian_R}, "EisensteinHu", m_redshift)[0]);

  return nu;
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_concentration_Duffy ()
{  
  const double fact = (m_cosmology->unit()) ? 1 : m_cosmology->little_h();
  const double Mpivot = 2.e12/fact; // in Msun/h
  
  return m_A*pow(m_mass/Mpivot, m_B)*pow(1.+m_redshift, m_C);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_rho_s_NFW ()
{
  // Compute rho_s
  const double conc = (this->*m_return_concentration)();
  const double rho_s = pow(conc, 3)*m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift)/(3.*(log(1.+conc)-conc/(1.+conc)));
  
  return rho_s;
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_rho_s_NFW_trunc ()
{
  // compute rho_s
  const double conc = (this->*m_return_concentration)();  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double r_t = m_trunc_fact*r_s*conc;
  const double tau = r_t/r_s;
  const double fact = tau*tau/(2.*pow(tau*tau+1.,3)*(1.+conc)*(tau*tau+conc*conc))*(conc*(tau*tau+1.)*(conc*(conc+1.) - tau*tau*(conc-1.)*(2.+3.*conc) - 2.*pow(tau,4)) + tau*(conc+1.)*(tau*tau+conc*conc)*(2.*(3.*tau*tau-1.)*std::atan(conc/tau) + tau*(tau*tau-3.)*log(tau*tau*(1.+conc)*(1.+conc)/(tau*tau+conc*conc)))); // Eq. 10 from Oguri & Hamana 2011
  
  const double rho_s = pow(conc, 3)*m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift)/(3.*fact);
  
  return rho_s;
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_rho_s_Einasto (const double nu)
{
  const double alpha = m_alpha_0 + m_alpha_nu*nu*nu;

  const double conc = (this->*m_return_concentration)();  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  
  // Compute rho_s
  std::function<double(double)> Integrand_in = [&] (double r) {
    const double rho_inner = exp((-2./alpha)*(pow(r/r_s, alpha)-1.)); // Einasto profile with rho_s=1
    return 4.*par::pi*r*r*rho_inner;
  };
  
  return m_mass/wrapper::gsl::GSL_integrate_qag(Integrand_in, 0., r_s*conc);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::m_rho_s_DK14 (const double nu)
{
  const double alpha = m_alpha_0 + m_alpha_nu*nu*nu;
  const double beta = m_beta;
  const double gamma = m_gamma_0*nu;

  const double conc = (this->*m_return_concentration)();  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double r_t = m_trunc_fact*r_s*conc;
  
  // Compute rho_s
  std::function<double(double)> Integrand_in = [&] (double r) {
    const double f_trans = pow(1. + pow(r/r_t, beta), -gamma/beta);
    const double rho_inner = exp((-2./alpha)*(pow(r/r_s, alpha)-1.)); // Einasto profile with rho_s=1
    const double densityInner = rho_inner*f_trans;
    return 4.*par::pi*r*r*densityInner;
  };

  std::function<double(double)> Integrand_out = [&] (double r) {
    const double rho_out = m_cosmology->rho_m(m_redshift)*(m_b_e*pow(r/(5*r_s*conc), -m_s_e) + 1.);
    return 4.*par::pi*r*r*rho_out;
  };

  const double enclosedMass_in = wrapper::gsl::GSL_integrate_qag(Integrand_in, 0., r_s*conc);
  const double enclosedMass_out = wrapper::gsl::GSL_integrate_qag(Integrand_out, 0., r_s*conc)*(1+m_AB_fact);
  const double rho_s = (m_mass - enclosedMass_out)/enclosedMass_in;
  
  return rho_s;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_rho_NFW (const std::vector<double> rad, const double nuisance)
{
  (void) nuisance;
  
  // compute rho_s and r_s
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double rho_s = m_rho_s_NFW();

  // compute rho
  std::vector<double> rho(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    rho[i] = rho_s/((rad[i]/r_s)*pow(1.+rad[i]/r_s, 2));
  
  return rho;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_rho_NFW_trunc (const std::vector<double> rad, const double nuisance)
{
  (void) nuisance;
  
  // compute rho_s and r_s
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double r_t = m_trunc_fact*r_s*conc;
  const double rho_s = m_rho_s_NFW_trunc();

  // compute rho
  std::vector<double> rho(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    rho[i] = rho_s/((rad[i]/r_s)*pow(1.+rad[i]/r_s, 2))*pow(r_t*r_t/(rad[i]*rad[i]+r_t*r_t),2);
  
  return rho;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_rho_Einasto (const std::vector<double> rad, const double nu_in)
{
  // profile parameters
  const double nu = (nu_in == -1) ? m_nu_vir() : nu_in;
  const double alpha = m_alpha_0 + m_alpha_nu*nu*nu;
  
  // compute rho_s and r_s
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double rho_s = m_rho_s_Einasto(nu);

  // compute rho
  std::vector<double> rho(rad.size());
  
  for (size_t i=0; i<rad.size(); i++)
    rho[i] = rho_s*exp( (-2./alpha)*(pow(rad[i]/r_s, alpha) - 1.) );
  
  return rho;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_rho_DK14 (const std::vector<double> rad, const double nu_in)
{
  // profile parameters
  const double nu = (nu_in == -1) ? m_nu_vir() : nu_in;
  const double alpha = m_alpha_0 + m_alpha_nu*nu*nu;
  const double beta = m_beta;
  const double gamma = m_gamma_0*nu;
  
  // compute rho_s and r_s
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double rho_s = m_rho_s_DK14(nu);
  const double r_t = m_trunc_fact*r_s*conc;

  // compute rho
  std::vector<double> rho(rad.size());
  
  for (size_t i=0; i<rad.size(); i++) {
    const double f_trans = pow(1. + pow(rad[i]/r_t, beta), -gamma/beta);
    const double rho_in  = rho_s*exp((-2./alpha)*(pow(rad[i]/r_s, alpha)-1.)); // Einasto profile
    const double rho_out = m_cosmology->rho_m(m_redshift)*(m_b_e*pow(rad[i]/(5*r_s*conc), -m_s_e) + 1.)*(1+m_AB_fact);
    rho[i] = rho_in*f_trans + rho_out;
    
  }
  
  return rho;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_NFW (const std::vector<double> rad)
{
  // Implementation of second part of eq. 4 from Golse et al. 2002 (https://ui.adsabs.harvard.edu/abs/2002A%26A...390..821G/abstract)
  std::function<double(double)> F = [] (double x) { 
    if (x < 1.)
      return (1 - std::acosh(1./x)/std::sqrt(1.-x*x))/(x*x-1.);
    else if (x == 1.)
      return 1./3.;
    else
      return (1 - std::acos(1./x)/std::sqrt(x*x-1.))/(x*x-1.);
  };

  // compute the density
  const double conc = (this->*m_return_concentration)();
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double rho_s = m_rho_s_NFW();

  std::vector<double> Sigma(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    Sigma[i] = 2.*rho_s*r_s*F(rad[i]/r_s)*1.e-12;

  // add the 2-halo term
  if (m_two_halo) {
    
    for (size_t i=0; i<rad.size(); i++)
      Sigma[i] += m_Sigma_2h_unbiased_interp->operator()(rad[i])*m_bias*(1+m_AB_fact);

  }
  
  return Sigma;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_NFW_trunc (const std::vector<double> rad)
{
  // Implementation of eq. A.5 from Baltz et al. 2009 (https://ui.adsabs.harvard.edu/abs/2009JCAP...01..015B/abstract)  
  const double conc = (this->*m_return_concentration)();
  
  const double overdensity = m_cosmology->rho_crit(m_redshift)*m_Delta_func(m_Delta, m_cosmology, m_redshift);
  const double r_s = pow( 3.*m_mass/(4.*cbl::par::pi*overdensity), 1./3. )/conc;
  const double r_t = m_trunc_fact*r_s*conc;
  const double tau = r_t/r_s;
  const double rho_s = m_rho_s_NFW_trunc();
  
  std::function<double(double)> F = [] (double x) {
    if (x < 1.)
      return std::acosh(1./x)/std::sqrt(1.-x*x);
    else if (x == 1.)
      return 1.;
    else
      return std::acos(1./x)/std::sqrt(x*x-1.);
  };

  // Implementation of eq. A.28 from Baltz et al. 2009
  std::function<double(double)> G = [&F] (double x) {
    if (x < 1.)
      return (F(x)-1.)/(1.-x*x);
    else if (x == 1.)
      return 1./3.;
    else
      return (1.-F(x))/(x*x-1.);
  };

  // Implementation of second part of eq. A.6 from Baltz et al. 2009
  std::function<double(double)> L = [&tau] (double x) {
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
    
    Sigma[i] = prefact*(a+b+c+d+e)*1.e-12;
  }

  // add the 2-halo term
  if (m_two_halo) {
    
    for (size_t i=0; i<rad.size(); i++)
      Sigma[i] += m_Sigma_2h_unbiased_interp->operator()(rad[i])*m_bias*(1+m_AB_fact);

  }
  
  return Sigma;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_generic (const std::vector<double> rad)
{
  double nu = -1;
  if (m_profile_author == "DK14" || m_profile_author == "Einasto")
    nu = m_nu_vir();
  
  // Interpolate DeltaRho
  std::vector<double> r_forInterp = cbl::logarithmic_bin_vector(300, 1.e-8, 150.);
  std::vector<double> rho_forInterp = (this->*m_rho_ptr)(r_forInterp, nu);
  
  cbl::glob::FuncGrid rho_interp (r_forInterp, rho_forInterp, "Spline");
  
  // Compute Sigma
  double dummy_R;
  std::function<double(double)> Integrand = [&] (double chi) {
    return rho_interp(sqrt(dummy_R*dummy_R + chi*chi));
  }; 
  
  std::vector<double> Sigma(rad.size());
  for (size_t i=0; i<rad.size(); i++) {
    
    dummy_R = rad[i];
    Sigma[i] = 2*wrapper::gsl::GSL_integrate_qag(Integrand, cbl::Min(r_forInterp), cbl::Max(r_forInterp))*1.e-12;

    if (m_two_halo && m_profile_author != "DK14")
      Sigma[i] += m_Sigma_2h_unbiased_interp->operator()(rad[i])*m_bias*(1+m_AB_fact);
    
  }
  
  return Sigma;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_mean_cen (const std::vector<double> rad)
{
  // Interpolate Sigma
  std::vector<double> r_forInterp = cbl::logarithmic_bin_vector(300, 1.e-8, 150.);
  std::vector<double> Sigma_forInterp = (this->*m_Sigma_cen_ptr)(r_forInterp);
  cbl::glob::FuncGrid Sigma_interp (r_forInterp, Sigma_forInterp, "Spline");
  
  // Compute Sigma_mean
  std::function<double(double)> mean_sigma_integrand = [&Sigma_interp] (double r) {
    return Sigma_interp(r)*r;
  };
  
  std::vector<double> Sigma_mean(rad.size());
  for (size_t i=0; i<rad.size(); i++)
    Sigma_mean[i] = cbl::wrapper::gsl::GSL_integrate_cquad(mean_sigma_integrand, 1.e-8, rad[i])*2./rad[i]/rad[i];
  
  return Sigma_mean;
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_mis (const std::vector<double> rad)
{
  // Interpolate Sigma_cen
  std::vector<double> r_forInterp = cbl::logarithmic_bin_vector(300, 1.e-8, 150.);
  std::vector<double> Sigma_forInterp = (this->*m_Sigma_cen_ptr)(r_forInterp);
  cbl::glob::FuncGrid Sigma_interp (r_forInterp, Sigma_forInterp, "Spline");
  
  // Define the integrands
  double dummy_R_s, dummy_R;
  cbl::glob::FuncGrid Sigma_R_Rs_interp, P_Rs_interp;
  
  std::function<double(double)> integrand = [&] (double theta) {
					      double new_rad2 = std::abs(dummy_R*dummy_R + dummy_R_s*dummy_R_s + 2*dummy_R*dummy_R_s*std::cos(theta));
					      return Sigma_interp(std::sqrt(new_rad2));
  };

  std::function<double(double)> Sigma_mis_integrand = [&] (double log_R_s) {
							double R_s = exp(log_R_s);
						      	return R_s * P_Rs_interp(R_s) * Sigma_R_Rs_interp(R_s);
						      };
  
  // Compute Sigma_mis
  std::vector<double> Sigma_mis(rad.size());

  if (m_single_profile) { // m_sigma_off is the miscentring scale 

    dummy_R_s = m_sigma_off;
    
    for (size_t i=0; i<rad.size(); i++) {
      dummy_R = rad[i];
      const double integral = cbl::wrapper::gsl::GSL_integrate_cquad(integrand, 0., cbl::par::pi);
      Sigma_mis[i] = integral / (cbl::par::pi);
    }
    
  }
  
  else { // m_sigma_off is the rms of the miscentring scale distribution

    // Define the integration range
    double Rs_min = 1.e-8;
    double Rs_max = m_sigma_off * sqrt(-2. * log(1. - 0.99)); // 99th quantile for a Rayleigh distribution

    // Define the vectors for the interpolation of Sigma_cen(R|R_s) and P(R_s)
    std::vector<double> Rs_forInterp = cbl::logarithmic_bin_vector(100, Rs_min, Rs_max);
    std::vector<double> Sigma_R_Rs_samples (Rs_forInterp.size(), 0.);
    std::vector<double> P_Rs_samples (Rs_forInterp.size(), 0.);

    // Interpolate P(R_s), which here is a Rayleigh distribution
    for (size_t j=0; j<Rs_forInterp.size(); j++)
      P_Rs_samples[j] = Rs_forInterp[j] / m_sigma_off / m_sigma_off * exp(-0.5 * pow(Rs_forInterp[j] / m_sigma_off, 2.));
    P_Rs_interp = cbl::glob::FuncGrid (Rs_forInterp, P_Rs_samples, "Spline");

    // Loop over the radial points
    for (size_t i=0; i<rad.size(); i++) {

      int departure_index = -1; // find the departure of Sigma_cen(R|R_s) from a constant
      dummy_R = rad[i];

      // Interpolate Sigma_cen(R|R_s)      
      for (size_t j=0; j<Rs_forInterp.size(); j++) {
	
	dummy_R_s = Rs_forInterp[j];
	Sigma_R_Rs_samples[j] = cbl::wrapper::gsl::GSL_integrate_cquad(integrand, 0., cbl::par::pi) / cbl::par::pi;

	if ( (rad[i] > Rs_forInterp[j]) && (std::abs( Sigma_R_Rs_samples[j] / Sigma_interp(rad[i]) - 1. ) > 0.01) && (departure_index == -1) )
	  departure_index = j;
	
      }
      Sigma_R_Rs_interp = cbl::glob::FuncGrid (Rs_forInterp, Sigma_R_Rs_samples, "Spline");

      // Compute the final quantity
      if (departure_index == -1)
	Sigma_mis[i] = cbl::wrapper::gsl::GSL_integrate_cquad(Sigma_mis_integrand, log(Rs_min), log(Rs_max));
      else {

	// the first integral part equals the constant Sigma_cen(R) times the cumulative Rayleigh distribution
	Sigma_mis[i] = Sigma_R_Rs_samples[0] * ( 1. - exp( - Rs_forInterp[departure_index] * Rs_forInterp[departure_index] / m_sigma_off / m_sigma_off / 2. ) );
	Sigma_mis[i] += cbl::wrapper::gsl::GSL_integrate_cquad(Sigma_mis_integrand, log(Rs_forInterp[departure_index]), log(Rs_max));

      }
      
    }

  }
  
  return Sigma_mis;
  
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_including_miscentering (const std::vector<double> rad)
{
  if (m_sigma_off > m_sigma_off_threshold && m_f_off > m_f_off_threshold) {
    const std::vector<double> Sigma_cen = (this->*m_Sigma_cen_ptr)(rad);
    const std::vector<double> Sigma_mis = m_Sigma_mis(rad);
  
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
  const std::vector<double> Sigma_mean = m_Sigma_mean_cen(rad);  

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
    const std::vector<double> rad_forInterp = cbl::logarithmic_bin_vector((size_t)(50), m_min_Smis, m_max_Smis);
    const std::vector<double> Sigma = m_Sigma_mis(rad_forInterp);
    
    cbl::glob::FuncGrid Sigma_interp (rad_forInterp, Sigma, "Spline");

    // Compute the mean surface density
    std::function<double(double)> mean_sigma_integrand = [&Sigma_interp] (double rad) {
      return Sigma_interp(rad)*rad;
    };
  
    std::vector<double> Sigma_mean(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      Sigma_mean[i] = cbl::wrapper::gsl::GSL_integrate_cquad(mean_sigma_integrand, 0., rad[i])*2./rad[i]/rad[i];

    // Compute DeltaSigma
    std::vector<double> DeltaSigma(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      DeltaSigma[i] = Sigma_mean[i] - Sigma_interp(rad[i]);
  
    return DeltaSigma;

  }
  else
    return m_DeltaSigma_cen(rad);
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_DeltaSigma_including_miscentering (const std::vector<double> rad)
{
  if (m_sigma_off > m_sigma_off_threshold && m_f_off > m_f_off_threshold) {
    const std::vector<double> DeltaSigma_cen = m_DeltaSigma_cen(rad);
    const std::vector<double> DeltaSigma_mis = m_DeltaSigma_mis(rad);
  
    std::vector<double> DeltaSigma_tot(rad.size());
    for (size_t i=0; i<rad.size(); i++)
      DeltaSigma_tot[i] = (1.-m_f_off)*DeltaSigma_cen[i] + m_f_off*DeltaSigma_mis[i];

    return DeltaSigma_tot;
  }
  else
    return m_DeltaSigma_cen(rad);
}


// =====================================================================================


std::vector<double> cbl::cosmology::HaloProfile::m_Sigma_2h_unbiased (const std::vector<double> rad)
{
  double Rho_m = m_cosmology->rho_m(m_redshift);
  double Dl = m_cosmology->D_A(m_redshift);
  
  double theta = 0; // in radians

  // Define the minimum and maximum l
  const double kl_min = 1.e-4;
  const double kl_max = 1.e2;

  // Interpolate P(k)
  const std::vector<double> kl_forInterp = cbl::logarithmic_bin_vector((size_t)(500), kl_min, kl_max);
  cosmology::PkXi PX(m_cosmology);
  std::vector<double> Pk = PX.Pk_matter(kl_forInterp, m_method_Pk, false, m_redshift, false, "test", -1, 1.e-4, 100.);
  cbl::glob::FuncGrid Pk_interp (kl_forInterp, Pk, "Linear");

  // Integrand function of the integral defining the 2-halo density profile
  std::function<double(double)> Func = [&] (double l)    
  {
    l = pow(10,l);
					 
    double X  = theta*l;
    double kl = l/((1+m_redshift)*Dl);
    double Pk_ = Pk_interp(kl);
    double J0 = gsl_sf_bessel_J0(X);
    return l*J0*Pk_*l;
  };

  // Compute quantities for linear extrapolations, if needed
  double Rmin = 1.e-2, Rmax = 30.;
  double x_max1=0, x_max2=0, y_max1=0, y_max2=0;
  
  if (cbl::Max(rad) > Rmax) {
    
    x_max1 = Rmax;
    x_max2 = Rmax-2.;

    theta = x_max1/Dl;
    const double integral1 = cbl::wrapper::gsl::GSL_integrate_cquad(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_max1 = abs( 1.e-12*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral1 );

    theta = x_max2/Dl;
    const double integral2 = cbl::wrapper::gsl::GSL_integrate_cquad(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_max2 = abs( 1.e-12*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral2 );
    
  }

  double x_min1=0, x_min2=0, y_min1=0, y_min2=0;
  
  if (cbl::Min(rad) < Rmin) {

    x_min1 = Rmin;
    x_min2 = 2.*Rmin;

    theta = x_min1/Dl;
    const double integral1 = cbl::wrapper::gsl::GSL_integrate_cquad(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_min1 = 1.e-12*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral1;

    theta = x_min2/Dl;
    const double integral2 = cbl::wrapper::gsl::GSL_integrate_cquad(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_min2 = 1.e-12*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral2;
    
  }

  // Compute the 2-halo term at all radial values
  std::vector<double> S_2h(rad.size());  
  for (size_t i=0; i<rad.size(); i++) {

    if (rad[i] > Rmax) { // Linear extrapolation at high radii

      S_2h[i] = pow( 10, std::max( log10(y_max2) + (log10(rad[i])-log10(x_max2)) / (log10(x_max1)-log10(x_max2)) * (log10(y_max1)-log10(y_max2)), log10(1.e-10)) );
      
    }

    else if (rad[i] < Rmin) { // Linear extrapolation at low radii

      S_2h[i] = pow( 10, std::max( log10(y_min2) + (log10(rad[i])-log10(x_min2)) / (log10(x_min1)-log10(x_min2)) * (log10(y_min1)-log10(y_min2)), log10(1.e-10)) );

    }

    else {
    
      theta = rad[i]/Dl;
      const double integral = cbl::wrapper::gsl::GSL_integrate_cquad(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
      S_2h[i] = 1.e-12*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral;

    }
    
  }

  return S_2h;
}


// =====================================================================================


void cbl::cosmology::HaloProfile::m_set_cM_relation (const std::string cM_author)
{
  m_isSet_cM_relation = true;

  if (cM_author=="Duffy") {
    if (m_redshift>2) ErrorCBL("the concentration-mass relation by Duffy et al. has been tested only at z<2", "set_cM_relation", "HaloProfile.cpp");

    if ( (m_halo_def == "critical" && m_Delta != 200.) || (m_halo_def == "mean" && 200.*m_cosmology->Omega_matter(m_redshift) != m_Delta_func(m_Delta, m_cosmology, m_redshift)) )
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
    else if (m_profile_author=="Einasto" || m_profile_author=="DK14") {
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


void cbl::cosmology::HaloProfile::m_set_profile (const bool two_halo, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact, const bool miscentering, const bool single_profile, const double sigma_off, const double f_off, const double b_e, const double s_e, const double AB_fact, const std::string bias, const std::string method_Pk, const double alpha_0, const double alpha_nu, const double beta, const double gamma_0)
{
  if (miscentering && single_profile && f_off != 1.)
    ErrorCBL("If a single profile is considered in the miscentering model, then f_off must be set equal to 1.", "m_set_profile", "HaloProfile.cpp");

  if (profile_author == "DK14" && two_halo)
    ErrorCBL("If the DK14 profile is assumed, the 2-halo term cannot be included.", "m_set_profile", "HaloProfile.cpp");

  m_sigma_off_threshold = 1.e-4;
  m_f_off_threshold = 1.e-4;
  m_min_Smis = 1.e-5;
  m_max_Smis = 300;

  m_two_halo = two_halo;
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

  m_b_e = b_e;
  m_s_e = s_e;
  m_AB_fact = AB_fact;

  m_bias_author = bias;
  m_method_Pk = method_Pk;

  m_alpha_0 = alpha_0;
  m_alpha_nu = alpha_nu;

  m_beta = beta;
  m_gamma_0 = gamma_0;

  // Set the overdensity factor Delta
  if (halo_def == "critical")
    m_Delta_func = [] (const double delta, std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const double redshift) {
      (void)cosmo; (void)redshift;
      return delta;};
  else if (halo_def == "mean")
    m_Delta_func = [] (const double delta, std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const double redshift) {
      return delta*cosmo->Omega_matter(redshift);};
  else if (halo_def == "vir")
    m_Delta_func = [] (const double delta, std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const double redshift) {
      (void)delta;
      return cosmo->Delta_c(redshift, "BryanNorman");};
  else
    ErrorCBL("wrong halo_def declaration!", "m_set_profile", "HaloProfile.cpp");

  // Set the function returning the concentration
  m_return_concentration = &cbl::cosmology::HaloProfile::m_return_set_concentration;

  // Set the profile  
  if (profile_author == "NFW") {
    
    m_rho_ptr = &cbl::cosmology::HaloProfile::m_rho_NFW;
    m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW;

  }
  
  else if (profile_author == "NFW_trunc") {
    
    m_rho_ptr = &cbl::cosmology::HaloProfile::m_rho_NFW_trunc;
    m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_NFW_trunc;
    
  }

  else if (profile_author == "Einasto") {
    
    m_rho_ptr = &cbl::cosmology::HaloProfile::m_rho_Einasto;
    m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_generic;
    
  }

  else if (profile_author == "DK14") {
    
    m_rho_ptr = &cbl::cosmology::HaloProfile::m_rho_DK14;
    m_Sigma_cen_ptr = &cbl::cosmology::HaloProfile::m_Sigma_generic;
    
  }
  
  else
    ErrorCBL("density profile author not recognised.", "m_set_profile", "HaloProfile.cpp");

  if (miscentering) {
    m_Sigma_ptr = &cbl::cosmology::HaloProfile::m_Sigma_including_miscentering;
    m_DeltaSigma_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_including_miscentering;
  }
  else {
    m_Sigma_ptr = m_Sigma_cen_ptr;
    m_DeltaSigma_ptr = &cbl::cosmology::HaloProfile::m_DeltaSigma_cen;
  }

  // Interpolate Sigma_2h
  if (m_two_halo)
    this->update_2halo();
  
}


// =====================================================================================


void cbl::cosmology::HaloProfile::update_2halo (const bool update_bias)
{
  m_two_halo = true;
  
  if (update_bias) {
    const double delta_bkg = m_Delta_func(m_Delta, m_cosmology, m_redshift)/m_cosmology->Omega_matter(m_redshift);

    Bias bias(m_cosmology);
    m_bias = bias.bias_halo(m_mass, m_redshift, m_bias_author, m_method_Pk, false, "test", "Linear", delta_bkg, -1., -1, 0.001, 100);
  }

  std::vector<double> r_forInterp = cbl::logarithmic_bin_vector(100, 1.e-8, 150.);
  std::vector<double> Sigma_forInterp = this->m_Sigma_2h_unbiased(r_forInterp);
  cbl::glob::FuncGrid Sigma_interp (r_forInterp, Sigma_forInterp, "Spline");

  m_Sigma_2h_unbiased_interp = make_shared<glob::FuncGrid>(Sigma_interp);
}


// =====================================================================================


void cbl::cosmology::HaloProfile::set_cosmology (const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const bool update_2halo)
{
  m_cosmology = move(cosmology);

  if (update_2halo && m_two_halo)
    this->update_2halo();
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::rho (const std::vector<double> rad)
{
  std::vector<double> rho = (this->*m_rho_ptr)(rad, -1.);

  if ( m_two_halo && (m_profile_author != "DK14") ) {

    double Rho_m = m_cosmology->rho_m(m_redshift);
    
    double delta_bkg = m_Delta_func(m_Delta, m_cosmology, m_redshift)/m_cosmology->Omega_matter(m_redshift);

    Bias bias_obj(m_cosmology);
    const double bias = bias_obj.bias_halo(m_mass, m_redshift, m_bias_author, m_method_Pk, false, "test", "Linear", delta_bkg, -1., -1, 0.001, 100);
  
    // Define the minimum and maximum l
    const double kl_min = 1.e-2;
    const double kl_max = 1.e2;

    // Interpolate P(k)  
    const std::vector<double> kl_forInterp = cbl::logarithmic_bin_vector((size_t)(500), kl_min, kl_max);
    cosmology::PkXi PX(m_cosmology);
    std::vector<double> Pk = PX.Pk_matter(kl_forInterp, m_method_Pk, false, m_redshift, false, "test", -1, 1.e-4, 100.);
    cbl::glob::FuncGrid Pk_interp (kl_forInterp, Pk, "Linear");

    // Compute rho_2halo for linear extrapolations, if needed    
    double R = 0;
    double Rmax = 15.;

    std::function<double(double)> Func = [&] (double k)    
    {					 
      double Pk_ = Pk_interp(k);
      double J0 = gsl_sf_bessel_J0(k*R);
      return pow(k,2) * J0 * Pk_;
    };
  
    double x_max1=0, x_max2=0, y_max1=0, y_max2=0;
  
    if (cbl::Max(rad) > Rmax) {
    
      x_max1 = Rmax;
      x_max2 = Rmax - 2.;

      R = x_max1;
      const double integral1 = cbl::wrapper::gsl::GSL_integrate_qag(Func, kl_min, kl_max);
      y_max1 = Rho_m * (1 + bias / (2.*cbl::par::pi*cbl::par::pi) * integral1) * (1+m_AB_fact);

      R = x_max2;
      const double integral2 = cbl::wrapper::gsl::GSL_integrate_qag(Func, kl_min, kl_max);
      y_max2 = Rho_m * (1 + bias / (2.*cbl::par::pi*cbl::par::pi) * integral2) * (1+m_AB_fact);
    
    }
  
    // Derive the 2-halo term    
    for (size_t i=0; i<rad.size(); i++) {

      if (rad[i] > 1.e-1) {

	if (rad[i] > Rmax) { // Linear extrapolation at high radii

	  rho[i] += pow( 10, std::max( log10(y_max2) + (log10(rad[i])-log10(x_max2)) / (log10(x_max1)-log10(x_max2)) * (log10(y_max1)-log10(y_max2)), log10(1.e-10)) );

	} else {
	
	  R = rad[i];
	  const double integral = cbl::wrapper::gsl::GSL_integrate_qag(Func, kl_min, kl_max);
	  rho[i] += Rho_m * (1 + bias / (2.*cbl::par::pi*cbl::par::pi) * integral) * (1+m_AB_fact);

	}

      }

    }
    
  }
  
  return rho;
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
    return m_Sigma_mis(rad);
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


std::vector<double> cbl::cosmology::HaloProfile::Sigma_2h (const std::vector<double> rad)
{
  if (m_two_halo)
    ErrorCBL("This is a safety error. If two_halo is true, the 2-halo contribution is already included in Sigma() and DeltaSigma().", "Sigma_2h", "HaloProfile.cpp");

  const double delta_bkg = m_Delta_func(m_Delta, m_cosmology, m_redshift)/m_cosmology->Omega_matter(m_redshift);

  Bias bias(m_cosmology);
  const double _bias = bias.bias_halo(m_mass, m_redshift, m_bias_author, m_method_Pk, false, "test", "Linear", delta_bkg, -1., -1, 0.001, 100);
  
  std::vector<double> Sigma = m_Sigma_2h_unbiased(rad);
  for (size_t i=0; i<rad.size(); i++)
    Sigma[i] = Sigma[i]*_bias*(1+m_AB_fact);
  
  return Sigma;
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

  return m_DeltaSigma_mis(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma_cen (const std::vector<double> rad)
{  
  if (m_miscentering == false)
    ErrorCBL("the miscentering is not set!", "DeltaSigma_cen", "HaloProfile.cpp");

  return m_DeltaSigma_cen(rad);
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::DeltaSigma_2h (const std::vector<double> rad)
{
  if (m_two_halo)
    ErrorCBL("This is a safety error. If two_halo is true, the 2-halo contribution is already included in Sigma() and DeltaSigma().", "DeltaSigma_2h", "HaloProfile.cpp");
  
  const double Rho_m = m_cosmology->rho_m(m_redshift);
  const double Dl = m_cosmology->D_A(m_redshift);

  const double delta_bkg = m_Delta_func(m_Delta, m_cosmology, m_redshift)/m_cosmology->Omega_matter(m_redshift);

  Bias bias(m_cosmology);
  const double _bias = bias.bias_halo(m_mass, m_redshift, m_bias_author, m_method_Pk, false, "test", "Linear", delta_bkg, -1., -1, 0.001, 100);
  
  double theta = 0.; // in radians

  // Define the minimum and maximum l
  const double kl_min = 1.e-4;
  const double kl_max = 1.e2;

  // Interpolate P(k)  
  const std::vector<double> kl_forInterp = cbl::logarithmic_bin_vector((size_t)(500), kl_min, kl_max);
  cosmology::PkXi PX(m_cosmology);
  std::vector<double> Pk = PX.Pk_matter(kl_forInterp, m_method_Pk, false, m_redshift, false, "test", -1, 1.e-4, 100.);
  cbl::glob::FuncGrid Pk_interp (kl_forInterp, Pk, "Linear");
  
  // Integrand function of the integral defining the 2-halo density profile
  std::function<double(double)> Func = [&] (double l)    
  {
    l = pow(10,l);
					 
    double X  = theta*l;
    double kl = l/((1+m_redshift)*Dl);         
    double Pk_ = Pk_interp(kl);
    double J2 = gsl_sf_bessel_Jn(2,X);
    return l*J2*Pk_*l;
  };
  
  // Compute quantities for linear extrapolations, if needed
  double Rmin = 1.e-2, Rmax = 30.;
  double x_max1=0, x_max2=0, y_max1=0, y_max2=0;
  
  if (cbl::Max(rad)>Rmax) {
    
    x_max1 = Rmax;
    x_max2 = Rmax-2.;

    theta = x_max1/Dl;
    const double integral1 = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_max1 = 1.e-12*_bias*Rho_m/(2.*cbl::par::pi*pow(Dl, 2)*pow(1+m_redshift, 3))*integral1*(1+m_AB_fact);

    theta = x_max2/Dl;
    const double integral2 = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_max2 = 1.e-12*_bias*Rho_m/(2.*cbl::par::pi*pow(Dl, 2)*pow(1+m_redshift, 3))*integral2*(1+m_AB_fact);
    
  }

  double x_min1 = 0., x_min2 = 0., y_min1 = 0., y_min2 = 0.;
  
  if (cbl::Min(rad)<Rmin) {

    x_min1 = Rmin;
    x_min2 = 2.*Rmin;

    theta = x_min1/Dl;
    const double integral1 = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_min1 = 1.e-12*_bias*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral1*(1+m_AB_fact);

    theta = x_min2/Dl;
    const double integral2 = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
    y_min2 = 1.e-12*_bias*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral2*(1+m_AB_fact);
    
  }

  // Compute the 2-halo term at all radial values
  std::vector<double> DS_2h(rad.size());  
  for (size_t i=0; i<rad.size(); i++) {

    if (rad[i] > Rmax) { // Linear extrapolation at high radii

      DS_2h[i] = pow( 10, std::max( log10(y_max2) + (log10(rad[i])-log10(x_max2)) / (log10(x_max1)-log10(x_max2)) * (log10(y_max1)-log10(y_max2)), log10(1.e-10)) );

    }

    else if (rad[i] < Rmin) { // Linear extrapolation at low radii

      DS_2h[i] = pow( 10, std::max( log10(y_min2) + (log10(rad[i])-log10(x_min2)) / (log10(x_min1)-log10(x_min2)) * (log10(y_min1)-log10(y_min2)), log10(1.e-10)) );

    }

    else {
      
      theta = rad[i]/Dl;
      const double integral = cbl::wrapper::gsl::GSL_integrate_qag(Func, log10(kl_min*(1.0+m_redshift)*Dl), log10(kl_max*(1.0+m_redshift)*Dl))*log(10);
      DS_2h[i] = 1.e-12*_bias*Rho_m/(2.*cbl::par::pi*pow(Dl,2)*pow(1+m_redshift,3))*integral*(1+m_AB_fact);

    }
  }

  return DS_2h;
}


// ============================================================================================


std::vector<double> cbl::cosmology::HaloProfile::g (const std::vector<double> rad, const std::vector<double> z_s, const std::vector<double> N_zs)
{
  // Inverse critical density  
  double inv1_Sigma_crit = 0, inv2_Sigma_crit = 0;
  for (size_t i=0; i<z_s.size(); i++)
    if (z_s[i] > m_redshift) {
      inv1_Sigma_crit += m_cosmology->D_A(m_redshift, z_s[i]) / m_cosmology->D_A(0., z_s[i]) * N_zs[i];
      inv2_Sigma_crit += pow(m_cosmology->D_A(m_redshift, z_s[i]) / m_cosmology->D_A(0., z_s[i]), 2) * N_zs[i];
    }
  
  double inv_Sigma_crit = inv2_Sigma_crit / inv1_Sigma_crit * 4 * cbl::par::pi * m_G_newt / (m_c_light * m_c_light) * m_cosmology->D_A(0., m_redshift) * 1.e12;

  // Convergence, shear, and reduced shear
  std::vector<double> convergence = (this->*m_Sigma_ptr)(rad);
  std::vector<double> shear = (this->*m_DeltaSigma_ptr)(rad);
  std::vector<double> reduced_shear (rad.size(), 0.);

  for (size_t i=0; i<rad.size(); i++) {
    convergence[i] *= inv_Sigma_crit;
    shear[i] *= inv_Sigma_crit;
    reduced_shear[i] = shear[i] / (1 - convergence[i]);
  }
  
  return reduced_shear;
}


// ============================================================================================


double cbl::cosmology::HaloProfile::density_profile_FourierSpace (const double kk)
{  
  const double conc = concentration();
  const double rho_s = m_cosmology->rho_crit(m_redshift)*m_cosmology->Delta_c(m_redshift)/3.*pow(conc, 3)/(log(1.+conc)-conc/(1.+conc));
  const double r_s = r_vir(m_mass, m_redshift)/conc;  
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


// =====================================================================================


double cbl::cosmology::HaloProfile::concentration_NFW_Duffy (const double Mass, const double redshift, const string halo_def) const
{
  double AA, BB, CC;

  if (redshift>2) ErrorCBL("The concentration-mass relation by Duffy et al. has been tested only at z<2", "concentration_NFW_Duffy", "Cosmology.cpp");

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
    else return ErrorCBL("halo_def not allowed!", "concentration_NFW_Duffy", "Cosmology.cpp");

  const double Mpivot = 2.e12; // in [Msun/h]
  
  return AA*pow(Mass/Mpivot, BB)*pow(1.+redshift, CC);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::M_vir (const double r_vir, const double redshift, const std::string author, const bool unit1) const
{
  return 4./3.*par::pi*pow(r_vir, 3)*m_cosmology->Delta_c(redshift, author)*m_cosmology->rho_crit(redshift, unit1);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::r_vir (const double M_vir, const double redshift, const std::string author, const bool unit1) const
{
  return pow(3*M_vir/(4.*par::pi*m_cosmology->Delta_c(redshift, author)*m_cosmology->rho_crit(redshift, unit1)), 1./3.);
}


// =====================================================================================


double cbl::cosmology::HaloProfile::c_vir (const double c_200, const double redshift, const std::string author) const
{
  const double a = -1.119*log10(m_cosmology->Delta_c(redshift, author))+3.537;
  const double b = -0.967*log10(m_cosmology->Delta_c(redshift, author))+2.181;
  return a*c_200+b;
}

