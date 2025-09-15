/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file
 *  Modelling/NumberCounts/ModelFunction_NumberCounts.cpp
 *
 *  @brief Functions to model the number counts
 *
 *  This file contains the implementation of the functions used to
 *  model the monopole of the number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo, Giorgio Lesci
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, giorgio.lesci2@unibo.it
 */


#include "ModelFunction_NumberCounts.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


double cbl::modelling::numbercounts::Filter_sigmaR (const double kk, const double radius)
{
  return pow(TopHat_WF(kk*radius),2);
}


// ===========================================================================================


double cbl::modelling::numbercounts::Filter_dsigmaR (const double kk, const double radius)
{
  return 2.*TopHat_WF(kk*radius)*TopHat_WF_D1(kk*radius)*kk;
}


// ===========================================================================================


void cbl::modelling::numbercounts::sigmaM_dlnsigmaM (double &sigmaM, double &dlnsigmaM, const double mass, const cbl::glob::FuncGrid interp_Pk, const double kmax, const double rho)
{
  double norm =  1./(2.*pow(par::pi, 2));
  double dRdM_fact = pow(3./(4.*par::pi*rho), 1./3.);

  double RR =  Radius(mass, rho);
  double dRdM =  dRdM_fact*pow(mass, -2./3.)/3.;

  auto integrand_sigmaR = [&] (const double kk)
  {
    return kk*kk*interp_Pk(kk)*Filter_sigmaR(kk, RR);
  };

  sigmaM = norm*wrapper::gsl::GSL_integrate_cquad(integrand_sigmaR, 1.e-4, kmax, 1.e-5);

  auto integrand_dsigmaR = [&] (const double kk)
  {
    return kk*kk*interp_Pk(kk)*Filter_dsigmaR(kk, RR);
  };

  dlnsigmaM = norm*wrapper::gsl::GSL_integrate_cquad(integrand_dsigmaR, 1.e-4, kmax, 1.e-5)*dRdM*(mass/(2*sigmaM));
  sigmaM = sqrt(sigmaM);
}


// ===========================================================================================


void cbl::modelling::numbercounts::sigmaM_dlnsigmaM (std::vector<double> &sigmaM, std::vector<double> &dlnsigmaM, const std::vector<double> mass, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax, const double rho)
{
  double norm =  1./(2.*pow(par::pi, 2));
  double dRdM_fact = pow(3./(4.*par::pi*rho), 1./3.);
  glob::FuncGrid Pk_interp(kk, Pk, interpType);
  sigmaM.resize(mass.size(), 0);
  dlnsigmaM.resize(mass.size(), 0);

  for (size_t i=0; i<mass.size(); i++) {

    double RR =  Radius(mass[i], rho);
    double dRdM =  dRdM_fact*pow(mass[i], -2./3.)/3.;

    auto integrand_sigmaR = [&] (const double kk)
    {
      return kk*kk*Pk_interp(kk)*Filter_sigmaR(kk, RR);
    };

    sigmaM[i] = norm*wrapper::gsl::GSL_integrate_cquad(integrand_sigmaR, 1.e-4, kmax, 1.e-5);

    auto integrand_dsigmaR = [&] (const double kk)
    {
      return kk*kk*Pk_interp(kk)*Filter_dsigmaR(kk, RR);
    };

    dlnsigmaM[i] = norm*wrapper::gsl::GSL_integrate_cquad(integrand_dsigmaR, 1.e-4, kmax, 1.e-5)*dRdM*(mass[i]/(2*sigmaM[i]));
    sigmaM[i] = sqrt(sigmaM[i]);
  }
}


// ===========================================================================================


std::vector<cbl::glob::FuncGrid> cbl::modelling::numbercounts::sigmaM_dlnsigmaM (const std::vector<double> mass, const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  const double rho = cosmology->rho_m(0., true);

  vector<double> sigmaM, dlnsigmaM; 

  modelling::numbercounts::sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, mass, kk, Pk, interpType, kmax, rho);

  vector<glob::FuncGrid> interp(2);
  interp[0] = glob::FuncGrid(mass, sigmaM, interpType);
  interp[1] = glob::FuncGrid(mass, dlnsigmaM, interpType);

  return interp;
}


// ===========================================================================================


double cbl::modelling::numbercounts::mass_function (const double mass, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_Pk, const double kmax)
{
  const double rho = cosmology->rho_m(0., true);

  double sigmaM, dlnsigmaM;

  sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, ((cosmology->unit()) ? mass : mass*cosmology->little_h()), interp_Pk, kmax, rho);
  double _Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;

  cosmology::MassFunction MF(cosmology);
  return MF.mass_function(mass, sigmaM, dlnsigmaM, redshift, model_MF, store_output, par::defaultString, _Delta);
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::mass_function (const std::vector<double> mass, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  vector<double> mass_function(mass.size());
  
  const double rho = cosmology->rho_m(0., true);

  vector<double> _mass = mass;

  if (!cosmology->unit())
    for(size_t i=0; i<mass.size(); i++)
      _mass[i] = mass[i]*cosmology->little_h();

  vector<double> sigmaM, dlnsigmaM;
  modelling::numbercounts::sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, _mass, kk, Pk, interpType, kmax, rho);

  double _Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;

  cosmology::MassFunction MF(cosmology);
  for (size_t i=0; i<mass.size(); i++) 
    mass_function[i] = MF.mass_function(mass[i], sigmaM[i], dlnsigmaM[i], redshift, model_MF, store_output, par::defaultString, _Delta);

  return mass_function;
}


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::mass_function (const std::vector<double> redshift, const std::vector<double> mass, const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  vector<vector<double>> mass_function(redshift.size(), vector<double>(mass.size()));
  
  const double rho = cosmology->rho_m(0., true);

  vector<double> _mass = mass;

  if (!cosmology->unit()) 
    for(size_t i=0; i<mass.size(); i++)
      _mass[i] = mass[i]*cosmology->little_h();

  vector<double> sigmaM, dlnsigmaM;
  modelling::numbercounts::sigmaM_dlnsigmaM(sigmaM, dlnsigmaM, _mass, kk, Pk, interpType, kmax, rho);

  cosmology::MassFunction MF(cosmology);
  for (size_t j=0; j<redshift.size(); j++) 
    for (size_t i=0; i<mass.size(); i++) 
      mass_function[j][i] = MF.mass_function(mass[i], sigmaM[i], dlnsigmaM[i], redshift[j], model_MF, store_output, par::defaultString, ((isDelta_critical) ? Delta/cosmology->Omega_matter(redshift[j]) : Delta));

  return mass_function;
}


// ===========================================================================================


double cbl::modelling::numbercounts::number_counts (const double redshift_min, const double redshift_max, const double Mass_min, const double Mass_max, const std::shared_ptr<cosmology::Cosmology> cosmology, const double Area, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const glob::FuncGrid interp_DN)
{
  double fact = (cosmology->unit()) ? 1 : cosmology->little_h();
  double nc = 0;
  cosmology::MassFunction MF(cosmology);

  auto integrand = [&] (const vector<double> x)
  {
    double Mass = pow(10, x[0])*pow(10, 14);
    return MF.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), x[1], interp_DN(x[1]), model_MF, store_output, par::defaultString, (isDelta_critical) ? Delta/cosmology->Omega_matter(x[1]) : Delta) * Area * cosmology->dV_dZdOmega(x[1], true) * pow(10,x[0]);
  };
  vector<vector<double>> integration_limits(2);
  integration_limits[0] = {log10(Mass_min/pow(10,14)), log10(Mass_max/pow(10,14))};
  integration_limits[1] = {redshift_min, redshift_max};

  wrapper::cuba::CUBAwrapper CW(integrand, 2);
  nc = CW.IntegrateVegas(integration_limits, false);
  
  return nc * pow(10,14) * log(10);
}


// ===========================================================================================


double cbl::modelling::numbercounts::counts_proxy (const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_bias, const double proxy_bias, const double z_err, const double proxy_err, const double Plambda_a, const double Plambda_b, const double Plambda_c, ProxyFunction fz, std::function<double(const double, const double)> z_error, std::function<double(const double, const double)> proxy_error, double (*response_fact)(const double, const double, const double, const double, const std::string, const double, const std::string, std::shared_ptr<void>), const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, const std::shared_ptr<cosmology::Cosmology> cosmology, const double Area, const std::string model_MF, const std::string model_bias, const bool store_output, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base, const double weight)
{  
  double fact = (cosmology->unit()) ? 1 : cosmology->little_h();
  shared_ptr<void> pp;
  auto cosmology_ptr = move(cosmology);

  // Declare the normalized mass and the redshift, used as constants in integrand_P_M__z (which is called in integrand)
  double normM=0; double the_redsh=0;

  // P(M|z) integrand
  auto integrand_P_M__z = [&] (const double x)
  {
    double log_lambda = x - log(proxy_pivot)/log(log_base);
    double log_f_z = log( fz(the_redsh, z_pivot, cosmology_ptr) )/log(log_base);
      
    double mean = alpha + beta*log_lambda + gamma*log_f_z;
    double sigma = abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
    double P_M__lambda_z = gaussian(normM, pp, {mean,sigma});      
    double P_lambda__z = Plambda_a * pow(pow(log_base,x),-Plambda_b) * exp(-Plambda_c*pow(log_base,x));
			    
    return P_M__lambda_z * P_lambda__z * pow(log_base,x);
  };

  // Total integrand
  cosmology::MassFunction MF(cosmology);
  
  auto integrand = [&] (const vector<double> x)
  {
    double Delta_ = (isDelta_critical) ? Delta/cosmology->Omega_matter(x[1]) : Delta;
    double Mass = pow(log_base,x[0])*mass_pivot;
    normM = x[0];
    the_redsh = x[1];

    // Compute P(M|lambda,z)
    double log_lambda = log(x[2]/proxy_pivot)/log(log_base);
    double log_f_z = log( fz(x[1], z_pivot, cosmology_ptr) )/log(log_base);
      
    double mean = alpha + beta*log_lambda + gamma*log_f_z;
    double sigma = abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
    double P_M__lambda_z = gaussian(normM, pp, {mean,sigma});

    // Compute P(lambda|z)
    double P_lambda__z = Plambda_a * pow(x[2],-Plambda_b) * exp(-Plambda_c*x[2]);

    // Compute P(M|z)
    double P_M__z=0;
    if (P_M__lambda_z*P_lambda__z > 0)
      P_M__z=wrapper::gsl::GSL_integrate_cquad(integrand_P_M__z,log(0.00001)/log(log_base),log(1000.)/log(log_base)) * log(log_base);
    else
      P_M__z = 1;

    // Compute the integrals of P(z|z) and P(lambda|lambda)
    double mean_Pz = x[1] + z_bias * (1+x[1]);
    double int_P_z = 0.5 * ( erf( (redshift_max - mean_Pz) / (sqrt(2)*z_error(z_err, redshift_max)) ) - erf( (redshift_min - mean_Pz) / (sqrt(2)*z_error(z_err, redshift_min)) ) );
    double mean_Plambda = x[2] + proxy_bias * (x[2]);
    double int_P_lambda = 0.5 * ( erf( (proxy_max - mean_Plambda) / (sqrt(2)*proxy_error(proxy_err, proxy_max)) ) - erf( (proxy_min - mean_Plambda) / (sqrt(2)*proxy_error(proxy_err, proxy_min)) ) );
    
    return response_fact(Mass, interp_sigmaM(Mass*fact), x[1], interp_DN(x[1]), model_bias, Delta_, "EisensteinHu", cosmology_ptr) * MF.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), x[1], interp_DN(x[1]), model_MF, store_output, par::defaultString, Delta_)*Area*cosmology->dV_dZdOmega(x[1], true) * pow(log_base,normM) * (P_M__lambda_z*P_lambda__z/P_M__z) * int_P_z * int_P_lambda;
  };
  
  // -------------------------------------------------------------

  // Find the minimum and maximum masses, given the parameters of the scaling relation
  double log_lambda_min = log((max(proxy_min - 3.5*proxy_error(proxy_err, proxy_min), 1.))/proxy_pivot)/log(log_base);
  double log_lambda_max = log((proxy_max + 3.5*proxy_error(proxy_err, proxy_max))/proxy_pivot)/log(log_base);
  double log_f_z_min = log( fz((max(redshift_min - 3.5*z_error(z_err, redshift_min), 0.)), z_pivot, cosmology_ptr) )/log(log_base);
  double log_f_z_max = log( fz((redshift_max + 3.5*z_error(z_err, redshift_min)), z_pivot, cosmology_ptr) )/log(log_base);

  double M1 = alpha + beta*log_lambda_min + gamma*log_f_z_min;
  double M2 = alpha + beta*log_lambda_max + gamma*log_f_z_min;
  double M3 = alpha + beta*log_lambda_min + gamma*log_f_z_max;
  double M4 = alpha + beta*log_lambda_max + gamma*log_f_z_max;

  double min1 = min(M1, M2);
  double min2 = min(min1, M3);
  double minM = min(min2, M4);
  double max1 = max(M1, M2);
  double max2 = max(max1, M3);
  double maxM = max(max2, M4);

  // Find the maximum value of the intrinsic scatter
  double s1 = abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
  double s2 = abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
  double s3 = abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );
  double s4 = abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );

  double maxs1 = max(s1, s2);
  double maxs2 = max(maxs1, s3);
  double max_intrinsic_scatter = max(maxs2, s4);

  // Define the integral limits
  int integral_dimension=3;
  vector<vector<double>> integration_limits(integral_dimension);
  integration_limits[0] = {max(minM-3.5*max_intrinsic_scatter,log(1.e10/mass_pivot)/log(log_base)), min(maxM+3.5*max_intrinsic_scatter,log(1.e16/mass_pivot)/log(log_base))};
  integration_limits[1] = {max(redshift_min - 3.5*z_error(z_err, redshift_min), 0.), redshift_max + 3.5*z_error(z_err, redshift_max)};
  integration_limits[2] = {max(proxy_min - 3.5*proxy_error(proxy_err, proxy_min), 0.00001), proxy_max + 3.5*proxy_error(proxy_err, proxy_max)};

  // Compute the integral
  wrapper::cuba::CUBAwrapper CW (integrand, integral_dimension);
  double nc;

  if (integration_limits[0][0] < integration_limits[0][1])
    nc = CW.IntegrateVegas(integration_limits,false);
  else
    nc = 0;

  return nc * mass_pivot * log(log_base) * weight;
}


// ===========================================================================================


double cbl::modelling::numbercounts::counts_proxy_classic(const std::shared_ptr<cosmology::Cosmology> cosmo, const double area, const double z_min, const double z_max, const double proxy_min, const double proxy_max, const double z_error_normalised, const std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc, const double Plambda_A_mu, const double Plambda_B_mu, const double Plambda_C_mu, const std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc, const double Plambda_A_sigma, const double Plambda_B_sigma, const double Plambda_C_sigma, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double purity, cbl::glob::FuncGrid2D completeness_interp, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double log_base, const bool isDelta_critical, const double Delta, const std::string model_MF, double (*response_fact)(const double, const double, const double, const double, const std::string, const double, const std::string, std::shared_ptr<void>), const std::string model_bias, const double q_HMF_correction, const double s_HMF_correction, double (*mass_conversion)(const double, const double, std::shared_ptr<cbl::cosmology::Cosmology>))
{
  if (*mass_conversion == NULL) {

    mass_conversion = [] (const double MDelta, const double z, std::shared_ptr<cbl::cosmology::Cosmology> cosmology) {
      (void)z; (void)cosmology; return MDelta;
    };
    
  }
  
  // true mass, proxy, and redshift ranges covered in the integrations
  // used for the following interpolations
  std::vector<std::vector<double>> int_limits = cbl::modelling::numbercounts::get_integration_limits
    (
     cosmo, z_min, z_max, proxy_min, proxy_max, z_error_normalised,
     Plambda_mean_fc, Plambda_A_mu, Plambda_B_mu, Plambda_C_mu, Plambda_std_fc, Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma,
     scaling_relation, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, z_pivot, proxy_pivot, mass_pivot, log_base
     ); // log(M/M_piv), proxy_tr, z_tr. In this order.

  // true redshift and log_mass_norm vectors used for interpolations
  const std::vector<double> redshift_vector = cbl::linear_bin_vector(20, int_limits[2][0], int_limits[2][1]);
  const std::vector<double> logM_norm_vector = cbl::linear_bin_vector(50, int_limits[0][0], int_limits[0][1]);
  
  // define and interpolate the integral over proxy_tr as a function of M and z_tr
  std::shared_ptr<void> ptr;
  double dummy_logM_norm = 0, dummy_z_tr = 0;
  
  std::function<double(double)> proxy_tr_dependent_integrand =
    [&](double proxy_tr) -> double {
      
      // P(proxy_tr|M,z)
      double log_proxy_tr = log(proxy_tr/proxy_pivot) / log(log_base);
      double mean = modelling::numbercounts::logProxy(dummy_z_tr, pow(log_base,dummy_logM_norm)*mass_pivot, alpha, beta, gamma, cosmo, scaling_relation, z_pivot, mass_pivot, log_base);
      double sigma = modelling::numbercounts::scatter_intr_proxy(dummy_z_tr, pow(log_base,dummy_logM_norm)*mass_pivot, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo, scaling_relation, z_pivot, mass_pivot, log_base);

      double P_proxy_M_z = 0.;
      if (sigma > 1.e-3)
	P_proxy_M_z = cbl::gaussian(log_proxy_tr, ptr, {mean, sigma}) / (log(log_base) * proxy_tr);
      else {
	if (std::abs(log_proxy_tr - mean) <= 1.e-3)
	  P_proxy_M_z = 1.;
      }

      // Integral of P(z_ob|z_tr)
      double mean_Pz = dummy_z_tr;
      double std_Pz = z_error_normalised * (1 + dummy_z_tr);
      double int_Pz = 0.5 * ( erf( (z_max - mean_Pz) / (sqrt(2) * std_Pz) ) - erf( (z_min - mean_Pz) / (sqrt(2) * std_Pz) ) );

      // Integral of P(proxy_ob|proxy_tr,z_tr)
      double mean_Pproxy = Plambda_mean_fc(proxy_tr, dummy_z_tr, Plambda_A_mu, Plambda_B_mu, Plambda_C_mu);
      double std_Pproxy = Plambda_std_fc(proxy_tr, dummy_z_tr, Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma);
      double int_Pproxy = 0.5 * ( erf( (proxy_max - mean_Pproxy) / (sqrt(2) * std_Pproxy) ) - erf( (proxy_min - mean_Pproxy) / (sqrt(2) * std_Pproxy) ) );

      double ztr2 = dummy_z_tr, proxytr2 = proxy_tr;
      if (ztr2 < completeness_interp.xmin())
	ztr2 = completeness_interp.xmin();
      if (ztr2 > completeness_interp.xmax())
	ztr2 = completeness_interp.xmax();
      if (proxytr2 < completeness_interp.ymin())
	proxytr2 = completeness_interp.ymin();
      if (proxytr2 > completeness_interp.ymax())
	proxytr2 = completeness_interp.ymax();
      
      return P_proxy_M_z * int_Pz * int_Pproxy * std::max(0., completeness_interp(ztr2, proxytr2));
            
    };

  std::vector<std::vector<double>> int_lambda_tr_samples(redshift_vector.size(), std::vector<double>(logM_norm_vector.size()));

  for (size_t i=0; i<redshift_vector.size(); i++) {
    for (size_t j=0; j<logM_norm_vector.size(); j++) {
      dummy_z_tr = redshift_vector[i];
      dummy_logM_norm = logM_norm_vector[j];
      int_lambda_tr_samples[i][j] = wrapper::gsl::GSL_integrate_cquad(proxy_tr_dependent_integrand, int_limits[1][0], int_limits[1][1]);
    }
  }

  cbl::glob::FuncGrid2D int_lambda_tr_interp (redshift_vector, logM_norm_vector, int_lambda_tr_samples, "Cubic");

  // define the full integrand and interpolate its integral over M, as a function of z
  cosmology::MassFunction MF(cosmo);
  double fact = (cosmo->unit()) ? 1 : cosmo->little_h();
  
  std::function<double(double)> counts_integrand =
    [&] (const double logM_norm) {

      double z_tr = dummy_z_tr;
      double Delta_new = (isDelta_critical) ? Delta / cosmo->Omega_matter(z_tr) : Delta;
      double mass = pow(log_base, logM_norm) * mass_pivot;

      double HMF = MF.mass_function(mass, interp_sigmaM(mass * fact), interp_DlnsigmaM(mass * fact), z_tr, interp_DN(z_tr), model_MF, false, par::defaultString, Delta_new);
      HMF *= s_HMF_correction * (log10(mass_conversion(mass * fact, z_tr, cosmo)) - 13.8) + q_HMF_correction;
      double counts_norm = response_fact(mass, interp_sigmaM(mass*fact), z_tr, interp_DN(z_tr), model_bias, Delta_new, "EisensteinHu", cosmo) * cosmo->dV_dZdOmega(z_tr, true) * HMF * (log(log_base) * mass);
      // the last factor above stems from the change of integration variable

      return counts_norm * int_lambda_tr_interp(z_tr, logM_norm);
            
    };
  
  std::vector<double> int_over_M_proxy_samples(redshift_vector.size());
  for (size_t i=0; i<redshift_vector.size(); i++) {
    dummy_z_tr = redshift_vector[i];
    int_over_M_proxy_samples[i] = wrapper::gsl::GSL_integrate_qag(counts_integrand, int_limits[0][0], int_limits[0][1]);
  }
  cbl::glob::FuncGrid int_over_M_proxy_interp(redshift_vector, int_over_M_proxy_samples, "Spline");

  // compute the counts
  std::function<double(double)> final_integrand =
    [&] (const double z_tr) {

      return int_over_M_proxy_interp(z_tr);
            
    };
  
  return area * wrapper::gsl::GSL_integrate_qag(final_integrand, int_limits[2][0], int_limits[2][1]) / purity;
}


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::get_integration_limits (const std::shared_ptr<cosmology::Cosmology> cosmology, const double z_min, const double z_max, const double mass_proxy_min, const double mass_proxy_max, const double z_error, const std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc, const double Plambda_A_mu, const double Plambda_B_mu, const double Plambda_C_mu, const std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc, const double Plambda_A_sigma, const double Plambda_B_sigma, const double Plambda_C_sigma, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base)
{
  (void) Plambda_mean_fc; (void) Plambda_A_mu; (void) Plambda_B_mu; (void) Plambda_C_mu;
  
  auto cosmo_ptr = move(cosmology);
  std::vector<std::vector<double>> int_limits(3, std::vector<double>(2)); // log(M/M_piv), proxy_tr, z_tr. In this order.

  // z_tr
  int_limits[2] = {std::max(0.0001, z_min - 3.5 * z_error * (1. + z_min)), z_max + 3.5 * z_error * (1. + z_max)};

  // proxy_tr  
  std::vector<double> proxies =
    {
      mass_proxy_min - 3.5 * Plambda_std_fc(mass_proxy_min, int_limits[2][0], Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma),
      mass_proxy_min - 3.5 * Plambda_std_fc(mass_proxy_min, int_limits[2][1], Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma),
      mass_proxy_max + 3.5 * Plambda_std_fc(mass_proxy_max, int_limits[2][0], Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma),
      mass_proxy_max + 3.5 * Plambda_std_fc(mass_proxy_max, int_limits[2][1], Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma),
    };
  
  int_limits[1] = {std::max( 0.0001, cbl::Min(proxies) ), cbl::Max(proxies)};

  // M
  std::vector<double> scatters =
    {
      cbl::modelling::numbercounts::scatter_intr_proxy(int_limits[2][0], 1.e16, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base),
      cbl::modelling::numbercounts::scatter_intr_proxy(int_limits[2][1], 1.e16, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base),
      cbl::modelling::numbercounts::scatter_intr_proxy(int_limits[2][0], 1.e6, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base),
      cbl::modelling::numbercounts::scatter_intr_proxy(int_limits[2][1], 1.e6, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base)
    };
  double scatter_intr_max = cbl::Max(scatters);
    
  std::vector<double> logz =
    {
      log( (scaling_relation)->data_model().fz(int_limits[2][0], z_pivot, cosmo_ptr) ) / log(logM_base),
      log( (scaling_relation)->data_model().fz(int_limits[2][1], z_pivot, cosmo_ptr) ) / log(logM_base)
    };

  std::vector<double> logproxy =
    {
      log(int_limits[1][0] / proxy_pivot) / log(logM_base) - 3.5 * scatter_intr_max,
      log(int_limits[1][1] / proxy_pivot) / log(logM_base) + 3.5 * scatter_intr_max
    };
  
  int_limits[0] =
    {
      std::min( logproxy[0] - gamma * cbl::Min(logz) - alpha, logproxy[0] - gamma * cbl::Max(logz) - alpha ) / beta,
      std::min(log(1.e16 / mass_pivot) / log(logM_base), std::max( logproxy[1] - gamma * cbl::Min(logz) - alpha, logproxy[1] - gamma * cbl::Max(logz) - alpha ) / beta)
    };

  if (int_limits[0][0] > int_limits[0][1])
    int_limits[0][0] = 0.99 * int_limits[0][1];

  return int_limits;
}


// ===========================================================================================


double cbl::modelling::numbercounts::scatter_intr_proxy (const double z_tr, const double M, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double M_pivot, const double logM_base)
{
  auto cosmo_ptr = move(cosmo);
  double log_M = log(M/M_pivot) / log(logM_base);
  double log_f_z_tr = log( (scaling_relation)->data_model().fz(z_tr, z_pivot, cosmo_ptr) ) / log(logM_base);
  return std::max( 0.01, std::abs(scatter0 + scatterM*pow(log_M, scatterM_exp) + scatterz*pow(log_f_z_tr, scatterz_exp)) );
}


// ===========================================================================================


double cbl::modelling::numbercounts::logProxy (const double z_tr, const double M, const double alpha, const double beta, const double gamma, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double M_pivot, const double logM_base)
{
  auto cosmo_ptr = move(cosmo);
  double log_M = log(M/M_pivot) / log(logM_base);
  double log_f_z_tr = log( (scaling_relation)->data_model().fz(z_tr, z_pivot, cosmo_ptr) ) / log(logM_base);
  return alpha + beta*log_M + gamma*log_f_z_tr;
}



// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::size_function (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> radii, const double redshift, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file)
{
  cosmology::SizeFunction SZ(cosmology);
  vector<double> size_function = SZ.size_function(radii, redshift, model, b_eff, slope, offset, deltav_NL, del_c, method_Pk, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  return size_function; 
}
