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
  return 2.*cbl::TopHat_WF(kk*radius)*cbl::TopHat_WF_D1(kk*radius)*kk;
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

  sigmaM = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_sigmaR, 1.e-4, kmax, 1.e-5);

  auto integrand_dsigmaR = [&] (const double kk)
			   {
			     return kk*kk*interp_Pk(kk)*Filter_dsigmaR(kk, RR);
			   };

  dlnsigmaM = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_dsigmaR, 1.e-4, kmax, 1.e-5)*dRdM*(mass/(2*sigmaM));
  sigmaM = sqrt(sigmaM);
}


// ===========================================================================================


void cbl::modelling::numbercounts::sigmaM_dlnsigmaM (std::vector<double> &sigmaM, std::vector<double> &dlnsigmaM, const std::vector<double> mass, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax, const double rho)
{
  double norm =  1./(2.*pow(par::pi, 2));
  double dRdM_fact = pow(3./(4.*par::pi*rho), 1./3.);
  cbl::glob::FuncGrid Pk_interp(kk, Pk, interpType);
  sigmaM.resize(mass.size(), 0);
  dlnsigmaM.resize(mass.size(), 0);

  for (size_t i=0; i<mass.size(); i++) {

    double RR =  Radius(mass[i], rho);
    double dRdM =  dRdM_fact*pow(mass[i], -2./3.)/3.;

    auto integrand_sigmaR = [&] (const double kk)
			    {
			      return kk*kk*Pk_interp(kk)*Filter_sigmaR(kk, RR);
			    };

    sigmaM[i] = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_sigmaR, 1.e-4, kmax, 1.e-5);

    auto integrand_dsigmaR = [&] (const double kk)
			     {
			       return kk*kk*Pk_interp(kk)*Filter_dsigmaR(kk, RR);
			     };

    dlnsigmaM[i] = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_dsigmaR, 1.e-4, kmax, 1.e-5)*dRdM*(mass[i]/(2*sigmaM[i]));
    sigmaM[i] = sqrt(sigmaM[i]);
  }
}


// ===========================================================================================


std::vector<cbl::glob::FuncGrid> cbl::modelling::numbercounts::sigmaM_dlnsigmaM (const std::vector<double> mass, cbl::cosmology::Cosmology cosmology, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  const double rho = cosmology.rho_m(0., true);

  vector<double> sigmaM, dlnsigmaM; 

  cbl::modelling::numbercounts::sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, mass, kk, Pk, interpType, kmax, rho);

  vector<cbl::glob::FuncGrid> interp(2);
  interp[0] = cbl::glob::FuncGrid(mass, sigmaM, interpType);
  interp[1] = cbl::glob::FuncGrid(mass, dlnsigmaM, interpType);

  return interp;
}


// ===========================================================================================


double cbl::modelling::numbercounts::mass_function (const double mass, cbl::cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_Pk, const double kmax)
{
  const double rho = cosmology.rho_m(0., true);

  double sigmaM, dlnsigmaM;

  sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, ((cosmology.unit()) ? mass : mass*cosmology.hh()), interp_Pk, kmax, rho);
  double _Delta = (isDelta_critical) ? Delta/cosmology.OmegaM(redshift) : Delta;

  return cosmology.mass_function(mass, sigmaM, dlnsigmaM, redshift, model_MF, store_output, cbl::par::defaultString, _Delta);
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::mass_function (const std::vector<double> mass, cbl::cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  vector<double> mass_function(mass.size());
  
  const double rho = cosmology.rho_m(0., true);

  vector<double> _mass = mass;

  if (!cosmology.unit())
    for(size_t i=0; i<mass.size(); i++)
      _mass[i] = mass[i]*cosmology.hh();

  vector<double> sigmaM, dlnsigmaM;
  cbl::modelling::numbercounts::sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, _mass, kk, Pk, interpType, kmax, rho);

  double _Delta = (isDelta_critical) ? Delta/cosmology.OmegaM(redshift) : Delta;

  for (size_t i=0; i<mass.size(); i++) 
    mass_function[i] = cosmology.mass_function(mass[i], sigmaM[i], dlnsigmaM[i], redshift, model_MF, store_output, cbl::par::defaultString, _Delta);

  return mass_function;
}


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::mass_function (const std::vector<double> redshift, const std::vector<double> mass, cbl::cosmology::Cosmology cosmology, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  vector<vector<double>> mass_function(redshift.size(), vector<double>(mass.size()));
  
  const double rho = cosmology.rho_m(0., true);

  vector<double> _mass = mass;

  if (!cosmology.unit()) 
    for(size_t i=0; i<mass.size(); i++)
      _mass[i] = mass[i]*cosmology.hh();

  vector<double> sigmaM, dlnsigmaM;
  cbl::modelling::numbercounts::sigmaM_dlnsigmaM(sigmaM, dlnsigmaM, _mass, kk, Pk, interpType, kmax, rho);

  for (size_t j=0; j<redshift.size(); j++) 
    for (size_t i=0; i<mass.size(); i++) 
      mass_function[j][i] = cosmology.mass_function(mass[i], sigmaM[i], dlnsigmaM[i], redshift[j], model_MF, store_output, par::defaultString, ((isDelta_critical) ? Delta/cosmology.OmegaM(redshift[j]) : Delta));

  return mass_function;
}


// ===========================================================================================


double cbl::modelling::numbercounts::number_counts (const double redshift_min, const double redshift_max, const double Mass_min, const double Mass_max, cbl::cosmology::Cosmology cosmology, const double Area, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const  cbl::glob::FuncGrid interp_DlnsigmaM)
{
  double fact = (cosmology.unit()) ? 1 : cosmology.hh();
  double nc = 0;

  auto integrand = [&cosmology,&fact,&model_MF,&isDelta_critical,&Delta,&Area,&interp_sigmaM,&interp_DlnsigmaM,&store_output] (const std::vector<double> x)
		   {
		     double Mass = pow(10,x[0])*pow(10,14);
		     return cosmology.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), x[1], model_MF, store_output, cbl::par::defaultString, (isDelta_critical) ? Delta/cosmology.OmegaM(x[1]) : Delta)*Area*cosmology.dV_dZdOmega(x[1], true) * pow(10,x[0]);
		   };
  std::vector<std::vector<double>> integration_limits(2);
  integration_limits[0] = {log10(Mass_min/pow(10,14)), log10(Mass_max/pow(10,14))};
  integration_limits[1] = {redshift_min, redshift_max};

  cbl::wrapper::cuba::CUBAwrapper CW(integrand, 2);
  nc = CW.IntegrateVegas(integration_limits, false);
  
  return nc * pow(10,14) * log(10);
}


// ===========================================================================================


double cbl::modelling::numbercounts::counts_proxy (const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_bias, const double proxy_bias, const double z_err, const double proxy_err, const double Plambda_a, const double Plambda_b, const double Plambda_c, std::function<double(const double, const double, const std::shared_ptr<void>)> fz, std::function<double(const double, const double)> z_error, std::function<double(const double, const double)> proxy_error, double (*response_fact)(const double, const double, const double, const double, const std::string, const double, const std::string, std::shared_ptr<void>), const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, cbl::cosmology::Cosmology cosmology, const double Area, const std::string model_MF, const std::string model_bias, const bool store_output, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base, const double weight)
{  
  double fact = (cosmology.unit()) ? 1 : cosmology.hh();
  std::shared_ptr<void> pp;
  auto cosmology_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmology);

  // Declare the normalized mass and the redshift, used as constants in integrand_P_M__z (which is called in integrand)
  double normM=0; double the_redsh=0;

  // P(M|z) integrand
  auto integrand_P_M__z = [&] (const double x)
			  {
			    double log_lambda = x - log(proxy_pivot)/log(log_base);
			    double log_f_z = log( fz(the_redsh, z_pivot, cosmology_ptr) )/log(log_base);
      
			    double mean = alpha + beta*log_lambda + gamma*log_f_z;
			    double sigma = std::abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
			    double P_M__lambda_z = cbl::gaussian(normM, pp, {mean,sigma});      
			    double P_lambda__z = Plambda_a * pow(pow(log_base,x),-Plambda_b) * exp(-Plambda_c*pow(log_base,x));
			    
			    return P_M__lambda_z * P_lambda__z * pow(log_base,x);
			  };

  // Total integrand
  auto integrand = [&] (const std::vector<double> x)
		   {
		     double Delta_ = (isDelta_critical) ? Delta/cosmology.OmegaM(x[1]) : Delta;
		     double Mass = pow(log_base,x[0])*mass_pivot;
		     normM = x[0];
		     the_redsh = x[1];

		     // Compute P(M|lambda,z)
		     double log_lambda = log(x[2]/proxy_pivot)/log(log_base);
		     double log_f_z = log( fz(x[1], z_pivot, cosmology_ptr) )/log(log_base);
      
		     double mean = alpha + beta*log_lambda + gamma*log_f_z;
		     double sigma = std::abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
		     double P_M__lambda_z = cbl::gaussian(normM, pp, {mean,sigma});

		     // Compute P(lambda|z)
		     double P_lambda__z = Plambda_a * pow(x[2],-Plambda_b) * exp(-Plambda_c*x[2]);

		     // Compute P(M|z)
		     double P_M__z=0;
		     if (P_M__lambda_z*P_lambda__z > 0)
		       P_M__z=cbl::wrapper::gsl::GSL_integrate_cquad(integrand_P_M__z,log(0.00001)/log(log_base),log(1000.)/log(log_base)) * log(log_base);
		     else
		       P_M__z = 1;

		     // Compute the integrals of P(z|z) and P(lambda|lambda)
		     double mean_Pz = x[1] + z_bias * (1+x[1]);
		     double int_P_z = 0.5 * ( erf( (redshift_max - mean_Pz) / (sqrt(2)*z_error(z_err, redshift_max)) ) - erf( (redshift_min - mean_Pz) / (sqrt(2)*z_error(z_err, redshift_min)) ) );
		     double mean_Plambda = x[2] + proxy_bias * (x[2]);
		     double int_P_lambda = 0.5 * ( erf( (proxy_max - mean_Plambda) / (sqrt(2)*proxy_error(proxy_err, proxy_max)) ) - erf( (proxy_min - mean_Plambda) / (sqrt(2)*proxy_error(proxy_err, proxy_min)) ) );
      
		     return response_fact(Mass, interp_sigmaM(Mass*fact), x[1], interp_DN(x[1]), model_bias, Delta_, "EisensteinHu", cosmology_ptr) * cosmology.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), x[1], interp_DN(x[1]), model_MF, store_output, cbl::par::defaultString, Delta_)*Area*cosmology.dV_dZdOmega(x[1], true) * pow(log_base,normM) * (P_M__lambda_z*P_lambda__z/P_M__z) * int_P_z * int_P_lambda;
		   };
  
  // -------------------------------------------------------------

  // Find the minimum and maximum masses, given the parameters of the scaling relation
  double log_lambda_min = log((std::max(proxy_min - 3.5*proxy_error(proxy_err, proxy_min), 1.))/proxy_pivot)/log(log_base);
  double log_lambda_max = log((proxy_max + 3.5*proxy_error(proxy_err, proxy_max))/proxy_pivot)/log(log_base);
  double log_f_z_min = log( fz((std::max(redshift_min - 3.5*z_error(z_err, redshift_min), 0.)), z_pivot, cosmology_ptr) )/log(log_base);
  double log_f_z_max = log( fz((redshift_max + 3.5*z_error(z_err, redshift_min)), z_pivot, cosmology_ptr) )/log(log_base);

  double M1 = alpha + beta*log_lambda_min + gamma*log_f_z_min;
  double M2 = alpha + beta*log_lambda_max + gamma*log_f_z_min;
  double M3 = alpha + beta*log_lambda_min + gamma*log_f_z_max;
  double M4 = alpha + beta*log_lambda_max + gamma*log_f_z_max;

  double min1 = std::min(M1, M2);
  double min2 = std::min(min1, M3);
  double minM = std::min(min2, M4);
  double max1 = std::max(M1, M2);
  double max2 = std::max(max1, M3);
  double maxM = std::max(max2, M4);

  // Find the maximum value of the intrinsic scatter
  double s1 = std::abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
  double s2 = std::abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
  double s3 = std::abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );
  double s4 = std::abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );

  double maxs1 = std::max(s1, s2);
  double maxs2 = std::max(maxs1, s3);
  double max_intrinsic_scatter = std::max(maxs2, s4);

  // Define the integral limits
  int integral_dimension=3;
  std::vector<std::vector<double>> integration_limits(integral_dimension);
  integration_limits[0] = {std::max(minM-3.5*max_intrinsic_scatter,log(1.e10/mass_pivot)/log(log_base)), std::min(maxM+3.5*max_intrinsic_scatter,log(1.e16/mass_pivot)/log(log_base))};
  integration_limits[1] = {std::max(redshift_min - 3.5*z_error(z_err, redshift_min), 0.), redshift_max + 3.5*z_error(z_err, redshift_max)};
  integration_limits[2] = {std::max(proxy_min - 3.5*proxy_error(proxy_err, proxy_min), 0.00001), proxy_max + 3.5*proxy_error(proxy_err, proxy_max)};

  // Compute the integral
  cbl::wrapper::cuba::CUBAwrapper CW (integrand, integral_dimension);
  double nc;

  if (integration_limits[0][0] < integration_limits[0][1])
    nc = CW.IntegrateVegas(integration_limits,false);
  else
    nc = 0;

  return nc * mass_pivot * log(log_base) * weight;
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::size_function (cbl::cosmology::Cosmology cosmology, const std::vector<double> radii, const double redshift, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file)
{
  
  vector<double> size_function = cosmology.size_function(radii, redshift, model, b_eff, slope, offset, deltav_NL, del_c, method_Pk, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  return size_function;
  
}
