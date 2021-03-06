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
 *  @file
 *  Modelling/DensityProfile/Modelling_DensityProfile.cpp
 *
 *  @brief Methods of the class Modelling_DensityProfile
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_DensityProfile, i.e. the common functions to model
 *  the galaxy cluster surface density profiles
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "Modelling_DensityProfile.h"

using namespace std;

using namespace cbl;


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const double redshift, const double contrast, const double trunc_fact, const double logM_base)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.cosmology->set_unit(true); // Force cosmological units

  m_data_model.cluster = make_shared<catalogue::Cluster>(cluster);
  m_data_model.redshift = redshift;
  m_data_model.contrast = contrast;
  m_data_model.trunc_fact = trunc_fact;
  m_data_model.logM_base = logM_base;
}


// ===========================================================================================

void cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior)
{
  m_data_model.Cpar = cosmo_param;

  const size_t nParams = cosmo_param.size()+4; // The total number of parameters is given by the cosmological ones + 4, since the density profile has 4 parameters (conc, logM, f_off, sigma_off)

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams);

  // Set the names and priors of the cosmological parameters
  for (size_t i=0; i<cosmo_param.size(); i++){
    Par_string[i] = CosmologicalParameter_name(cosmo_param[i]);
    param_prior[i] = cosmo_prior[i];
  }

  // Set the names and priors for the density profile parameters
  Par_string[cosmo_param.size()] = "concentration";
  param_prior[cosmo_param.size()] = concentration_prior;
  Par_string[cosmo_param.size()+1] = "logM";
  param_prior[cosmo_param.size()+1] = logM_prior;
  Par_string[cosmo_param.size()+2] = "f_off";
  param_prior[cosmo_param.size()+2] = f_off_prior;
  Par_string[cosmo_param.size()+3] = "sigma_off";
  param_prior[cosmo_param.size()+3] = sigma_off_prior;

  // input data used to construct the model
  auto inputs = make_shared<STR_Profile_data_model>(m_data_model);

  // set prior
  m_set_prior(param_prior);

  // construct the model
  if (m_2halo == false)
    m_model = make_shared<statistics::Model1D>(statistics::Model1D(&nfw_1halo_allBins, nParams, Par_type, Par_string, inputs));
  else
    ErrorCBL("The 2-halo term is not implemented yet!", "set_model_DensityProfile_cosmology", "Modelling_DensityProfile.cpp");
}
				  

// ===========================================================================================

double cbl::modelling::densityprofile::nfw_1halo (cbl::cosmology::Cosmology cosm, cbl::catalogue::Cluster cluster, const double radius, const double redshift, const double contrast, const double trunc_fact)
{
  // Calculate the radius enclosing a given mass, in the cosmology considered
  double H = cosm.HH(redshift)/cosm.HH(0.)*100./3.0857e+19; // in sec^-1
  double rho_crit = 3*H*H/8/cbl::par::pi/6.6732e-8/1.98892e+33*3.0857e+24*3.0857e+24*3.0857e+24; // in h^2*M_sun/Mpc^3
  double r_encl = pow( 3*cluster.mass()/(4*cbl::par::pi*contrast*rho_crit), 1./3. );

  double rscale = r_encl/cluster.concentration();
  double trunc = trunc_fact*r_encl;
  double tau = trunc/rscale;

  std::function<double(double)> F = [] (double x){
				      // both operations must be generalised to complex numbers when x<1
				      // in that case, both num and den are imaginary, and their ratio is real again
				      std::complex<double> num = std::acos(std::complex<double>(1./x,0.));
				      std::complex<double> den = std::sqrt(std::complex<double>(x*x-1.,0.));  
				      std::complex<double> res = num/den;						 
				      return std::abs(real(res));
				    };
  std::function<double(double)> L = [&tau] (double x){
				      return log(x/(sqrt(x*x+tau*tau)+tau));
				    };

  double rho_s = 1; 
  std::function<double(double)> integrate_mass = [&tau,&rscale,&rho_s] (double rad){
						   double x = rad/rscale;
						   double rho = rho_s/x/(1+x)/(1+x)*(tau*tau/(tau*tau+x*x))*(tau*tau/(tau*tau+x*x)); // in (M_sun/h)/((Mpc/h)^3)
						   return 4.*cbl::par::pi*rad*rad*rho;
						 };
  
  double a_meantau = 0, b_meantau = 0, inv_mis_scale = 1./cluster.sigma_off()/cluster.sigma_off(); 
  std::function<double(double)> tau_factor1 = [&a_meantau,&b_meantau,&inv_mis_scale] (double s){
						return (exp(-0.5*(s*s+a_meantau)*inv_mis_scale))/sqrt(b_meantau-a_meantau-s*s);
					      };
  std::function<double(double)> tau_factor2 = [&a_meantau,&b_meantau,&inv_mis_scale] (double s){
						return (exp(-0.5*(b_meantau-s*s)*inv_mis_scale))/sqrt(b_meantau-a_meantau-s*s);
					      };


  // Define the variable tau_for_misc
  double tau_for_misc = 0;
  
  std::function<double(double)> meantau_factor = [&tau_for_misc,&cluster,&a_meantau,&b_meantau,&inv_mis_scale,&tau_factor1,&tau_factor2] (double rad){
						   a_meantau = (tau_for_misc-rad)*(tau_for_misc-rad);
						   b_meantau = (tau_for_misc+rad)*(tau_for_misc+rad);
						   double c_meantau = (a_meantau+b_meantau)*0.5;
						   double upper_limit = std::min(sqrt(c_meantau-a_meantau),sqrt(25*cluster.sigma_off()*cluster.sigma_off()-a_meantau));
						   double meantau_factor = rad*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor1, 0, upper_limit);
						   if (c_meantau < 25*cluster.sigma_off()*cluster.sigma_off())
						     {
						       double int_left = std::max(0.,sqrt(b_meantau-25*cluster.sigma_off()*cluster.sigma_off()));
						       double int_right = sqrt(b_meantau-c_meantau);						     
						       meantau_factor += rad*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor2, int_left, int_right);
						     }
						   return meantau_factor;
						 };
  
  // Calculate rho_s given r_encl
  double mass_int = cbl::wrapper::gsl::GSL_integrate_cquad(integrate_mass,0.,r_encl);
  rho_s = cluster.mass()/mass_int; // scale density (=M0/(4 pi rs^3))

  std::function<double(double)> sigma_cen_fc = [&rscale,&tau,&rho_s,&F,&L] (double radius){
						 double x = radius/rscale;
						 double tausq = tau*tau; double tausq_xsq = tausq+x*x; double tau4th = tausq*tausq;
						 double prefact = rho_s*rscale*tau4th/(tausq+1.)/(tausq+1.)/(tausq+1.);
						 double a = 2*(tausq+1.)/(x*x-1.)*(1-F(x));
						 double b = 8*F(x);
						 double c = (tau4th-1.)/tausq/tausq_xsq;
						 double d = -cbl::par::pi*(4*tausq_xsq+tausq+1.)/pow(tausq_xsq,1.5);
						 double e = (tausq*(tau4th-1.)+(tausq_xsq)*(3*tau4th-6*tausq-1))*L(x)/(tau*tau*tau)/pow(tausq_xsq,1.5);    
						 return prefact*(a+b+c+d+e)/1.e+12; // in (M_sun/h)/((pc/h)^2) = h*M_sun/pc^2
					       };
  std::function<double(double)> sigma_cen_int = [&sigma_cen_fc] (double rad){
						  return sigma_cen_fc(rad)*2.*rad;
						};
  std::function<double(double)> meansigma_off_int = [&tau_for_misc,&radius,&cluster,&sigma_cen_fc,&meantau_factor] (double tau){
						      tau_for_misc = tau;
						      double int_left = std::max(0.,tau-5*cluster.sigma_off());
						      double int_right = std::min(radius,tau+5*cluster.sigma_off());
						      return tau * sigma_cen_fc(tau) * cbl::wrapper::gsl::GSL_integrate_cquad(meantau_factor, int_left, int_right);
						    };
  std::function<double(double)> sigma_off_int = [&radius,&cluster,&a_meantau,&b_meantau,&inv_mis_scale,&tau_factor1,&tau_factor2,&sigma_cen_fc] (double tau){
						  a_meantau = (tau-radius)*(tau-radius);
						  b_meantau = (tau+radius)*(tau+radius);
						  double c_meantau = (a_meantau+b_meantau)*0.5;
						  double upper_limit = std::min(sqrt(c_meantau-a_meantau),sqrt(25*cluster.sigma_off()*cluster.sigma_off()-a_meantau));
						  double sigma_off_int = radius*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor1, 0., upper_limit);
						  if (c_meantau < 25*cluster.sigma_off()*cluster.sigma_off())
						    {
						      double int_left = std::max(0.,sqrt(b_meantau-25*cluster.sigma_off()*cluster.sigma_off()));
						      double int_right = sqrt(b_meantau-c_meantau);						     
						      sigma_off_int += radius*2.*cbl::wrapper::gsl::GSL_integrate_cquad(tau_factor2, int_left, int_right);
						    }
						  return sigma_off_int * tau * sigma_cen_fc(tau);
						};
  
  double sigma_cen = sigma_cen_fc(radius);
  double meansigma_cen = cbl::wrapper::gsl::GSL_integrate_cquad(sigma_cen_int,0.,radius) /radius/radius;

  double deltasigma_cen = meansigma_cen - sigma_cen;
  
  double deltasigma_1h = 0, meansigma_off = 0, integrated_sigma_off = 0, deltasigma_off = 0;
  if (cluster.f_off() <= 1.e-5 || cluster.sigma_off() <= 1.e-5){
    deltasigma_1h = deltasigma_cen; 
  }else{
    meansigma_off = cbl::wrapper::gsl::GSL_integrate_cquad(meansigma_off_int,0.,radius+5.*cluster.sigma_off()) *2./cbl::par::pi/cluster.sigma_off()/cluster.sigma_off()/radius/radius;
    integrated_sigma_off = cbl::wrapper::gsl::GSL_integrate_cquad(sigma_off_int,std::max(0., radius-5*cluster.sigma_off()),radius+5.*cluster.sigma_off()) /radius/cbl::par::pi/cluster.sigma_off()/cluster.sigma_off();
    deltasigma_off = meansigma_off - integrated_sigma_off;
    deltasigma_1h = (1-cluster.f_off())*deltasigma_cen+cluster.f_off()*deltasigma_off;
  }
    
  return deltasigma_1h;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofile::nfw_1halo_allBins (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_Profile_data_model> pp = static_pointer_cast<STR_Profile_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // redefine the cluster object
  cbl::catalogue::Cluster cluster = *pp->cluster;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  // set the cluster parameters
  cluster.set_concentration(parameter[pp->Cpar.size()]);
  cluster.set_mass(pow(pp->logM_base, parameter[pp->Cpar.size()+1])*1.e14);
  cluster.set_f_off(parameter[pp->Cpar.size()+2]);
  cluster.set_sigma_off(parameter[pp->Cpar.size()+3]);

  std::vector<double> density_profile(radius.size());
  
  for (size_t j=0; j<radius.size(); j++)
    density_profile[j] = nfw_1halo(cosmo, cluster, radius[j], pp->redshift, pp->contrast, pp->trunc_fact);

  return density_profile;
}

// ===========================================================================================
