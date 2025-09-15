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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation1D_angular.cpp
 *
 *  @brief Functions to model the angular two-point correlation
 *  function
 *
 *  This file contains the implementation of the functions used to
 *  model the angular two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, massimilia.romanell2@unibo.it
 */

#include "Bias.h"
#include "ModelFunction_TwoPointCorrelation1D_angular.h"
#include "ModelFunction_TwoPointCorrelation.h"
#include "ModelFunction_NumberCounts.h"
#include "ModelFunction_TwoPointCorrelation_multipoles.h"
#include "ModelFunction_PowerSpectrum_Angular.h"
#include "Modelling_PowerSpectrum_Angular.h"
#include "Data1D.h"

using namespace cbl;
using namespace std;


// ===============================================================================================


std::vector<double> cbl::modelling::twopt::z_limits (const double zmin, const double zmax, const double z_abs_err)
{
  const double dz = (z_abs_err>=0) ? 3.*z_abs_err : 0.;

  vector<double> z_limits = { max(zmin-dz, 0.01), zmax+dz};

  return z_limits;
}


// ============================================================================================


double cbl::modelling::twopt::redshift_selection_function (double redshift, std::vector<double> parameter, std::string distribution_type, double zmin, double zmax)
{
  double distribution = 0.;
  if (distribution_type=="polynomial") {
    if (redshift<zmax && redshift>zmin)
      for (size_t i=0; i<parameter.size(); ++i)
	distribution += parameter[i]*pow(redshift, i);
  }
  else if (distribution_type=="skewnormal") {
    double normpdf=1./(parameter[1]*sqrt(2*par::pi))*exp(-pow((redshift-parameter[0]), 2)/(2*pow(parameter[1], 2)));
    double normcdf=0.5*(1+erf((parameter[2]*((redshift - parameter[0])/parameter[1]))/sqrt(2)));
    distribution=2*parameter[3]*normpdf*normcdf;
  }
  else if (distribution_type=="gaussian") {
    double normpdf=1./(parameter[1]*sqrt(2*par::pi))*exp(-pow((redshift-parameter[0]), 2)/(2*pow(parameter[1], 2)));
    distribution=parameter[2]*normpdf;
  }
  else
    ErrorCBL("distribution_type can be gaussian, polynomial or skewnormal", "redshift_selection_function", "ModelFunction_TwoPointCorrelation1D_angular.cpp", glob::ExitCode::_error_);
  
  return distribution;
  
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::wtheta (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
   
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  
  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  double bias = parameter[pp->Cpar.size()];

  //create zz and dn/dz vectors
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);
  //compute and interpolate growth factor and comoving distance
  vector<double> growth_factor, comoving_distance;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));    
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);

  //compute and interpolate dN_dz distribution
  vector<double> redshift_counts;
  for (size_t j=0; j<z_vector.size(); j++)
    redshift_counts.emplace_back(redshift_selection_function(z_vector[j],pp->dN_par, pp->distribution_type, pp->z_min, pp->z_max));   //pp->z_min and pp->z_min are here only to obtain zero outside limits, for polynomial distribution
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);
  
  // intepolate power spectrum
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);

  vector<double> rr, xi0;
  wrapper::fftlog::transform_FFTlog(rr, xi0, 1, kk, Pk0);
  if (pp->redshift_space) {  //include RSD
    double fact = xi_ratio(cosmo->linear_growth_rate((pp->z_max+pp->z_min)*0.5)/bias);
    for (size_t i=0; i<rr.size(); i++)
      xi0[i]*=fact;
  }
  glob::FuncGrid xi0_interp(rr, xi0, pp->interpType);
  
  //compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi0_interp, comoving_distance_interp, interp_DN, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs)/normalization_redshift_selection_function/normalization_redshift_selection_function;
  
  parameter[pp->Cpar.size()+1] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3);     //S8
  parameter[pp->Cpar.size()+2] = cosmo->Omega_baryon()/cosmo->Omega_matter();  //baryonic fraction

  return  w_theta;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::wtheta_MF_theory (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  double bias = parameter[pp->Cpar.size()];
  
  //create zz and dn/dz vectors
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);
  //compute and interpolate growth factor and comoving distance
  vector<double> growth_factor, comoving_distance;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));    
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);

  // intepolate power spectrum
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);

  vector<double> rr, xi0;
  wrapper::fftlog::transform_FFTlog(rr, xi0, 1, kk, Pk0);
  if (pp->redshift_space) {  //include RSD
    double fact = xi_ratio(cosmo->linear_growth_rate((pp->z_max+pp->z_min)*0.5)/bias);
    for (size_t i=0; i<rr.size(); i++)
      xi0[i]*=fact;
  }
  glob::FuncGrid xi0_interp(rr, xi0, pp->interpType);

  // generate mass function
  vector<double> Mass_vector=logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);
  vector<double> redshift_counts;
  double area_rad=pp->area_degrees/180/180*par::pi*par::pi;
  
  for (size_t i=0; i<z_vector.size(); ++i) {
    vector<double> MF = modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max); 
    glob::FuncGrid MF_interp(Mass_vector, MF, pp->interpType, BinType::_logarithmic_);
    
    auto integrand_MF = [&MF_interp, &cosmo, &z_vector, &i, &pp] (const double mass)
    {
      if (pp->z_abs_err>0) {
	//gaussian convolution, within the photometric bin
	double int_P_z=0.5 * ( erf( (pp->z_max - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) - erf( (pp->z_min - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) );
	return int_P_z*MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
      }
      return MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
    };
    redshift_counts.emplace_back(area_rad*wrapper::gsl::GSL_integrate_qag(integrand_MF, pp->Mass_min, pp->Mass_max));
  }
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);
  
  
  //compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi0_interp, comoving_distance_interp, interp_DN, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs)/normalization_redshift_selection_function/normalization_redshift_selection_function;
    
  parameter[pp->Cpar.size()+1] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3);     //S8
  parameter[pp->Cpar.size()+2] = cosmo->Omega_baryon()/cosmo->Omega_matter();  //baryonic fraction
  
  return  w_theta;
  
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::wtheta_bias_MF_theory (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  
  // intepolate power spectrum. For the mass function should be linear
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, false, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  double redshift = Average(pp->cluster_redshift);

  // create zz and dn/dz vectors
  // define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);

  //compute and interpolate comoving distance. Growth factore is 1, since it is already included in Pk
  vector<double> growth_factor, comoving_distance;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(1.);
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));    
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);

  // compute the redshift number counts
  double area_rad = pp->area_degrees/180/180*par::pi*par::pi;
  vector<double> redshift_counts;
  
  //generate mass function
  vector<double> Mass_vector = logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);
  
  for (size_t i=0; i<z_vector.size(); ++i) {
    vector<double> MF = modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max); 
    glob::FuncGrid MF_interp(Mass_vector, MF, pp->interpType, BinType::_logarithmic_);
    
    auto integrand_MF = [&MF_interp, &cosmo, &z_vector, &i, &pp] (const double mass)
    {
      if (pp->z_abs_err>0) {
	//gaussian convolution, within the photometric bin
	double int_P_z=0.5 * ( erf( (pp->z_max - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) - erf( (pp->z_min - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) );
	return int_P_z*MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
      }
      return MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
    };
    redshift_counts.emplace_back(area_rad*wrapper::gsl::GSL_integrate_qag(integrand_MF, pp->Mass_min, pp->Mass_max));

  }
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  // compute the bias
  double bias;
  cosmology::Bias bb(cosmo);
  //bias as average with mass and redshift of the cluster sample
  if (pp->bias_from_sample){
    vector<double> mass_grid = logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);

    bias = bb.bias_eff_mass(pp->mass, mass_grid, pp->cluster_redshift, pp->model_bias, pp->method_Pk, "mean_bias", false, pp->output_root, pp->Delta)[0];
  }   //theoretical bias with MF
  else {
    vector<double> selection_function;
    for (size_t i=0; i<Mass_vector.size(); i++)
      selection_function.emplace_back(1.); 
    glob::FuncGrid interp_SF(Mass_vector, selection_function, pp->interpType, BinType::_linear_);
    const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM (Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);
    double Delta_crit = (pp->isDelta_critical) ? pp->Delta_input : pp->Delta_input*cosmo->Omega_matter(redshift);
    bias=bb.bias_eff_selection_function(interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_SF, Min(pp->mass), Max(pp->mass), {redshift}, pp->model_bias, pp->model_MF, pp->method_Pk, 1., true, pp->output_root, Delta_crit, -1., pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, par::defaultString, true)[0];
  }

  //recompute the power spectrum
  vector<double> Pk = PX.Pk_matter(kk,  pp->method_Pk, pp->NL, redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  glob::FuncGrid pk_interp(kk, Pk, pp->interpType, BinType::_logarithmic_);
  auto pk_interp_ptr = make_shared<glob::FuncGrid>(pk_interp);
    
  //compute 3D correlation function
  vector<double> rr, xi;
  wrapper::fftlog::transform_FFTlog(rr, xi, 1, kk, Pk);
  glob::FuncGrid xi_interp2;

  if (pp->redshift_space) {  //include RSD
    double fact = xi_ratio(cosmo->linear_growth_rate(redshift)/bias);
    for (size_t i=0; i<rr.size(); i++)
      xi[i]*=fact;
    if (pp->include_quadrupole)
      xi_interp2=modelling::twopt::Xil_interp(kk, 2, "dispersion_Gauss", {cosmo->linear_growth_rate(redshift)/bias, 1., 0.}, {pk_interp_ptr});
  }
  glob::FuncGrid xi_interp(rr, xi, pp->interpType);

  //compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi_interp, comoving_distance_interp, interp_DN, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs, xi_interp2)/normalization_redshift_selection_function/normalization_redshift_selection_function;
    
  parameter[pp->Cpar.size()] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3);     //S8
  parameter[pp->Cpar.size()+1] = bias;  //bias

  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);

  return  w_theta;


}


//=======================================================================================================


vector<double> modelling::twopt::wtheta_bias_MF_theory_BAO (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  
  // parameter that control the BAO shape
  double sigmaNL = parameter[pp->Cpar.size()];
  if (pp->NL)
    ErrorCBL("This model works only with linear Pk!", "wtheta_bias_MF_theory_BAO", "ModelFunction_TwoPointCorrelation1D_angular.cpp",glob::ExitCode::_error_);

  //intepolate power spectrum. For the mass function should be linear
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, false, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
    
  double redshift= Average(pp->cluster_redshift);

  //create zz and dn/dz vectors
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);
  //compute and interpolate growth_factor and comoving distance
  vector<double> growth_factor, comoving_distance, ones;
  for (size_t i=0; i<z_vector.size(); ++i){
    ones.emplace_back(1.); 
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));    
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid ones_interp(z_vector, ones, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);


  //compute the redshift number counts
  double area_rad=pp->area_degrees/180/180*par::pi*par::pi;
  vector<double> redshift_counts;
  //generate mass function
  vector<double> Mass_vector=logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);

  for (size_t i=0; i<z_vector.size(); ++i){
  vector<double> MF=modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max); 
  glob::FuncGrid MF_interp(Mass_vector, MF, pp->interpType, BinType::_logarithmic_);
    
    auto integrand_MF = [&MF_interp, &cosmo, &z_vector, &i, &pp] (const double mass)
    {
      if (pp->z_abs_err>0){
	//gaussian convolution, within the photometric bin
	double int_P_z=0.5 * ( erf( (pp->z_max - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) - erf( (pp->z_min - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) );
	return int_P_z*MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
      }
      return MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
    };
    redshift_counts.emplace_back(area_rad*wrapper::gsl::GSL_integrate_qag(integrand_MF, pp->Mass_min, pp->Mass_max));
  }
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  //compute the bias
  double bias;
  cosmology::Bias bb(cosmo);
  //bias as average with mass and redshift of the cluster sample
  if (pp->bias_from_sample){
    vector<double> mass_grid = logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);
    bias = bb.bias_eff_mass(pp->mass, mass_grid, pp->cluster_redshift, pp->model_bias, pp->method_Pk, "mean_bias", true, pp->output_root, pp->Delta)[0];
  }   //theoretical bias with MF
  else{
    std::vector<double> selection_function;
    for (size_t i=0; i<Mass_vector.size(); i++)
      selection_function.emplace_back(1.); 
    glob::FuncGrid interp_SF(Mass_vector, selection_function, pp->interpType, BinType::_linear_);
    const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM (Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);
    
    double Delta_crit = (pp->isDelta_critical) ? pp->Delta_input : pp->Delta_input*cosmo->Omega_matter(redshift);

    bias=bb.bias_eff_selection_function(interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_SF, Min(pp->mass), Max(pp->mass), {redshift}, pp->model_bias, pp->model_MF, pp->method_Pk, 1., true, pp->output_root, Delta_crit, -1., pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, par::defaultString, true)[0];
  }

  parameter[parameter.size()-1] = bias;  //bias

  //compute linear Pk and PkNW (without BAO)
  vector<double> PkNW = PX.Pk_matter(kk, "EisensteinHu", false, redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  
  // compute the 2pcf signal
  vector<double> Pk(kk.size(), 0);
  
  for (size_t i =0; i<kk.size(); i++) {
    Pk[i] = PkNW[i]*(1.+(Pk0[i]*interp_DN(redshift)*interp_DN(redshift)/PkNW[i]-1.)*exp(-0.5*pow(kk[i]*sigmaNL, 2)));
  }
  glob::FuncGrid pk_interp(kk, Pk, pp->interpType, BinType::_logarithmic_);
  auto pk_interp_ptr = make_shared<glob::FuncGrid>(pk_interp);

  // compute the monopole of the two-point correlation function
  vector<double> rr, xi;
  wrapper::fftlog::transform_FFTlog(rr, xi, 1, kk, Pk);
  glob::FuncGrid xi_interp2, xi_interp4;
  //double fact = bias*bias;
  if (pp->redshift_space){
    double fact = xi_ratio(cosmo->linear_growth_rate(redshift)/bias);
    for (size_t i=0; i<rr.size(); i++)
      xi[i]*=fact;
    if (pp->include_quadrupole){  //interpolate the quadrupole and hexadecapole of the two-point correlation function
      xi_interp2=modelling::twopt::Xil_interp(kk, 2, "dispersion_Gauss", {cosmo->linear_growth_rate(redshift)/bias, 1., 0.}, {pk_interp_ptr});
      //xi_interp4=modelling::twopt::Xil_interp(kk, 4, "dispersion_Gauss", {cosmo->linear_growth_rate(redshift)/bias, 1., 0.}, {pk_interp_ptr});
    }
  }
    
  // interpolate the monopole of the two-point correlation function
  glob::FuncGrid xi_interp(rr, xi, pp->interpType);
 
  //compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi_interp, comoving_distance_interp, ones_interp, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs, xi_interp2, xi_interp4)/normalization_redshift_selection_function/normalization_redshift_selection_function;
    
  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);

  return  w_theta;
  
}


//=======================================================================================================


vector<double> modelling::twopt::wtheta_bias_MF_theory_BAO_poly (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  // parameter that control the BAO shape
  double sigmaNL = parameter[pp->Cpar.size()];

  if (pp->NL)
    ErrorCBL("This model works only with linear Pk!", "wtheta_bias_MF_theory_BAO", "ModelFunction_TwoPointCorrelation1D_angular.cpp",glob::ExitCode::_error_);

  //intepolate power spectrum. For the mass function should be linear
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, false, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  double redshift= Average(pp->cluster_redshift);

  
  //create zz and dn/dz vectors
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);
  //compute and interpolate growth_factor and comoving distance
  vector<double> growth_factor, comoving_distance, ones;
  for (size_t i=0; i<z_vector.size(); ++i){
    ones.emplace_back(1.); 
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));    
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid ones_interp(z_vector, ones, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);


  //compute the redshift number counts
  double area_rad=pp->area_degrees/180/180*par::pi*par::pi;
  vector<double> redshift_counts;
  //generate mass function
  vector<double> Mass_vector=logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);

  for (size_t i=0; i<z_vector.size(); ++i){
  vector<double> MF=modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max); 
  glob::FuncGrid MF_interp(Mass_vector, MF, pp->interpType, BinType::_logarithmic_);
    
    auto integrand_MF = [&MF_interp, &cosmo, &z_vector, &i, &pp] (const double mass)
    {
      if (pp->z_abs_err>0){
	//gaussian convolution, within the photometric bin
	double int_P_z=0.5 * ( erf( (pp->z_max - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) - erf( (pp->z_min - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) );
	return int_P_z*MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
      }
      return MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
    };
    redshift_counts.emplace_back(area_rad*wrapper::gsl::GSL_integrate_qag(integrand_MF, pp->Mass_min, pp->Mass_max));
  }
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  //compute the bias
  double bias;
  cosmology::Bias bb(cosmo);
  //bias as average with mass and redshift of the cluster sample
  if (pp->bias_from_sample){
    vector<double> mass_grid = logarithmic_bin_vector(pp->mass_step, pp->Mass_min, pp->Mass_max);
    bias = bb.bias_eff_mass(pp->mass, mass_grid, pp->cluster_redshift, pp->model_bias, pp->method_Pk, "mean_bias", true, pp->output_root, pp->Delta)[0];
  }   //theoretical bias with MF
  else{
    vector<double> selection_function;
    for (size_t i=0; i<Mass_vector.size(); i++)
      selection_function.emplace_back(1.); 
    glob::FuncGrid interp_SF(Mass_vector, selection_function, pp->interpType, BinType::_linear_);
    
    //vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
    //vector<double> Pk0=cosmo.Pk_matter(kk, pp->method_Pk, false, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
    const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM (Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);

    double Delta_crit = (pp->isDelta_critical) ? pp->Delta_input : pp->Delta_input*cosmo->Omega_matter(redshift);

    bias=bb.bias_eff_selection_function(interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_SF, Min(pp->mass), Max(pp->mass), {redshift}, pp->model_bias, pp->model_MF, pp->method_Pk, 1., true, pp->output_root, Delta_crit, -1., pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, par::defaultString, true)[0];
  }

  //compute linear Pk and PkNW (without BAO)
  vector<double> PkNW = PX.Pk_matter(kk, "EisensteinHu", false, redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  
  // compute the 2pcf signal
  vector<double> Pk(kk.size(), 0);
  
  for (size_t i =0; i<kk.size(); i++) {
    Pk[i] = PkNW[i]*(1.+(Pk0[i]*interp_DN(redshift)*interp_DN(redshift)/PkNW[i]-1.)*exp(-0.5*pow(kk[i]*sigmaNL, 2)));
  }
  vector<double> rr, xi;
  wrapper::fftlog::transform_FFTlog(rr, xi, 1, kk, Pk);

  // compute the monopole of the two-point correlation function
  double fact = bias*bias;
  if (pp->redshift_space)
    fact*=xi_ratio(cosmo->linear_growth_rate(redshift)/bias);
    
  for (size_t i =0; i<xi.size(); i++) {
    double poly = 0;  //A0+A1/r+A2/r/r
    for (int j = 0; j<3; j++){
      poly += parameter[pp->Cpar.size()+1+j]*pow(rr[i], -j);
    }
    xi[i] = fact*xi[i]+poly;
  }
  
  // return the monopole of the two-point correlation function

  glob::FuncGrid xi_interp(rr, xi, pp->interpType);
 
  //compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi_interp, comoving_distance_interp, ones_interp, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs)/normalization_redshift_selection_function/normalization_redshift_selection_function;
    
  parameter[pp->Cpar.size()+4] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3);     //S8
  parameter[pp->Cpar.size()+5] = bias;  //bias

  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);

  return  w_theta;
  
}


//=======================================================================================================


std::vector<double> cbl::modelling::twopt::wtheta_scaling_relation_cosmology (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  // scaling relation parameters
  vector<double> scalRel_pars;
  for (size_t i=0; i<parameter.size()-1; i++)
    scalRel_pars.emplace_back(parameter[i]);
  
  const double alpha = scalRel_pars[scalRel_pars.size()-8];
  const double beta = scalRel_pars[scalRel_pars.size()-7];
  const double gamma = scalRel_pars[scalRel_pars.size()-6];
  const double scatter0 = scalRel_pars[scalRel_pars.size()-5];
  const double scatterM = scalRel_pars[scalRel_pars.size()-4];
  const double scatterM_exp = scalRel_pars[scalRel_pars.size()-3];
  const double scatterz = scalRel_pars[scalRel_pars.size()-2];
  const double scatterz_exp = scalRel_pars[scalRel_pars.size()-1];

  double log_base = (pp->scaling_relation)->data_model().log_base;
  double mass_pivot = (pp->scaling_relation)->data_model().mass_pivot;
  double proxy_pivot = (pp->scaling_relation)->data_model().proxy_pivot;
  double redshift_pivot = (pp->scaling_relation)->data_model().redshift_pivot;

  // Interpolate the normalised amplitude of the growing mode and the comoving distance
  const vector<double> z_for_interp = linear_bin_vector(30, 0.0001, Max((pp->scaling_relation)->data_model().redshift)+5.*pp->z_abs_err); 
  vector<double> DN (z_for_interp.size(), 0.),  comoving_distance (z_for_interp.size(), 0.);
  for (size_t i=0; i<z_for_interp.size(); i++) {
    DN[i] = cosmo->DN(z_for_interp[i]);
    comoving_distance[i]=cosmo->D_C(z_for_interp[i]);    
  }
  
  glob::FuncGrid interp_DN(z_for_interp, DN, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_for_interp, comoving_distance, pp->interpType, BinType::_linear_);

  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the effective bias from the scaling relation
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  
  double bias = bias_from_scaling_relation(alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, log_base, mass_pivot, proxy_pivot, redshift_pivot, scalRel_pars, interp_DN, inputs);
  
  // set the value of the bias
  parameter[parameter.size()-1] = bias;

  //create zz and dn/dz vectors
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);
  //compute and interpolate dN_dz distribution
  vector<double> redshift_counts;
  for (size_t j=0; j<z_vector.size(); j++)
    redshift_counts.emplace_back(redshift_selection_function(z_vector[j],pp->dN_par, pp->distribution_type, pp->z_min, pp->z_max));   //pp->z_min and pp->z_min are here only to obtain zero outside limits, for polynomial distribution
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  // intepolate power spectrum
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);

  vector<double> rr, xi0;
  wrapper::fftlog::transform_FFTlog(rr, xi0, 1, kk, Pk0);
  if (pp->redshift_space) {  //include RSD
    double fact = xi_ratio(cosmo->linear_growth_rate((pp->z_max+pp->z_min)*0.5)/bias);
    for (size_t i=0; i<rr.size(); i++)
      xi0[i]*=fact;
  }
  glob::FuncGrid xi0_interp(rr, xi0, pp->interpType);
  
  //compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi0_interp, comoving_distance_interp, interp_DN, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs)/normalization_redshift_selection_function/normalization_redshift_selection_function;  
  
  return  w_theta;
}


//=========================================================================================


std::vector<double> cbl::modelling::twopt::wtheta_scaling_relation_cosmology_MF_theory (std::vector<double> theta, std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  // scaling relation parameters
  vector<double> scalRel_pars;
  for (size_t i=0; i<parameter.size()-4; i++)
    scalRel_pars.emplace_back(parameter[i]);
  
  const double alpha = scalRel_pars[pp->Cpar.size()];
  const double beta = scalRel_pars[pp->Cpar.size()+1];
  const double gamma = scalRel_pars[pp->Cpar.size()+2];
  const double scatter0 = scalRel_pars[pp->Cpar.size()+3];
  const double scatterM = scalRel_pars[pp->Cpar.size()+4];
  const double scatterM_exp = scalRel_pars[pp->Cpar.size()+5];
  const double scatterz = scalRel_pars[pp->Cpar.size()+6];
  const double scatterz_exp = scalRel_pars[pp->Cpar.size()+7];
  const double Plambda_a = parameter[pp->Cpar.size()+8];
  const double Plambda_b = parameter[pp->Cpar.size()+9];
  const double Plambda_c = parameter[pp->Cpar.size()+10];

  double log_base = (pp->scaling_relation)->data_model().log_base;
  double mass_pivot = (pp->scaling_relation)->data_model().mass_pivot;
  double proxy_pivot = (pp->scaling_relation)->data_model().proxy_pivot;
  double redshift_pivot = (pp->scaling_relation)->data_model().redshift_pivot;
   
  // Interpolate the normalised amplitude of the growing mode, the linear growth rate and comoving distance
  const vector<double> z_for_interp = linear_bin_vector(50, 0.0001, Max((pp->scaling_relation)->data_model().redshift)+5.*pp->z_abs_err);
  
  vector<double> DN (z_for_interp.size(), 0.), comoving_distance (z_for_interp.size(), 0.);
  for (size_t i=0; i<z_for_interp.size(); i++) {
    DN[i] = cosmo->DN(z_for_interp[i]);
    comoving_distance[i] = cosmo->D_C(z_for_interp[i]);    
  }
  
  glob::FuncGrid interp_DN(z_for_interp, DN, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_for_interp, comoving_distance, pp->interpType, BinType::_linear_);

  
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the effective bias from the scaling relation
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»


  double bias = bias_from_scaling_relation(alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, log_base, mass_pivot, proxy_pivot, redshift_pivot, scalRel_pars, interp_DN, inputs);
  // set the value of the bias
  parameter[parameter.size()-1] = bias;

  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Compute the angular correlation function
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  // define Pk
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  
  // interpolate sigmaM and its derivative
  const vector<double> Mass_vector = logarithmic_bin_vector(300, 1.e10, 1.e16);
  const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM(Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);

  //interpolate weights vector. Log binning
  double log_bin_size = (log10(pp->proxy_max)-log10(pp->proxy_min))/(pp->SF_weights.size());
  vector<double> lambda_star_vector=logarithmic_bin_vector(pp->SF_weights.size(), pow(10, log10(pp->proxy_min)+0.5*log_bin_size), pow(10, log10(pp->proxy_max)-0.5*log_bin_size));
  glob::FuncGrid weights_interp(lambda_star_vector, pp->SF_weights, pp->interpType, BinType::_logarithmic_);

  // create zz and dn/dz vectors
  // define z_vector, larger than photometric bin limits
  vector<double> z_limit = z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector = linear_bin_vector(12, z_limit[0], z_limit[1]);
  
  // interpolate and normalize redshift counts
  vector<double> redshift_counts;
  double area_rad = pp->area_degrees/180/180*par::pi*par::pi;
  
  for (size_t i=0; i<z_vector.size(); ++i)
    redshift_counts.emplace_back(redshift_distribution_proxy(z_vector[i], alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, pp->z_abs_err, pp->proxy_rel_err, Plambda_a, Plambda_b, Plambda_c, pp->fz, pp->z_min, pp->z_max, pp->proxy_min, pp->proxy_max, cosmo, area_rad, pp->model_MF, pp->Delta_input, pp->isDelta_critical, interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_DN, proxy_pivot, redshift_pivot, mass_pivot, log_base, weights_interp));

  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);
  vector<double> rr, xi0;
  wrapper::fftlog::transform_FFTlog(rr, xi0, 1, kk, Pk0);
  if (pp->redshift_space) { // include RSD
    double fact = xi_ratio(cosmo->linear_growth_rate((pp->z_max+pp->z_min)*0.5)/bias);
    for (size_t i=0; i<rr.size(); i++)
      xi0[i]*=fact;
  }
  
  glob::FuncGrid xi0_interp(rr, xi0, pp->interpType);
  
  // compute angular correlation function
  vector<double> w_theta(theta.size(), 0);
  for (size_t i=0; i<theta.size(); i++)
    w_theta[i] = integral_wtheta(theta[i], pp->theta_units, xi0_interp, comoving_distance_interp, interp_DN, redshift_counts_interp, bias, redshift_counts_interp, bias, true, parameter, inputs)/normalization_redshift_selection_function/normalization_redshift_selection_function;  

  return w_theta;
}


//==========================================================================================================================================


double cbl::modelling::twopt::redshift_distribution_proxy (const double redshift, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_abs_err, const double proxy_rel_err, const double Plambda_a, const double Plambda_b, const double Plambda_c, ProxyFunction fz, const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, std::shared_ptr<cosmology::Cosmology> cosmology, const double Area, const std::string model_MF, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base, cbl::glob::FuncGrid weights_interp)
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
  auto integrand = [&] (const vector<double> x)
  {
    double Delta_ = (isDelta_critical) ? Delta/cosmology_ptr->Omega_matter(redshift) : Delta;
    double Mass = pow(log_base,x[0])*mass_pivot;
    normM = x[0];
    the_redsh = redshift;

    // Compute P(M|lambda,z)
    double log_lambda = log(x[1]/proxy_pivot)/log(log_base);
    double log_f_z = log( fz(redshift, z_pivot, cosmology_ptr) )/log(log_base);
      
    double mean = alpha + beta*log_lambda + gamma*log_f_z;
    double sigma = abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
    double P_M__lambda_z = gaussian(normM, pp, {mean,sigma});
      
    // Compute P(lambda|z)
    double P_lambda__z = Plambda_a * pow(x[1],-Plambda_b) * exp(-Plambda_c*x[1]);
    
    // Compute P(M|z)
    double P_M__z = 0.;
    if (P_M__lambda_z*P_lambda__z > 0)
      P_M__z = wrapper::gsl::GSL_integrate_cquad(integrand_P_M__z,log(0.0001)/log(log_base),log(100.)/log(log_base)) * log(log_base);
    else
      P_M__z = 1;

    // Compute the integrals of P(z|z) and P(lambda|lambda)
    double mean_Pz = redshift;
    double int_P_z = 0.5 * ( erf( (redshift_max - mean_Pz) / (sqrt(2)*z_abs_err) ) - erf( (redshift_min - mean_Pz) / (sqrt(2)*z_abs_err) ) );
    double mean_Plambda = x[1];
    double int_P_lambda = 0.5 * ( erf( (proxy_max - mean_Plambda) / (sqrt(2)*proxy_rel_err*proxy_max) ) - erf( (proxy_min - mean_Plambda) / (sqrt(2)*proxy_rel_err*proxy_min) ) );
    double weights_int = weights_interp(x[1]);
    if (weights_int>1) weights_int = 1.;
    else if (weights_int<0) weights_int = 0.;
    cosmology::MassFunction MF(cosmology_ptr);
    return MF.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), redshift, interp_DN(redshift), model_MF, false, par::defaultString, Delta_) * pow(log_base,normM) * (P_M__lambda_z*P_lambda__z/P_M__z) * int_P_z * int_P_lambda*weights_int;
  };
  
  // -------------------------------------------------------------

  // Find the minimum and maximum masses, given the parameters of the scaling relation
  double log_lambda_min = log((max(proxy_min - 3.5*proxy_rel_err*proxy_min, 1.))/proxy_pivot)/log(log_base);
  double log_lambda_max = log((proxy_max + 3.5*proxy_rel_err)/proxy_pivot)/log(log_base);

  double log_f_z_min = log( fz((max(redshift_min - 3.5*z_abs_err, 0.)), z_pivot, cosmology_ptr) )/log(log_base);
  double log_f_z_max = log( fz((redshift_max + 3.5*z_abs_err), z_pivot, cosmology_ptr) )/log(log_base);

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
  int integral_dimension=2;
  vector<vector<double>> integration_limits(integral_dimension);
  integration_limits[0] = {max(minM-3.5*max_intrinsic_scatter,log(1.e10/mass_pivot)/log(log_base)), min(maxM+3.5*max_intrinsic_scatter,log(1.e16/mass_pivot)/log(log_base))};
  integration_limits[1] = {max(proxy_min - 2.*proxy_rel_err*proxy_min, 0.00001), proxy_max + 2.*proxy_rel_err*proxy_max};

  // Compute the integral
  wrapper::cuba::CUBAwrapper CW (integrand, integral_dimension);
  double dN_dz;
 
  if (integration_limits[0][0] < integration_limits[0][1])
    dN_dz = CW.IntegrateVegas(integration_limits, false)*Area*cosmology_ptr->dV_dZdOmega(redshift, true);
  else
    dN_dz = 0;
  return dN_dz * mass_pivot * log(log_base);
}


// =====================================================================================


double cbl::modelling::twopt::integral_wtheta (const double theta, const CoordinateUnits coordUnits, cbl::glob::FuncGrid xi0_interp, cbl::glob::FuncGrid comoving_distance_interp, cbl::glob::FuncGrid interp_DN, cbl::glob::FuncGrid redshift_counts_interp, const double bias, cbl::glob::FuncGrid redshift_counts_interp2, const double bias2, bool GSL, std::vector<double> par, std::shared_ptr<void> inputs, cbl::glob::FuncGrid xi2_interp, cbl::glob::FuncGrid xi4_interp)
{
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  
  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], par[i]);
  
  //set integration limits, more extended than photometric redshift bin
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);

  // relative measurement precision (it raises a warning message if not satisfied)
  const double prec = 1.e-1;
  const double abs_err = 5.;

  // default values (refer to the GSL documentation for different choiches)
  const int limit_size = 1000;
  const int rule = 6;
  
  double theta_rad = converted_angle (theta, coordUnits, CoordinateUnits::_radians_);

  double fact=1.; 
  if (GSL) {
    auto integrand = [&pp, &cosmo, &xi0_interp, &xi2_interp, &xi4_interp, &comoving_distance_interp, &interp_DN, &redshift_counts_interp, &redshift_counts_interp2, &theta_rad, &fact, &z_limit, &limit_size, &abs_err, &prec, &rule] (double z1)
    {
      double r1 = comoving_distance_interp(z1);
      auto integrand_z2 = [&pp, &cosmo, &xi0_interp, &xi2_interp, &xi4_interp, &comoving_distance_interp, &interp_DN, &redshift_counts_interp2, &r1, &theta_rad, &fact] (double z2) {
	double r2 = comoving_distance_interp(z2);
	double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
	if (pp->redshift_space && pp->include_quadrupole) {
	  double mu=(r2-r1)/ss*cos(theta_rad/2);
	  return interp_DN(z2)*(xi0_interp(ss)-P_2(mu)*xi2_interp(ss))*fact*redshift_counts_interp2(z2);
	  //return interp_DN(z2)*(xi0_interp(ss)-P_2(mu)*xi2_interp(ss)+P_4(mu)*xi4_interp(ss))*fact*redshift_counts_interp2(z2);

	}
	return interp_DN(z2)*xi0_interp(ss)*fact*redshift_counts_interp2(z2);
      };
		       
      return wrapper::gsl::GSL_integrate_cquad(integrand_z2, z_limit[0], z_limit[1], prec, abs_err)*redshift_counts_interp(z1)*interp_DN(z1);
    };
    return bias*bias2*wrapper::gsl::GSL_integrate_cquad(integrand, z_limit[0], z_limit[1], prec, abs_err);
  }
  else{
    auto integrand = [&] (vector<double> zz)
    {
      double r1 = comoving_distance_interp(zz[0]);
      double r2 = comoving_distance_interp(zz[1]);
      double ss = sqrt(pow(r1,2)+pow(r2,2)-2*r1*r2*cos(theta_rad));
      if (pp->redshift_space && pp->include_quadrupole) {
	double mu=(r2-r1)/ss*cos(theta_rad/2);
	return interp_DN(zz[1])*interp_DN(zz[0])*(xi0_interp(ss)-P_2(mu)*xi2_interp(ss))*redshift_counts_interp(zz[0])*redshift_counts_interp2(zz[1]);
      }
      return interp_DN(zz[1])*interp_DN(zz[0])*xi0_interp(ss)*redshift_counts_interp(zz[0])*redshift_counts_interp2(zz[1]);
    };

    wrapper::cuba::CUBAwrapper integrator(integrand, 2);

    return bias*bias2*integrator.IntegrateCuhre( {{z_limit[0], z_limit[1]}, {z_limit[0], z_limit[1]}});
  }
}


//=====================================================================================


std::vector<double> cbl::modelling::twopt::covariance_wtheta_bias_MF_theory(std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i = 0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
 
  //convert coordinates to compute covariance
  double fact=1.;
  if(pp->theta_units==cbl::CoordinateUnits::_degrees_)
    fact = par::pi/180; 
  if(pp->theta_units==cbl::CoordinateUnits::_arcminutes_)
    fact = par::pi/180/60; 
  if(pp->theta_units==cbl::CoordinateUnits::_arcseconds_)
    fact = par::pi/180/60/60; 

  for(size_t i = 0; i < theta.size(); i++)
    theta[i]*=fact; 

  double beta = parameter[pp->Cpar.size()]; 
  double alpha = parameter[pp->Cpar.size()+1];

  //compute cl, removing extra parameters (alpha, beta). Add shot noise (fixed to zero)
  std::vector<double> l = linear_bin_vector(200, 1., 200.);
  std::vector<double> parameter_cl = parameter;
  parameter_cl.pop_back();
  parameter_cl.pop_back();
  parameter_cl.emplace_back(0.);

  // set the cosmological parameters
  std::vector<double> Cl = (pp->angularpk)->likelihood()->get_m_model()->operator()(l, {parameter_cl});
  double area_rad=pp->area_degrees/180/180*par::pi*par::pi;
  double nObjects = pp->cluster_redshift.size();
  double shot_noise = area_rad/nObjects;
  double fsky = area_rad/4/par::pi;

  std::vector<double> diagonal_w_error(theta.size(), 0.0);
  std::vector<std::vector<double>> limits = bin_limits(theta);
  std::vector<double> lower_limits=limits[0];
  std::vector<double> upper_limits=limits[1]; 
   
  for (size_t t = 0; t < theta.size(); ++t) {
    diagonal_w_error[t] = 0.25 / (par::pi*par::pi) * std::pow(shot_noise * (1 + alpha), 2) / fsky / (std::cos(lower_limits[t]) - std::cos(upper_limits[t]));
    for (size_t i = 0; i < l.size(); ++i) {  
      double Pl = cbl::Legendre_bin_averaged(lower_limits[t], upper_limits[t], l[i]);
      diagonal_w_error[t] += (0.125 * (2 * l[i] + 1) / (par::pi*par::pi*fsky)* Pl * Pl* (beta * beta * Cl[i] * Cl[i] + 2 * beta * Cl[i] * shot_noise * (1 + alpha))); 
    }
  }
    
  return diagonal_w_error;
}


