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
 *  Modelling/AngularPowerSpectrum/ModelFunction_PowerSpectrum_Angular.cpp
 *
 *  @brief Functions to model the angular power spectrum
 *
 *  This file contains the implementation of the functions used to
 *  model the angular power spectrum
 *
 *  @authors Federico Marulli, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */


#include "ModelFunction_PowerSpectrum_Angular.h"

using namespace std;
using namespace cbl;


// ===============================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_mixed(std::vector<double> l_mixing, std::vector<std::vector<double>> mixing_matrix, std::vector<double> l, std::vector<double> Cl, double fsky){

  if(l_mixing[0]!=0)  cbl::ErrorCBL("the mixing matrix should start from l=0", "Cl_mixed", "ModelFunction_PowerSpectrum_Angular.cpp",cbl::glob::ExitCode::_error_); 
  
  double _fsky=1./fsky;
  std::vector<double> Cl_mixed(l.size(), 0);

  for (size_t j=0; j<l.size(); ++j)
    for (size_t i=0; i<l.size(); ++i)
      Cl_mixed[j]+=mixing_matrix[j+l[0]][i+l[0]]*Cl[i]*_fsky;
  
  return Cl_mixed;
  
}


// ============================================================================================


double cbl::modelling::angularpk::integral_limber_interp (const double l, std::vector<double> z_vector, std::vector<double> kk, cbl::glob::FuncGrid2D pk_interp, std::vector<double> par, std::shared_ptr<void> inputs){
    
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  //set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], par[i]);
  
  double lower_limit;
  double upper_limit;
  if(pp->z_min_bin2>0 && pp->z_max_bin2>0){
    lower_limit = min(pp->z_min, pp->z_min_bin2);
    upper_limit = max(pp->z_max, pp->z_max_bin2);
  }
  else{
    lower_limit = pp->z_min;
    upper_limit = pp->z_max;
  }

  // relative measurement precision (it raises a warning message if not satisfied)
  const double prec = 1.e-5;
  
  // default values (refer to the GSL documentation for different choiches)
  const int limit_size = 1000;
  const int rule = 6;
  
  std::function<double(double)> integrand_limber = [&l, &par, &pk_interp, &cosmo, &inputs, &z_vector, &kk] (double redshift)
  {
    shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
        
    double offset=0;
    double slope=0;
    if(par.size()>pp->Cpar.size()+4 && pp->dN_par_bin2.size()==0){ //Cpar.size()+2 because the first element is l
      offset = par[pp->Cpar.size()+1];
      slope = par[pp->Cpar.size()+2];
    }
    
    double distribution=0;
    double distribution_bin2=0;
    for (size_t i=0; i<pp->dN_par.size(); ++i){  //if dN_par.size=0 it does not enter into this cycle (offset and slope are defined outside)
      if(redshift<pp->z_max && redshift>pp->z_min)
	distribution += pp->dN_par[i]*pow(redshift, i);
      else distribution=0;
    }
    distribution+=par[pp->Cpar.size()]+par[pp->Cpar.size()+1]*redshift;
    if(pp->dN_par_bin2.size()>0){
      for (size_t i=0; i<pp->dN_par_bin2.size(); ++i)
	if(redshift<pp->z_max_bin2 && redshift>pp->z_min_bin2)
	  distribution_bin2 += pp->dN_par_bin2[i]*pow(redshift, i);
	else distribution_bin2=0;
      distribution_bin2+=offset+slope*redshift;
    }

    std::vector<double> pk_interpolated;
    std::vector<std::vector<double>> Pk_interp;
    double kk_exact=double(l+0.5)/cosmo.D_C(redshift);
    size_t index_k=10000;
    double kk_diff=10000;
    for(size_t i=0; i<kk.size(); ++i)
      if(abs(kk[i]-kk_exact)<kk_diff){
	kk_diff=abs(kk[i]-kk_exact);
	index_k=i;
      }
    if(index_k!=0) kk[index_k]=kk_exact;
    
    std::vector<std::vector<double>> interp_matrix={{redshift}, {kk[index_k]}};
    interp_matrix=cbl::transpose(interp_matrix);
    pk_interpolated=pk_interp.eval_func(interp_matrix);
    
    double _cc=1/cbl::par::cc;
    double inv_d2= 1/(cosmo.D_C(redshift)*cosmo.D_C(redshift));
    if(pp->dN_par_bin2.size()>0) return distribution*distribution_bin2*pk_interpolated[0]*inv_d2*cosmo.HH(redshift)*_cc;    
    return distribution*distribution*pk_interpolated[0]*inv_d2*cosmo.HH(redshift)*_cc;
    
  };
  
  std::vector<double> parameter_integrand;
  parameter_integrand.emplace_back(l);
  for (size_t i=0; i<par.size(); ++i) parameter_integrand.emplace_back(par[i]);  
  return cbl::wrapper::gsl::GSL_integrate_qag(integrand_limber, lower_limit, upper_limit, prec, limit_size, rule);

}


// ============================================================================================


double cbl::modelling::angularpk::integrand_limber_exact(double redshift, std::shared_ptr<void> inputs, std::vector<double> par){
  
  // structure contaning the required input data

  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;
  
  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], par[i+1]);  //i+1 because i=0 is l
  double offset=0;
  double slope=0;
  if(par.size()>pp->Cpar.size()+4 && pp->dN_par_bin2.size()==0){ //Cpar.size()+2 because the first element is l
    offset = par[pp->Cpar.size()+2];
    slope = par[pp->Cpar.size()+3];
  }
  
  // redefine the cosmology
  double distribution=0;
  double distribution_bin2=0;
  for (size_t i=0; i<pp->dN_par.size(); ++i){  //if dN_par.size=0 it does not enter into this cycle (offset and slope are definet outside)
    if(redshift<pp->z_max && redshift>pp->z_min){
      distribution += pp->dN_par[i]*pow(redshift, i);
    }
    else distribution=0;
  }
  distribution+=par[pp->Cpar.size()+1]+par[pp->Cpar.size()+2]*redshift;
  if(pp->dN_par_bin2.size()>0){
    for (size_t i=0; i<pp->dN_par_bin2.size(); ++i)
      if(redshift<pp->z_max_bin2 && redshift>pp->z_min_bin2){
	distribution_bin2 += pp->dN_par_bin2[i]*pow(redshift, i);
      }
      else distribution_bin2=0;
    distribution_bin2+=offset+slope*redshift;
  }
  
  double l=par[0];
  double _cc=1/cbl::par::cc;
  double Pk=cosmo.Pk_matter({double((l+0.5)/cosmo.D_C(redshift))}, pp->method_Pk, pp->NL, redshift, false, "test", pp->norm, pp->k_min, pp->k_max)[0];

  double inv_d2= 1/(cosmo.D_C(redshift)*cosmo.D_C(redshift));
  if(pp->dN_par_bin2.size()>0)   return distribution*distribution_bin2*Pk*inv_d2*cosmo.HH(redshift)*_cc;
 
  return distribution*distribution*Pk*inv_d2*cosmo.HH(redshift)*_cc;
  
}


// =====================================================================


double cbl::modelling::angularpk::integral_limber_exact (const double l, std::vector<double> parameter, std::shared_ptr<void> inputs){

  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  double lower_limit;
  double upper_limit;
  if(pp->z_min_bin2>0 && pp->z_max_bin2>0){
    lower_limit = min(pp->z_min, pp->z_min_bin2);
    upper_limit = max(pp->z_max, pp->z_max_bin2);
  }
  else{
    lower_limit = pp->z_min;
    upper_limit = pp->z_max;
  }

  // relative measurement precision (it raises a warning message if not satisfied)
  const double prec = 1.e-5;
  
  // default values (refer to the GSL documentation for different choiches)
  const int limit_size = 1000;
  const int rule = 6;

  std::vector<double> parameter_integrand;
  parameter_integrand.emplace_back(l);
  for (size_t i=0; i<parameter.size(); ++i) parameter_integrand.emplace_back(parameter[i]);
  return cbl::wrapper::gsl::GSL_integrate_qag(&integrand_limber_exact, inputs, parameter_integrand, lower_limit, upper_limit, prec, limit_size, rule);
   
}


// ===========================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_limber(const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);
  
  vector<double> Cl(l.size(), 0);
    
  // input parameters
  double bias = parameter[pp->Cpar.size()];
  
  if(pp->interpolate_power_spectrum==false)
    for (size_t i=0; i<l.size(); i++) Cl[i] = cbl::modelling::angularpk::integral_limber_exact(l[i], parameter, inputs);
  
  else{
    std::vector<double> kk, z_vector;
    double kk_step=0.0002, z_step=0.03;
    double kk_min=double((*std::min_element(std::begin(l), std::end(l))+0.5)/cosmo.D_C(pp->z_max)), kk_max=double((*max_element(std::begin(l), std::end(l))+0.5)/cosmo.D_C(pp->z_min));
    for(size_t i=0; i<(pp->z_max-pp->z_min)/z_step +1; ++i) z_vector.emplace_back(pp->z_min+i*z_step);
    for(size_t i=0; i<(kk_max-kk_min)/kk_step+1; ++i) kk.emplace_back(kk_min+i*kk_step);
    
    std::vector<std::vector<double>> Pk=cosmo.Pk_matter(kk, pp->method_Pk, pp->NL, z_vector, false, "test", pp->norm, pp->k_min, pp->k_max); 
    cbl::glob::FuncGrid2D pk_interp({z_vector},{kk}, Pk, "Linear");
    
    for (size_t i=0; i<l.size(); i++)
      Cl[i] = cbl::modelling::angularpk::integral_limber_interp(l[i], z_vector, kk, pk_interp, parameter, inputs);
  }

  if(pp->mixing_matrix.size()!=0)  Cl=Cl_mixed(pp->ll, pp->mixing_matrix, l, Cl, pp->fsky);
  
  if(pp->dN_par_bin2.size()>0){
    double bias_bin2 = parameter[pp->Cpar.size()+1];
    for (size_t i =0; i<Cl.size(); i++) Cl[i] *= bias*bias_bin2;
    parameter[pp->Cpar.size()+2] = cosmo.sigma8()*sqrt(cosmo.Omega_matter()/0.3);     //S8
    parameter[pp->Cpar.size()+3] = cosmo.Omega_baryon()/cosmo.Omega_matter();  //baryonic fraction
  }
  else {
    for (size_t i =0; i<Cl.size(); i++) Cl[i] *= bias*bias;
    parameter[pp->Cpar.size()+1] = cosmo.sigma8()*sqrt(cosmo.Omega_matter()/0.3);     //S8
    parameter[pp->Cpar.size()+2] = cosmo.Omega_baryon()/cosmo.Omega_matter();  //baryonic fraction
  }
  return Cl;
  
}


// ============================================================================================

