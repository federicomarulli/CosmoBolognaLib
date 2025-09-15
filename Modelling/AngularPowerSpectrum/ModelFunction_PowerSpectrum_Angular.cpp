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

#include "Bias.h"
#include "ModelFunction_PowerSpectrum_Angular.h"
#include "Modelling_PowerSpectrum_Angular.h"
#include "ModelFunction_NumberCounts.h"

using namespace std;
using namespace cbl;


// ===============================================================================================


std::vector<double> cbl::modelling::angularpk::z_limits (const double zmin, const double zmax, const double z_abs_err)
{
  const double dz = (z_abs_err>=0.) ? 3.*z_abs_err : 0.;
  
  vector<double> z_limits = { max(zmin-dz, 0.01), zmax+dz };

  return z_limits;
}

  
// ==================================================================================

  
std::vector<double> cbl::modelling::angularpk::ell_for_mixing (std::vector<double> l, double bin_width, bool limber_approximation)
{
  vector<double> new_l;
  int l_begin = 0;
  
  double delta_ell; // add some multipoles at the end, to avoid loss of power due to the truncation of Cl   //add few multipoles, to avoid the decrease due to mixing truncation (to cut later)
  
  if (bin_width<=1) delta_ell = 10;
  else delta_ell = 2*bin_width;
  if (!limber_approximation) l_begin = 2;
  for (int i=l_begin; i<l[l.size()-1]+delta_ell; i++) // mixing starts from 0, but the exact computetion of cl requires l>=2
    new_l.emplace_back(i);

  return new_l;
}


// ===============================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_mixed (std::vector<double> l_mixing, std::vector<std::vector<double>> mixing_matrix, std::vector<double> l, std::vector<double> Cl, double fsky, bool limber_approximation)
{

  bool weighted_average=true; 
  if (l_mixing[0]!=0 and limber_approximation)  ErrorCBL("the mixing matrix should start from l=0", "Cl_mixed", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);   
  if (l_mixing.size()<Cl.size())  ErrorCBL("the number of multipoles of the mixing matrix is lower than the number of multipoles of the angular power spectrum! Reduce the xrange or enlarge the mixing matrix.", "Cl_mixed", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);
  
  int lstart = 0; // in the limber approximation cl starts from l=0
  if (!limber_approximation) lstart = 2; // in exact case, cl starts from l=2, so we should reallign the mixing to skip l=0,1 of the matrix
  double _fsky = 1./fsky;
  vector<double> Cl_mixed(Cl.size(), 0);
  for (size_t j=0; j<Cl.size(); ++j)
    for (size_t i=0; i<Cl.size(); ++i)
      Cl_mixed[j] += mixing_matrix[j+lstart][i+lstart]*Cl[i]*_fsky;
  
  // in the exact case, extrapolate mixed cl to l=0,1
  if (!limber_approximation) { 
    glob::FuncGrid mixed_spectrum(l, Cl_mixed, "Spline", BinType::_linear_);
    Cl_mixed.insert(Cl_mixed.begin(), mixed_spectrum(1));
    Cl_mixed.insert(Cl_mixed.begin(), mixed_spectrum(0));
  }
  
  // Cl starts from l=0, but we need to consider the interval lmin-lmax
  // compute the average in l-bands
  vector<double> Cl_mixed_average(l.size()), l_average(l.size());
  double bin_width = l[1]-l[0];

  if (bin_width>1.) {
    int l_begin = int(max(l[0]-bin_width*0.5, 0.));
    int l_end = l[0]+bin_width*0.5;
    
    if (!weighted_average) { // if you want an arithmetic average of the power spectrum in Delta ell bands
      for (size_t i=0; i<l.size(); i++) {  
	int count = 0;
	for (int j=l_begin+1; j<l_end+1; j++) { // do not consider monopole in the average (also measurements should start from 1)
	  if (j+double(i)*bin_width>=double(Cl.size()))
	    break;
	  l_average[i]+=l_mixing[j+i*bin_width];
	  Cl_mixed_average[i]+=Cl_mixed[j+i*bin_width];
	  count=count+1;
	}
	l_average[i]/=count;	
	Cl_mixed_average[i]/=count;
      }  
      return Cl_mixed_average;
    }
    else{       //if you want to average the power spectrum in Delta ell bands, with weight (2*ell+1)
      for (size_t i=0; i<l.size(); i++) {  
	int count=0;
	double weight = 0;
	for (int j=l_begin+1; j<l_end+1; j++) { //do not consider monopole in the average (also measurements should start from 1)
	  if (j+double(i)*bin_width>=double(Cl.size()))
	    break;
	  l_average[i] += l_mixing[j+i*bin_width];
	  Cl_mixed_average[i] += Cl_mixed[j+i*bin_width]* (2 * l_mixing[j + i * bin_width] + 1);
	  weight += (2 * l_mixing[j + i * bin_width] + 1);
	  count ++;
	}
	l_average[i] /= count;	
	Cl_mixed_average[i] /= weight;
      }  
      return Cl_mixed_average;
    }    
  }
  else {
    vector<double> Cl_mixed_final;

    // Cl starts from l=0, but we need to consider the interval lmin-lmax
    for (size_t i=l[0]; i<l.size()+l[0]; i++)
      Cl_mixed_final.emplace_back(Cl_mixed[i]);
    
    return Cl_mixed_final;
  }
}


// ============================================================================================


double cbl::modelling::angularpk::redshift_selection_function (double redshift, std::vector<double> parameter, std::string distribution_type, double zmin, double zmax)
{
  double distribution=0.;
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
    ErrorCBL("distribution_type can be gaussian, polynomial or skewnormal", "redshift_selection_function", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);
  
  return distribution;
  
}


// ===========================================================================================


std::vector<double> cbl::modelling::angularpk::Cl (std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{ 
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  
  vector<double> l_mixing = l;
  
  if (pp->mixing_matrix.size()!=0)
    l = ell_for_mixing(l, (l[1]-l[0]), pp->limber_approximation);
  if (!pp->limber_approximation && l[0]<2) ErrorCBL("the exact computation of the power spectrum requires l>=2", "Cl", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);
  
  vector<double> Cl(l.size(), 0);

  
  // input parameters

  double bias = parameter[pp->Cpar.size()];
  double shot_noise= parameter[pp->Cpar.size()+1];
  
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);


  // define z_vector, larger than photometric bin limits

  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err); 
  if (pp->z_min_bin2>0 && pp->z_max_bin2>0) {
    vector<double> z_limit_bin2=z_limits(pp->z_min_bin2, pp->z_max_bin2, pp->z_abs_err); 
    z_limit[0] = min(z_limit[0], z_limit_bin2[0]);
    z_limit[1] = max(z_limit[1], z_limit_bin2[1]);
  }
  
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);


  // compute and interpolate growth factor, lgr and comoving distance

  vector<double> growth_factor, linear_growth_rate, comoving_distance;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    linear_growth_rate.emplace_back(cosmo->linear_growth_rate(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));    
  }
  
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid lgr(z_vector, linear_growth_rate, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);


  // compute and interpolate dN_dz distribution

  vector<double> redshift_counts;
  for (size_t j=0; j<z_vector.size(); j++)
    redshift_counts.emplace_back(redshift_selection_function(z_vector[j],pp->dN_par, pp->distribution_type, pp->z_min, pp->z_max));   //pp->z_min and pp->z_min are here only to obtain zero outside limits, for polynomial distribution 
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max); 
  glob::FuncGrid pk0_interp(kk, Pk0, pp->interpType);
  bool cross = false;
  
  if (pp->dN_par_bin2.size()>0) {
    cross=true;
    vector<double> redshift_counts_bin2;
    vector<double> z_limit_bin2 = z_limits(pp->z_min_bin2, pp->z_max_bin2, pp->z_abs_err); 
    for (size_t j=0; j<z_vector.size(); j++)
      redshift_counts_bin2.emplace_back(redshift_selection_function(z_vector[j],pp->dN_par_bin2, pp->distribution_type, pp->z_min_bin2, pp->z_max_bin2));//pp->z_min_bin2 and pp->z_min_bin2 are here only to obtain zero outside limits, for polynomial distribution 
    
    glob::FuncGrid redshift_counts_interp_bin2(z_vector, redshift_counts_bin2, pp->interpType, BinType::_linear_);
    double normalization_redshift_selection_function_bin2 = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp_bin2, z_limit_bin2[0], z_limit_bin2[1]);
    double bias_bin2 = parameter[pp->Cpar.size()+1];
    for (size_t i=0; i<Cl.size(); i++)
      Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp, comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp_bin2, bias_bin2, cross, parameter, inputs)/normalization_redshift_selection_function/normalization_redshift_selection_function_bin2;
    
    parameter[pp->Cpar.size()+2] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3); // S8
    parameter[pp->Cpar.size()+3] = cosmo->Omega_baryon()/cosmo->Omega_matter();     // baryonic fraction
  }
  
  else { // autocorrelation

    for (size_t i=0; i<l.size(); i++)
      Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp,comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp, 1., cross, parameter, inputs)/pow(normalization_redshift_selection_function, 2)+shot_noise*pp->shot_noise_normalization;
    parameter[pp->Cpar.size()+2] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3); // S8
    parameter[pp->Cpar.size()+3] = cosmo->Omega_baryon()/cosmo->Omega_matter();     // baryonic fraction

  }
     
  if (pp->mixing_matrix.size()!=0)
    Cl = Cl_mixed(pp->ll, pp->mixing_matrix, l_mixing, Cl, pp->fsky, pp->limber_approximation); // Cl starts from l=0, then we need to cut ls from 0 to lmin, and from lmax to ll max

  return Cl;
}


// ============================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_scaling_relation_cosmology (std::vector<double> l, std::shared_ptr<void> inputs, std::vector<double> &parameter)
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
  for (size_t i=0; i<parameter.size()-2; i++)
    scalRel_pars.emplace_back(parameter[i]);

  const double alpha = scalRel_pars[scalRel_pars.size()-8];
  const double beta = scalRel_pars[scalRel_pars.size()-7];
  const double gamma = scalRel_pars[scalRel_pars.size()-6];
  const double scatter0 = scalRel_pars[scalRel_pars.size()-5];
  const double scatterM = scalRel_pars[scalRel_pars.size()-4];
  const double scatterM_exp = scalRel_pars[scalRel_pars.size()-3];
  const double scatterz = scalRel_pars[scalRel_pars.size()-2];
  const double scatterz_exp = scalRel_pars[scalRel_pars.size()-1];

  //shot noise
  const double shot_noise=parameter[parameter.size()-2];    

  double log_base = (pp->scaling_relation)->data_model().log_base;
  double mass_pivot = (pp->scaling_relation)->data_model().mass_pivot;
  double proxy_pivot = (pp->scaling_relation)->data_model().proxy_pivot;
  double redshift_pivot = (pp->scaling_relation)->data_model().redshift_pivot;
  
  // Interpolate the normalised amplitude of the growing mode, the linear growth rate and comoving distance
  const vector<double> z_for_interp = linear_bin_vector(30, 0.0001, Max((pp->scaling_relation)->data_model().redshift)+5.*pp->z_abs_err);
  vector<double> DN (z_for_interp.size(), 0.), linear_growth_rate(z_for_interp.size(), 0.), comoving_distance (z_for_interp.size(), 0.);
  for (size_t i=0; i<z_for_interp.size(); i++) {
    DN[i] = cosmo->DN(z_for_interp[i]);
    linear_growth_rate[i] = cosmo->linear_growth_rate(z_for_interp[i]);
    comoving_distance[i] = cosmo->D_C(z_for_interp[i]);
  }
  glob::FuncGrid interp_DN(z_for_interp, DN, pp->interpType, BinType::_linear_);
  glob::FuncGrid lgr(z_for_interp, linear_growth_rate, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_for_interp, comoving_distance, pp->interpType, BinType::_linear_);
  
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the effective bias from the scaling relation
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  double bias = bias_from_scaling_relation(alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, log_base, mass_pivot, proxy_pivot, redshift_pivot, scalRel_pars, interp_DN, inputs);
  
  // set the value of the bias
  parameter[parameter.size()-1] = bias;

  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Compute the angular power spectrum
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  
  vector<double> l_mixing=l;
  
  if (pp->mixing_matrix.size()!=0)
    l=ell_for_mixing(l, (l[1]-l[0]), pp->limber_approximation);
  if (!pp->limber_approximation && l[0]<2) ErrorCBL("the exact computation of the power spectrum requires l>=2", "Cl_scaling_relation_cosmology", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);
  
  vector<double> Cl(l.size(), 0);

  //define and interpolate Pk
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);

  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max); 
  glob::FuncGrid pk0_interp(kk, Pk0, pp->interpType);
  
  //find minimum and maximum for true redshift, which is different from limits of the photometric redshift bins
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err); 
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);
  //double z_bin_width=(z_vector[1]-z_vector[0]);
  
  //interpolate and normalize redshift counts
  vector<double> redshift_counts;

  for (size_t j=0; j<z_vector.size(); j++)
    redshift_counts.emplace_back(redshift_selection_function(z_vector[j],pp->dN_par, pp->distribution_type, pp->z_min, pp->z_max));
  
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);
  
  //compute auto power spectrum
  bool cross=false; 
  for (size_t i=0; i<l.size(); i++)
    Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp, comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp, 1., cross, parameter, inputs)/pow(normalization_redshift_selection_function, 2)+shot_noise*pp->shot_noise_normalization;

  //convolve with mixing matrix
  if (pp->mixing_matrix.size()!=0)
    Cl=Cl_mixed(pp->ll, pp->mixing_matrix, l_mixing, Cl, pp->fsky, pp->limber_approximation);   //Cl starts from l=0, then we need to cut ls from 0 to lmin, and from lmax to ll max
  
  return Cl;
  
}


// ============================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_MF_theory (std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  vector<double> l_mixing=l;
  if (pp->mixing_matrix.size()!=0)
    l = ell_for_mixing(l, (l[1]-l[0]), pp->limber_approximation);
  if (!pp->limber_approximation && l[0]<2) ErrorCBL("the exact computation of the power spectrum requires l>=2", "Cl_MF_theory", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);

  vector<double> Cl(l.size(), 0);
  
  // input parameters
  double bias = parameter[pp->Cpar.size()];
  double shot_noise = parameter[pp->Cpar.size()+1];

  //create zz and dn/dz vectors
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);

  //compute and interpolate growth_factor, linear_growth_rate and comoving_distance
  vector<double> growth_factor, linear_growth_rate, comoving_distance;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    linear_growth_rate.emplace_back(cosmo->linear_growth_rate(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid lgr(z_vector, linear_growth_rate, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);

  //compute and interpolate the power spectrum
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);

  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max); 
  glob::FuncGrid pk0_interp(kk, Pk0, pp->interpType);
  
  //compute the redshift number counts
  double area_rad = pp->area_degrees/180/180*par::pi*par::pi;
  vector<double> redshift_counts;
  
  //generate mass function
  vector<double> Mass_vector = logarithmic_bin_vector(pp->Mass_step, pp->Mass_min, pp->Mass_max);
  
  for (size_t i=0; i<z_vector.size(); ++i) {
    vector<double> MF = modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max);

    glob::FuncGrid MF_interp(Mass_vector, MF, pp->interpType, BinType::_logarithmic_);

    auto integrand_MF = [&MF_interp, &cosmo, &z_vector, &i, &pp] (const double mass)
    {
      if (pp->z_abs_err>0) {
	// Gaussian convolution, within the photometric bin
	double int_P_z=0.5 * ( erf( (pp->z_max - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) - erf( (pp->z_min - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) );
	return int_P_z*MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
      }
      return MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
    };
    redshift_counts.emplace_back(area_rad*wrapper::gsl::GSL_integrate_qag(integrand_MF, pp->Mass_min, pp->Mass_max));
  }
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);    
  double redshift_counts_normalization=wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);
  
  bool cross=false;
  for (size_t i=0; i<l.size(); i++)
    Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp, comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp, 1., cross, parameter, inputs)/pow(redshift_counts_normalization, 2)+shot_noise*pp->shot_noise_normalization;
  parameter[pp->Cpar.size()+2] = cosmo->sigma8()*sqrt(cosmo->Omega_matter()/0.3);     //S8
  parameter[pp->Cpar.size()+3] = cosmo->Omega_baryon()/cosmo->Omega_matter();  //baryonic fraction
  
  if (pp->mixing_matrix.size()!=0)
    Cl=Cl_mixed(pp->ll, pp->mixing_matrix, l_mixing, Cl, pp->fsky, pp->limber_approximation);   //Cl starts from l=0, then we need to cut ls from 0 to lmin, and from lmax to ll max
  
  return Cl;
}


// ============================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_bias_MF_theory (std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  
  vector<double> l_mixing = l;
  if (pp->mixing_matrix.size()!=0)
    l = ell_for_mixing(l, (l[1]-l[0]), pp->limber_approximation);
  
  if (!pp->limber_approximation && l[0]<2) ErrorCBL("the exact computation of the power spectrum requires l>=2", "Cl_bias_MF_theory", "ModelFunction_PowerSpectrum_Angular.cpp", glob::ExitCode::_error_);

  vector<double> Cl(l.size(), 0);
  double redshift= cbl::Average(pp->cluster_redshift);

  // compute and interpolate the linear power spectrum, for mass function
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);

  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., true, pp->output_root, pp->norm, pp->k_min, pp->k_max); 
  glob::FuncGrid pk0_interp(kk, Pk0, pp->interpType);

  // define z_vector, larger than photometric bin limits
  vector<double> z_limit = z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector = linear_bin_vector(30, z_limit[0], z_limit[1]);

  // compute the redshift number counts
  double area_rad = pp->area_degrees/180/180*par::pi*par::pi;
  vector<double> redshift_counts;
  
  // generate mass function
  vector<double> Mass_vector=logarithmic_bin_vector(pp->Mass_step, pp->Mass_min, pp->Mass_max);
  
  for (size_t i=0; i<z_vector.size(); ++i) {
    vector<double> MF = modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max); 
    glob::FuncGrid MF_interp(Mass_vector, MF, pp->interpType, BinType::_logarithmic_);
    
    auto integrand_MF = [&MF_interp, &cosmo, &z_vector, &i, &pp] (const double mass)
    {
      if (pp->z_abs_err>0) {
	// Gaussian convolution, within the photometric bin
	double int_P_z=0.5 * ( erf( (pp->z_max - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) - erf( (pp->z_min - z_vector[i]) / (sqrt(2)*pp->z_abs_err) ) );
	return int_P_z*MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
      }
      return MF_interp(mass)*cosmo->dV_dZdOmega(z_vector[i], true);
    };
    redshift_counts.emplace_back(area_rad*wrapper::gsl::GSL_integrate_qag(integrand_MF, pp->Mass_min, pp->Mass_max));
  }
  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);    
  double redshift_counts_normalization = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  // compute and interpolate growth_factor, linear_growth_rate and comoving_distance
  vector<double> growth_factor, linear_growth_rate, comoving_distance, ones;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    linear_growth_rate.emplace_back(cosmo->linear_growth_rate(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));
    ones.emplace_back(1.); 
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid lgr(z_vector, linear_growth_rate, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);

  // compute the bias
  double bias;
  cosmology::Bias bb(cosmo);
  if (pp->bias_from_sample) {
    vector<double> mass_grid = logarithmic_bin_vector(pp->Mass_step, pp->Mass_min, pp->Mass_max);
    bias = bb.bias_eff_mass(pp->cluster_mass, mass_grid, pp->cluster_redshift, pp->model_bias, pp->method_Pk, "mean_bias", true, pp->output_root, pp->Delta)[0];  // Delta is already converted. If isDelta_critical=true, Delta=Delta_input/OmegaM 
  }
  else {
    vector<double> selection_function;
    for (size_t i=0; i<Mass_vector.size(); i++)
      selection_function.emplace_back(1.); 
    glob::FuncGrid interp_SF(Mass_vector, selection_function, pp->interpType, BinType::_linear_);
    const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM (Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);

    double Delta_crit = (pp->isDelta_critical) ? pp->Delta_input : pp->Delta_input*cosmo->Omega_matter(redshift);

    bias=bb.bias_eff_selection_function(interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_SF, Min(pp->cluster_mass), Max(pp->cluster_mass), {redshift}, pp->model_bias, pp->model_MF, pp->method_Pk, 1., false, pp->output_root, Delta_crit, -1., pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, par::defaultString, true)[0];
  }
  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);

  parameter[pp->Cpar.size()+1] = bias;
  double shot_noise=parameter[pp->Cpar.size()]; 

  //compute the angular power spectrum
  bool cross=false;
  for (size_t i=0; i<l.size(); i++)
    Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp, comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp, 1., cross, parameter, inputs)/pow(redshift_counts_normalization, 2);
  
  if (pp->mixing_matrix.size()!=0)
    Cl = Cl_mixed(pp->ll, pp->mixing_matrix, l_mixing, Cl, pp->fsky, pp->limber_approximation); // Cl starts from l=0, then we need to cut ls from 0 to lmin, and from lmax to ll max

  for (size_t i=0; i<Cl.size(); i++)
    Cl[i]+=shot_noise*pp->shot_noise_normalization;

  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);
      
  return Cl;
}


// ============================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_bias_MF_theory_BAO(std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
  
  vector<double> l_mixing=l;
  if (pp->mixing_matrix.size()!=0)
    l=ell_for_mixing(l, (l[1]-l[0]), pp->limber_approximation);

  if (!pp->limber_approximation && l[0]<2) ErrorCBL("the exact computation of the power spectrum requires l>=2", "Cl_MF_theory", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);

  vector<double> Cl(l.size(), 0);
  double redshift= cbl::Average(pp->cluster_redshift);

  //compute and interpolate the linear power spectrum, for mass function
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  
  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, false, 0., true, pp->output_root, pp->norm, pp->k_min, pp->k_max); 
    
  //define z_vector, larger than photometric bin limits
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err);
  vector<double> z_vector= linear_bin_vector(30, z_limit[0], z_limit[1]);

  //compute the redshift number counts
  double area_rad=pp->area_degrees/180/180*par::pi*par::pi;
  vector<double> redshift_counts;
  
  //generate mass function
  vector<double> Mass_vector=logarithmic_bin_vector(pp->Mass_step, pp->Mass_min, pp->Mass_max);

  for (size_t i=0; i<z_vector.size(); ++i) {
    vector<double> MF=modelling::numbercounts::mass_function(Mass_vector, cosmo, z_vector[i], pp->model_MF, false, pp->Delta_input, pp->isDelta_critical, kk, Pk0, pp->interpType, pp->k_max); 
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
  double redshift_counts_normalization=wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  //compute and interpolate growth_factor, linear_growth_rate and comoving_distance
  vector<double> growth_factor, linear_growth_rate, comoving_distance, ones;
  for (size_t i=0; i<z_vector.size(); ++i) {
    growth_factor.emplace_back(cosmo->DN(z_vector[i]));
    linear_growth_rate.emplace_back(cosmo->linear_growth_rate(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));
    ones.emplace_back(1.); 
  }
  glob::FuncGrid interp_DN(z_vector, growth_factor, pp->interpType, BinType::_linear_);
  glob::FuncGrid lgr(z_vector, linear_growth_rate, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_vector, comoving_distance, pp->interpType, BinType::_linear_);

  //compute the bias
  double bias;
  cosmology::Bias bb(cosmo);
  if (pp->bias_from_sample) {
    vector<double> mass_grid = logarithmic_bin_vector(pp->Mass_step, pp->Mass_min, pp->Mass_max);
    bias = bb.bias_eff_mass(pp->cluster_mass, mass_grid, pp->cluster_redshift, pp->model_bias, pp->method_Pk, "mean_bias", true, pp->output_root, pp->Delta)[0];  // Delta is already converted. If isDelta_critical=true, Delta=Delta_input/OmegaM 
  }
  else{
    vector<double> selection_function;
    for (size_t i=0; i<Mass_vector.size(); i++)
      selection_function.emplace_back(1.); 
    glob::FuncGrid interp_SF(Mass_vector, selection_function, pp->interpType, BinType::_linear_);
    const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM (Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);

    double Delta_crit = (pp->isDelta_critical) ? pp->Delta_input : pp->Delta_input*cosmo->Omega_matter(redshift);

    bias=bb.bias_eff_selection_function(interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_SF, Min(pp->cluster_mass), Max(pp->cluster_mass), {redshift}, pp->model_bias, pp->model_MF, pp->method_Pk, 1., false, pp->output_root, Delta_crit, -1., pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, par::defaultString, true)[0];
  }
  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);

  parameter[pp->Cpar.size()+2] = bias;
  double shot_noise=parameter[pp->Cpar.size()]; 
  double sigmaNL=parameter[pp->Cpar.size()+1];

  //BAO nonlinear correction
  //compute linear Pk and PkNW (without BAO)
  vector<double> Pk0NW = PX.Pk_matter(kk, "EisensteinHu", false, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  vector<double> Pk0_corr(kk.size(), 0);
  
  for (size_t i =0; i<kk.size(); i++)
    Pk0_corr[i] = Pk0NW[i]*(1.+(Pk0[i]/Pk0NW[i]-1.)*exp(-0.5*pow(kk[i]*sigmaNL, 2)));
  
  glob::FuncGrid pk0_interp(kk, Pk0_corr, pp->interpType, BinType::_logarithmic_);
  
  //compute the angular power spectrum
  bool cross=false;
  for (size_t i=0; i<l.size(); i++)
    Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp, comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp, 1., cross, parameter, inputs)/pow(redshift_counts_normalization, 2);
    
  if (pp->mixing_matrix.size()!=0)
    Cl=Cl_mixed(pp->ll, pp->mixing_matrix, l_mixing, Cl, pp->fsky, pp->limber_approximation);   //Cl starts from l=0, then we need to cut ls from 0 to lmin, and from lmax to ll max

  for (size_t i=0; i<Cl.size(); i++)
    Cl[i]+=shot_noise*pp->shot_noise_normalization;

  PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);
      
  return Cl;
}


//=========================================================================================


std::vector<double> cbl::modelling::angularpk::Cl_scaling_relation_cosmology_MF_theory (std::vector<double> l, std::shared_ptr<void> inputs, std::vector<double> &parameter)
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
  const double shot_noise = parameter[pp->Cpar.size()+11];

  double log_base = (pp->scaling_relation)->data_model().log_base;
  double mass_pivot = (pp->scaling_relation)->data_model().mass_pivot;
  double proxy_pivot = (pp->scaling_relation)->data_model().proxy_pivot;
  double redshift_pivot = (pp->scaling_relation)->data_model().redshift_pivot;

  // Interpolate the normalised amplitude of the growing mode, the linear growth rate and comoving distance
  const vector<double> z_for_interp = linear_bin_vector(30, std::max(Min((pp->scaling_relation)->data_model().redshift)-5.*pp->z_abs_err,0.001), Max((pp->scaling_relation)->data_model().redshift)+5.*pp->z_abs_err);
  vector<double> DN(z_for_interp.size(), 0.), linear_growth_rate(z_for_interp.size(), 0.), comoving_distance (z_for_interp.size(), 0.);
  for (size_t i=0; i<z_for_interp.size(); i++) {
    DN[i] = pp->cosmology->DN(z_for_interp[i]);
    linear_growth_rate[i] = pp->cosmology->linear_growth_rate(z_for_interp[i]);
    comoving_distance[i] = pp->cosmology->D_C(z_for_interp[i]);
  }
  
  glob::FuncGrid interp_DN (z_for_interp, DN, pp->interpType, BinType::_linear_);
  glob::FuncGrid lgr(z_for_interp, linear_growth_rate, pp->interpType, BinType::_linear_);
  glob::FuncGrid comoving_distance_interp(z_for_interp, comoving_distance, pp->interpType, BinType::_linear_);

  
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the effective bias from the scaling relation
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  double bias = bias_from_scaling_relation(alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, log_base, mass_pivot, proxy_pivot, redshift_pivot, scalRel_pars, interp_DN, inputs);

  // set the value of the bias
  parameter[parameter.size()-1] = bias;

  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Compute the angular power spectrum
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  
  vector<double> l_original=l;
  
  if (pp->mixing_matrix.size()!=0)
    l=ell_for_mixing(l, (l[1]-l[0]), pp->limber_approximation);
  
  if (!pp->limber_approximation && l[0]<2) ErrorCBL("the exact computation of the power spectrum requires l>=2", "Cl_scaling_relation_cosmology_MF_theory", "ModelFunction_PowerSpectrum_Angular.cpp",glob::ExitCode::_error_);
 
  vector<double> Cl(l.size(), 0);

  // define and interpolate Pk
  vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);

  cosmology::PkXi PX(cosmo);
  vector<double> Pk0 = PX.Pk_matter(kk, pp->method_Pk, pp->NL, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max); 
  glob::FuncGrid pk0_interp(kk, Pk0, pp->interpType);
  
  // interpolate sigmaM and its derivative
  const vector<double> Mass_vector = logarithmic_bin_vector(300, 1.e10, 1.e16);
  const vector<glob::FuncGrid> interp_sigmaM_dlnsigmaM = modelling::numbercounts::sigmaM_dlnsigmaM (Mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);

  // interpolate weights vector
  double log_bin_size = (log10(pp->proxy_max)-log10(pp->proxy_min))/(pp->SF_weights.size());
  vector<double> lambda_star_vector=logarithmic_bin_vector(pp->SF_weights.size(), pow(10, log10(pp->proxy_min)+0.5*log_bin_size), pow(10, log10(pp->proxy_max)-0.5*log_bin_size)); 
  glob::FuncGrid weights_interp(lambda_star_vector, pp->SF_weights, pp->interpType, BinType::_logarithmic_);
  
  // find minimum and maximum for true redshift, which is different from limits of the photometric redshift bins
  vector<double> z_limit = z_limits(pp->z_min, pp->z_max, pp->z_abs_err); 
  vector<double> z_vector = linear_bin_vector(12, z_limit[0], z_limit[1]);

  // interpolate and normalize redshift counts
  vector<double> redshift_counts;
  double area_rad=pp->area_degrees/180/180*par::pi*par::pi;
    
  for (size_t i=0; i<z_vector.size(); ++i)
    redshift_counts.emplace_back(redshift_distribution_proxy(z_vector[i], alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, pp->z_abs_err, pp->proxy_rel_err, Plambda_a, Plambda_b, Plambda_c, pp->fz, pp->z_min, pp->z_max, pp->proxy_min, pp->proxy_max, cosmo, area_rad, pp->model_MF, pp->Delta_input, pp->isDelta_critical, interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_DN, proxy_pivot, redshift_pivot, mass_pivot, log_base, weights_interp)*10);

  glob::FuncGrid redshift_counts_interp(z_vector, redshift_counts, pp->interpType, BinType::_linear_);
  double normalization_redshift_selection_function = wrapper::gsl::GSL_integrate_qag(redshift_counts_interp, z_limit[0], z_limit[1]);

  // compute auto power spectrum
  bool cross = false;
  for (size_t i=0; i<l.size(); i++)
    Cl[i] = modelling::angularpk::integral_cl(l[i], kk, pk0_interp, comoving_distance_interp, interp_DN, lgr, redshift_counts_interp, bias, redshift_counts_interp, 1., cross, parameter, inputs)/pow(normalization_redshift_selection_function, 2);
  
  // convolve with mixing matrix
  if (pp->mixing_matrix.size()!=0)
    Cl = Cl_mixed(pp->ll, pp->mixing_matrix, l_original, Cl, pp->fsky, pp->limber_approximation);

  for (size_t i=0; i<Cl.size(); i++)
    Cl[i]+=shot_noise*pp->shot_noise_normalization;

  return Cl;
   
}


//===========================================================================================


double cbl::modelling::angularpk::redshift_distribution_proxy (const double redshift, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_abs_err, const double proxy_rel_err, const double Plambda_a, const double Plambda_b, const double Plambda_c, ProxyFunction fz, const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double Area, const std::string model_MF, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base, cbl::glob::FuncGrid weights_interp)
{
  double fact = (cosmology->unit()) ? 1 : cosmology->little_h();
  shared_ptr<void> pp;
  auto cosmology_ptr = move(cosmology);

  // Declare the normalized mass and the redshift, used as constants in integrand_P_M__z (which is called in integrand)
  double normM = 0., the_redsh = 0.;

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
    double P_M__z = (P_M__lambda_z*P_lambda__z > 0) ? wrapper::gsl::GSL_integrate_cquad(integrand_P_M__z,log(0.0001)/log(log_base),log(100.)/log(log_base)) * log(log_base) : 1.;

    // Compute the integrals of P(z|z) and P(lambda|lambda)
    double mean_Pz = redshift;
    double int_P_z = 0.5 * ( erf( (redshift_max - mean_Pz) / (sqrt(2)*z_abs_err) ) - erf( (redshift_min - mean_Pz) / (sqrt(2)*z_abs_err) ) );
    double mean_Plambda = x[1];
    double int_P_lambda = 0.5 * ( erf( (proxy_max - mean_Plambda) / (sqrt(2)*proxy_rel_err*proxy_max) ) - erf( (proxy_min - mean_Plambda) / (sqrt(2)*proxy_rel_err*proxy_min) ) );
    double weights_int = weights_interp(x[1]);
    
    if (weights_int>1) weights_int = 1.;
    else if (weights_int<0) weights_int = 0.;

    cosmology::MassFunction mass_func(cosmology_ptr);

    return mass_func.mass_function(Mass, interp_sigmaM(Mass*fact), interp_DlnsigmaM(Mass*fact), redshift, interp_DN(redshift), model_MF, false, par::defaultString, Delta_)* pow(log_base,normM) * (P_M__lambda_z*P_lambda__z/P_M__z) * int_P_z * int_P_lambda*weights_int;
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
  integration_limits[1] = {max(proxy_min - 3.5*proxy_rel_err*proxy_min, 0.00001), proxy_max + 3.5*proxy_rel_err*proxy_max};
  // Compute the integral
  wrapper::cuba::CUBAwrapper CW (integrand, integral_dimension);
  double dN_dz;
 
  if (integration_limits[0][0] < integration_limits[0][1])
    dN_dz = CW.IntegrateVegas(integration_limits, false)*Area*cosmology_ptr->dV_dZdOmega(redshift, true);
  else
    dN_dz = 0;

  return dN_dz * mass_pivot * log(log_base);
}


// ===============================================================================================


double cbl::modelling::angularpk::integral_cl (const double l, std::vector<double> kk, cbl::glob::FuncGrid pk0_interp, cbl::glob::FuncGrid comoving_distance_interp, cbl::glob::FuncGrid interp_DN, cbl::glob::FuncGrid linear_growth_rate_interp, cbl::glob::FuncGrid redshift_counts_interp, const double bias, cbl::glob::FuncGrid redshift_counts_interp2, const double bias2, bool cross, std::vector<double> par, std::shared_ptr<void> inputs) {
  
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  
  //set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], par[i]);
  
  //set integration limits, more extended than photometric redshift bin
  vector<double> z_limit=z_limits(pp->z_min, pp->z_max, pp->z_abs_err); 
  if (pp->z_min_bin2>0 && pp->z_max_bin2>0) {
    vector<double> z_limit_bin2=z_limits(pp->z_min_bin2, pp->z_max_bin2, pp->z_abs_err); 
    z_limit[0] = min(z_limit[0], z_limit_bin2[0]);
    z_limit[1] = max(z_limit[1], z_limit_bin2[1]);
  }

  // relative measurement precision (it raises a warning message if not satisfied)
  const double prec = 1.e-2;
  const double abs_err = 5.;

  // default values (refer to the GSL documentation for different choiches)
  const int limit_size = 100;
  const int rule = 6;

  if (l>pp->ell_limber || pp->limber_approximation) {  // case l>ell_limber, multipole transition between exact power spectrum and Limber approximation

    const double prec = 1.e-3;
    const double abs_err = 0.;
    const int limit_size = 1000;

    function<double(double)> integrand_limber = [&l, &par, &pk0_interp, &interp_DN, &redshift_counts_interp, &redshift_counts_interp2, &cosmo, &inputs, &kk, &bias, &bias2] (double redshift)
    {
      shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

      double dc = cosmo->D_C(redshift), dn_interp = interp_DN(redshift), red_counts_interp = redshift_counts_interp(redshift);
      double _cc = 1/cbl::par::cc;
      double inv_d2 = 1/(dc*dc);
      
      if (pp->dN_par_bin2.size()>0) 
	return bias*bias2*red_counts_interp*redshift_counts_interp2(redshift)*pk0_interp((l+0.5)/dc)*dn_interp*dn_interp*inv_d2*cosmo->Hubble(redshift)*_cc;
      
      return bias*bias*red_counts_interp*red_counts_interp*pk0_interp((l+0.5)/dc)*dn_interp*dn_interp*inv_d2*cosmo->Hubble(redshift)*_cc;
      
    };

    return wrapper::gsl::GSL_integrate_qag(integrand_limber, z_limit[0], z_limit[1], prec, abs_err, limit_size, rule);
  }     //end if l>ell_limber
  
  else {   //if l<ell_limber exact power spectrum
    
    //-------------------------------
    //--------- Cl exact ------------
    //-------------------------------
   
    function<double(double)> integrand = [&pk0_interp, &comoving_distance_interp, &interp_DN, &linear_growth_rate_interp, &redshift_counts_interp, &bias, &redshift_counts_interp2, &bias2, &pp, &cosmo, &l, &cross, &z_limit, &prec, &abs_err, &limit_size, &rule] (const double k)   //integrand of the Cl exact
    {
      glob::FuncGrid redshift_distribution;
      if (cross)
	redshift_distribution=redshift_counts_interp2;
      else 
	redshift_distribution=redshift_counts_interp;

      auto integrand_Wk = [&] (const double redshift)   //integrand of the radial window function 
      {
	if ((l>7 and k*comoving_distance_interp(redshift)<l*0.7) or interp_DN(redshift)<0)
	  return 0.; 
	return redshift_distribution(redshift)*interp_DN(redshift)*jl(k*comoving_distance_interp(redshift), l);
      };
      
      auto integrand_Wkr = [&k, &comoving_distance_interp, &interp_DN, &linear_growth_rate_interp, &pp, &redshift_distribution, &cosmo, &l] (const double redshift)  //integrand of the radial window function, second term, accounting for RSD
      {
	if (l>7 and k*comoving_distance_interp(redshift)<l*0.7)
	  return 0.;
	return redshift_distribution(redshift)*interp_DN(redshift)*linear_growth_rate_interp(redshift)*((2.*l*l+2.*l-1)/(2.*l+3.)/(2.*l-1.)*jl(k*comoving_distance_interp(redshift), l)-l*(l-1.)/(2.*l-1.)/(2.*l+1.)*jl(k*comoving_distance_interp(redshift), l-2)-(l+1.)*(l+2.)/(2.*l+1.)/(2.*l+3.)*jl(k*comoving_distance_interp(redshift), l+2));
      };

      double integral_Wk = wrapper::gsl::GSL_integrate_cquad(integrand_Wk, z_limit[0], z_limit[1], 1.e-2); //, abs_err);//, limit_size, rule);  //integral of window function
      double integral_Wkr = 0;
	
      if (pp->redshift_space)  //in real space integral_Wkr=0
	integral_Wkr = cbl::wrapper::gsl::GSL_integrate_cquad(integrand_Wkr, z_limit[0], z_limit[1], 1.e-1); //, abs_err); //, limit_size, rule);

      if (cross) {  //cross correlation
	double integral_Wk_bin2 = cbl::wrapper::gsl::GSL_integrate_qag(integrand_Wk, z_limit[0], z_limit[1], prec, abs_err, limit_size, rule);  //integral of window function
	double integral_Wkr_bin2 = 0;
	  
	if (pp->redshift_space)  //in real space integral_Wkr=0
	  integral_Wkr_bin2 = wrapper::gsl::GSL_integrate_qag(integrand_Wkr, z_limit[0], z_limit[1], prec, abs_err, limit_size, rule);
	  
	return k*k*pk0_interp(k)*(integral_Wk*bias+integral_Wkr)*(integral_Wk_bin2*bias2+integral_Wkr_bin2);
      }
      else
	return k*k*pk0_interp(k)*pow((integral_Wk*bias+integral_Wkr), 2);
    };

    return 2/par::pi*wrapper::gsl::GSL_integrate_cquad(integrand, pp->k_min, pp->k_max, 1.e-2);
    
  } //close l<ell_limber
  
}


//================================================================================================================


double cbl::modelling::angularpk::bias_from_scaling_relation (const double alpha, const double beta, const double gamma,  const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, double log_base, double mass_pivot, double proxy_pivot, double redshift_pivot, std::vector<double> scalRel_pars, cbl::glob::FuncGrid interp_DN, std::shared_ptr<void> inputs)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters (contained in scalRel_pars)
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], scalRel_pars[i]);
  
  // interpolate sigmaM
  const vector<double> Mass_vector = logarithmic_bin_vector(300, 1.e10, 1.e16);

  cosmology::Sigma SG(cosmo);
  vector<double> sigma2M = SG.sigma2M(Mass_vector, pp->method_Pk, 0., false, pp->output_root, "Linear", 100.);
  vector<double> sigmaM(Mass_vector.size(), 0.);
  for (size_t i=0; i<Mass_vector.size(); i++)
    sigmaM[i] = sqrt(sigma2M[i]);

  //PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);
  glob::FuncGrid sigmaM_interp(Mass_vector, sigmaM, pp->interpType);

  // compute the bias
  cosmology::Bias bb(cosmo);
  double bias = 0.;
  
  if (pp->z_abs_err == -1 && pp->proxy_rel_err == -1) {
         
    vector<double> _bias(pp->scaling_relation->data()->xx().size());
  
    for (size_t i=0; i<pp->scaling_relation->data()->xx().size(); i++) {
       
      double log_lambda = log(pp->scaling_relation->data()->xx(i)/proxy_pivot)/log(log_base);
      double log_fz = log( (pp->scaling_relation)->data_model().fz((pp->scaling_relation)->data_model().redshift[i], redshift_pivot, cosmo) )/log(log_base);
     
      double scatter_intr = scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_fz, scatterz_exp);
     
      double log_mass = (pp->scaling_relation)->likelihood()->get_m_model()->operator()(pp->scaling_relation->data()->xx(i), scalRel_pars) + scatter_intr;
      double mass = pow(log_base, log_mass) * mass_pivot;
  
      double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo->Omega_matter((pp->scaling_relation)->data_model().redshift[i]) : pp->Delta_input;
      double z = (pp->scaling_relation)->data_model().redshift[i];
       
      _bias[i] = bb.bias_halo(mass, sigmaM_interp(mass), z, interp_DN(z), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
       
    }
  
    bias = Average(_bias);
    return bias;
  
  } else {
  
    // !!
    // WARNING: actually we only integrate over M, neglecting the integrals over z and proxy (such integrations are negligible).
    // !!
    // To verify yourself that integrating over z and proxy produces negligible differences, uncomment the commented "integrand" as well as the following line:
    //
    // _bias[i] = CW.IntegrateVegas(integration_limits,false);
    //
    // Beware that, in an example test, computing such full integral increases the time from 2 second to 3 minutes!
    // Definitely not worthy to compute.

    // define the integrand
    double dummy_proxy, dummy_z;
    shared_ptr<void> ptr;
  
    auto integrand = [&] (const double x)
    {
      double mass = pow(log_base,x)*mass_pivot;
                
      // compute P(M|lambda,z)
      double log_lambda = log(dummy_proxy/proxy_pivot)/log(log_base);
      double log_f_z = log( (pp->scaling_relation)->data_model().fz(dummy_z, redshift_pivot, cosmo) )/log(log_base);
                
      double mean = alpha + beta*log_lambda + gamma*log_f_z;
      double scatter_intr = abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
      double P_M__lambda_z = gaussian(x, ptr, {mean,scatter_intr});
  
      // compute the halo bias
      double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo->Omega_matter(dummy_z) : pp->Delta_input;
      double bias_halo = bb.bias_halo(mass, sigmaM_interp(mass), dummy_z, interp_DN(dummy_z), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
                
      return bias_halo*P_M__lambda_z;
    };
  
    /*
      auto integrand = [&] (const vector<double> x)
      {
  
      double mass = pow(log_base,x[0])*mass_pivot;
                
      // Compute P(M|lambda,z)
      double log_lambda = log(x[2]/proxy_pivot)/log(log_base);
      double log_f_z = log( (pp->scaling_relation)->data_model().fz(x[1], redshift_pivot, cosmo) )/log(log_base);
                
      double mean = alpha + beta*log_lambda + gamma*log_f_z;
      double scatter_intr = scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp);
      double P_M__lambda_z = (gaussian(x[0], ptr, {mean,scatter_intr}));
  
      // Compute P(z|z_ob)
      double Pz = gaussian(x[1], ptr, {dummy_z,pp->z_abs_err});
  
      // Compute P(proxy|proxy_ob)
      double Pproxy = gaussian(x[2], ptr, {dummy_proxy,pp->proxy_rel_err*dummy_proxy});
  
      // Compute the halo bias
      double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo->Omega_matter(x[1]) : pp->Delta_input;
      double bias_halo = cosmo->bias_halo(mass, sigmaM_interp(mass), x[1], interp_DN(x[1]), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
                
      return bias_halo * P_M__lambda_z * Pz * Pproxy;
              
      };
      vector<vector<double>> integration_limits(3);
      wrapper::cuba::CUBAwrapper CW (integrand, (int)(integration_limits.size()));
    */
     
     
    // Compute the bias    
    vector<double> _bias(pp->scaling_relation->data()->xx().size());
  
    for (size_t i=0; i<pp->scaling_relation->data()->xx().size(); i++) {
  
      dummy_proxy = pp->scaling_relation->data()->xx(i);
      dummy_z = (pp->scaling_relation)->data_model().redshift[i];
  
      // Find the minimum and maximum masses, given the parameters of the scaling relation
      double log_lambda_min = log(dummy_proxy*(1-pp->proxy_rel_err)/proxy_pivot)/log(log_base);
      double log_lambda_max = log(dummy_proxy*(1+pp->proxy_rel_err)/proxy_pivot)/log(log_base);
      double log_f_z_min = log( (pp->scaling_relation)->data_model().fz(dummy_z-pp->z_abs_err, redshift_pivot, cosmo) )/log(log_base);
      double log_f_z_max = log( (pp->scaling_relation)->data_model().fz(dummy_z+pp->z_abs_err, redshift_pivot, cosmo) )/log(log_base);
  
      double logM1 = alpha + beta*log_lambda_min + gamma*log_f_z_min;
      double logM2 = alpha + beta*log_lambda_max + gamma*log_f_z_min;
      double logM3 = alpha + beta*log_lambda_min + gamma*log_f_z_max;
      double logM4 = alpha + beta*log_lambda_max + gamma*log_f_z_max;
  
      double min1 = min(logM1, logM2);
      double min2 = min(min1, logM3);
      double min_logM = min(min2, logM4);
      double max1 = max(logM1, logM2);
      double max2 = max(max1, logM3);
      double max_logM = max(max2, logM4);
  
      // Find the maximum value of the intrinsic scatter
      double s1 = abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
      double s2 = abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
      double s3 = abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );
      double s4 = abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );
  
      double maxs1 = max(s1, s2);
      double maxs2 = max(maxs1, s3);
      double max_scatter_intr = max(maxs2, s4);
  
      // Integrate
      _bias[i] = wrapper::gsl::GSL_integrate_qag(integrand, max(min_logM-3.5*max_scatter_intr,log(Min(Mass_vector)/mass_pivot)/log(log_base)), min(max_logM+3.5*max_scatter_intr,log(Max(Mass_vector)/mass_pivot)/log(log_base)));
  
      /*
	integration_limits[0] = {max(min_logM-3.5*max_scatter_intr,log(Min(Mass_vector)/mass_pivot)/log(log_base)), min(max_logM+3.5*max_scatter_intr,log(Max(Mass_vector)/mass_pivot)/log(log_base))};
	integration_limits[1] = {Min(z_for_interp), Max(z_for_interp)};
	integration_limits[2] = {max(dummy_proxy - 3.5*pp->proxy_rel_err*dummy_proxy, 0.00001), dummy_proxy + 3.5*pp->proxy_rel_err*dummy_proxy};

	_bias[i] = CW.IntegrateVegas(integration_limits,false);
      */
       
    }
  
    bias = Average(_bias);
    return bias;
  
  }
}


// ============================================================================================


std::vector<double> cbl::modelling::angularpk::covariance_Cl_bias_MF_theory(std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  double binsize=l[1]-l[0];
  double nObjects=pp->cluster_redshift.size();
  double shot_noise=pp->fsky*4*cbl::par::pi/nObjects;

  vector<double> Cl = Cl_bias_MF_theory(l, inputs, parameter);
  vector<double> covariance; 

  double beta=parameter[pp->Cpar.size()+2]; 
  double alpha=parameter[pp->Cpar.size()+3]; 

  for (size_t i=0; i<l.size(); i++)
    covariance.emplace_back(2/(2*l[i]+1)/pp->fsky/binsize*pow(beta*Cl[i]+shot_noise*(1+alpha),2));
  
  return covariance; 
  
}


// ============================================================================================

double cbl::modelling::angularpk::window_function_den (const double zz, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz){

  double W_den_val; 
  W_den_val= cosmo->Hubble(zz)*nz->operator()(zz)*pow(cbl::par::cc,-1);     

  return W_den_val;
}
   
// ============================================================================================

double cbl::modelling::angularpk::window_function_rsd (const double zz, const double l, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz, cbl::glob::FuncGrid z_r, cbl::glob::FuncGrid bias){

  const std::vector<int> m{-1,0,1};
  const std::vector<double> Lm{-(l*(l-1))/((2*l-1)*sqrt((2*l-3)*(2*l+1))),(2*l*l+2*l-1)/((2*l-1)*(2*l+3)),-((l+1)*(l+2))/((2*l+3)*sqrt((2*l+5)*(2*l+1)))};
  const std::vector<double> zm{(2*l+1+4*m[0])/(2*l+1)*cosmo->D_C(zz),(2*l+1+4*m[1])/(2*l+1)*cosmo->D_C(zz),(2*l+1+4*m[2])/(2*l+1)*cosmo->D_C(zz)};
  const std::vector<double> z{z_r(zm[0]*cosmo->D_C(zz)),z_r(zm[1]*cosmo->D_C(zz)),z_r(zm[2]*cosmo->D_C(zz))};
   
  double W_rsd_val=0.; 
   
  for(size_t i=0; i< m.size();++i)
    W_rsd_val+= Lm[i]*cosmo->Hubble(z[i])*cosmo->linear_growth_rate(z[i])*nz->operator()(z[i])/bias(z[i]);     

  return W_rsd_val*pow(cbl::par::cc,-1);  
}
   
// ============================================================================================

double cbl::modelling::angularpk::window_function_mag (const double zz, const double l, const double z_max, cbl::glob::FuncGrid mag_bias, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz, cbl::glob::FuncGrid DC ){
   
  const double Lmu=l*(l+1)*pow(l+0.5,-2);
   
  auto integrand = [&zz,&nz,&mag_bias,&DC] (const double z){
    double N_z=nz->operator()(z);
    double bias=mag_bias(z);
    double D_=DC(zz)/DC(z);
    return bias*N_z*(1.-D_);
  };
  double integral = wrapper::gsl::GSL_integrate_cquad(integrand, zz, z_max);
     
  double W_mag_val = 1.5 *Lmu* pow(cosmo->Hubble(0.)*pow(cbl::par::cc,-1),2)* cosmo->Omega_matter(0.) * (1.0 + zz) *cosmo->D_C(zz)* integral;

   
  return W_mag_val;
}

// ============================================================================================


double cbl::modelling::angularpk::bias_inter (const double zz, const std::vector<double> z_center, const std::vector<double> bb)
{  
  double bias;
  auto b_z =make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(z_center,bb , "Linear"));
   
  if (zz <= z_center[0])
    bias=bb[0];
  else if (zz>= z_center[z_center.size()-1])
    bias=bb[bb.size()-1];
  else 
    bias=b_z->operator()(zz);
      
  return bias;
}
   
// ============================================================================================

double cbl::modelling::angularpk::C_ij_l (const double l, const double z_min, const double z_max, cbl::glob::FuncGrid W_i,  cbl::glob::FuncGrid W_j, cbl::glob::FuncGrid HH,cbl::glob::FuncGrid DC, cbl::glob::FuncGrid2D Pk,cbl::glob::FuncGrid func ){
        
  function<double(double)> integrand_limber = [&l, &W_i, &W_j,&HH,&DC,&Pk,&func] (double redshift)
  {
    double dc=DC(redshift), hh=HH(redshift),Wi=W_i(redshift),Wj=W_j(redshift),ff=func(redshift);
    double Pk_=Pk(redshift,(l+0.5)/dc);
    double inv_= 1/(dc*dc*hh);
    return Wi*Wj*ff*ff*Pk_*inv_;
          
  };

  return cbl::par::cc*wrapper::gsl::GSL_integrate_cquad(integrand_limber, z_min,z_max);
 
     
}  
   

// ============================================================================================


double cbl::modelling::angularpk::C_ij_xy (const double l,  const double z_min, const double z_max, cbl::glob::FuncGrid Wx_i,  cbl::glob::FuncGrid Wx_j,cbl::glob::FuncGrid Wy_i,cbl::glob::FuncGrid Wy_j, cbl::glob::FuncGrid HH,cbl::glob::FuncGrid DC, cbl::glob::FuncGrid2D Pk,cbl::glob::FuncGrid func){
     
  function<double(double)> integrand_limber_xy = [&l, &Wx_i, &Wx_j,&Wy_i, &Wy_j,&HH,&DC,&Pk,&func] (double redshift)
  {
    double dc=DC(redshift), hh=HH(redshift),Wxi=Wx_i(redshift),Wxj=Wx_j(redshift),Wyi=Wy_i(redshift),Wyj=Wy_j(redshift),ff=func(redshift);
    double Pk_=Pk(redshift,(l+0.5)/dc);
    double inv_= 1/(dc*dc*hh);
    double C_xx,C_xy,C_yy;
        
    C_xx=Wxi*Wxj*ff*ff*Pk_*inv_;
    C_xy=((Wxi*Wyj)+(Wxj*Wyi))*ff*Pk_*inv_;
    C_yy=Wyi*Wyj*Pk_*inv_;
        
    return C_xx+C_xy+C_yy;
          
  };

  return cbl::par::cc*wrapper::gsl::GSL_integrate_cquad(integrand_limber_xy, z_min,z_max);
      
}  

// ==============================================================================================

std::vector<double> cbl::modelling::angularpk::C_ij_gg (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{ 
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  
  // set the cosmological parameters 
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
    
  vector<double> z_vector= linear_bin_vector(150,pp->z_min,pp->z_max);
  vector<double> comoving_distance, hubble, bias_z;
  for (size_t i=0; i<z_vector.size(); ++i) {
    hubble.emplace_back(cosmo->Hubble(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));  
    bias_z.emplace_back(bias_inter(z_vector[i],pp->z_bias,pp->bias_all));  
  }

  glob::FuncGrid HH(z_vector, hubble, "Spline");
  glob::FuncGrid DC(z_vector, comoving_distance, "Spline");
  glob::FuncGrid bias(z_vector, bias_z, "Linear");
  
  cosmology::PkXi PX(cosmo);
  vector<double> kk = logarithmic_bin_vector(300, 1e-4, 1250.);
  vector<vector<double>> Pk = PX.Pk_matter(kk, "CAMB",false, z_vector, false); 
  glob::FuncGrid2D pk_interp(z_vector,kk, Pk, "Cubic");
  
  vector<vector<double>> WW (pp->nz_all.size());
  for (size_t j=0; j<WW.size(); ++j){
    for (size_t i=0; i<z_vector.size(); ++i) {
      WW[j].emplace_back(window_function_den(z_vector[i],cosmo,pp->nz_all[j]));
    }
  }
  
  vector<glob::FuncGrid> Wij (pp->nz_all.size());
  for (size_t i=0; i<Wij.size(); ++i)
    Wij[i]=glob::FuncGrid(z_vector, WW[i], "Spline");
  
  int n_size=0.5*pp->nz_all.size()*(pp->nz_all.size()-1)+pp->nz_all.size();
  int l_size=l.size()/n_size;
  vector<double> ells(l_size,0.);
  vector<double> Cl;

  for (int i=0; i<l_size;++i)
    ells[i]=l[i*n_size];    
 
  for (size_t i=0; i<ells.size(); ++i){
    for (size_t j=0; j<pp->nz_all.size(); ++j){
      for (size_t k=j; k<pp->nz_all.size();++k)
        Cl.push_back(C_ij_l(ells[i],pp->z_min,pp->z_max,Wij[j],Wij[k],HH, DC, pk_interp, bias));
    }
  }
    
  return Cl;
  
}

// ==============================================================================================

std::vector<double> cbl::modelling::angularpk::C_ij_GG (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{ 
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  
  // set the cosmological parameters 
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
    
  vector<double> z_vector= linear_bin_vector(150,pp->z_min,pp->z_max);
  vector<double> comoving_distance, hubble, bias_z, magbias_z;
  for (size_t i=0; i<z_vector.size(); ++i) {
    hubble.emplace_back(cosmo->Hubble(z_vector[i]));
    comoving_distance.emplace_back(cosmo->D_C(z_vector[i]));  
    bias_z.emplace_back(bias_inter(z_vector[i],pp->z_bias,pp->bias_all));
    magbias_z.emplace_back(bias_inter(z_vector[i],pp->z_bias,pp->mag_bias_all));  
  }

  glob::FuncGrid HH(z_vector, hubble, "Spline");
  glob::FuncGrid DC(z_vector, comoving_distance, "Spline");
  glob::FuncGrid z_DC(comoving_distance, z_vector, "Spline");
  glob::FuncGrid bias(z_vector, bias_z, "Linear");
  glob::FuncGrid magbias(z_vector, magbias_z, "Linear");

  
  int n_size=0.5*pp->nz_all.size()*(pp->nz_all.size()-1)+pp->nz_all.size();
  int l_size=l.size()/n_size;
  vector<double> ells(l_size,0.);
  vector<double> Cl;

  for (int i=0; i<l_size;++i)
    ells[i]=l[i*n_size];

  cosmology::PkXi PX(cosmo);
  vector<double> kk = logarithmic_bin_vector(300, 1e-4, 1250.);
  vector<vector<double>> Pk = PX.Pk_matter(kk, "CAMB",false, z_vector, false); 
  glob::FuncGrid2D pk_interp(z_vector,kk, Pk, "Cubic");
 
  vector<vector<vector<double>>> WW (l_size, vector<vector<double>> (pp->nz_all.size()));
  
  for (size_t i=0; i<ells.size(); ++i){
    for (size_t j=0; j<pp->nz_all.size(); ++j){
      for (size_t k=0; k<z_vector.size(); ++k){
	WW[i][j].emplace_back(window_function_den(z_vector[k],cosmo,pp->nz_all[j])+window_function_rsd (z_vector[k], ells[i], cosmo,pp->nz_all[j], z_DC, bias));  
      }
    }
  }
 
  vector<vector<glob::FuncGrid>> Wij (l_size, vector<glob::FuncGrid>(pp->nz_all.size()));
  for (size_t i=0; i<ells.size(); ++i){
    for (size_t j=0; j<pp->nz_all.size();++j)
      Wij[i][j]=glob::FuncGrid(z_vector, WW[i][j], "Spline");
  }

  vector<vector<vector<double>>> WW_mu (l_size, vector<vector<double>> (pp->nz_all.size()));
  
  for (size_t i=0; i<ells.size(); ++i){
    for (size_t j=0; j<pp->nz_all.size(); ++j){
      for (size_t k=0; k<z_vector.size(); ++k){
	WW_mu[i][j].emplace_back(window_function_mag(z_vector[k],ells[i], pp->z_max, magbias,cosmo,pp->nz_all[j],DC)); 
      }
    }
  }
 
  vector<vector<glob::FuncGrid>> Wij_mu (l_size, vector<glob::FuncGrid>(pp->nz_all.size()));
  for (size_t i=0; i<ells.size(); ++i){
    for (size_t j=0; j<pp->nz_all.size();++j)
      Wij_mu[i][j]=glob::FuncGrid(z_vector, WW_mu[i][j], "Spline");
  }
        
 
  for (size_t i=0; i<ells.size(); ++i){
    for (size_t j=0; j<pp->nz_all.size(); ++j){
      for (size_t k=j; k<pp->nz_all.size();++k)
        Cl.push_back(C_ij_xy(ells[i],pp->z_min,pp->z_max,Wij[i][j],Wij[i][k],Wij_mu[i][j],Wij_mu[i][k],HH, DC, pk_interp, bias));
    }
  }

  return Cl;
  
}



