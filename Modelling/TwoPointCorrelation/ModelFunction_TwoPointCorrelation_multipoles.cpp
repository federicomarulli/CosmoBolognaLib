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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_multipoles.cpp
 *
 *  @brief Functions to model the multipoles of the two-point
 *  correlation function
 *
 *  This file contains the implementation of the functions used to
 *  model the multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, massimilia.romanell2@unibo.it
 */


#include "ModelFunction_TwoPointCorrelation.h"
#include "ModelFunction_TwoPointCorrelation_multipoles.h"
#include "ModelFunction_NumberCounts.h"
#include "BAO.h"
#include "PkXiNonLinear.h"

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // input parameters
  vector<shared_ptr<glob::FuncGrid>> pk_interp(2);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_dispersion(1);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_Scoccimarro_CPT(3);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_TNS_CPT(17);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_eTNS_CPT(25);
  std::vector<double> Xi_ll;

  if (pp->Pk_mu_model=="dispersion_dewiggled") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_NW;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[2], parameter[3], parameter[4]/pp->sigma8_z, parameter[5]/pp->sigma8_z, parameter[6] }, pk_interp, pp->prec, parameter[0], parameter[1]);
  }
  
  else if (pp->Pk_mu_model=="dispersion_modecoupling") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk1loop;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[2]/pp->sigma8_z, parameter[3]/pp->sigma8_z, parameter[4], parameter[5] }, pk_interp, pp->prec, parameter[0], parameter[1]);
  }
  
  else if (pp->Pk_mu_model=="dispersion_Gauss" || pp->Pk_mu_model=="dispersion_Lorentz") {
    pk_interp_dispersion[0] = pp->func_Pk;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2] }, pk_interp_dispersion, pp->prec, parameter[3], parameter[4]);
  }
  
  else if (pp->Pk_mu_model=="Scoccimarro_Gauss" || pp->Pk_mu_model=="Scoccimarro_Lorentz") {
    pk_interp_Scoccimarro_CPT[0] = pp->func_Pk_DeltaDelta;
    pk_interp_Scoccimarro_CPT[1] = pp->func_Pk_DeltaTheta;
    pk_interp_Scoccimarro_CPT[2] = pp->func_Pk_ThetaTheta;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2] }, pk_interp_Scoccimarro_CPT, pp->prec, parameter[3], parameter[4]);
  }
  
  else if (pp->Pk_mu_model=="Scoccimarro_Pezzotta_Gauss" || pp->Pk_mu_model=="Scoccimarro_Pezzotta_Lorentz") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_nonlin;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2], parameter[3], parameter[4] }, pk_interp, pp->prec, parameter[5], parameter[6]);
  }
  
  else if (pp->Pk_mu_model=="Scoccimarro_Bel_Gauss" || pp->Pk_mu_model=="Scoccimarro_Bel_Lorentz") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_nonlin;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], parameter[7] }, pk_interp, pp->prec, parameter[8], parameter[9]);
  }
  
  else if (pp->Pk_mu_model=="TNS_Gauss" || pp->Pk_mu_model=="TNS_Lorentz") {
    pk_interp_TNS_CPT[0]  = pp->func_Pk_DeltaDelta;
    pk_interp_TNS_CPT[1]  = pp->func_Pk_DeltaTheta;
    pk_interp_TNS_CPT[2]  = pp->func_Pk_ThetaTheta;
    pk_interp_TNS_CPT[3]  = pp->func_Pk_A11;
    pk_interp_TNS_CPT[4]  = pp->func_Pk_A12;
    pk_interp_TNS_CPT[5]  = pp->func_Pk_A22;
    pk_interp_TNS_CPT[6]  = pp->func_Pk_A23;
    pk_interp_TNS_CPT[7]  = pp->func_Pk_A33;
    pk_interp_TNS_CPT[8]  = pp->func_Pk_B12;
    pk_interp_TNS_CPT[9]  = pp->func_Pk_B13;
    pk_interp_TNS_CPT[10] = pp->func_Pk_B14;
    pk_interp_TNS_CPT[11] = pp->func_Pk_B22;
    pk_interp_TNS_CPT[12] = pp->func_Pk_B23;
    pk_interp_TNS_CPT[13] = pp->func_Pk_B24;
    pk_interp_TNS_CPT[14] = pp->func_Pk_B33;
    pk_interp_TNS_CPT[15] = pp->func_Pk_B34;
    pk_interp_TNS_CPT[16] = pp->func_Pk_B44;

    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2] }, pk_interp_TNS_CPT, pp->prec, parameter[3], parameter[4]);
  }
  
  else if (pp->Pk_mu_model=="eTNS_Gauss" || pp->Pk_mu_model=="eTNS_Lorentz") {
    pk_interp_eTNS_CPT[0]  = pp->func_Pk_DeltaDelta;
    pk_interp_eTNS_CPT[1]  = pp->func_Pk_DeltaTheta;
    pk_interp_eTNS_CPT[2]  = pp->func_Pk_ThetaTheta;
    pk_interp_eTNS_CPT[3]  = pp->func_Pk_A11;
    pk_interp_eTNS_CPT[4]  = pp->func_Pk_A12;
    pk_interp_eTNS_CPT[5]  = pp->func_Pk_A22;
    pk_interp_eTNS_CPT[6]  = pp->func_Pk_A23;
    pk_interp_eTNS_CPT[7]  = pp->func_Pk_A33;
    pk_interp_eTNS_CPT[8]  = pp->func_Pk_B12;
    pk_interp_eTNS_CPT[9]  = pp->func_Pk_B13;
    pk_interp_eTNS_CPT[10] = pp->func_Pk_B14;
    pk_interp_eTNS_CPT[11] = pp->func_Pk_B22;
    pk_interp_eTNS_CPT[12] = pp->func_Pk_B23;
    pk_interp_eTNS_CPT[13] = pp->func_Pk_B24;
    pk_interp_eTNS_CPT[14] = pp->func_Pk_B33;
    pk_interp_eTNS_CPT[15] = pp->func_Pk_B34;
    pk_interp_eTNS_CPT[16] = pp->func_Pk_B44;
    pk_interp_eTNS_CPT[17]  = pp->func_Pk_b2d;
    pk_interp_eTNS_CPT[18]  = pp->func_Pk_b2v;
    pk_interp_eTNS_CPT[19] = pp->func_Pk_b22;
    pk_interp_eTNS_CPT[20] = pp->func_Pk_bs2d;
    pk_interp_eTNS_CPT[21] = pp->func_Pk_bs2v;
    pk_interp_eTNS_CPT[22] = pp->func_Pk_b2s2;
    pk_interp_eTNS_CPT[23] = pp->func_Pk_bs22;
    pk_interp_eTNS_CPT[24] = pp->func_sigma32Pklin;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2]/pp->sigma8_z, parameter[3], 0. }, pk_interp_eTNS_CPT, pp->prec, parameter[4], parameter[5]);
  }

  else ErrorCBL("the chosen model ("+pp->Pk_mu_model+") is not currently implemented!", "xiMultipoles", "ModelFunction_TwoPointCorrelation_multipoles.cpp");

  return Xi_ll;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles_sigma8_bias (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // input parameters

  // sigma8(z)
  double sigma8 = parameter[0];

  // bias
  double bias = parameter[1];

  // AP parameter that contains the distance information
  double alpha_perpendicular = 1.;

  // AP parameter that contains the distance information
  double alpha_parallel = 1.;

  // f(z)*sigma8(z)
  double fsigma8 = pp->linear_growth_rate_z*sigma8;

  // bias(z)*sigma8(z)
  double bsigma8 = bias*sigma8;

  // streaming scale
  double SigmaS = 0.;

  return Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {pp->sigmaNL_perp, pp->sigmaNL_par, fsigma8/pp->sigma8_z, bsigma8/pp->sigma8_z, SigmaS}, {pp->func_Pk, pp->func_Pk_NW}, pp->prec, alpha_perpendicular, alpha_parallel);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  vector<double> Xi(rad.size());

  // input parameters

  // AP parameter that contains the distance information
  double alpha_perpendicular = parameter[0];

  // AP parameter that contains the distance information
  double alpha_parallel = parameter[1];

  for (size_t i=0; i<rad.size(); i++) {
    const int mult = pp->dataset_order[i];
    const int ell = 2*mult;
    const double rr = rad[i];
    const double fact = (2*ell+1);

    const double B = parameter[2+mult];
    const double A0 = parameter[4+mult];
    const double A1 = parameter[6+mult];
    const double A2 = parameter[8+mult];

    auto integrand = [&] (const double mu_fid)
    {
      return twopt::Xi_polar(rr, mu_fid, alpha_perpendicular, alpha_parallel, pp->func_multipoles)*cbl::legendre_polynomial(mu_fid, ell);
    };

    Xi[i] = B*fact*wrapper::gsl::GSL_integrate_qag(integrand, 0, 1)+A0+A1/rr+A2/(rr*rr);
  }

  return Xi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles_linear_theoretical_bias_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
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

  if(pp->NL)
    cbl::ErrorCBL("This model works only with linear Pk!", "xi0_linear_theoretical_bias_BAO", "ModelFunction_TwoPointCorrelation1D_monopole.cpp",cbl::glob::ExitCode::_error_);
  
  cosmology::PkXi PX(cosmo);
  // compute linear Pk at redshift zero
  std::vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  std::vector<double> Pk0=PX.Pk_matter(kk, pp->method_Pk, false, 0., false, pp->output_root, pp->norm, pp->k_min, pp->k_max);

  double redshift= cbl::Average(pp->cluster_redshift);

  // compute the bias   //use this if you want to estimate bias directly from the sample
  double bias;
  cosmology::Bias bb(cosmo);
  if(pp->bias_from_sample){
    vector<double> mass_grid = logarithmic_bin_vector(pp->mass.size()/10, Min(pp->mass), Max(pp->mass));
    bias = bb.bias_eff_mass(pp->mass, mass_grid, pp->cluster_redshift, pp->model_bias, pp->method_Pk, "mean_bias", true, pp->output_root, pp->Delta)[0];
    PX.remove_output_Pk_tables(pp->method_Pk, false, 0., pp->output_root);
  }
  
  else{
    
    std::vector<double> mass_vector=logarithmic_bin_vector(30, Min(pp->mass), Max(pp->mass));
    std::vector<double> selection_function;
    for(size_t i=0; i<mass_vector.size(); i++)
      selection_function.emplace_back(1.); 
    cbl::glob::FuncGrid interp_SF(mass_vector, selection_function, pp->interpType, cbl::BinType::_linear_);
    
    const std::vector<cbl::glob::FuncGrid> interp_sigmaM_dlnsigmaM = cbl::modelling::numbercounts::sigmaM_dlnsigmaM (mass_vector, cosmo, kk, Pk0, pp->interpType, pp->k_max);

    double Delta_crit = (pp->isDelta_critical) ? pp->Delta_input : pp->Delta_input*cosmo->Omega_matter(redshift);

    bias=bb.bias_eff_selection_function(interp_sigmaM_dlnsigmaM[0], interp_sigmaM_dlnsigmaM[1], interp_SF, Min(pp->mass), Max(pp->mass), {redshift}, pp->model_bias, pp->model_MF, pp->method_Pk, 1., false, pp->output_root, Delta_crit, -1., pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, par::defaultString, true)[0];
    
  }
  parameter[parameter.size()-1] = bias;  //bias
  
  /// alpha and geometrical distorsion
  cbl::cosmology::BAO bao(cosmo);
  const double alpha_par = pp->HHfid/cosmo->Hubble(redshift)*pp->rsfid/bao.rs(pp->method_Pk);
  const double alpha_perp = cosmo->D_A(redshift)/pp->DAfid*pp->rsfid/bao.rs(pp->method_Pk);
  const double alpha = cosmo->D_V(redshift)/pp->DVfid*pp->rsfid/bao.rs(pp->method_Pk);
 
  std::vector<double> new_rad=rad;
  for(size_t i=0; i<new_rad.size(); i++)
    new_rad[i]*=alpha;
  
  std::vector<double> xi(rad.size(), 0);
  
  //compute PkNW (without BAO)
  std::vector<double> PkNW = PX.Pk_matter(kk, "EisensteinHu", false, redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
  //compute and interpolate Pk with IR-resummation
  vector<double> Pk(kk.size(), 0);
  double DN=cosmo->DN(redshift);
  for (size_t i=0; i<kk.size(); i++) 
    Pk[i] = PkNW[i]*(1.+(Pk0[i]*DN*DN/PkNW[i]-1.)*exp(-0.5*pow(kk[i]*sigmaNL, 2)));
  cbl::glob::FuncGrid pk_interp(kk, Pk, "Spline");
  auto pk_interp_ptr = std::make_shared<cbl::glob::FuncGrid>(pk_interp);

  //compute the 2pcf signal
  //damping due to photometric errors
  if(pp->z_abs_err!=0){
    // damping scale
    const double SigmaS = pp->z_abs_err*par::cc/cosmo->Hubble(redshift);
    //in redshift space, return also multipoles
    //we use dispersion model with f=1 and bias=bias/f, because we do not want linear growth rate in the damping exponential, but only in the Kaiser factor.
    if(pp->redshift_space){
      std::vector<double> xi_l=Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {cosmo->linear_growth_rate(redshift), bias, SigmaS/cosmo->linear_growth_rate(redshift)}, {pk_interp_ptr}, pp->prec, alpha_perp, alpha_par);

      return xi_l;
    }
    //in real space, return only monopole
    else    //linear growth rate = 0 in real
      return modelling::twopt::damped_Xi(new_rad, bias, 0., SigmaS, kk, pk_interp_ptr); 
  }
  //no photometric errors
  else{
    //in redshift space, return also multipoles
    //we put Sigma=0, because there is not exponential damping, but only Kaiser enhancement
    if(pp->redshift_space)
      return Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {cosmo->linear_growth_rate(redshift), bias, 0.}, {pk_interp_ptr}, pp->prec, alpha_perp, alpha_par);
    //in real space, return only monopole
    else{
      xi = wrapper::fftlog::transform_FFTlog(new_rad, 1, kk, Pk, 0);
      double fact = bias*bias;
      for (size_t i =0; i<xi.size(); i++)
	xi[i] = fact*xi[i]; 
      return xi; 
    }   
  }
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xil_damped_scaling_relation_sigmaz_cosmology (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);

  // damping scale
  const double SigmaS = par::cc*parameter[parameter.size()-2]/cosmo->Hubble(pp->redshift);
  
  // scaling relation parameters
  vector<double> scalRel_pars;
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    scalRel_pars.push_back(parameter[i]);
  for (size_t i=pp->Cpar.size(); i<parameter.size()-2; i++)
    scalRel_pars.emplace_back(parameter[i]);

  const double alpha = scalRel_pars[scalRel_pars.size()-8];
  const double beta = scalRel_pars[scalRel_pars.size()-7];
  const double gamma = scalRel_pars[scalRel_pars.size()-6];
  const double scatter0 = scalRel_pars[scalRel_pars.size()-5];
  const double scatterM = scalRel_pars[scalRel_pars.size()-4];
  const double scatterM_exp = scalRel_pars[scalRel_pars.size()-3];
  const double scatterz = scalRel_pars[scalRel_pars.size()-2];
  const double scatterz_exp = scalRel_pars[scalRel_pars.size()-1];

  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the effective bias from the scaling relation
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  const vector<double> z_for_DN = linear_bin_vector(100, 0.01, Max((pp->scaling_relation)->data_model().redshift));
  vector<double> DN (z_for_DN.size(), 0.);
  for (size_t i=0; i<z_for_DN.size(); i++)
    DN[i] = cosmo->DN(z_for_DN[i]);
  glob::FuncGrid DN_interp (z_for_DN, DN, "Spline");
  
  // Compute the bias
  double log_base = (pp->scaling_relation)->data_model().log_base;
  double mass_pivot = (pp->scaling_relation)->data_model().mass_pivot;
  double proxy_pivot = (pp->scaling_relation)->data_model().proxy_pivot;
  double redshift_pivot = (pp->scaling_relation)->data_model().redshift_pivot;

  double bias = bias_from_scaling_relation(alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, log_base, mass_pivot, proxy_pivot, redshift_pivot, scalRel_pars, DN_interp, inputs);

  // set the value of the bias
  parameter[parameter.size()-1] = bias;
    
  double alpha_perp=1.;
  double alpha_par=1.;

  
  // compute the power spectrum
  cosmology::PkXi PX(cosmo);
  // compute linear Pk at redshift zero
  std::vector<double> kk = logarithmic_bin_vector(300, pp->k_min, pp->k_max);
  std::vector<double> Xi_ll;

  if (pp->Pk_mu_model=="dispersion_Gauss" || pp->Pk_mu_model=="dispersion_Lorentz") {
    std::vector<double> Pk=PX.Pk_matter(kk, pp->method_Pk, false, pp->redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max);
    shared_ptr<glob::FuncGrid> Pk_interp = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, Pk, "Spline"));

    Xi_ll=Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {cosmo->linear_growth_rate(pp->redshift), bias, SigmaS/cosmo->linear_growth_rate(pp->redshift)}, {Pk_interp}, pp->prec, alpha_perp, alpha_par);
  }
  
  else if (pp->Pk_mu_model=="Scoccimarro_Gauss" || pp->Pk_mu_model=="Scoccimarro_Lorentz") {

    vector<shared_ptr<glob::FuncGrid>> pk_interp_Scoccimarro_CPT(3);
    cbl::cosmology::PkXiNonLinear PXN(cosmo);
    vector<vector<double>> Pk_terms = PXN.Pk_TNS_dd_dt_tt(kk, pp->method_Pk, pp->redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec);
    
    pk_interp_Scoccimarro_CPT[0] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_terms[0], "Spline"));
    pk_interp_Scoccimarro_CPT[1] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_terms[1], "Spline"));
    pk_interp_Scoccimarro_CPT[2] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_terms[2], "Spline"));
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {cosmo->linear_growth_rate(pp->redshift), bias, SigmaS/cosmo->linear_growth_rate(pp->redshift) }, pk_interp_Scoccimarro_CPT, pp->prec, alpha_perp, alpha_par);
  }

  else if (pp->Pk_mu_model=="TNS_Gauss" || pp->Pk_mu_model=="TNS_Lorentz") {
    vector<shared_ptr<glob::FuncGrid>> pk_interp_TNS_CPT(17);

    cbl::cosmology::PkXiNonLinear PXN(cosmo);
    vector<vector<double>> Pk_terms = PXN.Pk_TNS_dd_dt_tt(kk, pp->method_Pk, pp->redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec);
    vector<vector<double>> Pk_AB = PXN.Pk_TNS_AB_terms_1loop(kk, pp->method_Pk, pp->redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec);
    
    pk_interp_TNS_CPT[0]  = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_terms[0], "Spline")); 
    pk_interp_TNS_CPT[1] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_terms[1], "Spline"));
    pk_interp_TNS_CPT[2] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_terms[2], "Spline"));
    pk_interp_TNS_CPT[3] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[0], "Spline"));
    pk_interp_TNS_CPT[4] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[1], "Spline"));
    pk_interp_TNS_CPT[5] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[2], "Spline"));
    pk_interp_TNS_CPT[6] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[3], "Spline"));
    pk_interp_TNS_CPT[7] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[4], "Spline"));
    pk_interp_TNS_CPT[8] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[5], "Spline"));
    pk_interp_TNS_CPT[9] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[6], "Spline"));
    pk_interp_TNS_CPT[10] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[7], "Spline"));
    pk_interp_TNS_CPT[11] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[8], "Spline"));
    pk_interp_TNS_CPT[12] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[9], "Spline"));
    pk_interp_TNS_CPT[13] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[10], "Spline"));
    pk_interp_TNS_CPT[14] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[11], "Spline"));
    pk_interp_TNS_CPT[15] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[12], "Spline"));
    pk_interp_TNS_CPT[16] = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk, Pk_AB[13], "Spline"));

    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {cosmo->linear_growth_rate(pp->redshift), bias, SigmaS/cosmo->linear_growth_rate(pp->redshift) }, pk_interp_TNS_CPT, pp->prec, alpha_perp, alpha_par);
  }

    else ErrorCBL("the chosen model ("+pp->Pk_mu_model+") is not currently implemented!", "xil_damped_scaling_relation_sigmaz_cosmology", "ModelFunction_TwoPointCorrelation_multipoles.cpp");
    

  return Xi_ll;

}



