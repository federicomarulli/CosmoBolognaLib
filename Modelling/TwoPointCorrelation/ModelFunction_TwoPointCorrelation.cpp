/*******************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file
 *  ModelFunction/TwoPointCorrelation/ModelFunction_TwoPointCorrelation.cpp
 *
 *  @brief Global functions to model two-point correlation functions
 *  of any type
 *
 *  This file contains the implementation of the used to model
 *  two-point correlation functions of any type
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "ModelFunction_PowerSpectrum_multipoles.h"
#include "ModelFunction_TwoPointCorrelation.h"

using namespace std;

using namespace cbl;


// ===============================================================================

cbl::glob::FuncGrid cbl::modelling::twopt::Xil_interp (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par)
{
  vector<double> Pkl = modelling::powspec::Pk_l(kk, l, model, parameter, pk_interp, prec, alpha_perp, alpha_par);
  vector<double> rr, Xil;

  cbl::wrapper::fftlog::transform_FFTlog(rr, Xil, 1, kk, Pkl, l, 0, 2.*par::pi, 1);

  cbl::glob::FuncGrid interp(rr, Xil, "Spline");

  return interp;
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::Xi_l (const std::vector<double> rr, const int nmultipoles, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par)
{
  vector<vector<double>> Xil(3);

  for (int i=0; i<nmultipoles; i++) {
    double sign = (i%2==0) ? 1 : -1;
    Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec, alpha_perp, alpha_par).eval_func(rr);
    for (size_t j=0; j<rr.size(); j++)
      Xil[i][j] *= sign;
  }

  return Xil;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::Xi_l (const std::vector<double> rr, const std::vector<int> dataset_order, const std::vector<bool> use_pole, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par)
{
  vector<cbl::glob::FuncGrid> interp_Xil(3);
  vector<double> sign={1., -1., 1.};

  for (size_t i=0; i<3; i++) 
    if (use_pole[i])
      interp_Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec, alpha_perp, alpha_par);

  vector<double> Xil(rr.size());

  for (size_t i=0; i<rr.size(); i++)
    Xil[i] = sign[dataset_order[i]]*interp_Xil[dataset_order[i]](rr[i]);
  
  return Xil;
}


// ============================================================================================


double cbl::modelling::twopt::Xi_polar(const double rad_fid, const double mu_fid, const double alpha_perpendicular, const double alpha_parallel, const std::vector<std::shared_ptr<cbl::glob::FuncGrid>> xi_multipoles)
{
  const double apar2 = alpha_parallel*alpha_parallel;
  const double aperp2 = alpha_perpendicular*alpha_perpendicular;

  double mu_fid_sq = mu_fid*mu_fid;
  const double factor = sqrt(apar2*mu_fid_sq+(1.-mu_fid_sq)*aperp2);
  const double mu_true = mu_fid * alpha_parallel/factor;
  const double mu_true_sq = mu_true*mu_true;
  const double st = rad_fid*factor;

  return xi_multipoles[0]->operator()(st)+
         xi_multipoles[1]->operator()(st)*0.5*(3*mu_true_sq-1)+
	 xi_multipoles[2]->operator()(st)*0.125*(35*mu_true_sq*mu_true_sq-30*mu_true_sq+3);
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::Xi_rppi (const std::vector<double> rp, const std::vector<double> pi, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par)
{
  const int nr=200;
  const size_t nmultipoles = pk_interp.size();
  vector<double> sign = {1., -1., 1., 1, -1};
  vector<double> rr = linear_bin_vector(nr, min(Min(rp), Min(pi))*0.999, sqrt(Max(rp)*Max(rp)+Max(pi)*Max(pi))*1.001);

  vector<glob::FuncGrid> interp_Xil(nmultipoles);

  for (size_t i=0; i<nmultipoles; i++)
    interp_Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec, alpha_perp, alpha_par);

  vector<vector<double>> xi_rppi(rp.size(), vector<double>(pi.size(), 0));

  for (size_t i =0; i<rp.size(); i++)
    for (size_t j =0; j<pi.size(); j++) {
      double s = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/s;
      for (size_t l=0; l<nmultipoles; l++)
        xi_rppi[i][j] += sign[l]*interp_Xil[l](s)*legendre_polynomial (mu, l*2);
    }

  return xi_rppi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::wp_from_Xi_rppi (const std::vector<double> rp, const double pimax, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par) {

  vector<double> pi = linear_bin_vector(100, 1.e-4, pimax*1.001);

  vector<vector<double>> xi_rppi = modelling::twopt::Xi_rppi(rp, pi, model, parameter, pk_interp, prec, alpha_perp, alpha_par);
  vector<double> wp(rp.size());

  for (size_t i=0; i<rp.size(); i++) {
    glob::FuncGrid func(pi, xi_rppi[i], "Spline");
    wp[i] = func.integrate_qag(0., pimax);
  }

  return wp;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::damped_Xi (const std::vector<double> ss, const double bias, const double linear_growth_rate, const double SigmaS, const std::vector<double> kk, const std::shared_ptr<cbl::glob::FuncGrid> PkDM)
{
  vector<vector<double>> pk_terms = modelling::powspec::damped_Pk_terms(kk, linear_growth_rate, SigmaS, PkDM);

  vector<double> xi(ss.size(), 0);

  for (size_t i=0; i<pk_terms.size(); i++) {
    vector<double> xi_term = wrapper::fftlog::transform_FFTlog(ss, 1, kk, pk_terms[i], 0);
    for (size_t j=0; j<ss.size(); j++)
      xi[j] += pow(bias, 2-i)*xi_term[j];

  }

  return xi;
}


//================================================================================================================


double cbl::modelling::twopt::bias_from_scaling_relation(const double alpha, const double beta, const double gamma,  const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, double log_base, double mass_pivot, double proxy_pivot, double redshift_pivot, std::vector<double> scalRel_pars, cbl::glob::FuncGrid DN_interp, std::shared_ptr<void> inputs)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  
  //set the cosmological parameters (contained in scalRel_pars)
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], scalRel_pars[i]);

  // Interpolate sigmaM
  const vector<double> Mass_vector = logarithmic_bin_vector(300, 1.e10, 1.e16);
  cosmology::Sigma SG(cosmo);
  vector<double> sigma2M = SG.sigma2M(Mass_vector, pp->method_Pk, 0., false, pp->output_root, "Linear", 100.);
  vector<double> sigmaM(Mass_vector.size(), 0.);
  for (size_t i = 0; i < sigmaM.size(); i++)
    sigmaM[i] = sqrt(sigma2M[i]);

  glob::FuncGrid sigmaM_interp (Mass_vector, sigmaM, "Spline");
  
  // compute the bias
  double bias = 0;
  cosmology::Bias bb(cosmo);

  if (pp->z_abs_err == -1 && pp->proxy_rel_err == -1) {
         
    vector<double> _bias(pp->scaling_relation->data()->xx().size());
  
    for (size_t i=0; i<pp->scaling_relation->data()->xx().size(); i++) {
       
      double log_lambda = log(pp->scaling_relation->data()->xx(i)/proxy_pivot)/log(log_base);
      double log_fz = log( (pp->scaling_relation)->data_model().fz((pp->scaling_relation)->data_model().redshift[i],redshift_pivot,cosmo) )/log(log_base);
     
      double scatter_intr = scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_fz, scatterz_exp);
     
      double log_mass = (pp->scaling_relation)->likelihood()->get_m_model()->operator()(pp->scaling_relation->data()->xx(i), scalRel_pars) + scatter_intr;
      double mass = pow(log_base, log_mass) * mass_pivot;
  
      double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo->Omega_matter((pp->scaling_relation)->data_model().redshift[i]) : pp->Delta_input;
      double z = (pp->scaling_relation)->data_model().redshift[i];
       
      _bias[i] = bb.bias_halo(mass, sigmaM_interp(mass), z, DN_interp(z), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
       
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

    // Define the integrand
    double dummy_proxy, dummy_z;
    shared_ptr<void> ptr;
  
    auto integrand = [&] (const double x)
    {
      double mass = pow(log_base,x)*mass_pivot;
      // Compute P(M|lambda,z)
      double log_lambda = log(dummy_proxy/proxy_pivot)/log(log_base);
      double log_f_z = log( (pp->scaling_relation)->data_model().fz(dummy_z, redshift_pivot, cosmo) )/log(log_base);
                
      double mean = alpha + beta*log_lambda + gamma*log_f_z;
      double scatter_intr = abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
      double P_M__lambda_z = gaussian(x, ptr, {mean,scatter_intr});
  
      // Compute the halo bias
      double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo->Omega_matter(dummy_z) : pp->Delta_input;
      double bias_halo = bb.bias_halo(mass, sigmaM_interp(mass), dummy_z, DN_interp(dummy_z), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);

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
      double bias_halo = cosmo->bias_halo(mass, sigmaM_interp(mass), x[1], DN_interp(x[1]), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
                
      return bias_halo * P_M__lambda_z * Pz * Pproxy;
              
      };
      vector<vector<double>> integration_limits(3);
      wrapper::cuba::CUBAwrapper CW (integrand, (int)(integration_limits.size()));
    */
     
    // compute the bias    
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
