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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation1D_monopole.cpp
 *
 *  @brief Functions to model the monopole of the two-point
 *  correlation function
 *
 *  This file contains the implementation of the functions used to
 *  model the monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_TwoPointCorrelation1D_monopole.h"
#include "HaloProfile.h"

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_BAO_sigmaNL (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // input parameters

  // parameter that control the BAO shape
  double sigmaNL = parameter[0];

  // AP parameter that contains the distance information
  double alpha = parameter[1];

  vector<double> new_rad;

  for (size_t i=0; i<rad.size(); i++)
    new_rad.push_back(rad[i]*alpha);

  // compute the 2pcf signal
  
  vector<double> Pk(pp->kk.size(), 0);

  vector<double> Pklin = pp->func_Pk->y();
  vector<double> PkNW = pp->func_Pk_NW->y();

  for (size_t i =0; i<pp->kk.size(); i++) {
    Pk[i] = PkNW[i]*(1.+(Pklin[i]/PkNW[i]-1.)*exp(-0.5*pow(pp->kk[i]*sigmaNL, 2)));
  }

  vector<double> xi = wrapper::fftlog::transform_FFTlog(new_rad, 1, pp->kk, Pk, 0);
  
  // return the monopole of the two-point correlation function

  for (size_t i =0; i<xi.size(); i++) {
    
    double poly = 0;
    for (int j = 0;j<pp->poly_order; j++)
      poly += parameter[j+3]*pow(rad[i], -j);

    xi[i] = pow(parameter[2], 2)*xi[i]+poly;
  }

  return xi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // input parameters

  // AP parameter that contains the distance information
  double alpha = parameter[0];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[1];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[2];

  // return the redshift-space monopole of the two-point correlation function

  vector<double> xi(rad.size(), 0);

  for (size_t i =0; i<xi.size(); i++) {
    
    double poly = 0;
    for (int j = 0;j<pp->poly_order; j++)
      poly += parameter[j+3]*pow(rad[i], -j);

    xi[i] = xi_ratio(fsigma8, bsigma8)*pow(bsigma8/pp->sigma8_z, 2)*pp->func_xi->operator()(rad[i]*alpha)+poly;
  }

  return xi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_polynomial_LinearPoint (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  vector<double> model(rad.size(), 0);

  vector<double> poly_params(parameter.size()-3,0);
  for (size_t i=0; i<poly_params.size(); i++)
    poly_params[i] = parameter[i+3];

  for (size_t i=0; i<model.size(); i++)
    model[i] = wrapper::gsl::GSL_polynomial_eval(rad[i], NULL, poly_params);

  vector<double> deriv_coeff(pp->poly_order-1, 0);
  vector<vector<double>> roots;
  
  for (size_t i=1; i<poly_params.size(); i++)
    deriv_coeff[i-1] = i*poly_params[i];

  vector<double> D1_coeff(pp->poly_order-1, 0);
  for (size_t i=1; i<poly_params.size(); i++)
    D1_coeff[i-1] = i*poly_params[i];

  auto model_derivative = [&] (double rad)
    {
      double mm = 0.;
      for (size_t i=0; i<D1_coeff.size(); i++)
	mm += D1_coeff[i]*pow(rad, i);
      return mm;
    };
  
  double xmin=95, xmax=105;

  bool end=false;
  while (!end) {
    parameter[0] = wrapper::gsl::GSL_root_brent(model_derivative, 0., xmin, xmax, 1.e-10);
    if (fabs(parameter[0]-xmin)<0.1)
      xmin -=2;
    else if (fabs(parameter[0]-xmax)<0.1)
      xmax +=2;
    else
      end=true;
    if ((xmin<80) || (xmax>120)) {
      end=true;
      parameter[0]=0;
    }
  }

  parameter[1]=0;

  if ((parameter[0]>xmin) && (parameter[0]<xmax)) {
    xmin = parameter[0]-22; xmax = parameter[0]-2;
    end = false;
    while (!end) {
      parameter[1] = wrapper::gsl::GSL_root_brent(model_derivative, 0., xmin, xmax, 1.e-10);

      if (fabs(parameter[1]-xmin)<0.1)
	xmin -=2;
      else if (fabs(parameter[1]-xmax)<0.1) {
	end=true;
	parameter[0]=0;
	parameter[1]=0;
      }
      else{
	end=true;
      }

      if (xmin<60) {
	end=true;
	parameter[0]=0;
	parameter[1]=0;
      }
    }
  }

  parameter[2] = (parameter[0]+parameter[1])*0.5;

  return model;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_LinearPoint (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // input parameters

  // ap parameter that contains the distance information
  double alpha = parameter[3];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[4];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[5];

  // return the redshift-space monopole of the two-point correlation function

  vector<double> model(rad.size(), 0);

  for (size_t i =0; i<model.size(); i++) {
    
    double poly = 0;
    for (int j = 0;j<pp->poly_order; j++)
      poly += parameter[j+6]*pow(rad[i], -j);

    model[i] = xi_ratio(fsigma8, bsigma8)*pow(bsigma8/pp->sigma8_z, 2)*pp->func_xi->operator()(rad[i]*alpha)+poly;

  }

  // Linear point section
  auto model_derivative = [&] (double rr)
    { 
      double poly_der = 0.;
      for (int j = 1;j<pp->poly_order; j++) {
	poly_der += -j*parameter[j+6]*pow(rr, -j-1);
      }
      return xi_ratio(fsigma8, bsigma8)*pow(bsigma8/pp->sigma8_z, 2)*pp->func_xi->D1v(rr*alpha)+poly_der;
    };

  double xmin = 95., xmax = 105.;

  bool end = false;
  while (!end) {
    parameter[0] = wrapper::gsl::GSL_root_brent(model_derivative, 0., xmin, xmax, 1.e-10);
    if (fabs(parameter[0]-xmin)<0.1)
      xmin -= 2;
    else if (fabs(parameter[0]-xmax)<0.1)
      xmax += 2;
    else
      end=true;
    if ((xmin<70) || (xmax>160)) {
      end = true;
      parameter[0] = 0;
    }
  }

  parameter[1] = 0;
  if ((parameter[0]>xmin) && (parameter[0]<xmax)) {
    parameter[1] = wrapper::gsl::GSL_root_brent(model_derivative, 0., 40, parameter[0]-2, 1.e-10);
    if (parameter[0]-parameter[1]<2.1) {
      parameter[1] = 0;
      parameter[0] = 0;
    }
  }

  parameter[2] = (parameter[0]+parameter[1])*0.5;
  
  return model;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_sigma8_bias (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);


  // input likelihood parameters

  // sigma8(z)
  double sigma8 = parameter[0];

  // bias
  double bias = parameter[1];

  
  // fixed parameters 
  
  // f(z)*sigma8(z)
  double fsigma8 = pp->linear_growth_rate_z*sigma8;
  
  
  // return the redshift-space monopole of the two-point correlation function

  vector<double> xi(rad.size(), 0);

  const double fact = xi_ratio(fsigma8, bias*sigma8)*pow(bias, 2)*pow(sigma8/pp->sigma8_z, 2);
  
  for (size_t i =0; i<xi.size(); i++) 
    xi[i] = fact*pp->func_xi->operator()(rad[i]);

  return xi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_cosmology (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  
  // input parameters

  // AP parameter that contains the distance information
  double alpha = parameter[0];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[1];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[2];
  
  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<parameter.size(); ++i)
    pp->cosmology->set_parameter(pp->Cpar[i], parameter[i]);

  // return the redshift-space monopole of the two-point correlation function

  vector<double> xi(rad.size(), 0);

  for (size_t i =0; i<xi.size(); i++) {
    
    double poly = 0.;
    
    for (int j = 0;j<pp->poly_order; j++)
      poly += parameter[j+3]*pow(rad[i], -j);

    xi[i] = xi_ratio(fsigma8, bsigma8)*pp->cosmology->xi_matter(rad[i]*alpha, pp->method_Pk, pp->NL, pp->redshift, true, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par)/pow(pp->sigma8_z, 2)+poly;
  }
  if (!pp->store_output)
    pp->cosmology->remove_output_Pk_tables(pp->method_Pk, pp->NL, pp->redshift, pp->output_root);
  
  return xi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_bias_cosmology (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  cosmology::Cosmology cosmo = *pp->cosmology;

  
  // input likelihood parameters

  // bias
  double bias = parameter[0];

  
  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i+1]);

  
  // fixed parameters 
  
  /// alpha
  double alpha = cosmo.D_V(pp->redshift)/pp->DVfid;

  // return the redshift-space monopole of the two-point correlation function
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= alpha;

  return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->redshift, false, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_damped_scaling_relation_sigmaz (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // input likelihood parameters

  // scaling relation
  const double M0 = parameter[0];
  const double slope = parameter[1];
  const double scatter = parameter[2];

  // damping scale
  const double SigmaS = par::cc*parameter[3]/pp->HHfid;

  // bias computation

  vector<double> mass(pp->cluster_mass_proxy->ndata(), 0), _bias(pp->cluster_mass_proxy->ndata());

  for (int i=0; i<pp->cluster_mass_proxy->ndata(); i++) {

    bool isNan = true;
    double log10_proxy = 0;

    while (isNan) {
      log10_proxy = log10(pp->gau_ran->operator()()*pp->cluster_mass_proxy->error(i)+pp->cluster_mass_proxy->data(i));
      isNan = (log10_proxy!=log10_proxy);
    }
    double log10_mass = 14.+M0+slope*log10_proxy+pp->gau_ran->operator()()*(scatter);

    mass[i] = pow(10, log10_mass);

    _bias[i] = pp->cosmology->bias_halo(mass[i], pp->func_sigma->operator()(mass[i]), pp->cluster_mass_proxy->xx(i), pp->model_bias, false, par::defaultString, "Linear", pp->Delta);
  }

  const double bias = Average(_bias);

  // return the value of the bias
  parameter[4] = bias;

  return modelling::twopt::damped_Xi(rad, bias, pp->linear_growth_rate_z, SigmaS, pp->kk, pp->func_Pk);
}



// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_damped_scaling_relation_sigmaz_cosmology (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;

  // set the cosmological parameters
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  auto cosmo_ptr = std::make_shared<cbl::cosmology::Cosmology>(cosmo);

  // damping scale
  const double SigmaS = par::cc*parameter[parameter.size()-2]/cosmo.HH(pp->redshift);
  
  // scaling relation parameters
  std::vector<double> scalRel_pars;
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

  // compute the power spectrum
  vector<double> Pk = cosmo.Pk_matter(pp->kk, pp->method_Pk, pp->NL, pp->redshift, false, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par);

  std::shared_ptr<glob::FuncGrid> Pk_interp = make_shared<glob::FuncGrid>(glob::FuncGrid(pp->kk, Pk, "Spline"));
  
  
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Derive the effective bias from the scaling relation
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  // Interpolate sigmaM
  const std::vector<double> Mass_vector = cbl::logarithmic_bin_vector(300, 1.e10, 1.e16);
  std::vector<double> sigmaM (Mass_vector.size(), 0.);
  for (size_t i=0; i<sigmaM.size(); i++)
    sigmaM[i] = sqrt( cosmo.sigma2M(Mass_vector[i], pp->method_Pk, 0., false, "test", "Linear", 100.) );

  cbl::glob::FuncGrid sigmaM_interp (Mass_vector, sigmaM, "Spline");

  // Compute the bias
  double log_base = (pp->scaling_relation)->data_model().log_base;
  double mass_pivot = (pp->scaling_relation)->data_model().mass_pivot;
  double proxy_pivot = (pp->scaling_relation)->data_model().proxy_pivot;
  double redshift_pivot = (pp->scaling_relation)->data_model().redshift_pivot;
  
  double bias = 0;

  if (pp->z_abs_err == -1 && pp->proxy_rel_err == -1) {

    // Interpolate the normalised amplitude of the growing mode
    const std::vector<double> z_for_DN = cbl::linear_bin_vector(100, 0.0001, cbl::Max((pp->scaling_relation)->data_model().redshift));
    std::vector<double> DN (z_for_DN.size(), 0.);
    for (size_t i=0; i<z_for_DN.size(); i++)
      DN[i] = cosmo.DN(z_for_DN[i]);
    cbl::glob::FuncGrid DN_interp (z_for_DN, DN, "Spline");
    
    vector<double> _bias(pp->scaling_relation->data()->xx().size());

    for (size_t i=0; i<pp->scaling_relation->data()->xx().size(); i++) {
      
      double log_lambda = log(pp->scaling_relation->data()->xx(i)/proxy_pivot)/log(log_base);
      double log_fz = log( (pp->scaling_relation)->data_model().fz((pp->scaling_relation)->data_model().redshift[i],redshift_pivot,cosmo_ptr) )/log(log_base);
    
      double scatter_intr = scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_fz, scatterz_exp);
    
      double log_mass = (pp->scaling_relation)->likelihood()->get_m_model()->operator()(pp->scaling_relation->data()->xx(i), scalRel_pars) + scatter_intr;
      double mass = pow(log_base, log_mass) * mass_pivot;

      double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo.OmegaM((pp->scaling_relation)->data_model().redshift[i]) : pp->Delta_input;
      double z = (pp->scaling_relation)->data_model().redshift[i];
      
      _bias[i] = cosmo.bias_halo(mass, sigmaM_interp(mass), z, DN_interp(z), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
      
    }

    bias = Average(_bias);

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

    // Interpolate the normalised amplitude of the growing mode
    const std::vector<double> z_for_DN = cbl::linear_bin_vector(100, 0.0001, cbl::Max((pp->scaling_relation)->data_model().redshift)+5.*pp->z_abs_err);
    std::vector<double> DN (z_for_DN.size(), 0.);
    for (size_t i=0; i<z_for_DN.size(); i++)
      DN[i] = cosmo.DN(z_for_DN[i]);
    cbl::glob::FuncGrid DN_interp (z_for_DN, DN, "Spline");

    // Define the integrand
    double dummy_proxy, dummy_z;
    std::shared_ptr<void> ptr;

    auto integrand = [&] (const double x)
		     {
		       double mass = pow(log_base,x)*mass_pivot;
		       
		       // Compute P(M|lambda,z)
		       double log_lambda = log(dummy_proxy/proxy_pivot)/log(log_base);
		       double log_f_z = log( (pp->scaling_relation)->data_model().fz(dummy_z, redshift_pivot, cosmo_ptr) )/log(log_base);
		       
		       double mean = alpha + beta*log_lambda + gamma*log_f_z;
		       double scatter_intr = std::abs(scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp));
		       double P_M__lambda_z = cbl::gaussian(x, ptr, {mean,scatter_intr});

		       // Compute the halo bias
		       double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo.OmegaM(dummy_z) : pp->Delta_input;
		       double bias_halo = cosmo.bias_halo(mass, sigmaM_interp(mass), dummy_z, DN_interp(dummy_z), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
		       
		       return bias_halo * P_M__lambda_z;
		     }; 
    /*
    auto integrand = [&] (const std::vector<double> x)
		     {

		       double mass = pow(log_base,x[0])*mass_pivot;
		       
		       // Compute P(M|lambda,z)
		       double log_lambda = log(x[2]/proxy_pivot)/log(log_base);
		       double log_f_z = log( (pp->scaling_relation)->data_model().fz(x[1], redshift_pivot, cosmo_ptr) )/log(log_base);
		       
		       double mean = alpha + beta*log_lambda + gamma*log_f_z;
		       double scatter_intr = scatter0 + scatterM*pow(log_lambda, scatterM_exp) + scatterz*pow(log_f_z, scatterz_exp);
		       double P_M__lambda_z = (cbl::gaussian(x[0], ptr, {mean,scatter_intr}));

		       // Compute P(z|z_ob)
		       double Pz = cbl::gaussian(x[1], ptr, {dummy_z,pp->z_abs_err});

		       // Compute P(proxy|proxy_ob)
		       double Pproxy = cbl::gaussian(x[2], ptr, {dummy_proxy,pp->proxy_rel_err*dummy_proxy});

		       // Compute the halo bias
		       double Delta = (pp->isDelta_critical) ? pp->Delta_input/cosmo.OmegaM(x[1]) : pp->Delta_input;
		       double bias_halo = cosmo.bias_halo(mass, sigmaM_interp(mass), x[1], DN_interp(x[1]), pp->model_bias, false, par::defaultString, "Linear", Delta, -1, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
		       
		       return bias_halo * P_M__lambda_z * Pz * Pproxy;
		     
		     };
    std::vector<std::vector<double>> integration_limits(3);
    cbl::wrapper::cuba::CUBAwrapper CW (integrand, (int)(integration_limits.size()));
    */
    
    
    // Compute the bias    
    vector<double> _bias(pp->scaling_relation->data()->xx().size());

    for (size_t i=0; i<pp->scaling_relation->data()->xx().size(); i++) {

      dummy_proxy = pp->scaling_relation->data()->xx(i);
      dummy_z = (pp->scaling_relation)->data_model().redshift[i];

      // Find the minimum and maximum masses, given the parameters of the scaling relation
      double log_lambda_min = log(dummy_proxy*(1-pp->proxy_rel_err)/proxy_pivot)/log(log_base);
      double log_lambda_max = log(dummy_proxy*(1+pp->proxy_rel_err)/proxy_pivot)/log(log_base);
      double log_f_z_min = log( (pp->scaling_relation)->data_model().fz(dummy_z-pp->z_abs_err, redshift_pivot, cosmo_ptr) )/log(log_base);
      double log_f_z_max = log( (pp->scaling_relation)->data_model().fz(dummy_z+pp->z_abs_err, redshift_pivot, cosmo_ptr) )/log(log_base);

      double logM1 = alpha + beta*log_lambda_min + gamma*log_f_z_min;
      double logM2 = alpha + beta*log_lambda_max + gamma*log_f_z_min;
      double logM3 = alpha + beta*log_lambda_min + gamma*log_f_z_max;
      double logM4 = alpha + beta*log_lambda_max + gamma*log_f_z_max;

      double min1 = std::min(logM1, logM2);
      double min2 = std::min(min1, logM3);
      double min_logM = std::min(min2, logM4);
      double max1 = std::max(logM1, logM2);
      double max2 = std::max(max1, logM3);
      double max_logM = std::max(max2, logM4);

      // Find the maximum value of the intrinsic scatter
      double s1 = std::abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
      double s2 = std::abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_min, scatterz_exp) );
      double s3 = std::abs( scatter0 + scatterM*pow(log_lambda_min, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );
      double s4 = std::abs( scatter0 + scatterM*pow(log_lambda_max, scatterM_exp) + scatterz*pow(log_f_z_max, scatterz_exp) );

      double maxs1 = std::max(s1, s2);
      double maxs2 = std::max(maxs1, s3);
      double max_scatter_intr = std::max(maxs2, s4);

      // Integrate
      _bias[i] = wrapper::gsl::GSL_integrate_qag(integrand, std::max(min_logM-3.5*max_scatter_intr,log(cbl::Min(Mass_vector)/mass_pivot)/log(log_base)), std::min(max_logM+3.5*max_scatter_intr,log(cbl::Max(Mass_vector)/mass_pivot)/log(log_base)));

      /*
      integration_limits[0] = {std::max(min_logM-3.5*max_scatter_intr,log(cbl::Min(Mass_vector)/mass_pivot)/log(log_base)), std::min(max_logM+3.5*max_scatter_intr,log(cbl::Max(Mass_vector)/mass_pivot)/log(log_base))};
      integration_limits[1] = {cbl::Min(z_for_DN), cbl::Max(z_for_DN)};
      integration_limits[2] = {std::max(dummy_proxy - 3.5*pp->proxy_rel_err*dummy_proxy, 0.00001), dummy_proxy + 3.5*pp->proxy_rel_err*dummy_proxy};
      _bias[i] = CW.IntegrateVegas(integration_limits,false);
      */
      
    }

    bias = Average(_bias);
    
  }

  // set the value of the bias
  parameter[parameter.size()-1] = bias;

  return modelling::twopt::damped_Xi(rad, bias, cosmo.linear_growth_rate(pp->redshift), SigmaS, pp->kk, Pk_interp);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_damped_bias_sigmaz (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // input likelihood parameters

  // bias
  const double bias = parameter[0];

  // damping scale
  const double SigmaS = par::cc*parameter[1]/pp->HHfid;

  return modelling::twopt::damped_Xi(rad, bias, pp->linear_growth_rate_z, SigmaS, pp->kk, pp->func_Pk);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_sigma8_clusters (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  cosmology::Cosmology cosmo = *pp->cosmology;

  
  // input likelihood parameters
  
  // sigma8
  const double sigma8 = parameter[0];
  cosmo.set_sigma8(sigma8);
  
  // sigma8(z)
  const double sigma8_z = sigma8*cosmo.DN(pp->redshift);


  // output likelihood parameters
  
  // mean bias

  vector<double> bias(pp->cluster_mass_proxy->ndata());
  
  for (int k=0; k<pp->cluster_mass_proxy->ndata(); k++) {

    const double sigma8fid = pp->sigma8_z/cosmo.DN(pp->redshift);
    
    const double sigma = pp->cluster_mass_proxy->extra_info(0, k)*(sigma8/sigma8fid);
    
    bias[k] = cosmo.bias_halo(pp->cluster_mass_proxy->data(k), sigma, pp->cluster_mass_proxy->xx(k), pp->model_bias, false, pp->output_root, "Spline", pp->Delta, 1., pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
  }
  
  const double mean_bias = Average(bias);
  parameter[1] = mean_bias;


  // sigma8(z)
  parameter[2] = sigma8_z;

  
  // fixed parameters 
  
  // f(z)*sigma8(z)
  const double fsigma8 = pp->linear_growth_rate_z*sigma8_z;
  
  
  // return the redshift-space monopole of the two-point correlation function

  vector<double> xi(rad.size(), 0);

  const double fact = xi_ratio(fsigma8, mean_bias*sigma8_z)*pow(mean_bias, 2)*pow(sigma8_z/pp->sigma8_z, 2);
  
  for (size_t i =0; i<xi.size(); i++) 
    xi[i] = fact*pp->func_xi->operator()(rad[i]);

  return xi;

}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_one_cosmo_par_clusters (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  cosmology::Cosmology cosmo = *pp->cosmology;
  
  // input likelihood parameter: the cosmological parameter used to
  // compute the dark matter two-point correlation function in real
  // space
  cosmo.set_parameter(pp->Cpar[0], parameter[0]);

  // the linear bias 
  const double bias = pp->cosmopar_bias_interp_1D(parameter[0]);
  parameter[1] = bias;

  // the AP parameter
  const double alpha = cosmo.D_V(pp->redshift)/pp->DVfid;
  

  // return the redshift-space monopole of the two-point correlation function
  
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= alpha;

  return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->redshift, false, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_two_cosmo_pars_clusters (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
   // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  cosmology::Cosmology cosmo = *pp->cosmology;
  
  // input likelihood parameters: the cosmological parameters used to
  // compute the dark matter two-point correlation function in real
  // space
  cosmo.set_parameter(pp->Cpar[0], parameter[0]);
  cosmo.set_parameter(pp->Cpar[1], parameter[1]);

  // the linear bias
  const double bias = pp->cosmopar_bias_interp_2D(parameter[0], parameter[1]);
  parameter[2] = bias;
  
  /// the AP parameter
  const double alpha = cosmo.D_V(pp->redshift)/pp->DVfid;
  parameter[3] = alpha;

  
  // return the redshift-space monopole of the two-point correlation function
  
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= alpha;
  
  return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->redshift, false, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->file_par);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_cosmology_clusters (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  cosmology::Cosmology cosmo = *pp->cosmology;

  // input likelihood parameters
  
  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);


  // bias
 
  vector<double> mass_grid = logarithmic_bin_vector(pp->cluster_mass_proxy->ndata()/10, Min(pp->cluster_mass_proxy->data()), Max(pp->cluster_mass_proxy->data()));

  const double bias = cosmo.bias_eff_mass(pp->cluster_mass_proxy->data(), mass_grid, pp->cluster_mass_proxy->xx(), pp->model_bias, pp->method_Pk, pp->meanType, false, pp->output_root, pp->Delta)[0];

  
  // fixed parameters 
  
  /// alpha
  const double alpha = cosmo.D_V(pp->redshift)/pp->DVfid;

  // return the redshift-space monopole of the two-point correlation function
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= alpha;

  // return the redshift-space monopole of the two-point correlation function
  const double sigma8 = parameter[0];
  const double sigma8_z = sigma8*pp->cosmology->DN(pp->redshift);
  const double fsigma8 = pp->linear_growth_rate_z*sigma8_z;
  
  vector<double> xi(rad.size(), 0);

  const double fact = xi_ratio(fsigma8, bias*sigma8_z)*pow(bias, 2)*pow(sigma8_z/pp->sigma8_z, 2);
  
  for (size_t i =0; i<xi.size(); i++) 
    xi[i] = fact*pp->func_xi->operator()(new_rad[i]);
  
  return xi;
  
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi0_linear_cosmology_clusters_selection_function (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // set the test cosmology
  cosmology::Cosmology cosmo = *pp->test_cosmology;

  
  // ----- input likelihood parameters -----

  // set the cosmological parameters used to compute the dark matter
  // two-point correlation function in real space
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  // set the alpha parameter of the cluster mass scaling relation
  const double alpha = parameter[pp->Cpar.size()];


  // set the function to estimate the linear dark matter power spectrum at z=0, by interpolating it from a grid
  
  const vector<double> Pk_grid = cosmo.Pk_matter(pp->kk, pp->method_Pk, false, 0., false, pp->output_root, -1, pp->k_min, pp->k_max, pp->prec, pp->file_par);
  glob::FuncGrid interp_Pk(pp->kk, Pk_grid, "Spline");


  // compute sigma and dlnsigma, by interpolating them from a grid
  
  // set the cluster masses and rho_m
  const vector<double> mass = pp->mass;
  const double rho = cosmo.rho_m(0.);
  
  vector<double> sigma_grid, dnsigma_grid;

  for (size_t i=0; i<mass.size(); i++) {
    
    const double RR = Radius(mass[i], rho); 

    // compute sigma
    auto func_sigma = [&] (double kk)
      {
	return pow(TopHat_WF(kk*RR)*kk, 2)*interp_Pk(kk); 
      };
    sigma_grid.emplace_back(sqrt(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag(func_sigma, pp->k_min, pp->k_max, pp->prec)));

    // compute dlnsigma
    const double dRdM = pow(3./(4.*par::pi*rho), 1./3.)*pow(mass[i], -2./3.)/3.;
    auto func_dnsigma = [&] (const double kk)
      {
	double filter = 2*TopHat_WF(kk*RR)*TopHat_WF_D1(kk*RR)*kk*dRdM;
	return filter*pow(kk, 2)*interp_Pk(kk); 
      };
    dnsigma_grid.emplace_back(1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag(func_dnsigma, pp->k_min, pp->k_max, pp->prec));

  }
  
  glob::FuncGrid interp_sigma(mass, sigma_grid, "Spline");
  glob::FuncGrid interp_DnSigma(mass, dnsigma_grid, "Spline");

  
  // compute the bias
  const double bias = cosmo.bias_eff_selection_function(interp_sigma, interp_DnSigma, *pp->interp_SelectionFunction_cut, pp->Mass_min, pp->Mass_max, {pp->redshift}, pp->model_bias, pp->model_MF, "EisensteinHu", alpha, false, pp->output_root, pp->Delta, -1., "Spline", pp->norm, pp->k_min, pp->k_max, pp->prec)[0]; // check!!!
  parameter[pp->Cpar.size()+1] = bias;
  
  // set the AP factor
  const double AP_factor = cosmo.D_V(pp->redshift)/pp->DVfid;

  // rescale with the AP factor
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= AP_factor;
  
  // compute the real-space monopole of the two-point correlation function at z=0, by Fourier transforming the P(k)
  vector<double> xi = wrapper::fftlog::transform_FFTlog(new_rad, 1, pp->kk, Pk_grid, 0);
  
  // compute the redshift-space monopole at z=pp->redshift (scaling by D(z)/D(0) the monopole at z=0) 
  const double fact = bias*bias*xi_ratio(cosmo.linear_growth_rate(pp->redshift)/bias)*pow(cosmo.DN(pp->redshift), 2);
  for (size_t i=0; i<xi.size(); i++)
    xi[i] *= fact;

  return xi;
}


// ============================================================================================


double cbl::modelling::twopt::Ncen (const double Mass, const double lgMmin, const double sigmalgM) 
{
  const double Nc = 0.5*(1.+gsl_sf_erf(log10(Mass)-lgMmin)/(sqrt(2.)*sigmalgM)); 
  return (Nc<0 || std::isnan(Nc)) ? 0. : Nc;
}


// ============================================================================================


double cbl::modelling::twopt::Nsat (const double Mass, const double lgMmin, const double sigmalgM, const double lgM0, const double lgM1, const double alpha) 
{
  const double Ns = Ncen(Mass, lgMmin, sigmalgM)*pow((Mass-pow(10.,lgM0))/pow(10.,lgM1),alpha);
  return (Ns<0 || std::isnan(Ns)) ? 0. : Ns;
}


// ============================================================================================


double cbl::modelling::twopt::Navg (const double Mass, const double lgMmin, const double sigmalgM, const double lgM0, const double lgM1, const double alpha) 
{
  return Ncen(Mass, lgMmin, sigmalgM) + Nsat(Mass, lgMmin, sigmalgM, lgM0, lgM1, alpha);
}


// ============================================================================================


double cbl::modelling::twopt::ng (const double lgMmin, const double sigmalgM, const double lgM0, const double lgM1, const double alpha, const std::shared_ptr<void> inputs)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  auto ng_integrand = [&] (const double lgmass)
    {
      return pp->interpMF(pow(10., lgmass))*Navg(pow(10., lgmass), lgMmin, sigmalgM, lgM0, lgM1, alpha)*pow(10., lgmass)*log(10.);
    };
  
  return wrapper::gsl::GSL_integrate_qag(ng_integrand, log10(1.e10), log10(1.e16));
}


// ============================================================================================


double cbl::modelling::twopt::bias (const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const std::shared_ptr<void> inputs)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);
  
  auto func = [&] (double mass)
    {
      // halo mass function -> it depends on cosmology
      const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, false, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
      
      // halo bias -> it depends on cosmology
      const double bias_halo = pp->cosmology->bias_halo(mass, pp->func_sigma->operator()(mass), pp->redshift, pp->model_bias, false, pp->output_root, pp->interpType, pp->Delta, pp->kk, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
      
      // average number of galaxies in each halo -> it depends on
      // galaxy evolution
      const double NN = Navg(mass, Mmin, sigmalgM, M0, M1, alpha);
      
      return dndM*NN*bias_halo;
    };

  return 1./ng(Mmin, sigmalgM, M0, M1, alpha, inputs)*wrapper::gsl::GSL_integrate_qag(func, pp->Mh_min, pp->Mh_max);
}


// ============================================================================================


double cbl::modelling::twopt::NcNs (const double Mass, const double lgMmin, const double sigmalgM, const double lgM0, const double lgM1, const double alpha)
{
  return Ncen(Mass, lgMmin, sigmalgM)*Nsat(Mass, lgMmin, sigmalgM, lgM0, lgM1, alpha);
}


// ============================================================================================


double cbl::modelling::twopt::NsNs1 (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha)
{
  return pow(Nsat(Mass, Mmin, sigmalgM, M0, M1, alpha), 2);
}


// ============================================================================================


double cbl::modelling::twopt::Pk_cs (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter) {
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double lgMmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double lgM0 = parameter[2];
  const double lgM1 = parameter[3];
  const double alpha = parameter[4];

  // build a HaloProfile object for computing the density profile in Fourier space
  cbl::cosmology::HaloProfile halo_profile (*pp->cosmology, pp->redshift, "Duffy", 0., 200., pp->profile, pp->halo_def);  
  
  auto Pk_cs_integrand = [&] (const double lgmass) {
			   halo_profile.set_mass(pow(10., lgmass));
			   return pp->interpMF(pow(10., lgmass))*NcNs(pow(10., lgmass), lgMmin, sigmalgM, lgM0, lgM1, alpha)*halo_profile.density_profile_FourierSpace(kk)*pow(10., lgmass)*log(10.);
			 };

  return 2./pow(ng(lgMmin, sigmalgM, lgM0, lgM1, alpha, inputs), 2)*wrapper::gsl::GSL_integrate_qag(Pk_cs_integrand, log10(pp->Mh_min), log10(pp->Mh_max));

}


// ============================================================================================


double cbl::modelling::twopt::Pk_ss (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter) {

  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double lgMmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double lgM0 = parameter[2];
  const double lgM1 = parameter[3];
  const double alpha = parameter[4];

  // build a HaloProfile object for computing the density profile in Fourier space
  cbl::cosmology::HaloProfile halo_profile (*pp->cosmology, pp->redshift, "Duffy", 0., 200., pp->profile, pp->halo_def);   
    
  auto Pk_ss_integrand = [&] (const double lgmass) {
			   halo_profile.set_mass(pow(10., lgmass));
			   return pp->interpMF(pow(10., lgmass))*NsNs1(pow(10., lgmass), parameter[0], parameter[1], parameter[2], parameter[3], parameter[4])*pow(halo_profile.density_profile_FourierSpace(kk), 2)*pow(10., lgmass)*log(10.);
			 };
  
  return 1./pow(ng(lgMmin, sigmalgM, lgM0, lgM1, alpha, inputs), 2)*wrapper::gsl::GSL_integrate_qag(Pk_ss_integrand, log10(pp->Mh_min), log10(pp->Mh_max));
}


// ============================================================================================


double cbl::modelling::twopt::Pk_1halo (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return Pk_cs(kk, inputs, parameter)+Pk_ss(kk, inputs, parameter);   
}


// ============================================================================================


double cbl::modelling::twopt::Pk_2halo (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter) {

  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double lgMmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double lgM0 = parameter[2];
  const double lgM1 = parameter[3];
  const double alpha = parameter[4];

  // build a HaloProfile object for computing the density profile in Fourier space
  cbl::cosmology::HaloProfile halo_profile (*pp->cosmology, pp->redshift, "Duffy", 0., 200., pp->profile, pp->halo_def); 

  auto Pk_2halo_integrand = [&] (const double lgmass) {
			      halo_profile.set_mass(pow(10., lgmass));
			      return pp->interpMF(pow(10., lgmass))*Navg(pow(10., lgmass), parameter[0], parameter[1], parameter[2], parameter[3], parameter[4])*pp->interpBias(pow(10., lgmass))*halo_profile.density_profile_FourierSpace(kk)*pow(10., lgmass)*log(10.);
  };
  return pp->interpPk(kk)*1./pow(ng(lgMmin, sigmalgM, lgM0, lgM1, alpha, inputs),2.)*pow(wrapper::gsl::GSL_integrate_qag(Pk_2halo_integrand, log10(pp->Mh_min), log10(pp->Mh_max)), 2.);
}


// ============================================================================================


double cbl::modelling::twopt::Pk_HOD (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return Pk_1halo(kk, inputs, parameter)+Pk_2halo(kk, inputs, parameter);  
}


// ============================================================================================


vector<double> cbl::modelling::twopt::xi_1halo (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  //structure contaning the required input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  vector<double> Pk(pp->kkvec.size(), 0);

  #pragma omp parallel num_threads(omp_get_max_threads())
  {  
    #pragma omp for schedule(static, 2)
      for (size_t i=0; i<pp->kkvec.size(); i++)
        Pk[i] = Pk_1halo(pp->kkvec[i], inputs, parameter);
  }

  return wrapper::fftlog::transform_FFTlog(rad, 1, pp->kkvec, Pk, 0);
}


// ============================================================================================


vector<double> cbl::modelling::twopt::xi_2halo (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  //structure contaning the required input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  vector<double> Pk(pp->kkvec.size(), 0);

  #pragma omp parallel num_threads(omp_get_max_threads())
  {  
    #pragma omp for schedule(static, 2)
      for (size_t i=0; i<pp->kkvec.size(); i++)
        Pk[i] = Pk_2halo(pp->kkvec[i], inputs, parameter);
  }

  return wrapper::fftlog::transform_FFTlog(rad, 1, pp->kkvec, Pk, 0);

}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi_HOD (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  const vector<double> xi1h = xi_1halo(rad, inputs, parameter);
  const vector<double> xi2h = xi_2halo(rad, inputs, parameter);
  vector<double> xi(rad.size());
  
#pragma omp parallel num_threads(omp_get_max_threads())
  {
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<rad.size(); i++) 
      xi[i] = xi1h[i]+xi2h[i];
  }
  
  return xi;
}


// ============================================================================================


double cbl::modelling::twopt::xi_zspace (FunctionVectorVectorPtrVectorRef func, const double rp, const double pi, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data

  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  
  // input parameters

  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];

  
  const double rad = sqrt(pow(rp, 2)+pow(pi, 2));
  const double mu = rp/rad;

  const double beta2 = pp->cosmology->linear_growth_rate(pp->redshift)/bias(Mmin, sigmalgM, M0, M1, alpha, inputs);

  const double fact_xi0 = 1.+2./3.*beta2+1./5.*beta2*beta2;
  const double fact_xi2 = (4./3.*beta2+4./7.*beta2*beta2);
  const double fact_xi4 = 8./35.*beta2*beta2;

  auto func3 = [&] (double yy) { return func({yy}, inputs, parameter)[0]*pow(yy, 2); };
  const double J3 = pow(rad, -3)*wrapper::gsl::GSL_integrate_qag(func3, 0., rad);
 
  auto func5 = [&] (double yy) { return func({yy}, inputs, parameter)[0]*pow(yy, 4); };
  const double J5 = pow(rad, -5)*wrapper::gsl::GSL_integrate_qag(func5, 0., rad);
  
  const double xi0 = fact_xi0*func({rad}, inputs, parameter)[0];
  const double xi2 = fact_xi2*(func({rad}, inputs, parameter)[0]-3.*J3);
  const double xi4 = fact_xi4*(func({rad}, inputs, parameter)[0]+15./2.*J3-35./2.*J5);
  
  return xi0 + xi2*legendre_polynomial(mu, 2) + xi4*legendre_polynomial(mu, 4);
}


// ============================================================================================


double cbl::modelling::twopt::xi_1halo_zspace (const double rp, const double pi, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return xi_zspace(xi_1halo, rp, pi, inputs, parameter);
}


// ============================================================================================


double cbl::modelling::twopt::xi_2halo_zspace (const double rp, const double pi, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return xi_zspace(xi_2halo, rp, pi, inputs, parameter);
}


// ============================================================================================


double cbl::modelling::twopt::xi_HOD_zspace (const double rp, const double pi, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return xi_zspace(xi_HOD, rp, pi, inputs, parameter);
}

