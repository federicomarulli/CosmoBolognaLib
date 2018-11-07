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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_TwoPointCorrelation1D_monopole.h"

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

  vector<double> xi = fftlog::transform_FFTlog(new_rad, 1, pp->kk, Pk, 0);
  
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
    model[i] = gsl::GSL_polynomial_eval(rad[i], NULL, poly_params);

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
    parameter[0] = gsl::GSL_root_brent(model_derivative, 0., xmin, xmax, 1.e-10);
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
      parameter[1] = gsl::GSL_root_brent(model_derivative, 0., xmin, xmax, 1.e-10);

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

    model[i] =  xi_ratio(fsigma8, bsigma8)*pow(bsigma8/pp->sigma8_z, 2)*pp->func_xi->operator()(rad[i]*alpha)+poly;

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
    parameter[0] = gsl::GSL_root_brent(model_derivative, 0., xmin, xmax, 1.e-10);
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
    parameter[1] = gsl::GSL_root_brent(model_derivative, 0., 40, parameter[0]-2, 1.e-10);
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

    xi[i] = xi_ratio(fsigma8, bsigma8)*pp->cosmology->xi_DM(rad[i]*alpha, pp->method_Pk, pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par)/pow(pp->sigma8_z, 2)+poly;
  }

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

  return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->NL, pp->redshift, pp->output_dir, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->step, pp->prec, pp->file_par);
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

    _bias[i] = pp->cosmology->bias_halo(mass[i], pp->func_sigma->operator()(mass[i]), pp->cluster_mass_proxy->xx(i), pp->model_bias, par::defaultString, "Linear", pp->Delta);
  }

  const double bias = Average(_bias);

  // return the value of the bias
  parameter[4] = bias;

  return modelling::twopt::damped_Xi(rad, bias, pp->linear_growth_rate_z, SigmaS, pp->kk, pp->func_Pk);
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
  cosmo.set_sigma8(parameter[0]);

  
  // input likelihood parameters
  
  // sigma8
  const double sigma8 = parameter[0];

  // sigma8(z)
  const double sigma8_z = sigma8*cosmo.DD(pp->redshift)/cosmo.DD(0.);


  // output likelihood parameters
  
  // mean bias

  vector<double> bias(pp->cluster_mass_proxy->ndata());
  
  for (int k=0; k<pp->cluster_mass_proxy->ndata(); k++) {

    const double sigma8fid = pp->sigma8_z*cosmo.DD(0.)/cosmo.DD(pp->redshift);
    
    const double sigma = pp->cluster_mass_proxy->extra_info(0, k)*(sigma8/sigma8fid);
    
    bias[k] = cosmo.bias_halo(pp->cluster_mass_proxy->data(k), sigma, pp->cluster_mass_proxy->xx(k), pp->model_bias, pp->output_root, "Spline", pp->Delta, 1., pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk);
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
  
  return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->NL, pp->redshift, pp->output_dir, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->step, pp->prec, pp->file_par);
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
  
  return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->NL, pp->redshift, pp->output_dir, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->step, pp->prec, pp->file_par);
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

  const double bias = cosmo.bias_eff_mass(pp->cluster_mass_proxy->data(), mass_grid, pp->cluster_mass_proxy->xx(), pp->model_bias, pp->method_Pk, pp->meanType, pp->output_root, pp->Delta)[0];

  
  // fixed parameters 
  
  /// alpha
  const double alpha = cosmo.D_V(pp->redshift)/pp->DVfid;

  // return the redshift-space monopole of the two-point correlation function
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= alpha;
  
  //return cosmo.xi0_Kaiser(new_rad, bias, pp->method_Pk, pp->NL, pp->redshift, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->step, pp->prec, pp->file_par);


  // return the redshift-space monopole of the two-point correlation function
  const double sigma8 = parameter[0];
  const double sigma8_z = sigma8*pp->cosmology->DD(pp->redshift)/pp->cosmology->DD(0.);
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
  
  const vector<double> Pk_grid = cosmo.Pk(pp->kk, pp->method_Pk, false, 0., pp->output_dir, pp->output_root, -1, pp->k_min, pp->k_max, pp->prec, pp->file_par);
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
    sigma_grid.emplace_back(sqrt(1./(2.*pow(par::pi, 2))*gsl::GSL_integrate_qag(func_sigma, pp->k_min, pp->k_max, pp->prec)));

    // compute dlnsigma
    const double dRdM = pow(3./(4.*par::pi*rho), 1./3.)*pow(mass[i], -2./3.)/3.;
    auto func_dnsigma = [&] (const double kk)
      {
	double filter = 2*TopHat_WF(kk*RR)*TopHat_WF_D1(kk*RR)*kk*dRdM;
	return filter*pow(kk, 2)*interp_Pk(kk); 
      };
    dnsigma_grid.emplace_back(1./(2.*pow(par::pi, 2))*gsl::GSL_integrate_qag(func_dnsigma, pp->k_min, pp->k_max, pp->prec));

  }
  
  glob::FuncGrid interp_sigma(mass, sigma_grid, "Spline");
  glob::FuncGrid interp_DnSigma(mass, dnsigma_grid, "Spline");

  
  // compute the bias
  const double bias = cosmo.bias_eff_selection_function(interp_sigma, interp_DnSigma, *pp->interp_SelectionFunction_cut, pp->Mass_min, pp->Mass_max, {pp->redshift}, pp->model_bias, pp->model_MF, "EisensteinHu", alpha, pp->output_root, pp->Delta, -1., "Spline", pp->norm, pp->k_min, pp->k_max, pp->prec)[0]; // check!!!
  parameter[pp->Cpar.size()+1] = bias;
  
  // set the AP factor
  const double AP_factor = cosmo.D_V(pp->redshift)/pp->DVfid;

  // rescale with the AP factor
  vector<double> new_rad = rad;
  for (size_t i=0; i<rad.size(); i++)
    new_rad[i] *= AP_factor;
  
  // compute the real-space monopole of the two-point correlation function at z=0, by Fourier transforming the P(k)
  vector<double> xi = fftlog::transform_FFTlog(new_rad, 1, pp->kk, Pk_grid, 0);

  // compute the redshift-space monopole at z=pp->redshift
  const double fact = pow(bias, 2)*xi_ratio(cosmo.linear_growth_rate(pp->redshift, 1), bias)*pow(cosmo.DD(pp->redshift)/cosmo.DD(0), 2);
  for (size_t i=0; i<xi.size(); i++)
    xi[i] *= fact;

  return xi;
}


// ============================================================================================


double cbl::modelling::twopt::Ncen (const double Mass, const double Mmin, const double sigmalgM) 
{
  const double Nc = 0.5*(1.+gsl_sf_erf((log10(Mass)-log10(Mmin))/(sqrt(2.)*sigmalgM))); 
  return (Nc<0 || std::isnan(Nc)) ? 0. : Nc;
}


// ============================================================================================


double cbl::modelling::twopt::Nsat (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha) 
{
  const double Ns = Ncen(Mass, Mmin, sigmalgM)*pow((Mass-M0)/M1, alpha);
  return (Ns<0 || std::isnan(Ns)) ? 0. : Ns;
}


// ============================================================================================


double cbl::modelling::twopt::Navg (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha) 
{
  return Ncen(Mass, Mmin, sigmalgM) + Nsat(Mass, Mmin, sigmalgM, M0, M1, alpha);
}


// ============================================================================================


double cbl::modelling::twopt::ng_integrand (const double mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const std::shared_ptr<void> inputs)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);
  
  // halo mass function -> it depends on cosmology
  const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
  
  // average number of galaxies in each halo -> it depends on
  // galaxy evolution
  const double NN = Navg(mass, Mmin, sigmalgM, M0, M1, alpha);
  
  return dndM*NN;
}


// ============================================================================================


double cbl::modelling::twopt::ng (const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const std::shared_ptr<void> inputs)
{
  return gsl::GSL_integrate_qag(bind(&modelling::twopt::ng_integrand, std::placeholders::_1, Mmin, sigmalgM, M0, M1, alpha, inputs), 1.e10, 1.e16);
}


// ============================================================================================


double cbl::modelling::twopt::bias (const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha, const std::shared_ptr<void> inputs)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);
  
  auto func = [&] (double mass)
    {
      // halo mass function -> it depends on cosmology
      const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
      
      // halo bias -> it depends on cosmology
      const double bias_halo = pp->cosmology->bias_halo(mass, pp->func_sigma->operator()(mass), pp->redshift, pp->model_bias, pp->output_root, pp->interpType, pp->Delta, pp->kk, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
      
      // average number of galaxies in each halo -> it depends on
      // galaxy evolution
      const double NN = Navg(mass, Mmin, sigmalgM, M0, M1, alpha);
      
      return dndM*NN*bias_halo;
    };

  return 1./ng(Mmin, sigmalgM, M0, M1, alpha, inputs)*gsl::GSL_integrate_qag(func, pp->Mh_min, pp->Mh_max);
}


// ============================================================================================


double cbl::modelling::twopt::NcNs (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha)
{
  return Ncen(Mass, Mmin, sigmalgM)*Nsat(Mass, Mmin, sigmalgM, M0, M1, alpha);
}


// ============================================================================================


double cbl::modelling::twopt::NsNs1 (const double Mass, const double Mmin, const double sigmalgM, const double M0, const double M1, const double alpha)
{
  return pow(Nsat(Mass, Mmin, sigmalgM, M0, M1, alpha), 2);
}


// ============================================================================================


double cbl::modelling::twopt::Pk_cs_numerator_integrand (const double mass, const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];

  // halo mass function -> it depends on cosmology
  const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
  
  // density profile -> it depends on cosmology
  const double uk = pp->cosmology->density_profile_FourierSpace(kk, mass, pp->redshift, pp->model_cM, pp->profile, pp->halo_def);
  
  // mean number of central-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = NcNs(mass, Mmin, sigmalgM, M0, M1, alpha);
  
  return NN*dndM*uk;
}


// ============================================================================================


double cbl::modelling::twopt::Pk_cs (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];
  
  return 2./pow(ng(Mmin, sigmalgM, M0, M1, alpha, inputs), 2)*gsl::GSL_integrate_qag(bind(&modelling::twopt::Pk_cs_numerator_integrand, std::placeholders::_1, kk, inputs, parameter), pp->Mh_min, pp->Mh_max);
}


// ============================================================================================


double cbl::modelling::twopt::Pk_ss_numerator_integrand (const double mass, const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];

  // halo mass function -> it depends on cosmology
  const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
      
  // density profile -> it depends on cosmology
  const double uk2 = pow(pp->cosmology->density_profile_FourierSpace(kk, mass, pp->redshift, pp->model_cM, pp->profile, pp->halo_def), 2);
  
  // mean number of satellite-satellite galaxy pairs -> it depends
  // on galaxy evolution
  const double NN = NsNs1(mass, Mmin, sigmalgM, M0, M1, alpha);
      
  return NN*dndM*uk2;
}


// ============================================================================================


double cbl::modelling::twopt::Pk_ss (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];
  
  return 1./pow(ng(Mmin, sigmalgM, M0, M1, alpha, inputs), 2)*gsl::GSL_integrate_qag(bind(&modelling::twopt::Pk_ss_numerator_integrand, std::placeholders::_1, kk, inputs, parameter), pp->Mh_min, pp->Mh_max);
}


// ============================================================================================


double cbl::modelling::twopt::Pk_1halo (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return Pk_cs(kk, inputs, parameter)+Pk_ss(kk, inputs, parameter);   
}

  
// ============================================================================================


double cbl::modelling::twopt::Pk_2halo (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];
  
  auto func = [&] (double mass)
    {
      // halo mass function -> it depends on cosmology
      const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);

      // halo bias -> it depends on cosmology
      const double bias = pp->cosmology->bias_halo(mass, pp->func_sigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->interpType, pp->Delta, kk, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
      
      // density profile -> it depends on cosmology
      const double uk = pp->cosmology->density_profile_FourierSpace(kk, mass, pp->redshift, pp->model_cM, pp->profile, pp->halo_def);
      
      // mean number of galaxy pairs -> it depends on galaxy evolution
      const double NN = Navg(mass, Mmin, sigmalgM, M0, M1, alpha);

      return NN*dndM*bias*uk;
    };
  
  return pp->func_Pk->operator()(kk)*pow(1./ng(Mmin, sigmalgM, M0, M1, alpha, inputs)*gsl::GSL_integrate_qag(func, pp->Mh_min, pp->Mh_max), 2);

}


// ============================================================================================


double cbl::modelling::twopt::Pk_HOD (const double kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  return Pk_1halo(kk, inputs, parameter)+Pk_2halo(kk, inputs, parameter);  
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi_1halo (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];
  
  // multiplicative factor
  const double fact = 1./(2.*pow(par::pi*ng(Mmin, sigmalgM, M0, M1, alpha, inputs), 2));
  
  // vector that will be return in output
  vector<double> xi(rad.size(), fact);

  // limits of the integral
  vector<vector<double>> integration_limits(2);
  integration_limits[0] = {0., pp->k_max};
  integration_limits[1] = {1.e10, 1.e16};
  
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<xi.size(); i++) {
      
      // integrand function
      auto func = [&] (const vector<double> kk_mass)
	{
	  return (2.*Pk_cs_numerator_integrand(kk_mass[1], kk_mass[0], inputs, parameter)+Pk_ss_numerator_integrand(kk_mass[1], kk_mass[0], inputs, parameter))*kk_mass[0]*sin(kk_mass[0]*rad[i])/rad[i];
	};
      
      // wrapper to CUBA libraries
      cuba::CUBAwrapper CW(func, 2);
      
      // integrate with Cuhre
      xi[i] *= CW.IntegrateCuhre(integration_limits);
      
    }
    
  }

  return xi;
}


// ============================================================================================


std::shared_ptr<glob::FuncGrid> cbl::modelling::twopt::func_2halo (const std::vector<double> kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the HOD input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];

  vector<double> fact(kk.size());
  
  for (size_t i=0; i<kk.size(); i++) {
  
    auto func = [&] (double mass)
      {
	// halo mass function -> it depends on cosmology
	const double dndM = pp->cosmology->mass_function(mass, pp->func_sigma->operator()(mass), pp->func_dlnsigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->Delta, pp->interpType, pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
	
	// halo bias -> it depends on cosmology
	const double bias = pp->cosmology->bias_halo(mass, pp->func_sigma->operator()(mass), pp->redshift, pp->model_MF, pp->output_root, pp->interpType, pp->Delta, kk[i], pp->norm, pp->k_min, pp->k_max, pp->prec, pp->method_Pk, pp->input_file, pp->is_parameter_file);
	
	// density profile -> it depends on cosmology
	const double uk = pp->cosmology->density_profile_FourierSpace(kk[i], mass, pp->redshift, pp->model_cM, pp->profile, pp->halo_def);
	
	// mean number of galaxy pairs -> it depends on galaxy evolution
	const double NN = Navg(mass, Mmin, sigmalgM, M0, M1, alpha);
	
	return NN*dndM*bias*uk;
      };
    
    fact[i] = pow(gsl::GSL_integrate_qag(func, pp->Mh_min, pp->Mh_max), 2);
  }
  
  return move(unique_ptr<glob::FuncGrid>(new glob::FuncGrid(kk, fact, "Spline")));
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xi_2halo (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_HOD> pp = static_pointer_cast<STR_data_HOD>(inputs);

  // input parameters
  const double Mmin = parameter[0];
  const double sigmalgM = parameter[1];
  const double M0 = parameter[2];
  const double M1 = parameter[3];
  const double alpha = parameter[4];
  
  // vector that will be return in output
  vector<double> xi(rad.size());

  // grid used to speed up the computation
  auto func = func_2halo(logarithmic_bin_vector(50, 1.e-4, 500.), inputs, parameter);
  
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<xi.size(); i++) {
      
      auto integrand = [&] (double kk) { return pp->func_Pk->operator()(kk)*func->operator()(kk)*kk*sin(kk*rad[i])/rad[i]; };

      xi[i] = 1./(2.*pow(par::pi*ng(Mmin, sigmalgM, M0, M1, alpha, inputs),2))*gsl::GSL_integrate_qag(integrand, 0., pp->k_max);
      
    }

  }

  return xi;
  
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
  const double J3 = pow(rad, -3)*gsl::GSL_integrate_qag(func3, 0., rad);
 
  auto func5 = [&] (double yy) { return func({yy}, inputs, parameter)[0]*pow(yy, 4); };
  const double J5 = pow(rad, -5)*gsl::GSL_integrate_qag(func5, 0., rad);
  
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

