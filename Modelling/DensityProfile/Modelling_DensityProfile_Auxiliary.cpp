/********************************************************************
 *  Copyright (C) 2024 by Federico Marulli and Giorgio Lesci        *
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
 *  Modelling/DensityProfile/Modelling_DensityProfile_Auxiliary.cpp
 *
 *  @brief Methods of the class Modelling_DensityProfile_Auxiliary
 *
 *  This file contains the implementation of the methods of
 *  Modelling_DensityProfile_Auxiliary, used by the class Modelling_DensityProfile
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */


#include "Modelling_DensityProfile_Auxiliary.h"

using namespace std;

using namespace cbl;


// ===========================================================================================

double cbl::modelling::densityprofileaux::splashback_radius (const std::vector<double> rho, const std::vector<double> radius)
{
  // Interpolate the logarithm of rho
  std::vector<double> log_rho (rho.size(), 0.);
  std::vector<double> log_R (radius.size(), 0.);
  for (size_t i=0; i<log_rho.size(); i++) {
    log_rho[i] = log10(rho[i]);
    log_R[i] = log10(radius[i]);
  }
  cbl::glob::FuncGrid log_rho_interp (log_R, log_rho, "Spline");

  // rho derivative
  auto fc = [&] (const double x)
	    {		       
	      return log_rho_interp(x);
	    };
  
  std::vector<double> log_R_forDeriv = cbl::logarithmic_bin_vector( 100, cbl::Min(radius), std::min(cbl::Max(radius), 8.) ); // At most up to 8 Mpc
  for (size_t i=0; i<log_R_forDeriv.size(); i++) { log_R_forDeriv[i] = log10(log_R_forDeriv[i]); }
  
  std::vector<double> log_rho_deriv (log_R_forDeriv.size(), 0.);
  for (size_t i=0; i<log_R_forDeriv.size(); i++)
    log_rho_deriv[i] = cbl::wrapper::gsl::GSL_derivative(fc, log_R_forDeriv[i], 0.01);

  // Interpolate the derivative
  cbl::glob::FuncGrid log_rho_deriv_interp (log_R_forDeriv, log_rho_deriv, "Spline");

  // Find its minimum
  auto fc_deriv = [&] (const double x)
		  {		       
		    return log_rho_deriv_interp(x);
		  };

  const double log_Rsp = cbl::wrapper::gsl::GSL_minimize_1D(fc_deriv, 0.5 * (cbl::Min(log_R_forDeriv) + cbl::Max(log_R_forDeriv)), cbl::Min(log_R_forDeriv), cbl::Max(log_R_forDeriv));

  return pow(10, log_Rsp);
}

// ===========================================================================================

double cbl::modelling::densityprofileaux::logMass (const double z_tr, const double proxy_tr, const double alpha, const double beta, const double gamma, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double logM_base)
{
  auto cosmo_ptr = move(cosmo);
  double log_proxy_tr = log(proxy_tr/proxy_pivot) / log(logM_base);
  double log_f_z_tr = log( (scaling_relation)->data_model().fz(z_tr, z_pivot, cosmo_ptr) ) / log(logM_base);
  return alpha + beta*log_proxy_tr + gamma*log_f_z_tr;
}

// ===========================================================================================

double cbl::modelling::densityprofileaux::scatter_intr_M (const double z_tr, const double proxy_tr, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double logM_base)
{
  auto cosmo_ptr = move(cosmo);
  double log_proxy_tr = log(proxy_tr/proxy_pivot) / log(logM_base);
  double log_f_z_tr = log( (scaling_relation)->data_model().fz(z_tr, z_pivot, cosmo_ptr) ) / log(logM_base);
  return std::abs(scatter0 + scatterM*pow(log_proxy_tr, scatterM_exp) + scatterz*pow(log_f_z_tr, scatterz_exp));
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofileaux::min_max_logMass (const double z_min, const double z_max, const double proxy_min, const double proxy_max, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double logM_base)
{
  auto cosmo_ptr = move(cosmo);
  
  // Find the minimum and maximum masses, given the parameters of the scaling relation and the intrinsic scatter
  double logM1 = logMass(z_min, proxy_min, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);
  double logM2 = logMass(z_min, proxy_max, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);
  double logM3 = logMass(z_max, proxy_min, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);
  double logM4 = logMass(z_max, proxy_max, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);

  double min1 = min(logM1, logM2);
  double min2 = min(min1, logM3);
  double min_logM = min(min2, logM4);
  double max1 = max(logM1, logM2);
  double max2 = max(max1, logM3);
  double max_logM = max(max2, logM4);

  // Find the maximum value of the intrinsic scatter
  double maxs1 = max(
		     scatter_intr_M(z_min, proxy_min, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base),
		     scatter_intr_M(z_min, proxy_max, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base)
		     );
  double maxs2 = max(
		     maxs1,
		     scatter_intr_M(z_max, proxy_min, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base)
		     );
  double scatter = max(
		       maxs2,
		       scatter_intr_M(z_max, proxy_max, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base)
		       );

  std::vector<double> out = {min_logM - 3.5*scatter, max_logM + 3.5*scatter};
			   
  return out;
}

// ===========================================================================================

std::vector<std::vector<cbl::glob::FuncGrid2D>> cbl::modelling::densityprofileaux::interpolated_profiles (const std::vector<double> radius, const std::vector<double> redshift_vector, const std::vector<double> Mass_vector, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const double AB_fact, const std::string bias_author, const std::string method_Pk, std::function<double(const double, const double, const double, const double, const double)> conc_scaling_relation_func, const double c0, const double cM, const double cz, const double sigma_off0, const double sigma_offM, const double sigma_offz, const double Rt, const double z_pivot, const double mass_pivot, const std::shared_ptr<cosmology::Cosmology> cosmo_fid, const double b_e, const double s_e, const double alpha_0, const double alpha_nu, const double beta, const double gamma_0)
{
  auto cosmo_ptr = move(cosmo);
  auto cosmo_fid_ptr = (cosmo_fid == NULL) ? move(cosmo) : move(cosmo_fid);
  
  // Compute the profile at the given redshift and mass values
  std::vector<std::vector<std::vector<double>>> DeltaSigma_cen_samples(redshift_vector.size(), std::vector<std::vector<double>>(Mass_vector.size(), std::vector<double>(radius.size())));
  std::vector<std::vector<std::vector<double>>> Sigma_cen_samples(redshift_vector.size(), std::vector<std::vector<double>>(Mass_vector.size(), std::vector<double>(radius.size())));
  std::vector<std::vector<std::vector<double>>> DeltaSigma_off_samples(redshift_vector.size(), std::vector<std::vector<double>>(Mass_vector.size(), std::vector<double>(radius.size())));
  std::vector<std::vector<std::vector<double>>> Sigma_off_samples(redshift_vector.size(), std::vector<std::vector<double>>(Mass_vector.size(), std::vector<double>(radius.size())));

  for (size_t i=0; i<redshift_vector.size(); i++) {

    const double redshift = redshift_vector[i];
    std::vector<double> radius_rescaled (radius.size());
    for (size_t i_R=0; i_R<radius.size(); i_R++) radius_rescaled[i_R] = radius[i_R] * cosmo_ptr->D_A(0., redshift) / cosmo_fid_ptr->D_A(0., redshift);
    
    cosmology::HaloProfile halo_profile_cen (cosmo_ptr, false, redshift, 0., 0., Delta, profile_author, halo_def, Rt, false, false, 0., 0., b_e, s_e, AB_fact, bias_author, method_Pk, alpha_0, alpha_nu, beta, gamma_0);
    cosmology::HaloProfile halo_profile_off (cosmo_ptr, false, redshift, 0., 0., Delta, profile_author, halo_def, Rt, true, false, 0., 1., b_e, s_e, AB_fact, bias_author, method_Pk, alpha_0, alpha_nu, beta, gamma_0);

    if (two_halo) {

      halo_profile_cen.update_2halo(false);
      halo_profile_off.update_2halo(false);
      
    }

    const double delta_bkg = halo_profile_cen.Delta() / cosmo_ptr->Omega_matter(redshift);
    cosmology::Bias bias(cosmo_ptr);
    
    for (size_t j=0; j<Mass_vector.size(); j++) {

      const double mass = Mass_vector[j];
      halo_profile_cen.set_mass(mass);
      halo_profile_off.set_mass(mass);
      
      const double conc = conc_scaling_relation_func( c0, cM, cz, log10(mass/mass_pivot), log10((1+redshift)/(1+z_pivot)) );
      halo_profile_cen.set_concentration(conc);
      halo_profile_off.set_concentration(conc);
      
      const double sigma_off = sigma_off0 * pow(mass/mass_pivot, sigma_offM) * pow((1+redshift)/(1+z_pivot), sigma_offz);
      halo_profile_cen.set_sigma_off(sigma_off);
      halo_profile_off.set_sigma_off(sigma_off);

      if (two_halo) {
	
	const double bias_halo = bias.bias_halo(mass, redshift, bias_author, method_Pk, false, "test", "Linear", delta_bkg, -1., -1, 0.001, 100);
	halo_profile_cen.set_bias(bias_halo);
	halo_profile_off.set_bias(bias_halo);
	
      }

      DeltaSigma_cen_samples[i][j] = halo_profile_cen.DeltaSigma(radius_rescaled);
      DeltaSigma_off_samples[i][j] = halo_profile_off.DeltaSigma(radius_rescaled);

      Sigma_cen_samples[i][j] = halo_profile_cen.Sigma(radius_rescaled);
      Sigma_off_samples[i][j] = halo_profile_off.Sigma(radius_rescaled);

    }
  }

  // Interpolate the profile. A 2D interpolation is performed for each radial point.
  std::vector<cbl::glob::FuncGrid2D> DeltaSigma_cen_interp(radius.size());
  std::vector<cbl::glob::FuncGrid2D> Sigma_cen_interp(radius.size());
  std::vector<cbl::glob::FuncGrid2D> DeltaSigma_off_interp(radius.size());
  std::vector<cbl::glob::FuncGrid2D> Sigma_off_interp(radius.size());
  
  for (size_t k=0; k<radius.size(); k++) {

    std::vector<std::vector<double>> DeltaSigma_cen_subsamples(redshift_vector.size(), std::vector<double>(Mass_vector.size()));
    std::vector<std::vector<double>> Sigma_cen_subsamples(redshift_vector.size(), std::vector<double>(Mass_vector.size()));
    std::vector<std::vector<double>> DeltaSigma_off_subsamples(redshift_vector.size(), std::vector<double>(Mass_vector.size()));
    std::vector<std::vector<double>> Sigma_off_subsamples(redshift_vector.size(), std::vector<double>(Mass_vector.size()));
    
    for (size_t j=0; j<Mass_vector.size(); j++) {
      for (size_t i=0; i<redshift_vector.size(); i++) {

	DeltaSigma_cen_subsamples[i][j] = DeltaSigma_cen_samples[i][j][k];
	DeltaSigma_off_subsamples[i][j] = DeltaSigma_off_samples[i][j][k];

	Sigma_cen_subsamples[i][j] = Sigma_cen_samples[i][j][k];
	Sigma_off_subsamples[i][j] = Sigma_off_samples[i][j][k];
	
      }
    }

    cbl::glob::FuncGrid2D DeltaSigma_interp_cen (redshift_vector, Mass_vector, DeltaSigma_cen_subsamples, "Cubic");
    cbl::glob::FuncGrid2D DeltaSigma_interp_off (redshift_vector, Mass_vector, DeltaSigma_off_subsamples, "Cubic");
    DeltaSigma_cen_interp[k] = DeltaSigma_interp_cen;
    DeltaSigma_off_interp[k] = DeltaSigma_interp_off;

    cbl::glob::FuncGrid2D Sigma_interp_cen (redshift_vector, Mass_vector, Sigma_cen_subsamples, "Cubic");
    cbl::glob::FuncGrid2D Sigma_interp_off (redshift_vector, Mass_vector, Sigma_off_subsamples, "Cubic");
    Sigma_cen_interp[k] = Sigma_interp_cen;
    Sigma_off_interp[k] = Sigma_interp_off;
        
  }

  std::vector<std::vector<cbl::glob::FuncGrid2D>> out = {DeltaSigma_cen_interp, DeltaSigma_off_interp, Sigma_cen_interp, Sigma_off_interp};

  return out;
}

// ===========================================================================================

cbl::glob::FuncGrid cbl::modelling::densityprofileaux::Sigma_crit_factor_interp(const std::vector<double> redshift_vector, const std::vector<double> zl_bin_edges_for_N_zs, const std::vector<double> zs, const std::vector<std::vector<double>> N_zs, const std::shared_ptr<cosmology::Cosmology> cosmo)
{
  auto cosmo_ptr = move(cosmo);
  
  const double G_newt = cbl::par::GN / pow(3.086e+16 * 1.e6, 3.) * 1.98855e30;
  const double c_light = cbl::par::cc * 1.e3 / (3.086e+16 * 1.e6);

  cbl::glob::FuncGrid out;

  if (zl_bin_edges_for_N_zs.size() > 2) {

    std::vector<double> inv_Sigma_crit_samples (zl_bin_edges_for_N_zs.size() - 1, 0.);
    std::vector<double> zl_bin_centres (zl_bin_edges_for_N_zs.size() - 1, 0.);
      
    for (size_t i=0; i<zl_bin_centres.size(); i++) {

      zl_bin_centres[i] = 0.5 * (zl_bin_edges_for_N_zs[i] + zl_bin_edges_for_N_zs[i+1]);
      double inv1_Sigma_crit = 0, inv2_Sigma_crit = 0;

      for (size_t i_zs=0; i_zs<zs.size(); i_zs++) {
	
	if (zs[i_zs] > zl_bin_centres[i]) {
	  inv1_Sigma_crit += cosmo_ptr->D_A(zl_bin_centres[i], zs[i_zs]) / cosmo_ptr->D_A(0., zs[i_zs]) * N_zs[i][i_zs];
	  inv2_Sigma_crit += pow(cosmo_ptr->D_A(zl_bin_centres[i], zs[i_zs]) / cosmo_ptr->D_A(0., zs[i_zs]), 2) * N_zs[i][i_zs];
	}

      }
  
      inv_Sigma_crit_samples[i] = inv2_Sigma_crit / inv1_Sigma_crit * 4 * cbl::par::pi * G_newt / (c_light * c_light) * cosmo_ptr->D_A(0., zl_bin_centres[i]) * 1.e12;

    }

    out = cbl::glob::FuncGrid(zl_bin_centres, inv_Sigma_crit_samples, "Spline");

  } else {

    std::vector<double> inv_Sigma_crit_samples (redshift_vector.size(), 0.);
    std::vector<double> zl_bin_centres = redshift_vector;

    for (size_t i=0; i<zl_bin_centres.size(); i++) {

      double inv1_Sigma_crit = 0, inv2_Sigma_crit = 0;

      for (size_t i_zs=0; i_zs<zs.size(); i_zs++) {
	
	if (zs[i_zs] > zl_bin_centres[i]) {
	  inv1_Sigma_crit += cosmo_ptr->D_A(zl_bin_centres[i], zs[i_zs]) / cosmo_ptr->D_A(0., zs[i_zs]) * N_zs[0][i_zs];
	  inv2_Sigma_crit += pow(cosmo_ptr->D_A(zl_bin_centres[i], zs[i_zs]) / cosmo_ptr->D_A(0., zs[i_zs]), 2) * N_zs[0][i_zs];
	}

      }
  
      inv_Sigma_crit_samples[i] = inv2_Sigma_crit / inv1_Sigma_crit * 4 * cbl::par::pi * G_newt / (c_light * c_light) * cosmo_ptr->D_A(0., zl_bin_centres[i]) * 1.e12;

    }

    out = cbl::glob::FuncGrid(zl_bin_centres, inv_Sigma_crit_samples, "Spline");

  }
  return out;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofileaux::compute_model_density_scaling_relation (const std::vector<double> redshifts, const std::vector<double> mass_proxies, const std::vector<double> radius, const std::string probe, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_off, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_off, const cbl::glob::FuncGrid inv_Sigma_crit_interp, const std::shared_ptr<cosmology::Cosmology> cosmology, const double f_off0, const double f_offM, const double f_offz, const double OB_fact, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double purity, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base, const std::vector<double> purity_background, const std::vector<int> N_zs_index)
{
  auto cosmo_ptr = move(cosmology);
  
  // »»»»»»»»»»»»»»»»»»»»
  // Define the integrand
  // »»»»»»»»»»»»»»»»»»»»

  // function computing P(M|proxy,z)
  std::shared_ptr<void> ptr;
  
  std::function<double(double, double, double)> P_M_proxy_z = 
    [&](double logM_norm, double z_tr, double proxy_tr) -> double {
      
      double mean = modelling::densityprofileaux::logMass(z_tr, proxy_tr, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);
      double sigma = modelling::densityprofileaux::scatter_intr_M(z_tr, proxy_tr, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);

      double result = 0.;
      if (sigma > 1.e-3)
	result = cbl::gaussian(logM_norm, ptr, {mean, sigma});
      else {
	if (std::abs(logM_norm - mean) <= 1.e-3)
	  result = 1.;
      }

      return result;
	
    };

  // integrand function
  cbl::glob::FuncGrid2D DeltaSigma_Rk_interp, Sigma_Rk_interp;
  double proxy_tr = 0, z_tr = 0;
  
  std::function<double(double)> integrand;

  if (probe == "DeltaSigma" && redshifts.size() > 1) {
    
    integrand = [&] (const double x)
		{   
		  // Halo profile
		  double mass = pow(logM_base,x)*mass_pivot;

		  double z_tr2 = z_tr;
		  if (z_tr2 < DeltaSigma_Rk_interp.xmin())
		    z_tr2 = DeltaSigma_Rk_interp.xmin();
		  if (z_tr2 > DeltaSigma_Rk_interp.xmax())
		    z_tr2 = DeltaSigma_Rk_interp.xmax();
		  
		  double profile = DeltaSigma_Rk_interp(z_tr2, mass);
		       
		  return profile * P_M_proxy_z(x, z_tr2, proxy_tr);
		};
    
  } else if (probe == "g" && redshifts.size() > 1) {

    integrand = [&] (const double x)
		{
		  // Halo profile
		  double mass = pow(logM_base,x)*mass_pivot;

		  double z_tr2 = z_tr;
		  if (z_tr2 < DeltaSigma_Rk_interp.xmin())
		    z_tr2 = DeltaSigma_Rk_interp.xmin();
		  if (z_tr2 > DeltaSigma_Rk_interp.xmax())
		    z_tr2 = DeltaSigma_Rk_interp.xmax();
		  
		  double profile = DeltaSigma_Rk_interp(z_tr2, mass) * inv_Sigma_crit_interp(z_tr2) / std::max(1 - Sigma_Rk_interp(z_tr2, mass) * inv_Sigma_crit_interp(z_tr2), 0.9999);
		  
		  return profile * P_M_proxy_z(x, z_tr2, proxy_tr);
		};

  } else
    ErrorCBL("Wrong declaration of the probe, or low number of objects (i.e., <= 1) in the sample!", "compute_model_density_scaling_relation", "Modelling_DensityProfile_Auxiliary.cpp");


  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
  // Compute the profile for each object in the sample.
  // If the number of objects is not too large (see below),
  // compute the exact expected profile for each object.
  // Otherwise, interpolate.
  // »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

  // mass proxy and redshift vectors for the interpolation
  const std::vector<double> proxy_vector = cbl::linear_bin_vector(40, cbl::Min(mass_proxies), cbl::Max(mass_proxies));

  const double z_min = cbl::Min(redshifts);
  const double z_max = cbl::Max(redshifts);
  const std::vector<double> redshift_vector2 = cbl::linear_bin_vector(std::max(4, (int)((z_max-z_min)/0.1)), z_min, z_max);

  // final profile vectors
  std::vector<double> profile_cen_expected(radius.size(), 0.);
  std::vector<double> profile_off_expected(radius.size(), 0.);

  if (redshifts.size() > proxy_vector.size() * redshift_vector2.size()) { // interpolate

    // Compute the expected profile in each point
    std::vector<std::vector<std::vector<double>>> profile_cen_expected_samples(redshift_vector2.size(), std::vector<std::vector<double>>(proxy_vector.size(), std::vector<double>(radius.size())));
    std::vector<std::vector<std::vector<double>>> profile_off_expected_samples(redshift_vector2.size(), std::vector<std::vector<double>>(proxy_vector.size(), std::vector<double>(radius.size())));

    for (size_t i=0; i<redshift_vector2.size(); i++) {
      for (size_t j=0; j<proxy_vector.size(); j++) {

	z_tr = redshift_vector2[i];
	proxy_tr = proxy_vector[j];
	std::vector<double> int_limits = modelling::densityprofileaux::min_max_logMass(z_tr, z_tr, proxy_tr, proxy_tr, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);
	
	// Expected profiles
	for (size_t k=0; k<radius.size(); k++) {

	  // centred profile
	  DeltaSigma_Rk_interp = DeltaSigma_interp_cen[k];
	  Sigma_Rk_interp = Sigma_interp_cen[k];

	  profile_cen_expected_samples[i][j][k] += wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0], int_limits[1]);

	  // miscentred profile
	  DeltaSigma_Rk_interp = DeltaSigma_interp_off[k];
	  Sigma_Rk_interp = Sigma_interp_off[k];

	  profile_off_expected_samples[i][j][k] += wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0], int_limits[1]);
    
	}
	
      }

    }

    // Interpolation of the expected profiles
    std::vector<cbl::glob::FuncGrid2D> profile_cen_expected_interp(radius.size());
    std::vector<cbl::glob::FuncGrid2D> profile_off_expected_interp(radius.size());
    
    for (size_t k=0; k<radius.size(); k++) {

      std::vector<std::vector<double>> profile_cen_expected_subsamples(redshift_vector2.size(), std::vector<double>(proxy_vector.size()));
      std::vector<std::vector<double>> profile_off_expected_subsamples(redshift_vector2.size(), std::vector<double>(proxy_vector.size()));
      
      for (size_t j=0; j<proxy_vector.size(); j++) {
	for (size_t i=0; i<redshift_vector2.size(); i++) {

	  profile_cen_expected_subsamples[i][j] = profile_cen_expected_samples[i][j][k];
	  profile_off_expected_subsamples[i][j] = profile_off_expected_samples[i][j][k];

	}
      }

      cbl::glob::FuncGrid2D interp_cen (redshift_vector2, proxy_vector, profile_cen_expected_subsamples, "Cubic");
      cbl::glob::FuncGrid2D interp_off (redshift_vector2, proxy_vector, profile_off_expected_subsamples, "Cubic");
      profile_cen_expected_interp[k] = interp_cen;
      profile_off_expected_interp[k] = interp_off;

    }

    // Compute the profile for each object
    for (size_t i_obj=0; i_obj<mass_proxies.size(); i_obj++) {

      // avoid roundoff issues
      double zz=0, mass_proxy=0;

      if (redshifts[i_obj] < cbl::Min(redshift_vector2))
	zz = cbl::Min(redshift_vector2);
      else if (redshifts[i_obj] > cbl::Max(redshift_vector2))
	zz = cbl::Max(redshift_vector2);
      else
	zz = redshifts[i_obj];

      if (mass_proxies[i_obj] < cbl::Min(proxy_vector))
	mass_proxy = cbl::Min(proxy_vector);
      else if (mass_proxies[i_obj] > cbl::Max(proxy_vector))
	mass_proxy = cbl::Max(proxy_vector);
      else
	mass_proxy = mass_proxies[i_obj];
      
      for (size_t k=0; k<radius.size(); k++) {

	// compute the models
	if (probe == "DeltaSigma") {
	  profile_cen_expected[k] += profile_cen_expected_interp[k](zz, mass_proxy);
	  profile_off_expected[k] += profile_off_expected_interp[k](zz, mass_proxy);
	} else {
	  profile_cen_expected[k] += profile_cen_expected_interp[k](zz, mass_proxy) * purity_background[N_zs_index[i_obj]];
	  profile_off_expected[k] += profile_off_expected_interp[k](zz, mass_proxy) * purity_background[N_zs_index[i_obj]];
	}
	
      }

    }

  }

  else { // do not interpolate

    // Loop over the objects in the cluster sample
    for (size_t i_obj=0; i_obj<mass_proxies.size(); i_obj++) {

      z_tr = redshifts[i_obj];
      proxy_tr = mass_proxies[i_obj];
      std::vector<double> int_limits = modelling::densityprofileaux::min_max_logMass(z_tr, z_tr, proxy_tr, proxy_tr, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, proxy_pivot, logM_base);
      
      // Expected profiles    
      for (size_t k=0; k<radius.size(); k++) {

	// centred profile
	DeltaSigma_Rk_interp = DeltaSigma_interp_cen[k];
	Sigma_Rk_interp = Sigma_interp_cen[k];

	if (probe == "DeltaSigma")
	  profile_cen_expected[k] += wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0], int_limits[1]);
	else
	  profile_cen_expected[k] += wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0], int_limits[1]) * purity_background[N_zs_index[i_obj]];

	// miscentred profile
	DeltaSigma_Rk_interp = DeltaSigma_interp_off[k];
	Sigma_Rk_interp = Sigma_interp_off[k];

	if (probe == "DeltaSigma")
	  profile_off_expected[k] += wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0], int_limits[1]);
	else
	  profile_off_expected[k] += wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0], int_limits[1]) * purity_background[N_zs_index[i_obj]];
    
      }

    }

  }
  
  // Average profile
  const double f_off = f_off0 * pow(cbl::Average(mass_proxies) / proxy_pivot, f_offM) * pow((1+cbl::Average(redshifts))/(1+z_pivot), f_offz);

  std::vector<double> profile_expected(radius.size(), 0.);
  for (size_t k=0; k<radius.size(); k++)
    profile_expected[k] = (1 + OB_fact) * purity / redshifts.size() * ((1 - f_off) * profile_cen_expected[k] + f_off * profile_off_expected[k]);
  
  return profile_expected;
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofileaux::compute_model_density_scaling_relation (std::string MF_author, std::string halo_def, double overdensity, const double z_min, const double z_max, const double mass_proxy_min, const double mass_proxy_max, const double z_error, const std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc, const double Plambda_A_mu, const double Plambda_B_mu, const double Plambda_C_mu, const std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc, const double Plambda_A_sigma, const double Plambda_B_sigma, const double Plambda_C_sigma, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const std::vector<double> radius, const std::string probe, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_off, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_off, const cbl::glob::FuncGrid inv_Sigma_crit_interp, const std::shared_ptr<cosmology::Cosmology> cosmology, const double f_off0, const double f_offM, const double f_offz, const double OB_fact, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double purity, cbl::glob::FuncGrid2D completeness_interp, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base, cbl::glob::FuncGrid purity_background_interp, const std::vector<double> purity_background_interp_range)
{  
  if (f_offM != 0 || f_offz != 0)
    ErrorCBL("In the current implementation, f_off must be constant!", "compute_model_density_scaling_relation", "Modelling_DensityProfile_Auxiliary.cpp");
  
  auto cosmo_ptr = move(cosmology);
  cosmology::MassFunction MF(cosmo_ptr);
  double fact = (cosmo_ptr->unit()) ? 1 : cosmo_ptr->little_h();
  
  bool isDelta_critical = false;
  if (halo_def == "critical")
    isDelta_critical = true;
  else if (halo_def == "mean")
    isDelta_critical = false;
  else
    ErrorCBL("Wrong declaration of the overdensity definition!", "compute_model_density_scaling_relation", "Modelling_DensityProfile_Auxiliary.cpp");


  // »»»»»»»»»»»»»»»»»»»»
  //  Integration limits
  // »»»»»»»»»»»»»»»»»»»»

  std::vector<std::vector<double>> int_limits = cbl::modelling::numbercounts::get_integration_limits
    (
     cosmology, z_min, z_max, mass_proxy_min, mass_proxy_max, z_error,
     Plambda_mean_fc, Plambda_A_mu, Plambda_B_mu, Plambda_C_mu, Plambda_std_fc, Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma,
     scaling_relation, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, z_pivot, proxy_pivot, mass_pivot, logM_base
     ); // log(M/M_piv), proxy_tr, z_tr. In this order.
  
  
  // »»»»»»»»»»»»»»»»»»»»
  // Define the integrand
  // »»»»»»»»»»»»»»»»»»»»

  std::shared_ptr<void> ptr;


  // -------------- interpolate the integral over proxy_tr as a function of M and z_tr --------------

  double dummy_logM_norm = 0, dummy_z_tr = 0;
  std::function<double(double)> proxy_tr_dependent_integrand =
    [&](double proxy_tr) -> double {
      
      // P(proxy_tr|M,z)
      double log_proxy_tr = log(proxy_tr/proxy_pivot) / log(logM_base);
      double mean = modelling::numbercounts::logProxy(dummy_z_tr, pow(logM_base,dummy_logM_norm)*mass_pivot, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base);
      double sigma = modelling::numbercounts::scatter_intr_proxy(dummy_z_tr, pow(logM_base,dummy_logM_norm)*mass_pivot, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base);

      double P_proxy_M_z = 0.;
      if (sigma > 1.e-3)
	P_proxy_M_z = cbl::gaussian(log_proxy_tr, ptr, {mean, sigma}) / (log(logM_base) * proxy_tr);
      else {
	if (std::abs(log_proxy_tr - mean) <= 1.e-3)
	  P_proxy_M_z = 1.;
      }

      // Integral of P(z_ob|z_tr)
      double mean_Pz = dummy_z_tr;
      double std_Pz = z_error * (1 + dummy_z_tr);
      double int_Pz = 0.5 * ( erf( (z_max - mean_Pz) / (sqrt(2) * std_Pz) ) - erf( (z_min - mean_Pz) / (sqrt(2) * std_Pz) ) );

      // Integral of P(proxy_ob|proxy_tr,z_tr)
      double mean_Pproxy = Plambda_mean_fc(proxy_tr, dummy_z_tr, Plambda_A_mu, Plambda_B_mu, Plambda_C_mu);
      double std_Pproxy = Plambda_std_fc(proxy_tr, dummy_z_tr, Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma);
      double int_Pproxy = 0.5 * ( erf( (mass_proxy_max - mean_Pproxy) / (sqrt(2) * std_Pproxy) ) - erf( (mass_proxy_min - mean_Pproxy) / (sqrt(2) * std_Pproxy) ) );

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

  const std::vector<double> redshift_vector = cbl::linear_bin_vector(20, int_limits[2][0], int_limits[2][1]);
  const std::vector<double> logM_norm_vector = cbl::linear_bin_vector(50, int_limits[0][0], int_limits[0][1]);
  std::vector<std::vector<double>> int_lambda_tr_samples(redshift_vector.size(), std::vector<double>(logM_norm_vector.size()));

  for (size_t i=0; i<redshift_vector.size(); i++) {
    for (size_t j=0; j<logM_norm_vector.size(); j++) {
      dummy_z_tr = redshift_vector[i];
      dummy_logM_norm = logM_norm_vector[j];
      int_lambda_tr_samples[i][j] = wrapper::gsl::GSL_integrate_cquad(proxy_tr_dependent_integrand, int_limits[1][0], int_limits[1][1]);
    }
  }

  cbl::glob::FuncGrid2D int_lambda_tr_interp (redshift_vector, logM_norm_vector, int_lambda_tr_samples, "Cubic");

  // -------------- function returning the integrand function for the counts --------------
  std::function<double(double)> counts_integrand =
    [&] (const double logM_norm) {

      double z_tr = dummy_z_tr;
      double Delta = (isDelta_critical) ? overdensity/cosmo_ptr->Omega_matter(z_tr) : overdensity;
      double mass = pow(logM_base, logM_norm) * mass_pivot;

      double counts_norm = cosmo_ptr->dV_dZdOmega(z_tr, true) * MF.mass_function(mass, interp_sigmaM(mass*fact), interp_DlnsigmaM(mass*fact), z_tr, interp_DN(z_tr), MF_author, false, par::defaultString, Delta) * (log(logM_base) * mass); // the last factor stems from the change of integration variable

      return counts_norm * int_lambda_tr_interp(z_tr, logM_norm);
            
  };

  // -------------- define the base integrand function for the numerator, that will be used to interpolate the integral over M --------------
  cbl::glob::FuncGrid2D DeltaSigma_cen_Rk_interp, DeltaSigma_off_Rk_interp, Sigma_cen_Rk_interp, Sigma_off_Rk_interp;

  std::function<double(double)> integrand;

  if (probe == "DeltaSigma") {
    
    integrand = [&] (const double logM_norm)
		{
		  double M = pow(logM_base, logM_norm) * mass_pivot;
		  double profile = (1. - f_off0) * DeltaSigma_cen_Rk_interp(dummy_z_tr, M) + f_off0 * DeltaSigma_off_Rk_interp(dummy_z_tr, M);
		  
		  return profile * counts_integrand(logM_norm);
		};
    
  } else if (probe == "g") {

    integrand = [&] (const double logM_norm)
		{
		  double M = pow(logM_base, logM_norm) * mass_pivot;
		  double z_tr = dummy_z_tr;
		  double profile_cen = DeltaSigma_cen_Rk_interp(z_tr, M) * inv_Sigma_crit_interp(z_tr) / std::max(1 - Sigma_cen_Rk_interp(z_tr, M) * inv_Sigma_crit_interp(z_tr), 0.9999);
		  double profile_off = DeltaSigma_off_Rk_interp(z_tr, M) * inv_Sigma_crit_interp(z_tr) / std::max(1 - Sigma_off_Rk_interp(z_tr, M) * inv_Sigma_crit_interp(z_tr), 0.9999);
		  double profile = (1. - f_off0) * profile_cen + f_off0 * profile_off;

		  if (z_tr < purity_background_interp_range[0])
		    z_tr = purity_background_interp_range[0];
		  else if (z_tr > purity_background_interp_range[1])
		    z_tr = purity_background_interp_range[1];
		       
		  return purity_background_interp(z_tr) * profile * counts_integrand(logM_norm);
		};

  } else
    ErrorCBL("Wrong declaration of the probe!", "compute_model_density_scaling_relation", "Modelling_DensityProfile_Auxiliary.cpp");


  // -------------- define the integrands used in the final computation --------------
  cbl::glob::FuncGrid int_over_M_proxy_numerator_interp, int_over_M_proxy_denominator_interp;

  std::function<double(double)> integrand_numerator =
    [&] (const double z_tr) {

      return int_over_M_proxy_numerator_interp(z_tr);
            
  };

  std::function<double(double)> integrand_denominator =
    [&] (const double z_tr) {

      return int_over_M_proxy_denominator_interp(z_tr);
            
  };
  
  
  // »»»»»»»»»»»»»»»»»»»»
  //      Integrate
  // »»»»»»»»»»»»»»»»»»»»

  std::vector<double> profile_expected(radius.size(), 0.);
  for (size_t k=0; k<radius.size(); k++) {

    // 
    DeltaSigma_cen_Rk_interp = DeltaSigma_interp_cen[k];
    Sigma_cen_Rk_interp = Sigma_interp_cen[k];
    DeltaSigma_off_Rk_interp = DeltaSigma_interp_off[k];
    Sigma_off_Rk_interp = Sigma_interp_off[k];

    // interpolate the numerator's integral over M, as a function of z
    std::vector<double> int_over_M_proxy_numerator_samples(redshift_vector.size());
    for (size_t i=0; i<redshift_vector.size(); i++) {
      dummy_z_tr = redshift_vector[i];
      int_over_M_proxy_numerator_samples[i] = wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0][0], int_limits[0][1]);
    }
    int_over_M_proxy_numerator_interp = cbl::glob::FuncGrid(redshift_vector, int_over_M_proxy_numerator_samples, "Spline");

    // final integral for the numerator
    if (probe == "DeltaSigma")
      profile_expected[k] += wrapper::gsl::GSL_integrate_qag(integrand_numerator, int_limits[2][0], int_limits[2][1]);
    else
      profile_expected[k] += wrapper::gsl::GSL_integrate_qag(integrand_numerator, int_limits[2][0], int_limits[2][1]);
    
  }

  // interpolate the denominator's integral over M, as a function of z
  std::vector<double> int_over_M_proxy_denominator_samples(redshift_vector.size());
  for (size_t i=0; i<redshift_vector.size(); i++) {
    dummy_z_tr = redshift_vector[i];
    int_over_M_proxy_denominator_samples[i] = wrapper::gsl::GSL_integrate_qag(counts_integrand, int_limits[0][0], int_limits[0][1]);
  }
  int_over_M_proxy_denominator_interp = cbl::glob::FuncGrid(redshift_vector, int_over_M_proxy_denominator_samples, "Spline");

  const double denominator = wrapper::gsl::GSL_integrate_qag(integrand_denominator, int_limits[2][0], int_limits[2][1]);
  
  // Average profile  
  for (size_t k=0; k<radius.size(); k++)
    profile_expected[k] *= (1 + OB_fact) * purity / denominator;
  
  return profile_expected;  
}

// ===========================================================================================

std::vector<double> cbl::modelling::densityprofileaux::compute_mass_expected_value (const std::vector<double> z_ob, const std::vector<double> mass_proxy_ob, std::string MF_author, std::string halo_def, double overdensity, const double z_error, const std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc, const double Plambda_A_mu, const double Plambda_B_mu, const double Plambda_C_mu, const std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc, const double Plambda_A_sigma, const double Plambda_B_sigma, const double Plambda_C_sigma, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const std::shared_ptr<cosmology::Cosmology> cosmology, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, cbl::glob::FuncGrid2D completeness_interp, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base)
{
  if (z_ob.size() != mass_proxy_ob.size())
    ErrorCBL("The sizes of observed redshift and mass proxy vectors must be identical!", "mass_expected_value", "Modelling_DensityProfile_Auxiliary.cpp");
  
  auto cosmo_ptr = move(cosmology);
  cosmology::MassFunction MF(cosmo_ptr);
  double fact = (cosmo_ptr->unit()) ? 1 : cosmo_ptr->little_h();
  
  bool isDelta_critical = false;
  if (halo_def == "critical")
    isDelta_critical = true;
  else if (halo_def == "mean")
    isDelta_critical = false;
  else
    ErrorCBL("Wrong declaration of the overdensity definition!", "mass_expected_value", "Modelling_DensityProfile_Auxiliary.cpp");


  // -------------------------------------------------------------------------------------------
  // ------------------- First, interpolate mass values vs z_ob and proxy_ob -------------------
  // -------------------------------------------------------------------------------------------
  
  // Define z_ob and proxy_ob vectors used for interpolation
  const std::vector<double> proxy_vec = cbl::logarithmic_bin_vector(50, 0.99 * cbl::Min(mass_proxy_ob), 1.01 * cbl::Max(mass_proxy_ob));
  const std::vector<double> redshift_vec = cbl::linear_bin_vector(20, 0.99 * cbl::Min(z_ob), 1.01 * cbl::Max(z_ob));
  
  std::vector<std::vector<double>> MASS_POINTS (redshift_vec.size(), std::vector<double> (proxy_vec));

  for (size_t i_z=0; i_z<redshift_vec.size(); i_z++) {
    for (size_t i_p=0; i_p<proxy_vec.size(); i_p++) {

      // »»»»»»»»»»»»»»»»»»»»
      //  Integration limits
      // »»»»»»»»»»»»»»»»»»»»

      std::vector<std::vector<double>> int_limits = cbl::modelling::numbercounts::get_integration_limits
	(
	 cosmology, redshift_vec[i_z], redshift_vec[i_z], proxy_vec[i_p], proxy_vec[i_p], z_error,
	 Plambda_mean_fc, Plambda_A_mu, Plambda_B_mu, Plambda_C_mu, Plambda_std_fc, Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma,
	 scaling_relation, alpha, beta, gamma, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, z_pivot, proxy_pivot, mass_pivot, logM_base
	 ); // log(M/M_piv), proxy_tr, z_tr. In this order.
  
  
      // »»»»»»»»»»»»»»»»»»»»
      // Define the integrand
      // »»»»»»»»»»»»»»»»»»»»

      std::shared_ptr<void> ptr;


      // -------------- interpolate the integral over proxy_tr as a function of M and z_tr --------------

      double dummy_logM_norm = 0, dummy_z_tr = 0;
      std::function<double(double)> proxy_tr_dependent_integrand =
	[&](double proxy_tr) -> double {
      
	  // P(proxy_tr|M,z)
	  double log_proxy_tr = log(proxy_tr/proxy_pivot) / log(logM_base);
	  double mean = modelling::numbercounts::logProxy(dummy_z_tr, pow(logM_base,dummy_logM_norm)*mass_pivot, alpha, beta, gamma, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base);
	  double sigma = modelling::numbercounts::scatter_intr_proxy(dummy_z_tr, pow(logM_base,dummy_logM_norm)*mass_pivot, scatter0, scatterM, scatterM_exp, scatterz, scatterz_exp, cosmo_ptr, scaling_relation, z_pivot, mass_pivot, logM_base);

	  double P_proxy_M_z = 0.;
	  if (sigma > 1.e-3)
	    P_proxy_M_z = cbl::gaussian(log_proxy_tr, ptr, {mean, sigma}) / (log(logM_base) * proxy_tr);
	  else {
	    if (std::abs(log_proxy_tr - mean) <= 1.e-3)
	      P_proxy_M_z = 1.;
	  }

	  // P(z_ob|z_tr)
	  double mean_Pz = dummy_z_tr;
	  double std_Pz = z_error * (1 + dummy_z_tr);
	  double Pz = cbl::gaussian(redshift_vec[i_z], ptr, {mean_Pz, std_Pz});

	  // P(proxy_ob|proxy_tr,z_tr)
	  double mean_Pproxy = Plambda_mean_fc(proxy_tr, dummy_z_tr, Plambda_A_mu, Plambda_B_mu, Plambda_C_mu);
	  double std_Pproxy = Plambda_std_fc(proxy_tr, dummy_z_tr, Plambda_A_sigma, Plambda_B_sigma, Plambda_C_sigma);
	  double Pproxy = cbl::gaussian(proxy_vec[i_p], ptr, {mean_Pproxy, std_Pproxy});

	  double ztr2 = dummy_z_tr, proxytr2 = proxy_tr;
	  if (ztr2 < completeness_interp.xmin())
	    ztr2 = completeness_interp.xmin();
	  if (ztr2 > completeness_interp.xmax())
	    ztr2 = completeness_interp.xmax();
	  if (proxytr2 < completeness_interp.ymin())
	    proxytr2 = completeness_interp.ymin();
	  if (proxytr2 > completeness_interp.ymax())
	    proxytr2 = completeness_interp.ymax();

	  return P_proxy_M_z * Pz * Pproxy * std::max(0., completeness_interp(ztr2, proxytr2));
            
	};

      const std::vector<double> redshift_vector = cbl::linear_bin_vector(20, int_limits[2][0], int_limits[2][1]);
      const std::vector<double> logM_norm_vector = cbl::linear_bin_vector(50, int_limits[0][0], int_limits[0][1]);
      std::vector<std::vector<double>> int_lambda_tr_samples(redshift_vector.size(), std::vector<double>(logM_norm_vector.size()));

      for (size_t i=0; i<redshift_vector.size(); i++) {
	for (size_t j=0; j<logM_norm_vector.size(); j++) {
	  dummy_z_tr = redshift_vector[i];
	  dummy_logM_norm = logM_norm_vector[j];
	  int_lambda_tr_samples[i][j] = wrapper::gsl::GSL_integrate_cquad(proxy_tr_dependent_integrand, int_limits[1][0], int_limits[1][1]);
	}
      }

      cbl::glob::FuncGrid2D int_lambda_tr_interp (redshift_vector, logM_norm_vector, int_lambda_tr_samples, "Cubic");

      // -------------- function returning the integrand function for the counts --------------
      std::function<double(double)> counts_integrand =
	[&] (const double logM_norm) {

	  double z_tr = dummy_z_tr;
	  double Delta = (isDelta_critical) ? overdensity/cosmo_ptr->Omega_matter(z_tr) : overdensity;
	  double mass = pow(logM_base, logM_norm) * mass_pivot;

	  double counts_norm = cosmo_ptr->dV_dZdOmega(z_tr, true) * MF.mass_function(mass, interp_sigmaM(mass*fact), interp_DlnsigmaM(mass*fact), z_tr, interp_DN(z_tr), MF_author, false, par::defaultString, Delta) * (log(logM_base) * mass); // the last factor stems from the change of integration variable

	  return counts_norm * int_lambda_tr_interp(z_tr, logM_norm);
            
	};

      // -------------- define the base integrand function for the numerator, that will be used to interpolate the integral over M --------------
    
      std::function<double(double)> integrand =
	[&] (const double logM_norm)
	{
	  return logM_norm * counts_integrand(logM_norm);
	};

      // -------------- define the integrands used in the final computation --------------
      cbl::glob::FuncGrid int_over_M_proxy_numerator_interp, int_over_M_proxy_denominator_interp;

      std::function<double(double)> integrand_numerator =
	[&] (const double z_tr) {

	  return int_over_M_proxy_numerator_interp(z_tr);
            
	};

      std::function<double(double)> integrand_denominator =
	[&] (const double z_tr) {

	  return int_over_M_proxy_denominator_interp(z_tr);
            
	};
  
  
      // »»»»»»»»»»»»»»»»»»»»
      //      Integrate
      // »»»»»»»»»»»»»»»»»»»»

      // interpolate the numerator's integral over M, as a function of z
      std::vector<double> int_over_M_proxy_numerator_samples(redshift_vector.size());
      for (size_t i=0; i<redshift_vector.size(); i++) {
	dummy_z_tr = redshift_vector[i];
	int_over_M_proxy_numerator_samples[i] = wrapper::gsl::GSL_integrate_qag(integrand, int_limits[0][0], int_limits[0][1]);
      }
      int_over_M_proxy_numerator_interp = cbl::glob::FuncGrid(redshift_vector, int_over_M_proxy_numerator_samples, "Spline");

      MASS_POINTS[i_z][i_p] = wrapper::gsl::GSL_integrate_qag(integrand_numerator, int_limits[2][0], int_limits[2][1]);
    

      // interpolate the denominator's integral over M, as a function of z
      std::vector<double> int_over_M_proxy_denominator_samples(redshift_vector.size());
      for (size_t i=0; i<redshift_vector.size(); i++) {
	dummy_z_tr = redshift_vector[i];
	int_over_M_proxy_denominator_samples[i] = wrapper::gsl::GSL_integrate_qag(counts_integrand, int_limits[0][0], int_limits[0][1]);
      }
      int_over_M_proxy_denominator_interp = cbl::glob::FuncGrid(redshift_vector, int_over_M_proxy_denominator_samples, "Spline");

      const double denominator = wrapper::gsl::GSL_integrate_qag(integrand_denominator, int_limits[2][0], int_limits[2][1]);
  
      // Mass expected value  
      MASS_POINTS[i_z][i_p] = pow(logM_base, MASS_POINTS[i_z][i_p] / denominator);

    }
  }

  cbl::glob::FuncGrid2D MASS_INTERP (redshift_vec, proxy_vec, MASS_POINTS, "Cubic");


  // -------------------------------------------------------------------------------------------
  // -------------------------- Finally derive the mass of each object -------------------------
  // -------------------------------------------------------------------------------------------
  
  std::vector<double> MASS (z_ob.size(), 0.);
  for (size_t i=0; i<z_ob.size(); i++) MASS[i] = MASS_INTERP(z_ob[i], mass_proxy_ob[i]);
  
  return MASS;  
}

