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
 *  @file Headers/Modelling_DensityProfile_Auxiliary.h
 *
 *  @brief Collection of auxiliary functions used in the class Modelling_DensityProfile
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */

#ifndef __MODELLINGDPROFILEAUX__
#define __MODELLINGDPROFILEAUX__

#include "Bias.h"
#include "HaloProfile.h"
#include "Modelling_MassObservableRelation.h"
#include "ModelFunction_NumberCounts.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace densityprofileaux {
      
      /**
       * @brief Function computing the splashback radius as the minimum of 
       * \f$ {\rm d}\log\rho \, / \, {\rm d}\log R \f$
       *
       * @param rho the 3D mass density
       *
       * @param radius 
       *
       * @return the splashback radius
       *
       */
      double splashback_radius (const std::vector<double> rho, const std::vector<double> radius);
      
      /**
       * @brief function computing the logarithm of the 
       * normalised mass derived from the scaling relation
       *
       * @param z_tr true redshift
       *
       * @param proxy_tr true mass proxy
       *
       * @param alpha normalisation of the scaling relation
       *
       * @param beta slope of the scaling relation
       *
       * @param gamma scaling relation redshift evolution
       * parameter
       *
       * @param cosmo Cosmology pointer
       *
       * @param scaling_relation mass-observable relation pointer
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param proxy_pivot proxy pivot in the scaling relation
       *
       * @param logM_base base of the mass logarithm
       *
       * @return logarithm of the normalised mass derived from the scaling relation
       *
       */
      double logMass (const double z_tr, const double proxy_tr, const double alpha, const double beta, const double gamma, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double logM_base);
      
      /**
       * @brief function computing the intrinsic scatter
       * of P(M|proxy,z)
       *
       * @param z_tr true redshift
       *
       * @param proxy_tr true mass proxy
       *
       * @param scatter0 normalisation
       *
       * @param scatterM mass evolution
       *
       * @param scatterM_exp exponent of the mass evolution
       *
       * @param scatterz z evolution
       *
       * @param scatterz_exp exponent of the z evolution
       *
       * @param cosmo Cosmology pointer
       *
       * @param scaling_relation mass-observable relation pointer
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param proxy_pivot proxy pivot in the scaling relation
       *
       * @param logM_base base of the mass logarithm
       *
       * @return intrinsic scatter
       *
       */
      double scatter_intr_M (const double z_tr, const double proxy_tr, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double logM_base);
      
      /**
       * @brief function computing the minimum and maximum logarithm
       * of the normalised mass from the mean of P(M|proxy,z)
       *
       * @param z_min minimum redshift
       *
       * @param z_max maximum redshift
       *
       * @param proxy_min minimum mass proxy
       *
       * @param proxy_max maximum mass proxy
       *
       * @param alpha normalisation of the scaling relation
       *
       * @param beta slope of the scaling relation
       *
       * @param gamma scaling relation redshift evolution
       * parameter
       *
       * @param scatter0 normalisation
       *
       * @param scatterM mass evolution
       *
       * @param scatterM_exp exponent of the mass evolution
       *
       * @param scatterz z evolution
       *
       * @param scatterz_exp exponent of the z evolution
       *
       * @param cosmo Cosmology pointer
       *
       * @param scaling_relation mass-observable relation pointer
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param proxy_pivot proxy pivot in the scaling relation
       *
       * @param logM_base base of the mass logarithm
       *
       * @return log-mass range
       *
       */
      std::vector<double> min_max_logMass (const double z_min, const double z_max, const double proxy_min, const double proxy_max, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double logM_base);
      
      /**
       * @brief function interpolating the halo profiles
       *
       * @param radius radial points
       *
       * @param redshift_vector redshift points used for the interpolation
       *
       * @param Mass_vector mass points used for interpolation
       *
       * @param cosmo Cosmology pointer
       *
       * @param profile_author author(s) of the cluster density profile.
       * See available options in cbl::cosmology::HaloProfile
       *
       * @param two_halo if true, include the 2-halo contribution
       *
       * @param halo_def the halo definition; available options are:
       * "critical", "vir", "mean"
       *
       * @param Delta overdensity factor which needs to be multiplied 
       * to the critical density in order to define an overdensity
       *
       * @param AB_fact \f$\sigma_{\rm AB}\f$
       *
       * @param bias_author author(s) of the halo bias function
       *
       * @param method_Pk method used for the computation of the
       * power spectrum
       *
       * @param conc_scaling_relation_func function returning the concentration
       *
       * @param c0 \f$c_0\f$
       *
       * @param cM \f$c_M\f$
       *
       * @param cz \f$c_z\f$
       *
       * @param sigma_off0 \f$ \sigma_{\rm off,0}\f$
       *
       * @param sigma_offM \f$ \sigma_{\rm off,M} \f$
       *
       * @param sigma_offz \f$ \sigma_{\rm off,z}\f$
       *
       * @param Rt \f$F_t\f$
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param cosmo_fid pointer to the fiducial cosmology. If it
       * is not NULL, the radii are rescaled to account for geometric distortions
       *
       *  @param b_e normalisation of the DK14 outer profile
       *
       *  @param s_e slope of the DK14 outer profile
       *
       *  @param alpha_0 \f$\alpha_0\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_0=0.155\f$ by default.
       *
       *  @param alpha_nu \f$\alpha_\nu\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_\nu=0.0095\f$ by default.
       *
       *  @param beta \f$\beta\f$ parameter in the DK14 profile. Following
       *  DK14, \f$\beta=4\f$ by default.
       *
       *  @param gamma_0 \f$\gamma_0\f$ parameter in the expression
       *  \f$\gamma=\gamma_0\nu\f$ appearing in the DK14 profile. Following
       *  DK14, \f$\gamma_0=4\f$ by default.
       *
       * @return in this order: fully-centred DeltaSigma, fully-miscentred
       * DeltaSigma, fully-centred Sigma, fully-miscentred Sigma
       *
       */
       std::vector<std::vector<cbl::glob::FuncGrid2D>> interpolated_profiles (const std::vector<double> radius, const std::vector<double> redshift_vector, const std::vector<double> Mass_vector, const std::shared_ptr<cosmology::Cosmology> cosmo, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const double AB_fact, const std::string bias_author, const std::string method_Pk, std::function<double(const double, const double, const double, const double, const double)> conc_scaling_relation_func, const double c0, const double cM, const double cz, const double sigma_off0, const double sigma_offM, const double sigma_offz, const double Rt, const double z_pivot, const double mass_pivot, const std::shared_ptr<cosmology::Cosmology> cosmo_fid=NULL, const double b_e=1., const double s_e=1.5, const double alpha_0=0.155, const double alpha_nu=0.0095, const double beta=4., const double gamma_0=4.);
       
       /**
       * @brief function computing the interpolation of
       *
       * \f$\langle\Sigma_{\rm crit}^{-1}\rangle^{-1}\langle\Sigma_{\rm crit}^{-2}\rangle\f$
       *
       * @param redshift_vector cluster redshift points where the function
       * is interpolated
       *
       * @param zl_bin_edges_for_N_zs lens redshift bin edges where
       * the \f$N(z_{\rm s})\f$ distributions are defined
       *
       * @param zs central bin values where the binned \f$N(z_{\rm s})\f$
       * distribution is computed
       *
       * @param N_zs binned \f$N(z_{\rm s})\f$ distributions
       *
       * @param cosmo Cosmology pointer
       *
       * @return interpolated inverse critical surface density
       *
       */
       cbl::glob::FuncGrid Sigma_crit_factor_interp(const std::vector<double> redshift_vector, const std::vector<double> zl_bin_edges_for_N_zs, const std::vector<double> zs, const std::vector<std::vector<double>> N_zs, const std::shared_ptr<cosmology::Cosmology> cosmo);
       
       /**
       * @brief Compute the excess density profile model
       * in all the radial bins, with mass derived from the 
       * scaling relation. The average profile is weighted over the
       * observed distribution of objects
       *
       * @param redshifts vector containing the redshift values of the clusters
       * in the sample
       *
       * @param mass_proxies vector containing the mass proxy values of the clusters
       * in the sample
       *
       * @param radius the radius array
       *
       * @param probe "DeltaSigma" \f$\rightarrow\f$ excess surface density; 
       *  "g" \f$\rightarrow\f$ reduced shear
       *
       * @param DeltaSigma_interp_cen interpolated centred DeltaSigma
       *
       * @param DeltaSigma_interp_off interpolated centred DeltaSigma
       *
       * @param Sigma_interp_cen interpolated centred Sigma
       *
       * @param Sigma_interp_off interpolated centred Sigma
       *
       * @param inv_Sigma_crit_interp interpolated inverse critical density
       *
       * @param cosmology cosmological model
       *
       * @param f_off0 \f$ f_{\rm off,0} \f$
       *
       * @param f_offM \f$ f_{\rm off,M} \f$
       *
       * @param f_offz \f$ f_{\rm off,z} \f$
       *
       * @param OB_fact orientation boost factor, \f$\mathcal{B}_O\f$.
       * The average profile is multiplied by \f$(1+\mathcal{B}_O)\f$
       *
       * @param alpha normalisation of the scaling relation
       *
       * @param beta slope of the scaling relation
       *
       * @param gamma scaling relation redshift evolution
       * parameter
       *
       * @param scatter0 scatter normalisation
       *
       * @param scatterM scatter mass proxy dependence parameter
       *
       * @param scatterM_exp scatter mass proxy exponent parameter
       *
       * @param scatterz scatter redshift dependence parameter
       *
       * @param scatterz_exp scatter redshift exponent parameter
       *
       * @param purity the sample purity
       *
       * @param scaling_relation mass-observable relation pointer
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param proxy_pivot proxy pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param logM_base base of the mass logarithm
       *
       * @param purity_background the purity of the background 
       * source catalogue as a function of the lens redshift
       *
       * @param N_zs_index in which bin the cluster redshifts lie
       *
       * @return halo model in each radial bin
       *
       */
      std::vector<double> compute_model_density_scaling_relation (const std::vector<double> redshifts, const std::vector<double> mass_proxies, const std::vector<double> radius, const std::string probe, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_off, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_off, const cbl::glob::FuncGrid inv_Sigma_crit_interp, const std::shared_ptr<cosmology::Cosmology> cosmology, const double f_off0, const double f_offM, const double f_offz, const double OB_fact, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double purity, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base, const std::vector<double> purity_background, const std::vector<int> N_zs_index);
      
      /**
       * @brief Compute the excess density profile model
       * in all the radial bins, with mass derived from a 
       * scaling relation. The average profile is weighted over the
       * theoretical halo mass function
       *
       * @param MF_author author(s) of the mass function model
       *
       * @param halo_def overdensity definition
       *
       * @param overdensity overdensity
       *
       * @param z_min minimum redshift
       *
       * @param z_max maximum redshift
       *
       * @param mass_proxy_min minimum mass proxy
       *
       * @param mass_proxy_max maximum mass proxy
       *
       * @param z_error \f$\sigma_{z,0}\f$ in the equation 
       * \f$\sigma_z=\sigma_{z,0}(1+z)\f$, where \f$z\f$ is the true redshift and 
       * \f$\sigma_z\f$ represents the total statistical uncertainty on redshift.
       * In particular, \f$\sigma_z\f$ is the rms of the Gaussian PDF 
       * \f$P(z_{\rm ob}|z)\f$
       *
       * @param Plambda_mean_fc Function computing the mean of \f$P(\lambda_{\rm ob}|\lambda,z)\f$
       *
       * @param Plambda_A_mu \f$A_\mu\f$
       *
       * @param Plambda_B_mu \f$B_\mu\f$
       *
       * @param Plambda_C_mu \f$C_\mu\f$
       *
       * @param Plambda_std_fc Function computing the rms of \f$P(\lambda_{\rm ob}|\lambda,z)\f$
       *
       * @param Plambda_A_sigma \f$A_\sigma\f$
       *
       * @param Plambda_B_sigma \f$B_\sigma\f$
       *
       * @param Plambda_C_sigma \f$C_\sigma\f$
       *
       * @param interp_sigmaM interpolating function of sigma(M)
       *
       * @param interp_DlnsigmaM interpolating function of dln(sigma(M))/dM
       *
       * @param interp_DN interpolated amplitude of the growing mode
       *
       * @param radius the radius array
       *
       * @param probe "DeltaSigma" \f$\rightarrow\f$ excess surface density; 
       *  "g" \f$\rightarrow\f$ reduced shear
       *
       * @param DeltaSigma_interp_cen interpolated centred DeltaSigma
       *
       * @param DeltaSigma_interp_off interpolated centred DeltaSigma
       *
       * @param Sigma_interp_cen interpolated centred Sigma
       *
       * @param Sigma_interp_off interpolated centred Sigma
       *
       * @param inv_Sigma_crit_interp interpolated inverse critical density
       *
       * @param cosmology cosmological model
       *
       * @param f_off0 \f$ f_{\rm off,0} \f$
       *
       * @param f_offM \f$ f_{\rm off,M} \f$
       *
       * @param f_offz \f$ f_{\rm off,z} \f$
       *
       * @param OB_fact orientation boost factor, \f$\mathcal{B}_O\f$.
       * The average profile is multiplied by \f$(1+\mathcal{B}_O)\f$
       *
       * @param alpha normalisation of the scaling relation
       *
       * @param beta slope of the scaling relation
       *
       * @param gamma scaling relation redshift evolution
       * parameter
       *
       * @param scatter0 scatter normalisation
       *
       * @param scatterM scatter mass proxy dependence parameter
       *
       * @param scatterM_exp scatter mass proxy exponent parameter
       *
       * @param scatterz scatter redshift dependence parameter
       *
       * @param scatterz_exp scatter redshift exponent parameter
       *
       * @param purity the sample purity
       *
       * @param completeness_interp completeness interpolant
       *
       * @param scaling_relation mass-observable relation pointer
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param proxy_pivot proxy pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param logM_base base of the mass logarithm
       *
       * @param purity_background the purity of the background 
       * source catalogue interpolated as a function of the lens redshift
       *
       * @param purity_background_interp_range vector containing minum and
       * maximum of the cluster redshift range used for interpolating the 
       * background purity
       *
       * @return halo model in each radial bin
       *
       */
      std::vector<double> compute_model_density_scaling_relation (std::string MF_author, std::string halo_def, double overdensity, const double z_min, const double z_max, const double mass_proxy_min, const double mass_proxy_max, const double z_error, const std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc, const double Plambda_A_mu, const double Plambda_B_mu, const double Plambda_C_mu, const std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc, const double Plambda_A_sigma, const double Plambda_B_sigma, const double Plambda_C_sigma, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const std::vector<double> radius, const std::string probe, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> DeltaSigma_interp_off, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_cen, const std::vector<cbl::glob::FuncGrid2D> Sigma_interp_off, const cbl::glob::FuncGrid inv_Sigma_crit_interp, const std::shared_ptr<cosmology::Cosmology> cosmology, const double f_off0, const double f_offM, const double f_offz, const double OB_fact, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double purity, cbl::glob::FuncGrid2D completeness_interp, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base, const cbl::glob::FuncGrid purity_background, const std::vector<double> purity_background_interp_range);
      
      /**
       * @brief Compute the mass expected values from the following formula
       *
       * \f$\langle M(\lambda^*_{\text{ob},i},z_{\text{ob},j})\rangle = 
       * \frac{1}{N(\lambda^*_{\text{ob},i},z_{\text{ob},j})}
       * \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,\frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}
       * \int_{0}^{\infty} {\rm d} M_{\rm tr}\,\,M_{\rm tr}\,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}
       * \int_{0}^{\infty}{\rm d}\lambda^*_{\rm tr}\,\,P(\lambda^*_{\rm tr}| M_{\rm tr},z_{\rm tr})
       * \,\,P(z_{{\rm ob},j}|z_{\rm tr})
       * \,\,P(\lambda^*_{{\rm ob},i}|\lambda^*_{\rm tr})\,.\f$
       *
       * To enhance computational speed, the masses are derived from an
       * interpolation on a grid of observed redshift and mass proxy values.
       *
       * @param z_ob vector of cluster observed redshifts
       *
       * @param mass_proxy_ob vector of cluster observed mass proxies
       *
       * @param MF_author author(s) of the mass function model
       *
       * @param halo_def overdensity definition
       *
       * @param overdensity overdensity
       *
       * @param z_error \f$\sigma_{z,0}\f$ in the equation 
       * \f$\sigma_z=\sigma_{z,0}(1+z)\f$, where \f$z\f$ is the true redshift and 
       * \f$\sigma_z\f$ represents the total statistical uncertainty on redshift.
       * In particular, \f$\sigma_z\f$ is the rms of the Gaussian PDF 
       * \f$P(z_{\rm ob}|z)\f$
       *
       * @param Plambda_mean_fc Function computing the mean of \f$P(\lambda_{\rm ob}|\lambda,z)\f$
       *
       * @param Plambda_A_mu \f$A_\mu\f$
       *
       * @param Plambda_B_mu \f$B_\mu\f$
       *
       * @param Plambda_C_mu \f$C_\mu\f$
       *
       * @param Plambda_std_fc Function computing the rms of \f$P(\lambda_{\rm ob}|\lambda,z)\f$
       *
       * @param Plambda_A_sigma \f$A_\sigma\f$
       *
       * @param Plambda_B_sigma \f$B_\sigma\f$
       *
       * @param Plambda_C_sigma \f$C_\sigma\f$
       *
       * @param interp_sigmaM interpolating function of sigma(M)
       *
       * @param interp_DlnsigmaM interpolating function of dln(sigma(M))/dM
       *
       * @param interp_DN interpolated amplitude of the growing mode
       *
       * @param cosmology cosmological model
       *
       * @param alpha normalisation of the scaling relation
       *
       * @param beta slope of the scaling relation
       *
       * @param gamma scaling relation redshift evolution
       * parameter
       *
       * @param scatter0 scatter normalisation
       *
       * @param scatterM scatter mass proxy dependence parameter
       *
       * @param scatterM_exp scatter mass proxy exponent parameter
       *
       * @param scatterz scatter redshift dependence parameter
       *
       * @param scatterz_exp scatter redshift exponent parameter
       *
       * @param completeness_interp completeness interpolant
       *
       * @param scaling_relation mass-observable relation pointer
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param proxy_pivot proxy pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param logM_base base of the mass logarithm
       *
       * @return the mass expected values divided by mass_pivot
       *
       */
      std::vector<double> compute_mass_expected_value (const std::vector<double> z_ob, const std::vector<double> mass_proxy_ob, std::string MF_author, std::string halo_def, double overdensity, const double z_error, const std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc, const double Plambda_A_mu, const double Plambda_B_mu, const double Plambda_C_mu, const std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc, const double Plambda_A_sigma, const double Plambda_B_sigma, const double Plambda_C_sigma, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const std::shared_ptr<cosmology::Cosmology> cosmology, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, cbl::glob::FuncGrid2D completeness_interp, const std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double logM_base);
      
    }
  }
}

#endif
