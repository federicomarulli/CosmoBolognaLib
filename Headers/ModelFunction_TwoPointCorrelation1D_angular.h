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
 *  @file Headers/ModelFunction_TwoPointCorrelation1D_angular.h
 *
 *  @brief Functions to model the angular two-point correlation function
 *
 *  This file contains all the prototypes of the functions used to
 *  model the angular two-point correlation function
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, massimilia.romanell2@unibo.it
 */

#ifndef __MODFUNCTWOP1DANG__
#define __MODFUNCTWOP1DANG__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {
      
      using ProxyFunction = std::function<double(double, double, std::shared_ptr<void>)>;
      
      /**
       *  @brief The angular correlation function. The function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser
       *  limit, that is neglecting non-linearities in dynamics
       *  and bias. Specifically:
       * 
       *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
       *  f \cdot b + \frac{1}{5}f^2
       *  \right] \cdot \xi_{\rm DM}(s) \f]
       *
       *  the parameters of the true dN/dz distribution ("gaussian",
       *  "skewnormal" or "polynomial") in a given photometric
       *  redshift bin are fixed by the user. The angular correlation
       *  function is computed using the input cosmological parameters
       *
       *  @param theta the angular vector at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @param parameter 1D vector containing the parameter
       *  
       *  @return the angular correlation function 
       */
      std::vector<double> wtheta (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief The angular correlation function. The function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and
       *  bias. Specifically:
       * 
       *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
       *  f \cdot b + \frac{1}{5}f^2
       *  \right] \cdot \xi_{\rm DM}(s) \f]
       *
       *  the dN/dz distribution is computed from the theoretical mass
       *  function:
       *
       *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
       *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
       *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
       *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
       *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
       *  , \f]
       *
       *  @param theta the angular scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological parameters used to compute the
       *  two-point correlation function
       *
       *  @param parameter 1D vector
       *
       *  @return the angular correlation function
       */
      std::vector<double> wtheta_MF_theory (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       *  @brief The angular correlation function. The function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and bias,
       *  only if the user choose redshift space in
       *  set_data_model_bias. Specifically:
       * 
       *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
       *  f \cdot b + \frac{1}{5}f^2
       *  \right] \cdot \xi_{\rm DM}(s) \f]
       *
       *  The model has N cosmological parameters. 
       *  The linear effective bias is computed theoretically for each cosmology in two ways: 
       *  using the provided halo masses with cbl::cosmology::Bias::bias_eff_mass, 
       *  or with cbl::cosmology::Bias::bias_eff_selection_function. The dN/dz distribution is 
       *  computed from the theoretical mass function:
       *
       *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
       *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
       *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
       *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
       *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
       *  , \f]
       *
       *  @param theta the angular scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector
       *
       *  @return the angular correlation function
       */
      std::vector<double> wtheta_bias_MF_theory (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief The angular correlation function, including BAO. The function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and bias,
       *  only if the user choose redshift space in
       *  set_data_model_bias. Specifically:
       * 
       *  \f[\xi_0(s) = b^2 \left[ 1 + \frac{2 \beta}{3}
       *  \cdot + \frac{\beta^2}{5} \right] \cdot
       *  \xi_{\rm DM}(s, \Sigma_{NL})\f]
       *
       *  The model has N cosmological parameters, plus: 
       *    - \f$\Sigma_{NL}\f$ damping at BAO
       *  The bias is considered a derived parameter, from cosmology
       *
       *  where \f$ b \f$ is the linear effective bias, computed theoretically for each cosmology in two ways: 
       *  using the provided halo masses with cbl::cosmology::Bias::bias_eff_mass, 
       *  or with cbl::cosmology::Bias::bias_eff_selection_function. The dN/dz distribution is 
       *  computed from the theoretical mass function:
       *
       *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
       *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
       *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
       *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
       *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
       *  , \f]
       *
       *  @param theta the angular scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector
       *
       *  @return the angular correlation function
       */
      std::vector<double> wtheta_bias_MF_theory_BAO (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief The angular correlation function, including BAO. The function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and bias,
       *  only if the user choose redshift space in
       *  set_data_model_bias. Specifically:
       * 
       *  \f[\xi_0(s) = b^2 \left[ 1 + \frac{2 \beta}{3}
       *  \cdot + \frac{\beta^2}{5} \right] \cdot
       *  \xi_{\rm DM}(s, \Sigma_{NL})\ + A_0 + A_1/s+A_2/s^2\f]
       *
       *  The model has N cosmological parameters, plus: 
       *    - \f$A_0, A_1, A2\f$ parameters of the polynomial
       *    - \f$\Sigma_{NL}\f$ damping at BAO
       *  The bias is considered a derived parameter, from cosmology
       *
       *  where \f$ b \f$ is the linear effective bias, computed theoretically for each cosmology in two ways: 
       *  using the provided halo masses with cbl::cosmology::Bias::bias_eff_mass, 
       *  or with cbl::cosmology::Bias::bias_eff_selection_function. The dN/dz distribution is 
       *  computed from the theoretical mass function:
       *
       *  \f[ \frac{{\rm d}N }{{\rm d} z_{\rm tr}}(\Delta z_{\text{ob}})=\,\,\Omega 
       *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
       *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
       *  \int_{\Delta z_{\text{ob}}}{\rm d} z_{\rm ob} 
       *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
       *  , \f]
       *
       *  @param theta the angular scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector
       *
       *  @return the angular correlation function
       */
      std::vector<double> wtheta_bias_MF_theory_BAO_poly (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief The angular correlation function. The function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and
       *  bias. Specifically:
       * 
       *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
       *  f \cdot b + \frac{1}{5}f^2
       *  \right] \cdot \xi_{\rm DM}(s) \f]
       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the damped two-point correlation monopole.
       */
      std::vector<double> wtheta_scaling_relation_cosmology (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief the angular correlation function, with mass/intrinsic richness scaling relation, using a theoretical mass function to compute dN/dz
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser limit,
       *  that is neglecting non-linearities in dynamics and
       *  bias. Specifically:
       * 
       *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
       *  f \cdot b + \frac{1}{5}f^2
       *  \right] \cdot \xi_{\rm DM}(s) \f]       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the angular power spectrum computed with mass/intrinsic richness scaling relation.
       */
      std::vector<double> wtheta_scaling_relation_cosmology_MF_theory (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief limits of the true redshift interval, given the
       * photometric bin width and the photometric error
       *
       * @param zmin zmin
       *
       * @param zmax zmax
       *
       * @param z_abs_err the absolute redshift error, alone or
       * multiplied for (1+z_mean), with z_mean is the average of the
       * photometric redshift bin
       *
       * @return limits of the true redshift interval
       *
       */
      std::vector<double> z_limits (const double zmin, const double zmax, const double z_abs_err);

      /**
       *  @brief compute the redshift selection function at given z 
       *
       *  @param redshift redshift
       *
       *  @param parameter vector containing the fit parameter of the
       *  corresponding redshift bin
       *
       *  @param distribution_type the shape of the dN/dz
       *  distribution.  Can be "polynomial", for a fit of the
       *  photometric redshift distribution, "gaussian" or
       *  "skewnormal", in the correct case of convolution with P(z
       *  true|z_phot)
       *
       *  @param redshift_min minimum redshift of photometric bin (for
       *  the normalization of polynomial distribution)
       *
       *  @param redshift_max maximum redshift of photometric bin (for
       *  the normalization of polynomial distribution)
       *
       *  @return the redshift selection function at given z
       */
      double redshift_selection_function (double redshift, std::vector<double> parameter, std::string distribution_type, double redshift_min, double redshift_max);

      /**
       * @brief compute the redshift distribution as function of mass
       * proxy
       *
       * @param redshift redshift
       *
       * @param alpha alpha
       *
       * @param beta beta
       *
       * @param gamma gamma
       *
       * @param scatter0 scatter0
       *
       * @param scatterM scatterM
       *
       * @param scatterM_exp scatterM_exp
       *
       * @param scatterz scatterz
       *
       * @param scatterz_exp scatterz_exp
       *
       * @param z_abs_err abs. error on z
       *
       * @param proxy_rel_err rel. error on proxy
       *
       * @param Plambda_a Plambda_a
       *
       * @param Plambda_b Plambda_b
       *
       * @param Plambda_c Plambda_c
       *
       * @param fz redshift evolution function in the scaling relation
       *
       * @param redshift_min minimum redshift
       *
       * @param redshift_max maximum redshift
       *
       * @param proxy_min minimum mass proxy
       *
       * @param proxy_max maximum mass proxy
       *
       * @param cosmology the cosmology 
       *
       * @param Area the area in degrees
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       * Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
       * (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
       * Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
       * MF by Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       * @param Delta \f$\Delta\f$, the overdensity
       *
       * @param isDelta_critical \f$\rightarrow\f$ \f$\Delta\f$ is the
       * overdensity defined with respect to the critical density
       *
       * @param interp_sigmaM interpolating function of \f$
       * \sigma(M)\f$
       *
       * @param interp_DlnsigmaM interpolating function of \f$
       * \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       *
       * @param interp_DN interpolated amplitude of the growing mode
       *
       * @param proxy_pivot mass proxy pivot in the scaling relation
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param log_base logarithmic base used in the scaling relation
       *
       * @param weights_interp weight derived from the selection
       * function, interpolated to include different richness bins
       *
       * @return values of redshift distribution as a function of
       * mass proxy
       */
      double redshift_distribution_proxy (const double redshift, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp, const double z_abs_err, const double proxy_rel_err, const double Plambda_a, const double Plambda_b, const double Plambda_c, ProxyFunction fz, const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, std::shared_ptr<cosmology::Cosmology> cosmology, const double Area, const std::string model_MF, const double Delta, const bool isDelta_critical, const cbl::glob::FuncGrid interp_sigmaM, const cbl::glob::FuncGrid interp_DlnsigmaM, const cbl::glob::FuncGrid interp_DN, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base, cbl::glob::FuncGrid weights_interp);

      /**
       *  @brief the integrand function for
       *  the calculus of angular correlation function 
       *
       *  the function computes:
       *
       *  \f[w(\theta) = \int\int dz_1 dz_2 \phi(z_1) \phi(z_2) \xi(s), \f]
       *
       *  where,
       *
       *  \f$ s=\sqrt{r^2(z_1)+r^2(z_2)-2r_1r_2\cos(\theta)} \f$
       *  and: 
       *  \f[ \phi(z)=\frac{\frac{dN(z)}{dz}}{\int dz \frac{dN(z)}{dz}}\f]
       *  
       *  Redshift-space distorsions are modelled in the Kaiser
       *  limit, that is neglecting non-linearities in dynamics
       *  and bias. Specifically:
       * 
       *  \f[\xi(s) = \left[ b^2 + \frac{2}{3}
       *  f \cdot b + \frac{1}{5}f^2
       *  \right] \cdot \xi_{\rm DM}(s) \f]
       *
       *  where \f$ b \f$ is the bias.
       *
       *  @param theta the multipole at which angular power spectrum is computed
       *
       *  @param coordUnits the coordinate units
       *
       *  @param xi0_interp the correlation function power spectrum at z=0
       *
       *  @param comoving_distance_interp the interpolated comoving distance
       *
       *  @param interp_DN the interpolated normalised amplitude of the growing mode
       *
       *  @param redshift_counts_interp the interpolated dN/dz distribution
       *
       *  @param bias the bias
       *
       *  @param redshift_counts_interp_bin2 the interpolated dN/dz distribution of the second redshift bin
       *
       *  @param bias_bin2 the bias of the second redshift bin
       *
       *  @param GSL true->use GSL integration, false->use CUBA
       *
       *  @param parameter the list of parameters, containing the
       *  cosmological parameters
       *
       *  @param inputs pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @param xi2_interp the interpolated correlation function
       *  quadrupole at z=0
       *
       *  @param xi4_interp the interpolated correlation function
       *  hexadecapole at z=0
       *
       *  @return the integral of the angular correlation function  
       */
      double integral_wtheta (const double theta, const CoordinateUnits coordUnits, cbl::glob::FuncGrid xi0_interp, cbl::glob::FuncGrid comoving_distance_interp, cbl::glob::FuncGrid interp_DN, cbl::glob::FuncGrid redshift_counts_interp, const double bias, cbl::glob::FuncGrid redshift_counts_interp_bin2, const double bias_bin2, bool GSL, std::vector<double> parameter, std::shared_ptr<void> inputs, cbl::glob::FuncGrid xi2_interp=cbl::glob::FuncGrid(), cbl::glob::FuncGrid xi4_interp=cbl::glob::FuncGrid());

      /**
       *  @brief analytic covariance matrix for the angular correlation function
       *
       *  The function computes the analytic covariance matrix for angular correlation function, 
       *  using a theoretical mass function to compute dN/dz and theoretical bias
       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the angular power spectrum
       *
       *  @param parameter 1D vector containing the model parameters 
       *
       *  @return the covariance matrix for the angular correlation function
       */
      std::vector<double> covariance_wtheta_bias_MF_theory (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    }
  }
}

#endif
