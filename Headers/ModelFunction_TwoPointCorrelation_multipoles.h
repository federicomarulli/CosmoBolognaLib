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
 *  @file Headers/ModelFunction_TwoPointCorrelation_multipoles.h
 *
 *  @brief Functions to model the multipoles of the two-point
 *  correlation function
 *
 *  This file contains all the prototypes of the functions used to
 *  model the multipoles of the two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCTWOPMULT__
#define __MODFUNCTWOPMULT__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {

      /**
       * @brief return multipoles of the two-point correlation
       * function
       *
       * The functions computes the multipoles of the two-point
       * correlation function
       *
       * \f[ \xi_l(s) = i^l \int \frac{\mathrm{d} k}{2\pi^2} k^2
       * P_l(k) j_l(ks); \f]
       *
       * where \f$j_l(ks)\f$ are the Bessel functions.
       *
       * The function takes as inputs four fundamental parameters
       *    - \f$\alpha_{\perp}\f$
       *    - \f$\alpha_{\parallel}\f$
       *    - \f$f(z)\sigma_8(z)\f$
       *    - \f$b(z)\sigma_8(z)\f$
       *    - \f$ \Sigma_S \f$
       *
       * @param rad the scale at which the model is computed
       *
       * @param inputs pointer to the structure that contains the data
       * required to construct the model
       *
       * @param parameter 4D vector containing the input parameters
       *
       * @return the multipoles of the two-point correlation function
       */
      std::vector<double> xiMultipoles (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       * @brief return multipoles of the two-point 
       * correlation function, intended for anisotropic BAO
       * measurements (Ross et al. 2017). In only works with
       * monopole and quadrupole.
       *
       * The functions computes the multipoles of the 
       * two-point correlation function:
       *
       * \f[ \xi_0(s) = B_0\xi_0(s, \alpha_{\perp},
       * \alpha_{\parallel})+A_0^0+\frac{A_0^1}{s}+\frac{A_0^2}{s^2}
       * \f]
       *	 
       * \f[ \xi_2(s) = \frac{5}{2}\left[B_2\xi_{\mu2}(s,
       * \alpha_{\perp}, \alpha_{\parallel})-B_0\xi_0(s,
       * \alpha_{\perp}, \alpha_{\parallel})\right]
       * +A_2^0+\frac{A_2^1}{s}+\frac{A_2^2}{s^2} \f]
       *
       * where \f$\xi_0(s, \alpha_{\perp}, \alpha_{\parallel})\f$ is
       * the monopole computed at the fiducial cosmology,
       * \f$\xi_{\mu2}(s, \alpha_{\perp}, \alpha_{\parallel}) =
       * 3\int_0^1\mathrm{d}\mu\mu^2\xi(s, \mu, \alpha_{\perp},
       * \alpha_{\parallel})\f$.
       *
       * The function takes as inputs ten parameters
       *    - \f$\alpha_{\perp}\f$
       *    - \f$\alpha_{\parallel}\f$
       *    - \f$B_0\f$
       *    - \f$B_2\f$
       *    - \f$A^0_0\f$
       *    - \f$A^0_1\f$
       *    - \f$A^0_2\f$
       *    - \f$A^2_0\f$
       *    - \f$A^2_1\f$
       *    - \f$A^2_2\f$
       *
       * @param rad the scale at which the model is computed
       *
       * @param inputs pointer to the structure that contains the data
       * required to construct the model
       *
       * @param parameter 10D vector containing the input parameters
       *
       * @return the multipoles of the two-point correlation function
       */
      std::vector<double> xiMultipoles_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief return multipoles of the two-point correlation
       *  function
       *
       *  The functions computes the multipoles of the two-point
       *  correlation function
       *
       *  \f[ \xi_l(s) = i^l \int \frac{\mathrm{d} k}{2\pi^2} k^2
       *  P_l(k) j_l(ks); \f]
       *
       *  where \f$j_l(ks)\f$ are the Bessel functions.
       *
       *  The function takes as inputs 2 parameters
       *    - \f$\sigma_8(z)\f$
       *    - \f$b(z)\f$
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  data required to construct the model
       *
       *  @param parameter vector containing the input parameters
       *
       *  @return the multipoles of the two-point correlation function
       */
      std::vector<double> xiMultipoles_sigma8_bias (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief Model for the two-point correlation function, with damping at the BAO peak.
       *
       *  The function computes the multipoles of the two-point
       *  correlation function:
       *
       *  \f[ \xi_l(s) = i^l \int \frac{\mathrm{d} k}{2\pi^2} k^2
       *  P_l(k) j_l(ks), \f]
       *
       *  where \f$j_l(ks)\f$ are the Bessel functions
       *  and \f$P_l(k)\f$ is given by:
       *
       *  \f[ P_l(k)=\frac{2l+1}{2} \int^{+1}_{-1}
       *  \mathrm{d} \mu P\left(k, \mu\right) L_l\left(\mu\right) , \f]
       *
       *  where \f$ L_l\left(\mu\right) \f$ are the Legendre polynomials.
       *  Redshift space distorsion are modelled in the Kaiser limit, 
       *  if the user chooses redshift space in set_data_model_bias.
       *  Otherwise, they are neglected. Photometric errors are modelled
       *  with a Gaussian damping factor:
       *
       *   \f[ P(k, \mu)=P(k)\, b_\mathrm{eff}^2
       *   \left(1+\frac{f}{b_\mathrm{eff}} \mu^2\right)^2 \mathrm{exp}(-k^2\mu^2\sigma^2) , \f]
       *
       *  where \f$\sigma=\frac{c\sigma_z}{H(z)}\f$.
       *  We account for the nonlinear damping at the BAO peak:
       *
       *  \f[ P(k)=[P_\mathrm{lin}(k)-P_\mathrm{nw}(k)]e^{-k^2\Sigma^2_{\mathrm{NL}}/2}+P_\mathrm{nw}(k), \f]
       *
       *  where \f$P_\mathrm{lin}\f$ is the linear power spectrum
       *  and \f$P_\mathrm{nw}\f$ is the de-wiggled power spectrum,
       *  computed with the parametrisation of Eiseinstein_Hu.
       *
       *  The model has N cosmological parameters, plus: 
       *    - \f$\Sigma_{NL}\f$ damping at BAO scale
       *
       *  The linear effective bias is considered a derived parameter,
       *  from cosmology, and computed from the input cluster masses
       *  in two ways: 
       *  using the provided halo masses with cbl::cosmology::bias_eff_mass 
       *  or with cbl::cosmology::bias_eff_selection_function 
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the monopole of the two-point correlation function
       */
      std::vector<double> xiMultipoles_linear_theoretical_bias_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

       /**
       *  @brief the damped two-point correlation;
       *
       *  @param rad the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  cosmological paramters used to compute the dark matter
       *  two-point correlation function
       *
       *  @param parameter 1D vector containing the linear bias
       *
       *  @return the damped two-point correlation.
       */
      std::vector<double> xil_damped_scaling_relation_sigmaz_cosmology (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    }
  }
}

#endif
