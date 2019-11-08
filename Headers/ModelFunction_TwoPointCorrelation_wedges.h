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
 *  @file Headers/ModelFunction_TwoPointCorrelation_wedges.h
 *
 *  @brief Functions to model the wedges of the two-point correlation
 *  function
 *
 *  This file contains all the prototypes of the functions used to
 *  model the wedges of the two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCTWOPWED__
#define __MODFUNCTWOPWED__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {

      /**
       *  @brief the model wedges of the two-point correlation
       *  function
       *
       *  this function computes the wedges of the two-point
       *  correlation function either with the de-wiggled model or
       *  with the mode-coupling model (Kazin et al. 2012):
       *
       *  \f[ \xi(\Delta \mu, s) = \xi_0(s)+\frac{1}{2} \left( \frac{
       *  \mu_{max}^3- \mu_{min}^3}{\mu_{max}-\mu_{min}} -1
       *  \right)\xi_2(s)+ \frac{1}{8} \left( \frac{ 7
       *  \left(\mu_{max}^5-\mu_{min}^5\right)-10 \left(\mu_{max}^3-
       *  \mu_{min}^3 \right)}{\mu_{max}-\mu_{min}}+3 \right)
       *  \xi_4(s).  \f]
       *
       *  where \f$\xi_0(s), \xi_2(s), \xi_4(s)\f$ are the two-point
       *  correlation function multipoles computed by
       *  cbl::modelling::twopt::Xi_l
       *
       *  @param rr vector of scales to compute wedges 
       * 
       *  @param nWedges the number of wedges
       *
       *  @param mu_integral_limits the \f$\mu\f$ integral limits
       *  used to measure the wedges
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the wedges of the two-point correlation function.
       */
      std::vector<std::vector<double>> xi_Wedges (const std::vector<double> rr, const int nWedges, const std::vector<std::vector<double>> mu_integral_limits, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp=1, const double alpha_par=1.);

      /**
       *  @brief the wedge of the two-point correlation function
       *
       *  The function computes the wedges of the two-point
       *  correlation function (Kazin et al. 2012):
       *
       *  \f[ \xi(\Delta \mu, s) = \xi_0(s)+\frac{1}{2} \left( \frac{
       *  \mu_{max}^3- \mu_{min}^3}{\mu_{max}-\mu_{min}} -1
       *  \right)\xi_2(s)+ \frac{1}{8} \left( \frac{ 7
       *  \left(\mu_{max}^5-\mu_{min}^5\right)-10 \left(\mu_{max}^3-
       *  \mu_{min}^3 \right)}{\mu_{max}-\mu_{min}}+3 \right)
       *  \xi_4(s).  \f]
       *
       *  where \f$\xi_0(s), \xi_2(s), \xi_4(s)\f$ are the two-point
       *  correlation function multipoles
       *
       *  @param rr vector of scales to compute wedges 
       *
       *  @param dataset_order vector that specify the wedges
       *  to be computed for each scale 
       *	
       *  @param mu_integral_limits the \f$\mu\f$ integral limits
       *  used to measure the wedges
       *
       *  @param model the \f$P(k,\mu)\f$ model; the possible options
       *  are: dispersion_dewiggled, dispersion_modecoupling
       *
       *  @param parameter vector containing parameter values
       *
       *  @param pk_interp vector containing power spectrum
       *  interpolating functions
       *
       *  @param prec the integral precision
       *
       *  @param alpha_perp the shift transverse to the l.o.s.
       *
       *  @param alpha_par the shift parallel to the l.o.s.
       *
       *  @return the wedges of the two-point correlation function.
       */
      std::vector<double> xi_Wedges (const std::vector<double> rr, const std::vector<int> dataset_order, const std::vector<std::vector<double>> mu_integral_limits, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec=1.e-5, const double alpha_perp=1., const double alpha_par=1.);
      
      /**
       *  @brief the model wedges of the two-point correlation
       *  function
       *
       *  this functions computes the wedges of the two-point either
       *  with the de-wiggled model or with the mode-coupling model;
       *  specifically, the wedges are computed by
       *  cbl::modelling::twopt::xi_Wedges
       *
       *  the wedges of the two-point correlation functions are defined as follows:
       *
       *  \f[ \xi_w(s) = \frac{1}{\mu_{max}-\mu_{min}}
       *  \int_{\mu_{min}}^{\mu_{max}} \mathrm{\mu} \xi(s', \mu'); \f]
       *
       *  where \f$ \xi(s', \mu')\f$ is the polar two-point
       *  correlation function computed at shifted positions \f$s' =
       *  s\sqrt{\mu^2\alpha^2_{\parallel}+(1-\mu^2)\alpha^2_{\perp}}\f$,
       *  \f$\mu' =
       *  \mu\alpha_{\parallel}/\sqrt{\mu^2\alpha^2_{\parallel}+(1-\mu^2)\alpha^2_{\perp}}\f$
       *
       * @param rad the scale at which the model is computed
       *
       * @param inputs pointer to the structure that contains the
       * de-wiggled power spectrum, the number of wedges and
       * \f$\sigma_8(z)\f$, computed at a given (fixed) cosmology
       *
       * @param parameter vector containing the input parameters
       *
       * @return the wedges of the two-point correlation function
       */
      std::vector<double> xiWedges (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       * @brief return the wedges of the two-point correlation
       * function, intended for anisotropic BAO measurements
       *
       * the wedges of the two-point correlation function are computed
       * as follows (Ross et al. 2017):
       *
       *  \f[ \xi_{\perp}(s) = B_{\perp}\xi_{\perp}(s,
       *  \alpha_{\perp},
       *  \alpha_{\parallel})+A_{\perp}^0+\frac{A_{\perp}^1}{s}+\frac{A_{\perp}^2}{s^2};
       *  \\ \xi_{\parallel}(s) = B_{\parallel}\xi_{\parallel}(s,
       *  \alpha_{\parallel},
       *  \alpha_{\parallel})+A_{\parallel}^0+\frac{A_{\parallel}^1}{s}+\frac{A_{\parallel}^2}{s^2};
       *  \\ \f]
       *
       *  where \f$\xi_{\perp}\f$, \f$\xi_{\parallel}\f$ are the two
       *  wedges of the two-point correlation function.
       *
       *  The function takes as inputs ten parameters
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
       * @param inputs pointer to the structure that contains the
       * de-wiggled power spectrum, the number of wedges and
       * \f$\sigma_8(z)\f$, computed at a given (fixed) cosmology
       *
       * @param parameter 10D vector containing the input parameters
       *
       * @return the wedges of the two-point correlation function
       */
      std::vector<double> xiWedges_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

    }
  }
}

#endif
