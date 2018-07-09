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
       *  @brief the wedges of the two-point correlation function
       *
       *  this functions computes the wedges of the two-point
       *  correlation function
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
       *  The function takes as inputs four fundamental parameters -
       *  \f$\alpha_{\perp}\f$ - \f$\alpha_{\parallel}\f$ -
       *  \f$f(z)\sigma_8(z)\f$ - \f$b(z)\sigma_8(z)\f$ -
       *  \f$\Sigma_S\f$
       *
       * @param rad the scale at which the model is computed
       *
       * @param inputs pointer to the structure that contains the
       * de-wiggled power spectrum, the number of multipoles and
       * \f$\sigma_8(z)\f$, computed at a given (fixed) cosmology
       *
       * @param parameter 4D vector containing the input parameters
       *
       * @return the wedges of the two-point correlation function
       */
      std::vector<double> xiWedges (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       * @brief return the wedges of the two-point 
       * correlation function, intended for anisotropic BAO
       * measurements (Ross et al. 2017)
       *
       * The functions computes the wedges of the 
       * two-point correlation function
       *
       * where \f$j_l(ks)\f$ are the bessel functions.
       *
       * The function takes as inputs ten parameters
       *    - \f$\alpha_{\perp}\f$
       *    - \f$\alpha_{\parallel}\f$
       *    - \f$B_0\f$
       *    - \f$B_2\f$
       *    - \f$A^{\perp}_0\f$
       *    - \f$A^{\perp}_1\f$
       *    - \f$A^{\perp}_2\f$
       *    - \f$A^{\parallel}_0\f$
       *    - \f$A^{\parallel}_1\f$
       *    - \f$A^{\parallel}_2\f$
       *
       *
       * @param rad the scale at which the model is computed
       * @param inputs pointer to the structure that contains the
       * de-wiggled power spectrum, the number of multipoles and
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
