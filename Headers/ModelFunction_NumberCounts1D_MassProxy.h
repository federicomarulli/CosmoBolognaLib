/********************************************************************
 *  Copyright (C) 2021 by Giorgio Lesci and Federico Marulli        *
 *  giorgio.lesci2@unibo.it, federico.marulli3@unibo.it             *
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
 *  @file Headers/ModelFunction_NumberCounts1D_MassProxy.h
 *
 *  @brief Global functions to model number counts as a function of a mass proxy
 *
 *  This file contains all the prototypes of the global functions used
 *  to model number counts as a funciton of mass proxy
 *  
 *  @author Giorgio Lesci, Federico Marulli
 *
 *  @author giorgio.lesci2@unibo.it, federico.marulli3@unibo.it
 */

#ifndef __MODFUNCNCMP__
#define __MODFUNCNCMP__

#include "ModelFunction_NumberCounts.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace numbercounts {
    
      /**
       * @brief compute the number counts as a function
       * of the mass proxy, with the following model:
       *
       * \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle 
       *  = w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\,\,\Omega 
       *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
       *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
       *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
       *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
       *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
       *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
       *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
       *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
       *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}). \f$
       *
       *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is expressed as:
       *
       *  \f$P(\lambda_{\rm tr}|M_{\rm tr},z_{\rm tr})= 
       *  P(M_{\rm tr}|\lambda_{\rm tr},z_{\rm tr})\,
       *  P(\lambda_{\rm tr}|z_{\rm tr})\,/\,P( M_{\rm tr}|z_{\rm tr}),\f$
       *
       * @param proxy mass proxy bin centers
       *
       * @param inputs inputs to compute the predicted counts
       *
       * @param parameter vector containing cosmological parameters
       *
       * @return the number counts as a function of mass proxy
       */
      std::vector<double> number_counts_proxy (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief compute the number counts as a function
       * of the mass proxy, with the following model:
       *
       * \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle 
       *  = w(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\,\,\Omega 
       *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
       *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
       *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
       *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
       *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
       *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
       *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
       *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
       *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}). \f$
       *
       *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is 
       *  a log-normal whose mean is given by the mass-mass proxy
       *  relation, i.e. 
       *
       *  \f$\log (\lambda/\lambda_{\rm piv}) = \alpha + \beta 
       *  \log (M/M_{\rm piv}) + \gamma \log (f(z)),\f$
       *
       *  and whose standard deviation is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as:
       *
       *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
       *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
       *
       * @param proxy mass proxy bin centers
       *
       * @param inputs inputs to compute the predicted counts
       *
       * @param parameter vector containing cosmological parameters
       *
       * @return the number counts as a function of mass proxy
       */
      std::vector<double> number_counts_proxy_classic (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief make the counts model a simple model describing the observed counts, 
       * i.e. this function returns 1
       *
       * @param Mass halo mass
       *
       * @param Sigma mass variance at z=0
       *
       * @param redshift redshift
       *
       * @param model_bias bias model
       *
       * @param Delta overdensity
       *
       * @param method_SS method used to compute the power spectrum
       *
       * @param cosmo cosmological model
       *
       * @return return 1
       *
       */
      double no_response (const double Mass, const double Sigma, const double redshift, const std::string model_bias, const double Delta, const std::string method_SS, std::shared_ptr<void> cosmo);
      
      /**
       * @brief introduce the halo bias in the model, making the counts model
       * a response function for computing the super-sample covariance
       *
       * @param Mass halo mass
       *
       * @param Sigma mass variance at z=0
       *
       * @param redshift redshift
       *
       * @param model_bias bias model
       *
       * @param Delta overdensity
       *
       * @param method_SS method used to compute the power spectrum
       *
       * @param cosmo cosmological model
       *
       * @return the response function factor, i.e. the halo bias
       *
       */
      double bias_response (const double Mass, const double Sigma, const double redshift, const std::string model_bias, const double Delta, const std::string method_SS, std::shared_ptr<void> cosmo);
      
    }
  }
}

#endif
