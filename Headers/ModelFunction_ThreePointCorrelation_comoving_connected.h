/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  Headers/ModelFunction_ThreePointCorrelatio_comoving_connected.h
 *
 *  @brief Functions to model the connected three-point correlation
 *  function in comoving coordinates
 *
 *  This file contains the implementation of the functions used to
 *  model the connected three-point correlation function in comoving
 *  coordinates
 *  
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __MODFUNCTHREEPCOMCON__
#define __MODFUNCTHREEPCOMCON__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace threept {

      /**
       *  @brief model for the connected three-point
       *  correlation function
       *
       *  Model for the connected three-point correlation function
       *  described in Slepian & Eisenstein (2017) and implemented in
       *  cbl::cosmology::Cosmology::zeta_RSD
       *
       *  the model has 4 parameters: 
       *    - \f$ b_1\f$
       *    - \f$ b_2\f$
       *    - \f$ b_t\f$
       *    - \f$ \beta = f/b_1\f$
       *
       *  the dark matter two-point correlation function is fixed and
       *  provided in input
       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  dark matter power spectrum
       *
       *  @param parameter 4D vector containing the input parameters
       *
       *  @return the redshift space monopole of the connected
       *  three-point correlation function
       */
      std::vector<double> zeta_RSD (const std::vector<double> theta, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

    }
  }
}

#endif
