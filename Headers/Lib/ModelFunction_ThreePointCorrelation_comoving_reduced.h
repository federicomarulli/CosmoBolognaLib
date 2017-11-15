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
 *  Headers/Lib/ModelFunction_ThreePointCorrelation.h
 *
 *  @brief Functions to model the reduced three-point correlation
 *  function in comoving coordinates
 *
 *  This file contains the implementation of the functions used to
 *  model the reduced three-point correlation function in comoving
 *  coordinates
 *  
 *  @author Federico Marulli, Michele Moresco
 *
 *  @author federico.marulli3@unbo.it, michele.moresco@unibo.it
 */

#ifndef __MODFUNCTHREEPCOMRED__
#define __MODFUNCTHREEPCOMRED__

#include "ModelFunction_ThreePointCorrelation.h"


// ============================================================================


namespace cosmobl {

  namespace modelling {

    namespace threept {

      /**
       *  @brief model for the reduced three-point correlation
       *  function
       *
       *  The function computes:
       *  \f[ Q_{NL, lb} = \frac{1}{b_1}\left( Q_{DM}+\frac{b_2}{b_1} \right) \f]
       *
       *  the model has 6 parameters: 
       *     - \f$b_1\f$ the linear bias
       *     - \f$b_2\f$ the non-linear bias
       *
       *  the dark matter reduced three-point correlation function is fixed and
       *  provided in input
       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  dark matter two-point correlation function
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the reduced three-point correlation function
       */
      vector<double> Q_nonlinear_localbias (const vector<double> theta, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the reduced three-point correlation
       *  function with non-local contributions
       *
       *  The function computes:
       *  \f[ Q_{NL, lb} = \frac{1}{b_1}\left( Q_{DM}+\frac{b_2}{b_1} + g_2 Q_{nl}\right) \f]
       *
       *  the model has 6 parameters: 
       *     - \f$b_1\f$ the linear bias
       *     - \f$b_2\f$ the non-linear bias
       *     - \f$g_2\f$ the non-local parameter
       *
       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  galaxy matter three-point correlation function
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the reduced three-point correlation function
       */
      vector<double> Q_nonlinear_nonlocalbias (const vector<double> theta, const shared_ptr<void> inputs, vector<double> &parameter);

      /**
       *  @brief model for the reduced three-point correlation
       *  function with non-local contributions with alpha shift
       *  on Pk
       *
       *  The function computes:
       *  \f[ Q_{NL, lb} = \frac{1}{b_1}\left( Q_{DM}+\frac{b_2}{b_1} + g_2 Q_{nl}\right) \f]
       *
       *  the model has 6 parameters: 
       *     - \f$b_1\f$ the linear bias
       *     - \f$b_2\f$ the non-linear bias
       *     - \f$g_2\f$ the non-local parameter
       *     - \f$\alpha\f$ the shift
       *
       *
       *  @param theta the scale at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the
       *  galaxy matter three-point correlation function
       *
       *  @param parameter 6D vector containing the input parameters
       *
       *  @return the reduced three-point correlation function
       */
      vector<double> Q_nonlinear_nonlocalbias_alpha (const vector<double> theta, const shared_ptr<void> inputs, vector<double> &parameter);

    }
  }
}

#endif
