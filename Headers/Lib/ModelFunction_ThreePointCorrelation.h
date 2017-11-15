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
 *  @file Headers/Lib/ModelFunction_ThreePointCorrelation.h
 *
 *  @brief Functions to model the three-point correlation function
 *
 *  This file contains the prototypes of the functions used to model
 *  the three-point correlation function
 *  
 *  @author Federico Marulli, Michele Moresco
 *
 *  @author federico.marulli3@unbo.it, michele.moresco@unibo.it
 */

#ifndef __MODFUNCTHREEP__
#define __MODFUNCTHREEP__

#include "Cosmology.h"


// ============================================================================


namespace cosmobl {

  namespace modelling {

    namespace threept {

      /**
       *  @struct STR_data_model
       *  @brief the structure STR_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of the two-point correlation function
       */
      struct STR_data_model_threept {

	/// Q dark matter
	vector<double> Q_DM;

  /// cosmology
  shared_ptr<cosmology::Cosmology> cosmology;

  /// 1st side of the triangle
  double r1;

  /// 2nd side of the triangle
  double r2;

  /// theta
  vector<double> theta;

  /// model for the 3PCF
  string model;

  /// k vector
  vector<double> kk;

  /// Dark matter power spectrum
  vector<double> Pk_DM;

	/**
	 *  @brief default constructor
	 *  @return object of type STR_data_model
	 */
	STR_data_model_threept () = default;
      };

    }
  }
}

#endif
