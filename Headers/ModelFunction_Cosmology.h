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
 *  @file Headers/ModelFunction_Cosmology.h
 *
 *  @brief Global functions to model two-point correlation functions
 *  of any type
 *
 *  This file contains all the prototypes of the global functions used
 *  to model two-point correlation functions of any type
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCCOSM__
#define __MODFUNCCOSM__

#include "Cosmology.h"
#include "Modelling_Cosmology_DistancePrior.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace cosmo {

      struct STR_data_model_cosmology{

	/// cosmology
	std::shared_ptr<cbl::cosmology::Cosmology> cosmology;

	/// cosmological parameters
	std::vector<cbl::cosmology::CosmologicalParameter> Cpar;

	/// data type
	std::vector<std::string> data_type;

	std::shared_ptr<CMB_DistancePrior> distance_prior;

      };

      /**
       *  @brief generic cosmological measurements
       *
       *  the function provides theoretical predictions for the following
       *  quantities:
       *
       *  - DV
       *  - DV/rs
       *  - rs/DV
       *  - F_AP
       *  - DA
       *  - DA/rs
       *  - rs/DA
       *  - HH
       *  - HH*rs
       *  - DH
       *  - DH/rs
       *  - rs/DH
       *  - DL
       *
       *  @param redshift the redshift
       *  @param data_type the data type
       *  @param cosmology the cosmology
       *
       *  @return cosmological measurements
       */     
      double cosmological_measurements(const double redshift, const std::string data_type, const cbl::cosmology::Cosmology cosmology);

      /**
       *  @brief model for the generic cosmological measurements
       *
       *  the function provides theoretical predictions for the following
       *  quantities:
       *
       *  - \f$D_V\f$
       *  - \f$D_V/r_s\f$
       *  - \f$r_s/D_V\f$
       *  - \f$F_{AP}\f$
       *  - \f$D_A\f$
       *  - \f$D_A/r_s\f$
       *  - \f$r_s/D_A\f$
       *  - \f$D_M\f$
       *  - \f$D_M/r_s\f$
       *  - \f$r_s/D_M\f$
       *  - \f$HH\f$
       *  - \f$HH*r_s\f$
       *  - \f$D_H\f$
       *  - \f$D_H/r_s\f$
       *  - \f$r_s/D_H\f$
       *  - \f$D_L\f$
       *
       *  @param redshift the redshift at which the model is computed
       *
       *  @param inputs pointer to the structure that contains inputs
       *  quantities
       *
       *  @param parameter cosmological parameters
       *
       *  @return cosmological measurements
       */
      std::vector<double> cosmological_measurements_model(const std::vector<double> redshift, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       *  @brief model for the generic cosmological measurements plus distance
       *  prior from CMB experiments
       *
       *  the function provides theoretical predictions for the following
       *  quantities:
       *
       *  - \f$D_V\f$
       *  - \f$D_V/r_s\f$
       *  - \f$r_s/D_V\f$
       *  - \f$D_A\f$
       *  - \f$D_A/r_s\f$
       *  - \f$r_s/D_A\f$
       *  - \f$D_M\f$
       *  - \f$D_M/r_s\f$
       *  - \f$r_s/D_M\f$
       *  - \f$HH\f$
       *  - \f$HH*r_s\f$
       *  - \f$D_H\f$
       *  - \f$D_H/r_s\f$
       *  - \f$r_s/D_H\f$
       *  - \f$D_L\f$
       *
       *  @param redshift the redshift at which the model is computed
       *
       *  @param inputs pointer to the structure that contains inputs
       *  quantities
       *
       *  @param parameter cosmological parameters
       *  
       *  @return cosmological measurements
       */
      std::vector<double> cosmological_measurements_model_CMB_DistancePrior (const std::vector<double> redshift, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    }
  }
}

#endif
