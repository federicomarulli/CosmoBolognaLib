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
 *  @file Headers/ModelFunction_NumberCounts1D_Redshift.h
 *
 *  @brief Global functions to model number counts
 *
 *  This file contains all the prototypes of the global functions used
 *  to model redshift number counts
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCNCR__
#define __MODFUNCNCR__

#include "ModelFunction_NumberCounts.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace numbercounts {

      /**
       * @brief compute the number density as a function
       * of the redshift
       *
       * @param redshift redshift bins
       *
       * @param inputs inputs to compute the mass function
       *
       * @param parameter vector containing cosmological parameters
       *
       * @return the number density as a function of redshift
       */
      std::vector<double> number_density_redshift (const std::vector<double> redshift, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

      /**
       * @brief compute the number counts as a function
       * of redshift
       *
       * @param redshift redshift bins
       *
       * @param inputs inputs to compute the mass function
       *
       * @param parameter vector containing cosmological parameters
       *
       * @return the number counts as a function of redshift
       */
      std::vector<double> number_counts_redshift (const std::vector<double> redshift, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

    }
  }
}

#endif
