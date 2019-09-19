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
 *  @file Headers/ModelFunction_NumberCounts1D_Size.h
 *
 *  @brief Global functions to model number counts
 *
 *  This file contains all the prototypes of the global functions used
 *  to model redshift number counts
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCNCS__
#define __MODFUNCNCS__

#include "ModelFunction_NumberCounts.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace numbercounts {
   
      /**
       * @brief compute the size function as a function
       * of the void effective radii
       *
       * @param radii radius bins
       *
       * @param inputs inputs to compute the size function
       *
       * @param parameter vector containing the cosmological
       * parameters and eventually the bias used to re-parameterise
       * the void size function
       *
       * @return the size function as a function of void size
       */
      std::vector<double> size_function_model (const std::vector<double> radii, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
    }
  }
}

#endif
