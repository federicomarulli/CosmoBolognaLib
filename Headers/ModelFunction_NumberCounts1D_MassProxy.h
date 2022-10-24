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
       * of the mass proxy
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
      
    }
  }
}

#endif
