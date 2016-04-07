/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_2D.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_2D
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_2D, used to model the 2D two-point
 *  correlation function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLING2D__
#define __MODELLING2D__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_2D
     *  Modelling_TwoPointCorrelation_2D.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_2D.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_2D
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation_2D, used to model the 2D
     *  two-point correlation function
     *
     */
    class Modelling_TwoPointCorrelation_2D : public Modelling_TwoPointCorrelation {
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class
       *  Modelling_TwoPointCorrelation_2D
       */
      Modelling_TwoPointCorrelation_2D () {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation_2D () {}

      ///@}
      
    };
  }
}

#endif
