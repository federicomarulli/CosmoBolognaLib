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
 *  @file Headers/Modelling_TwoPointCorrelation2D_polar.h
 *
 *  @brief The class Modelling_TwoPointCorrelation2D_polar
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation2D_polar, used for modelling the 2D
 *  two-point correlation function in Polar coordinates,
 *  &xi;(r<SUB>p</SUB>,&pi;)
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINT2DPOL__
#define __MODELLINGTWOPOINT2DPOL__


#include "Modelling_TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cbl {
  
  namespace modelling {

    namespace twopt {
      
      /**
       *  @class Modelling_TwoPointCorrelation2D_polar
       *  Modelling_TwoPointCorrelation2D_polar.h
       *  "Headers/Modelling_TwoPointCorrelation2D_polar.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation2D_polar
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation2D_polar, used for modelling
       *  the 2D two-point correlation function in polar coordinates
       *
       */
      class Modelling_TwoPointCorrelation2D_polar : public Modelling_TwoPointCorrelation2D {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  Modelling_TwoPointCorrelation2D_polar
	 */
	Modelling_TwoPointCorrelation2D_polar () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation2D_polar
	 */
	Modelling_TwoPointCorrelation2D_polar (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
	  : Modelling_TwoPointCorrelation2D(twop) {}

	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_monopole
	 */
	Modelling_TwoPointCorrelation2D_polar (const std::shared_ptr<data::Data> twop_dataset)
	  : Modelling_TwoPointCorrelation2D() { set_data(twop_dataset); }
      
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation2D_polar () = default;
	
	///@}

      
      };
    }
  }
}

#endif
