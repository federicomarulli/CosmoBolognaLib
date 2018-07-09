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
 *  @file Headers/Modelling_TwoPointCorrelation1D_angular.h
 *
 *  @brief The class Modelling_TwoPointCorrelation1D_angular
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation1D_angular, used to model the angular
 *  of two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINTANG__
#define __MODELLINGTWOPOINTANG__


#include "Modelling_TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cbl {
  
  namespace modelling {

    namespace twopt {
    
      /**
       *  @class Modelling_TwoPointCorrelation1D_angular
       *  Modelling_TwoPointCorrelation1D_angular.h
       *  "Headers/Modelling_TwoPointCorrelation1D_angular.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation1D_angular
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation1D_angular, used for modelling
       *  the angular two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation1D_angular : public Modelling_TwoPointCorrelation1D {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  Modelling_TwoPointCorrelation1D_angular
	 */
	Modelling_TwoPointCorrelation1D_angular () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation1D_angular
	 */
	Modelling_TwoPointCorrelation1D_angular (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
	  : Modelling_TwoPointCorrelation1D(twop) {}
	
	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation1D_angular
	 */
	Modelling_TwoPointCorrelation1D_angular (const std::shared_ptr<data::Data> twop_dataset)
	  : Modelling_TwoPointCorrelation1D(twop_dataset, cbl::measure::twopt::TwoPType::_1D_angular_) {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation1D_angular () = default;
	
	///@}

      };
    }
  }
}

#endif
