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
 *  Headers/Lib/Modelling_ThreePointCorrelation_comoving_connected.h
 *
 *  @brief The class Modelling_ThreePointCorrelation_comoving_connected
 *
 *  This file defines the interface of the class
 *  Modelling_ThreePointCorrelation_comoving_connected, that contains
 *  all the methods to model the connected three-point correlation
 *  function in comoving coordinates
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGTHREEPOINTCOMCON__
#define __MODELLINGTHREEPOINTCOMCON__


#include "Modelling_ThreePointCorrelation.h"
#include "ModelFunction_ThreePointCorrelation_comoving_connected.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    namespace threept {
      
      /**
       *  @class Modelling_ThreePointCorrelation_comoving_connected
       *  Modelling_ThreePointCorrelation_comoving_connected.h
       *  "Headers/Lib/Modelling_ThreePointCorrelation_comoving_connected.h"
       *
       *  @brief The class
       *  Modelling_ThreePointCorrelation_comoving_connected
       *
       *  This file defines the interface of the class
       *  Modelling_ThreePointCorrelation_comoving_connected, that
       *  contains all the methods to model the connected three-point
       *  correlation function in comoving coordinates
       *
       */
      class Modelling_ThreePointCorrelation_comoving_connected : public Modelling_ThreePointCorrelation {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  Modelling_ThreePointCorrelation_comoving_connected
	 */
	Modelling_ThreePointCorrelation_comoving_connected () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param threep the three-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_ThreePointCorrelation_comoving_connected
	 */
	Modelling_ThreePointCorrelation_comoving_connected (const shared_ptr<cosmobl::measure::threept::ThreePointCorrelation> threep)
	  : Modelling_ThreePointCorrelation(threep) {}

	/**
	 *  @brief constructor
	 *  
	 *  @param threep_dataset the dataset containing the
	 *  three-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_ThreePointCorrelation_comoving_connected
	 */
	Modelling_ThreePointCorrelation_comoving_connected (const shared_ptr<data::Data> threep_dataset)
	  : Modelling_ThreePointCorrelation() { set_data(threep_dataset); }
	
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_ThreePointCorrelation_comoving_connected () = default;
	
	///@}
      
      };
    }
  }
}

#endif
