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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_deprojected.h
 *
 *  @brief The class Modelling_TwoPointCorrelatoin_deprojected
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_deprojected, used to model the
 *  deprojected two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINTDEPROJ__
#define __MODELLINGTWOPOINTDEPROJ__


#include "Modelling_TwoPointCorrelation_projected.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {

    namespace twopt {
    
      /**
       *  @class Modelling_TwoPointCorrelation_deprojected
       *  Modelling_TwoPointCorrelation_deprojected.h
       *  "Headers/Lib/Modelling_TwoPointCorrelation_deprojected.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation_deprojected
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation_deprojected, used for modelling
       *  the deprojected two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation_deprojected : public Modelling_TwoPointCorrelation_projected {
      
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_TwoPointCorrelation_deprojected
	 */
	Modelling_TwoPointCorrelation_deprojected () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_deprojected
	 */
	Modelling_TwoPointCorrelation_deprojected (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop)
	  : Modelling_TwoPointCorrelation_projected(twop) {}

	/**
	 *  @brief constructor
	 *  
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_TwoPointCorrelation_deprojected
	 */
	Modelling_TwoPointCorrelation_deprojected (const shared_ptr<data::Data> twop_dataset)
	  : Modelling_TwoPointCorrelation_projected() { set_data(twop_dataset); }
      
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation_deprojected () = default;
	
	///@}

      };
    }
  }
}

#endif
