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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_projected.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_projected
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_projected, used to model the
 *  projected two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGPROJ__
#define __MODELLINGPROJ__


#include "Modelling_TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_projected
     *  Modelling_TwoPointCorrelation_projected.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_projected.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_projected
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation_projected, used to model the
     *  projected of the two-point correlation function
     *
     */
    class Modelling_TwoPointCorrelation_projected : public Modelling_TwoPointCorrelation1D {

    public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingTwoPointCorrelation_projected
	 */
      Modelling_TwoPointCorrelation_projected () = default;
      
      /**
       *  @brief constructor
       *  
       *  @param twop the two-point correlation function to model
       *
       *  @return object of type
       *  Modelling_TwoPointCorrelation_projected
       */
      Modelling_TwoPointCorrelation_projected (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
	: Modelling_TwoPointCorrelation1D(twop) {}

      /**
       *  @brief constructor
       *  
       *  @param twop_dataset the dataset containing the two-point
       *  correlation function to model
       *
       *  @return object of type
       *  Modelling_TwoPointCorrelation_projected
       */
      Modelling_TwoPointCorrelation_projected (const shared_ptr<data::Data> twop_dataset)
	: Modelling_TwoPointCorrelation1D() { set_data(twop_dataset); }
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Modelling_TwoPointCorrelation_projected () = default;
	
      ///@}

      /**
       * @brief set the fiducial model for dark matter 
       * two point correlation function
       *
       *  @return none
       */
      void set_fiducial_xiDM () override;

    };
  }
}

#endif
