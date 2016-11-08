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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_cartesian.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_cartesian
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_cartesian, used for modelling the 2D
 *  two-point correlation function in Cartesian coordinates,
 *  &xi;(r<SUB>p</SUB>,&pi;)
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGCART__
#define __MODELLINGCART__


#include "Modelling_TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_cartesian
     *  Modelling_TwoPointCorrelation_cartesian.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_cartesian.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_cartesian
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation_cartesian, used for modelling
     *  the 2D two-point correlation function in cartesian coordinates
     *
     */
    class Modelling_TwoPointCorrelation_cartesian : public Modelling_TwoPointCorrelation2D {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Modelling_TwoPointCorrelation_cartesian
       */
      Modelling_TwoPointCorrelation_cartesian () = default;

      /**
       *  @brief constructor
       *  
       *  @param twop the two-point correlation function to model
       *
       *  @return object of type Modelling_TwoPointCorrelation_cartesian
       */
      Modelling_TwoPointCorrelation_cartesian (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
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
      Modelling_TwoPointCorrelation_cartesian (const shared_ptr<data::Data> twop_dataset)
	: Modelling_TwoPointCorrelation2D() { set_data(twop_dataset); }
      
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation_cartesian () = default;
	
      ///@}

      /**
       *  @brief set the fiducial model for dark matter 
       *  two point correlation function
       *
       *  @return none
       */
      void set_fiducial_xiDM () override;

      
    };
  }
}

#endif
