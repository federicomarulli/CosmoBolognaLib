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
 *  @file Headers/Lib/Modelling.h
 *
 *  @brief The class Modelling
 *
 *  This file defines the interface of the class
 *  Modelling, used for modelling any kind of measurements
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLING__
#define __MODELLING__


#include "Likelihood.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of functions and classes used for modelling
   *  
   * The \e modelling namespace contains all the functions and classes
   * used to model any kind of measurements
   */
  namespace modelling {
    
    /**
     *  @class Modelling Modelling.h
     *  "Headers/Lib/Modelling.h"
     *
     *  @brief The class Modelling
     *
     *  This file defines the interface of the base class Modelling,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling {

    protected:
    
      /// input data to be modelled
      shared_ptr<Data> m_data;

      /// input model
      shared_ptr<statistics::Model> m_model;

      /// likelihood
      statistics::Likelihood m_likelihood;
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Modelling
       */
      Modelling () {}

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling () {}

      ///@}


    };
  }
}

#endif
