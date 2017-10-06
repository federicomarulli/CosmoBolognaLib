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
 *  @file Headers/Lib/Measure.h
 *
 *  @brief The class Measure
 *
 *  This file defines the interface of the class Measure, used to
 *  manage any kind of cosmological measures
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MEASURE__
#define __MEASURE__


#include "Catalogue.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of the functions and classes used to manage
   *  any kind of <B> measure </B>
   *  
   *  The \e measure namespace contains all the functions and classes
   *  used for any kind of cosmological measure
   */
  namespace measure {

      /**
       *  @enum ErrorType
       *  @brief the two-point correlation function error type
       */
      enum ErrorType { 

	/// Poissonian error
	_Poisson_,

	/// Jackknife resampling
	_Jackknife_,

	/// Bootstrap resampling
	_Bootstrap_,
  
	/// Jackknife resampling, test
	_JackknifeTest_,

	/// Bootstrap resampling, test
	_BootstrapTest_,

	/// No error computed
	_None_	
      };

    /**
     *  @class Measure Measure.h
     *  "Headers/Lib/Measure.h"
     *
     *  @brief The class Measure
     *
     *  This is the base class used to measure the two-point
     *  correlation function
     *
     */
    class Measure {
    
    protected:
      
      /// the dataset of the measure
      shared_ptr<data::Data> m_dataset;
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Measure
       */
      Measure () = default;

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Measure () = default;

      ///@}
	
	
      /**
       *  @name Member functions to get the private/protected members
       */
      ///@{
      
      /**
       *  @brief get the protected member dataset
       *  @return a shared pointer to the dataset
       */
      virtual shared_ptr<data::Data> dataset () const { return m_dataset; }
      
      ///@}
      
    };
  }
}
  
#endif
