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
 *  @file Headers/Measure.h
 *
 *  @brief The class Measure
 *
 *  This file defines the interface of the class Measure, used to
 *  manage any kind of cosmological measures
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __MEASURE__
#define __MEASURE__


#include "Catalogue.h"


// ===================================================================================================


namespace cbl {
  
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
      enum class ErrorType { 

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
  
	/// Jackknife resampling, by objects
	_JackknifeObjects_,

	/// Bootstrap resampling, by objects
	_BootstrapObjects_,

      /// No error computed
      _None_	
    };

    /**
     * @brief return a vector containing the
     * ErrorType names
     * @return a vector containing the
     * ErrorType names
     */
    inline std::vector<std::string> ErrorTypeNames ()
    { return {"Poisson", "Jackknife", "Bootstrap", "JackknifeTest", "BootstrapTest", "JackknifeObject", "BootstrapObject", "None"}; }

    /**
     *
     * @brief cast an enum of type ErrorType
     * from its index
     * @param errorTypeIndex the errorType index
     * @return object of class ErrorType
     */
    inline ErrorType ErrorTypeCast (const int errorTypeIndex)
    { return castFromValue<ErrorType>(errorTypeIndex); }

    /**
     * @brief cast an enum of type ErrorType
     * from its name
     * @param errorTypeName the errorType name
     * @return object of class ErrorType
     */
    inline ErrorType ErrorTypeCast (const std::string errorTypeName)
    { return castFromName<ErrorType>(errorTypeName, ErrorTypeNames()); }

    /**
     * @brief cast an enum of type ErrorType
     * from indeces
     * @param errorTypeIndeces the errorType indeces
     * @return vector of objects of class ErrorType
     */
    inline std::vector<ErrorType> ErrorTypeCast (const std::vector<int> errorTypeIndeces)
    { return castFromValues<ErrorType>(errorTypeIndeces); }  

    /**
     * @brief cast an enum of type ErrorType
     * from thier names
     * @param errorTypeNames the errorType names
     * @return vector of objects of class ErrorType
     */
    inline std::vector<ErrorType> ErrorTypeCast (const std::vector<std::string> errorTypeNames)
    { return castFromNames<ErrorType>(errorTypeNames, ErrorTypeNames()); }

    /**
     *  @class Measure Measure.h
     *  "Headers/Measure.h"
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
      std::shared_ptr<data::Data> m_dataset;
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      Measure () = default;

      /**
       *  @brief default destructor
       *  
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
      virtual std::shared_ptr<data::Data> dataset () const { return m_dataset; }
      
      ///@}
      
    };
  }
}
  
#endif
