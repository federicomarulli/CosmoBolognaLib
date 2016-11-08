/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Tommaso Ronconi      *
 *  federico.marulli3@unibo.it, tommaso.ronconi@outlook.it          *
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
 *  @file Headers/Lib/ReadParameters.h
 *
 *  @brief The class ReadParameters
 *
 *  This file defines the methods of the class ReadParameters,
 *  used to read generic parameter files (*.ini)
 *
 *  @authors Federico Marulli, Tommaso Ronconi
 *
 *  @authors federico.marulli3@unibo.it, tommaso.ronconi@outlook.it
 */

#ifndef __READPARA__
#define __READPARA__

#include "Func.h"


// ============================================================================


namespace cosmobl {

  namespace glob {
    
    /**
     * @class ReadParameters ReadParameters.h
     * "Headers/Lib/ReadParameters.h"
     *
     * @brief The class ReadParameters
     *
     */
    class ReadParameters {

    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       * @brief
       *    Default empty constructor
       * @return
       *    Empty object of class ReadParameters
       */
      ReadParameters () = default;

      /**
       * @brief
       *    Constructor
       * @param parameter_file
       *    string with path and name of the parameter file to read
       * @return
       *    object of class ReadParameters
       */
      ReadParameters (const string &parameter_file);

      /**
       * @brief
       *   Default destructor
       */
      ~ReadParameters () = default;

      ///@}

  
      /**
       *  @name Member functions used to get private/protected parameters
       */
      ///@{

      /**
       * @brief
       *    Template method to get parameter value
       * @param key
       *    string, parameter name
       * @return
       *    value of the requested parameter of type T, 
       *    error message if parameter not found
       */
      template <class T> T find (const string & key) const;

      /**
       * @brief
       *    Template method to get parameter value with default value
       * @param key
       *    string, parameter name
       * @param default_value
       *    typeT parameter default value to be returned if parameter not found
       * @return
       *    value associated to key, or default_value if parameter not found
       */
      template <class T> T find (const string & key, T & default_value) const;
      
      ///@}

      
    private:

      /// map with all the parameter name/value couples
      std::unordered_map<string, string> m_parameters;

      /**
       * @brief
       *    Remove white spaces treading and leading each string (private function)
       * @param inStr
       *    string read in the parameter file (after '=')
       * @return
       *    input string with treading and leading white spaces removed
       */
      string m_trim (string & inStr);

    }; // class ReadParameters

    
    /// template method to get the parameter value in the requested T type
    template <class T> T ReadParameters::find (const string & key) const {
      T value;
      
      if (m_parameters.find(key) != m_parameters.end()) {
	stringstream tmpVal(m_parameters.at(key));
	tmpVal >> std::boolalpha >> value ;

      }
      else {
	string errorMsg = "[ReadParameters] Parameter "+key+" not found";
	ErrorCBL(errorMsg); return 0;
      }

      return value;
    }
  
    /// template method to get the parameter value in the requested T type with default
    template <class T> T ReadParameters::find (const string & key, T & default_value) const {

      T value;
      if (m_parameters.find(key) != m_parameters.end()) {
	stringstream tmpVal(m_parameters.at(key));
	tmpVal >> std::boolalpha >> value;

      } else {
	std::cout <<"Parameter " << key << " not found. Using default " << default_value << ::std::endl;
	value = default_value;
      }

      return value;
    }
    
  } // namespace glob
} // namespace cosmobl


#endif //__READPARA__
