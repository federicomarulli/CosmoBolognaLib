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
 *  @file Headers/ParameterFile.h
 *
 *  @brief The class ParameterFile
 *
 *  This file defines the methods of the class ParameterFile,
 *  used to manage generic parameter files
 *
 *  @authors Alfonso Veropalumbo, Federico Marulli
 *
 *  @authors alfonso.veropalumbo@uniroma3.it, federico.marulli3@unibo.it
 */

#ifndef __PARAM__
#define __PARAM__

#include "Kernel.h"


// ============================================================================


namespace cbl {

  namespace glob {

    /// Map type used in ParameterFile class
    typedef std::unordered_map<std::string, std::vector<std::string>> parameter_map;

    /**
     * @class ParameterFile ParameterFile.h
     * "Headers/ParameterFile.h"
     *
     * @brief The class ParameterFile
     *
     */
    class ParameterFile
    {

      private:

	/// map with all the vector type parameter name/value couples
	parameter_map m_parameters; 

	/**
	 * @brief
	 *    Remove white spaces treading and leading each std::string (private function)
	 * @param inStr
	 *    std::string read in the parameter file (after '=')
	 * @return
	 *    input std::string with treading and leading white spaces removed
	 */
	std::string m_trim (const std::string inStr);

	/**
	 * @brief
	 *    Stores values contained in between curly brackets in a vector of std::string (private function)
	 * @param inStr
	 *    std::string read in the parameter file (after '=', in between curly brackets '{', '}');
	 *    values inside brackets must be separated by comas ','
	 * @return
	 *    a vector of std::string values
	 */
	std::vector<std::string> m_trim_vect (const std::string inStr);

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 * @brief
	 *    Default empty constructor
	 * @return
	 *    Empty object of class ParameterFile
	 */
	ParameterFile () = default;

	/**
	 * @brief
	 *    Constructor
	 * @param parameter_file
	 *    std::string with path and name of the parameter file to read
	 * @return
	 *    object of class ParameterFile
	 */
	ParameterFile (const std::string parameter_file);

	/**
	 * @brief
	 *   Default destructor
	 */
	~ParameterFile () = default;

	///@}

	/**
	 *  @name Member functions used to get private/protected parameters
	 */
	///@{
	/**
	 * @brief member to read a parameter file
	 *
	 * @param parameter_file string with path and name of 
	 * the parameter file to read
	 *
	 * @return none
	 */

	void read(const std::string parameter_file);

	/**
	 * @brief member to write a parameter file
	 *
	 * @param parameter_file string with path and name of 
	 * the parameter file to write
	 *
	 * @return none
	 */
	void write(const std::string parameter_file);

	///@}

	/**
	 *  @name Member functions used to get private/protected parameters
	 */
	///@{
	
	/**
	 * @brief Method to get entries of the parameter file
	 *
	 * @param key std::string, parameter name
	 * @return reference to the values of the requested parameter of type T, 
	 * 	   if the key is not found is automatically created.
	 */
	std::vector<std::string> & operator[] (const std::string key); 

	/** 
	 * @brief Method to get entries of the parameter file
	 *
	 * @param key std::string, parameter name
	 * 
	 * @return value of the requested parameter of type T, 
	 *    if the key is not found returns the default values
	 */
	std::vector<std::string> get (const std::string key, const std::vector<std::string> default_values) const;

	///@}

    }; // class ParameterFile

  } // namespace glob
} // namespace cbl


#endif //__PARAM__
