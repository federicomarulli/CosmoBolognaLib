/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/DerivedParameter.h
 *
 *  @brief The class DerivedParameter
 *
 *  This file defines the interface of the class DerivedParameter,
 *  used to manage derived model parameters in statistical analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __OUTPARAM__
#define __OUTPARAM__

#include "Parameter.h"


// ===================================================================================================


namespace cosmobl {

  namespace statistics {


    /**
     * @class DerivedParameter DerivedParameter.h
     * "Headers/Lib/DerivedParameter.h"
     *
     *  @brief The class DerivedParameter
     *
     *  This class is used to define the output parameters of models
     */
    class DerivedParameter : public Parameter {
      
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class Parameter
	 */
	DerivedParameter () = default;

	/**
	 *  @brief constructor for DerivedParameter	 
	 *
	 *  @param name parameter name
	 *
	 *  @return object of class Parameter
	 */
	DerivedParameter (const string name) : Parameter(ParameterType::_DerivedParameter_, name) {}

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~DerivedParameter () = default;

	/**
	 *  @name Member functions used to interact with posterior distribution
	 */
	///@{ 

	/**
	 * @brief set the posterior distribution
	 * from the chain
	 *
	 * @param start the chain values
	 *
	 * @param thin the number of parallel walkers
	 *
	 * @param seed the distribution seed
	 *
	 * @return none
	 */
	void set_posterior (const int start, const int thin, const int seed);

	/**
	 * @brief set the posterior distribution
	 * from the chain
	 *
	 * @param chain_values the chain values
	 *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @param seed the distribution seed
	 *
	 * @return none
	 */
	void set_posterior (const vector<double> chain_values, const int nwalkers, const int seed);

	/**
	 * @brief set the posterior distribution
	 * from the chain
	 *
	 * @param chain_values the chain values
	 *
	 * @param seed the distribution seed
	 *
	 * @return none
	 */
	void set_posterior (const vector<vector<double>> chain_values, const int seed);

	///@}

    };
  }
}

#endif
