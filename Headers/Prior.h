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
 *  @file Headers/Prior.h
 *
 *  @brief The class Prior 
 *
 *  This file defines the interface of the class Prior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PRIOR__
#define __PRIOR__

#include "Distribution.h"

// ===================================================================================================


namespace cbl {

  namespace statistics {

    /**
     * @var typedf Prior_function
     * @brief definition of a function for computation of 
     * the Prior
     */
    typedef std::function<double (const std::vector<double> , const std::shared_ptr<void>)> Prior_function;

    /**
     *  @class Prior Prior.h "Headers/Prior.h"
     *
     *  @brief The class Prior
     *
     *  This class is used to define the distribution
     */
    class Prior {

    protected:

      /// prior function
      Prior_function m_prior_function;

      /// inputs of the prior function
      std::shared_ptr<void> m_prior_function_inputs;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  @return object of class Prior
       */
      Prior () {}

      /**
       * @brief constructor
       *
       * @param prior_function the prior function
       *
       * @param prior_function_inputs inputs for the prior function
       *
       * @return object of class Prior
       */
      Prior (const Prior_function prior_function, const std::shared_ptr<void> prior_function_inputs);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Prior () = default;

      ///@}
	
      /**
       * compute the prior for input
       * parameters
       *
       * @param parameters the input parameters
       *
       * @return the prior function
       */
      double operator () (const std::vector<double> parameters);
	
      /**
       * @brief compute the logarithm of the prior for the input
       * parameters
       *
       * @param parameter vector containing the input parameters
       *
       * @return the prior function
       */
      double log (const std::vector<double> parameter);
	
    };
  }
}

#endif
