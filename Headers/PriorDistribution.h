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
 *  @file Headers/PriorDistribution.h
 *
 *  @brief The class PriorDistribution 
 *
 *  This file defines the interface of the class PriorDistribution
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PRIORDIST__
#define __PRIORDIST__

#include "Distribution.h"


// ===================================================================================================


namespace cbl {

  namespace statistics {

    /**
     *  @class PriorDistribution PriorDistribution.h "Headers/PriorDistribution.h"
     *
     *  @brief The class PriorDistribution
     *
     *  This class is used to define the distribution
     */
    class PriorDistribution : public glob::Distribution {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  
       */
      PriorDistribution () : Distribution () {}

      /**
       *  @brief constructor of a constant distribution
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param value the value to be returned
       *
       *  
       */
      PriorDistribution (const glob::DistributionType priorType, const double value) : Distribution(priorType, value) {}

      /**
       *  @brief constructor of a flat distribution
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  
       */
      PriorDistribution (const glob::DistributionType priorType, const double xmin, const double xmax, const int seed=1) : Distribution(priorType, xmin, xmax, seed) {}

      /**
       *  @brief constructor
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param prior_params parameters of the distribution function
       *  or discrete list of values for discrete distribution
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  
       */
      PriorDistribution (const glob::DistributionType priorType, const std::vector<double> prior_params, const double xmin, const double xmax, const int seed=1) : 
	Distribution(priorType, prior_params, xmin, xmax, seed) {}

      /**
       *  @brief constructor
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param prior_func the functional form of the distribution
       *
       *  @param prior_fixed_pars the fixed parameters
       *
       *  @param prior_pars the distribution parameters
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  
       */
      PriorDistribution (const glob::DistributionType priorType, const distribution_func prior_func, const std::shared_ptr<void> prior_fixed_pars, const std::vector<double> prior_pars, const double xmin, const double xmax, const int seed=1)
	: Distribution(priorType, prior_func, prior_fixed_pars, prior_pars, xmin, xmax, seed) {}

      /**
       *  @brief constructor
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param discrete_values list of discrete values 
       *
       *  @param weights list of weights for discrete values
       *
       *  @param seed the distribution seed for random sampling
       *
       *  
       */
      PriorDistribution (const glob::DistributionType priorType, const std::vector<double> discrete_values, const std::vector<double> weights, const int seed=1)
	: Distribution(priorType, discrete_values, weights, seed) {}

      /**
       * @brief constructor
       *
       * @param priorType the type of distribution to be created
       *
       * @param var vector containing binned values
       *
       * @param dist list of distribution values for each bin
       *
       * @param nbin the number of bins
       *
       * @param interpolationType the kind of interpolation
       *
       * @param seed the distribution seed for random sampling
       *
       * 
       */
      PriorDistribution (const glob::DistributionType priorType, const std::vector<double> var, const std::vector<double> dist, const int nbin, const std::string interpolationType, const int seed=1)
	: Distribution(priorType, var, dist, nbin, interpolationType, seed) {}

      /**
       *  @brief default destructor
       *
       *  
       */
      ~PriorDistribution () = default;

      ///@}

    };
  }
}

#endif
