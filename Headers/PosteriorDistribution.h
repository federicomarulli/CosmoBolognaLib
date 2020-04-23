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
 *  @file Headers/PosteriorDistribution.h
 *
 *  @brief The class PosteriorDistribution 
 *
 *  This file defines the interface of the class PosteriorDistribution
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __POSTERIORDIST__
#define __POSTERIORDIST__

#include "Distribution.h"


// ===================================================================================================


namespace cbl {

  namespace statistics {

    /**
     *  @class PosteriorDistribution PosteriorDistribution.h "Headers/PosteriorDistribution.h"
     *
     *  @brief The class PosteriorDistribution
     *
     *  This class is used to define the distribution
     */
    class PosteriorDistribution : public glob::Distribution {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  @return object of class PosteriorDistribution
       */
      PosteriorDistribution () : Distribution() {}

      /**
       *  @brief constructor of a constant distribution
       *
       *  @param posteriorType the type of distribution to be created
       *
       *  @param value the value to be returned
       *
       *  @return object of class PosteriorDistribution
       */
      PosteriorDistribution (const glob::DistributionType posteriorType, const double value) : Distribution(posteriorType, value) {}

      /**
       *  @brief constructor
       *
       *  @param distributionType the type of distribution to be created
       *
       *  @param discrete_values list of discrete values 
       *
       *  @param weights list of weights for discrete values
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return object of class PosteriorDistribution
       */
      PosteriorDistribution (const glob::DistributionType distributionType, const std::vector<double> discrete_values, const std::vector<double> weights, const int seed=1) : Distribution(distributionType, discrete_values, weights, seed) {}

      /**
       * @brief constructor
       *
       * @param distributionType the type of distribution to be created
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
       * @return object of class PosteriorDistribution
       */
      PosteriorDistribution (const glob::DistributionType distributionType, const std::vector<double> var, const std::vector<double> dist, const int nbin, const std::string interpolationType, const int seed=1) : Distribution(distributionType, var, dist, nbin, interpolationType, seed) {}

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~PosteriorDistribution () = default;

      ///@}

    };
  }
}

#endif
