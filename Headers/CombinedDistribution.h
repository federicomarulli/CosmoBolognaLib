/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Sofia Contarini      *
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
 *  @file Headers/CombinedDistribution.h
 *
 *  @brief The class CombinedDistribution 
 *
 *  This file defines the interface of the class CombinedDistribution
 *
 *  @authors Federico Marulli, Sofia Contarini
 *
 *  @authors federico.marulli3@unibo.it, sofia.contarini3@unibo.it
 */

#ifndef __COMBDISTR__
#define __COMBDISTR__

#include "Distribution.h"
#include "EigenWrapper.h"
#include "ChainMesh.h"

// ===================================================================================================


namespace cbl {

  namespace glob {

    /**
     *  @class Distribution CombinedDistribution.h "Headers/Distribution.h"
     *
     *  @brief The class CombinedDistribution
     *
     *  This class is used to define the N-dimensional distribution
     */
    class CombinedDistribution : public Distribution {
          
    protected:

      /// vector of objects Distribution
      std::vector<std::shared_ptr<Distribution>> m_distributionVec;

      /// the distribution lower limits
      std::vector<double> m_xMinVec;

      /// the distribution upper limits
      std::vector<double> m_xMaxVec;

      /// the distribution inputs
      std::shared_ptr<void> m_inputs;

      /// the n-dimensional probability distribution function
      nDim_distribution_func m_Func;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       */
      CombinedDistribution () : Distribution () {}

      /**
       *  @brief constructor of a constant distribution
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param value the value to be returned
       *  
       */
      CombinedDistribution (const glob::DistributionType priorType, const double value) : Distribution(priorType, value) {}

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
       */
      CombinedDistribution (const glob::DistributionType priorType, const double xmin, const double xmax, const int seed=1)
	: Distribution(priorType, xmin, xmax, seed) {}
      
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
       */
      CombinedDistribution (const glob::DistributionType priorType, const std::vector<double> prior_params, const double xmin, const double xmax, const int seed=1)
	: Distribution(priorType, prior_params, xmin, xmax, seed) {}
      
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
       */
      CombinedDistribution (const glob::DistributionType priorType, const distribution_func prior_func, const std::shared_ptr<void> prior_fixed_pars, const std::vector<double> prior_pars, const double xmin, const double xmax, const int seed=1)
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
       */
      CombinedDistribution (const glob::DistributionType priorType, const std::vector<double> discrete_values, const std::vector<double> weights, const int seed=1)
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
       */
      CombinedDistribution (const glob::DistributionType priorType, const std::vector<double> var, const std::vector<double> dist, const int nbin, const std::string interpolationType, const int seed=1)
	: Distribution(priorType, var, dist, nbin, interpolationType, seed) {}

      /**
       * @brief constructor of multidimensional distributions
       *
       * @param distributionType the type of combined distribution to
       * be created (only Gaussian type available for the moment)
       *
       * @param meanVec vector containing the mean of the distributions
       *
       * @param covMat the covariance matrix of the multidimensional distribution 
       *
       * @param xMinVec vector containing the minima of the distributions
       *
       * @param xMaxVec vector containing the maxima of the distributions
       *
       * @param seed the distribution seed for random sampling
       *
       */
      CombinedDistribution (const DistributionType distributionType, const std::vector<double> meanVec, const std::vector<std::vector<double>> covMat, const std::vector<double> xMinVec, const std::vector<double> xMaxVec, const int seed=3213);

      /**
       * @brief constructor of multidimensional distributions from
       * external chains
       *
       * @param filename the name of the file to read
       *
       * @param path the path where the file is stored
       *
       * @param columns_to_read vector of integers correspondent to
       * the columns to read. The first columns has index 1. The last
       * number of the vector represents the column with the values of
       * the Posterior distribution
       *
       * @param skip_nlines the number of lines to skip at the
       * beginning of the file
       *
       * @param type_data the type of data read as last column: 0
       * \f$\rightarrow\f$ log(posterior), 1 \f$\rightarrow\f$
       * posterior, 2 \f$\rightarrow\f$ chi2
       *
       * @param normalize if true the posterior distribution is
       * normalized between 0 and 1
       *
       * @param distNum the maximum number of points used to
       * interpolate
       *
       * @param rMAX the radius used to search for close points
       *
       * @param cell_size the size of cells of the normalised chain
       * mesh (each side has lenght 100)
       *
       */
      CombinedDistribution (const std::string filename, const std::string path, const std::vector<int> columns_to_read, const int skip_nlines=0, const int type_data=0, const bool normalize=true, const int distNum=200, const double rMAX=2, const double cell_size=2);

      /**
       *  @brief default destructor
       */
      ~CombinedDistribution () = default;

      ///@}

      /**
       * @brief evaluate distribution 
       *
       * @param xx the n-dimensional coordinate for distribution calculation
       *
       * @return the distribution value
       */
      double operator [] (std::vector<double> xx);

      /**
       * @brief the i-th distribution 
       *
       * @param index the value for distribution calculation
       *
       * @return shared pointer to the object of the class
       * Distribution correspondent to a certain index
       */
      std::shared_ptr<Distribution> get_distribution (int i) {return m_distributionVec[i];}

      /**
       * @brief the size of the distribution vector
       *
       * @return the size of the distribution vector
       */
      size_t get_size_distribution () {return m_distributionVec.size();}

    };
  }
}

#endif
