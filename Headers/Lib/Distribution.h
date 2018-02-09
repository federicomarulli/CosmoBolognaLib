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
 *  @file Headers/Lib/Distribution.h
 *
 *  @brief The class Distribution 
 *
 *  This file defines the interface of the class Distribution
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DISTR__
#define __DISTR__

#include "Data2D_extra.h"


// ===================================================================================================


namespace cosmobl {

  namespace glob {

    /**
     *  @enum DistributionType
     *  @brief the distribution type
     */
    enum DistributionType {

      /// Constant function
      _ConstantDistribution_,

      /// Identity function
      _UniformDistribution_,

      /// Gaussian function
      _GaussianDistribution_,

      /// Poisson function
      _PoissonDistribution_,

      /// Custom function
      _CustomDistribution_,
      
      /// User defined Tabulated Function
      _InterpolatedDistribution_,

      /// Discrete distribution
      _DiscreteDistribution_

    };
    

    /**
     *  @class Distribution Distribution.h "Headers/Lib/Distribution.h"
     *
     *  @brief The class Distribution
     *
     *  This class is used to define the distribution
     */
    class Distribution {
          
    protected:

      /// the type of distribution
      DistributionType m_distributionType;

      /// the probability distribution function
      distribution_func m_func;

      /// the distribution random generator
      shared_ptr<random::RandomNumbers> m_distribution_random;

      /// the distribution lower limit
      double m_xmin;

      /// the distribution upper limit
      double m_xmax;

      /// parameters of the distribution function
      vector<double> m_distribution_func_pars;

      /// void pointer for the distribution function
      shared_ptr<void> m_distribution_func_fixed_pars;

      /// distribution normalization
      double m_distribution_normalization;

      /// distribution mean
      double m_mean;

      ///distribution variance
      double m_variance;

      /**
       * @brief set distribution normalization 
       * @return none
       */
      void m_set_distribution_normalization ();

      /**
       * @brief integrand of the moments distribution
       *
       * @param xx the integration variable
       * @param order moment order
       *
       * @return integrand
       */
      double m_moments_integrator (const double xx, const unsigned int order)
      { return this->operator()(xx)*pow(xx, order); }

      /**
       * @brief integrand of the central moments 
       * of the distribution
       *
       * @param xx the integration variable
       * @param order moment order
       *
       * @return integrand
       */
      double m_central_moments_integrator (const double xx, const unsigned int order)
      { return this->operator()(xx)*pow(xx-m_mean, order); }

      /**
       * @brief integrand of the percentile 
       * of the distribution
       *
       * @param xx the integration variable
       *
       * @return integrand
       */
      double m_percentile_integrator (const double xx);

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  @return object of class Distribution
       */
      Distribution () : Distribution(DistributionType::_UniformDistribution_, 0., 0.) {}

      /**
       *  @brief constructor of a constant distribution
       *
       *  @param distributionType the type of distribution to be created
       *
       *  @param value the value to be returned
       *
       *  @return object of class Distribution
       */
      Distribution (const DistributionType distributionType, const double value);

      /**
       *  @brief constructor of a flat distribution
       *
       *  @param distributionType the type of distribution to be created
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return object of class Distribution
       */
      Distribution (const DistributionType distributionType, const double xmin, const double xmax, const int seed=3213);

      /**
       *  @brief constructor
       *
       *  @param distributionType the type of distribution to be created
       *
       *  @param distribution_params parameters of the distribution function or discrete
       *  list of values for discrete distribution
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return object of class Distribution
       */
      Distribution (const DistributionType distributionType, const vector<double> distribution_params, const double xmin, const double xmax, const int seed=1);

      /**
       *  @brief constructor
       *
       *  @param distributionType the type of distribution to be created
       *
       *  @param func the functional form of the distribution
       *
       *  @param distribution_fixed_pars the fixed parameters
       *
       *  @param distribution_pars the distribution parameters
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return object of class Distribution
       */
      Distribution (const DistributionType distributionType, const distribution_func func, const shared_ptr<void> distribution_fixed_pars, const vector<double> distribution_pars, const double xmin, const double xmax, const int seed=1);

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
       *  @return object of class Distribution
       */
      Distribution (const DistributionType distributionType, const vector<double> discrete_values, const vector<double> weights, const int seed=1);

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
       * @return none
       */
      Distribution (const DistributionType distributionType, const vector<double> var, const vector<double> dist, const int nbin, const string interpolationType, const int seed=1);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Distribution () = default;

      ///@}

      /**
       * @brief return the distribution type 
       *
       * @return the distribution type
       */
      double distributionType () const
      {
	return m_distributionType;
      }

      /**
       * @brief evaluate distribution 
       *
       * @param xx the value for distribution calculation
       *
       * @return the distribution value
       */
      double operator() (double xx) 
      {
	if (xx<m_xmin || xx>m_xmax) return 0;
	else return m_func(xx, m_distribution_func_fixed_pars, m_distribution_func_pars)/m_distribution_normalization;
      }

      /**
       * @brief evaluate log-distribution 
       *
       * @param xx the value for distribution calculation
       *
       * @return the log-distribution value
       */
      double log_distribution (double xx) 
      {
	if (xx<m_xmin || xx>m_xmax) return par::defaultDouble;
	else return log(m_func(xx, m_distribution_func_fixed_pars, m_distribution_func_pars)/m_distribution_normalization);
      }

      /**
       * @brief set distribution seed
       * @param seed the distribution seed
       * return none
       */
      void set_seed (const int seed) {m_distribution_random->set_seed(seed);}

      /**
       * @brief set the distribution limits 
       *
       * @param xmin lower limit of the distribution
       *
       * @param xmax upper limit of the distribution
       *
       * @return none
       */
      void set_limits (const double xmin, const double xmax);

      /**
       * @brief set a constant distribution
       *
       * @param value the value to be returned
       *
       * @return none
       */
      void set_constant_distribution (const double value);

      /**
       * @brief set an uniform distribution with input limits and seed
       *
       * @param xmin lower limit of the distribution
       *
       * @param xmax upper limit of the distribution
       *
       * @param seed the distribution seed for random sampling
       *
       * @return none
       */
      void set_uniform_distribution (const double xmin, const double xmax, const int seed=1);

      /**
       * @brief set normal distribution 
       *
       * @param mean the normal distribution mean
       *
       * @param sigma the normal distribution standard deviation
       *
       * @param seed the distribution seed for random sampling
       *
       * @return none
       */
      void set_gaussian_distribution (const double mean, const double sigma, const int seed=1);

      /**
       * @brief set poisson distribution  
       *
       * @param mean the poisson distribution mean
       *
       * @param seed the distribution seed for random sampling
       *
       * @return none
       */       
      void set_poisson_distribution (const double mean, const int seed=1);

      /**
       *  @brief set a custom distribution
       *
       *  @param func the functional form of the distribution
       *
       *  @param distribution_fixed_pars the fixed parameters
       *
       *  @param distribution_pars the distribution parameters
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return object of class Distribution
       */
      void set_custom_distribution (const distribution_func func, const shared_ptr<void> distribution_fixed_pars, const vector<double> distribution_pars, const int seed=1);

      /**
       * @brief set discrete distribution values and weights 
       *
       * @param discrete_values vector containing discrete values
       *
       * @param weights list of weights for discrete values
       *
       * @param seed the distribution seed for random sampling
       *
       * @return none
       */       
      void set_discrete_values (const vector<double> discrete_values, const vector<double> weights, const int seed=1);

      /**
       * @brief set discrete distribution values and weights 
       *
       * @param var vector containing binned values
       *
       * @param dist list of distribution values for each bin
       *
       * @param interpolationType the kind of interpolation
       *
       * @param seed the distribution seed for random sampling
       *
       * @return none
       */    
      void set_binned_distribution (const vector<double> var, const vector<double> dist, const string interpolationType="Spline", const int seed=1);

      /**
       * @brief return the private member m_xmin
       *
       * @return the distribution lower limit
       */
      double xmin () const { return m_xmin; }

      /**
       * @brief return the private member m_xmax
       *
       * @return the distribution upper limit
       */
      double xmax () const { return m_xmax; }

      /**
       * @brief the distribution range
       *
       * @return the prior range defined as 
       * xmax-xmin
       */
      double distribution_range () const {return (xmax()-xmin());}

      /**
       * @brief check if a value is included in the distribution limits
       * @param value the value to be checked
       * @return 0 &rarr; not included in distribution range; 1 &rarr; included in distribution range
       */
      bool isIncluded (const double value) const;

      /**
       * @brief sample a value from the distribution
       *
       * @return value of the parameter
       */
      double sample () const;

      /**
       * @brief sample a value from the distribution
       *
       * @param seed the seed for random number generation
       *
       * @return value of the parameter
       */
      double sample (const int seed);

      /**
       * @brief sample values from the distribution
       *
       * @param nvalues the number of points to be generated
       *
       * @return values of the parameter
       */
      vector<double> sample_vector (const int nvalues);

      /**
       * @brief return the mean value 
       * of the distribution
       *
       * @return the mean value of the distribution
       */
      double mean ();

      /**
       * @brief return the standard deviation 
       * of the distribution
       *
       * @return the distribution std
       */
      double variance ();

      /**
       * @brief return the standard deviation 
       * of the distribution
       *
       * @return the distribution std
       */
      double std ();

      /**
       * @brief return the skewness 
       * of the distribution
       *
       * @return the distribution std
       */
      double skewness ();

      /**
       * @brief return the kurtosis 
       * of the distribution
       *
       * @return the distribution std
       */
      double kurtosis ();

      /**
       * @brief return the moments of the 
       * distribution distribution 
       *
       * @return the distribution moments
       */
      vector<double> moments ();

      /**
       * @brief return the median 
       * of the distibution
       *
       * @return the distribution median
       */
      double median () { return percentile(50); }

      /**
       * @brief return the i-th percentile
       * of the distribution
       *
       * @param i the percentile
       *
       * @return the distribution median
       */
      double percentile (const unsigned int i);

      /**
       * @brief return the distribution mode
       *
       * @return the distribution mode
       */
      double mode ();

    };
  }
}

#endif
