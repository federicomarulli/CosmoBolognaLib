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
 *  @file Headers/Distribution.h
 *
 *  @brief The class Distribution 
 *
 *  This file defines the interface of the class Distribution
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __DISTR__
#define __DISTR__

#include "Data2D_extra.h"
#include "Func.h"


// ===================================================================================================


namespace cbl {

  namespace glob {

    /**
     *  @enum DistributionType
     *  @brief the distribution type
     */
    enum class DistributionType {

      /// Constant function
      _Constant_,

      /// Identity function
      _Uniform_,

      /// Gaussian function
      _Gaussian_,

      /// Poisson function
      _Poisson_,

      /// Custom function
      _Custom_,

      /// User defined Tabulated Function
      _Interpolated_,

      /// Discrete distribution
      _Discrete_

    };

    /**
     * @brief return a vector containing the
     * DistributionType names
     * @return a vector containing the
     * DistributionType names
     */
    inline std::vector<std::string> DistributionTypeNames () { return {"Constant", "Uniform", "Gaussian", "Poisson", "Custom", "Interpolated", "Discrete"}; }

    /**
     * @brief cast an enum of type DistributionType
     * from its index
     * @param distributionTypeIndex the distributionType index
     * @return object of class DistributionType
     */
    inline DistributionType DistributionTypeCast (const int distributionTypeIndex) { return castFromValue<DistributionType>(distributionTypeIndex); }

    /**
     * @brief cast an enum of type DistributionType
     * from its name
     * @param distributionTypeName the distributionType name
     * @return object of class DistributionType
     */
    inline DistributionType DistributionTypeCast (const std::string distributionTypeName) { return castFromName<DistributionType>(distributionTypeName, DistributionTypeNames()); }

    /**
     * @brief cast an enum of type DistributionType
     * from indeces
     * @param distributionTypeIndeces the distributionType indeces
     * @return vector of objects of class DistributionType
     */
    inline std::vector<DistributionType> DistributionTypeCast (const std::vector<int> distributionTypeIndeces) { return castFromValues<DistributionType>(distributionTypeIndeces); } 

    /**
     * @brief cast an enum of type DistributionType
     * from thier names
     * @param distributionTypeNames the distributionType names
     * @return vector of objects of class DistributionType
     */
    inline std::vector<DistributionType> DistributionTypeCast (const std::vector<std::string> distributionTypeNames) { return castFromNames<DistributionType>(distributionTypeNames, DistributionTypeNames()); }

    
    /**
     *  @class Distribution Distribution.h "Headers/Distribution.h"
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
      std::shared_ptr<random::RandomNumbers> m_distribution_random;

      /// the distribution lower limit
      double m_xmin;

      /// the distribution upper limit
      double m_xmax;

      /// parameters of the distribution function
      std::vector<double> m_distribution_func_pars;

      /// void pointer for the distribution function
      std::shared_ptr<void> m_distribution_func_fixed_pars;

      /// distribution normalization
      double m_distribution_normalization;
      
      /// natural log of distribution normalization
      double m_log_distribution_normalization;

      /// distribution mean
      double m_mean;

      ///distribution variance
      double m_variance;

      /**
       * @brief set distribution normalization 
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
       */
      Distribution () : Distribution (DistributionType::_Uniform_, 0., 0.) {}

      /**
       *  @brief constructor of a constant distribution
       *
       *  @param distributionType the type of distribution to be created
       *
       *  @param value the value to be returned
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
       */
      Distribution (const DistributionType distributionType, const std::vector<double> distribution_params, const double xmin, const double xmax, const int seed=1);

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
       */
      Distribution (const DistributionType distributionType, const distribution_func func, const std::shared_ptr<void> distribution_fixed_pars, const std::vector<double> distribution_pars, const double xmin, const double xmax, const int seed=1);

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
       */
      Distribution (const DistributionType distributionType, const std::vector<double> discrete_values, const std::vector<double> weights, const int seed=1);

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
       */
      Distribution (const DistributionType distributionType, const std::vector<double> var, const std::vector<double> dist, const int nbin, const std::string interpolationType, const int seed=1);

      /**
       *  @brief default destructor
       */
      ~Distribution () = default;

      ///@}

      /**
       * @brief return the distribution type 
       *
       * @return the distribution type
       */
      DistributionType distributionType () const
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
      double operator () (double xx);

      /**
       * @brief evaluate log-distribution 
       *
       * @param xx the value for distribution calculation
       *
       * @return the log-distribution value
       */
      double log_distribution (double xx);

      /**
       * @brief set distribution seed
       * @param seed the distribution seed
       */
      void set_seed (const int seed) {m_distribution_random->set_seed(seed);}

      /**
       * @brief set the distribution limits 
       *
       * @param xmin lower limit of the distribution
       *
       * @param xmax upper limit of the distribution
       */
      void set_limits (const double xmin, const double xmax);

      /**
       * @brief set a constant distribution
       *
       * @param value the value to be returned
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
       */
      void set_gaussian_distribution (const double mean, const double sigma, const int seed=1);

      /**
       * @brief set poisson distribution  
       *
       * @param mean the poisson distribution mean
       *
       * @param seed the distribution seed for random sampling
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
       */
      void set_custom_distribution (const distribution_func func, const std::shared_ptr<void> distribution_fixed_pars, const std::vector<double> distribution_pars, const int seed=1);

      /**
       * @brief set discrete distribution values and weights 
       *
       * @param discrete_values vector containing discrete values
       *
       * @param weights list of weights for discrete values
       *
       * @param seed the distribution seed for random sampling
       */       
      void set_discrete_values (const std::vector<double> discrete_values, const std::vector<double> weights, const int seed=1);

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
       */    
      void set_binned_distribution (const std::vector<double> var, const std::vector<double> dist, const std::string interpolationType="Spline", const int seed=1);

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
      std::vector<double> sample_vector (const int nvalues);

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
      std::vector<double> moments ();

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

      /**
       *  @brief derive and store the number distribution of a given
       *  std::vector 
       *  @param [out] xx std::vector containing the binned values of the
       *  variable 
       *  @param [out] fx std::vector containing the binned values of the
       *  distribution
       *  @param [out] err std::vector containing the binned Poisson errors
       *  @param [in] FF std::vector containing the given set of data
       *  @param [in] WW std::vector containing the weights
       *  @param [in] nbin the number of bins
       *  @param [in] linear true &rarr; linear binning; false &rarr; logarithmic
       *  binning
       *  @param [in] file_out the output file where the distribution is
       *  stored
       *  @param [in] fact factor used to normalized the distribution
       *  @param [in] V1 the minimum limit of the distribution
       *  @param [in] V2 the maximum limit of the distribution
       *  @param [in] bin_type true &rarr; dn/dvar; false &rarr; dn/dlogvar
       *  @param [in] conv true &rarr; compute the Gaussian convolvolution of
       *  the distribution; false &rarr; do not convolve
       *  @param [in] sigma &sigma; of the Gaussian kernel
       */
      void get_distribution (std::vector<double> &xx, std::vector<double> &fx, std::vector<double> &err, const std::vector<double> FF, const std::vector<double> WW, const int nbin, const bool linear=true, const std::string file_out=par::defaultString, const double fact=1., const double V1=par::defaultDouble, const double V2=par::defaultDouble, const bool bin_type=true, const bool conv=false, const double sigma=0.);

    };
  }
}

#endif
