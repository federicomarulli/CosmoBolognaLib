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
 *  @file Headers/RandomNumbers.h
 *
 *  @brief Class functions used to generate random numbers
 * 
 *  This file contains the classes to generate random numbers from
 *  various distributions
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __RANNUM__
#define __RANNUM__ 


#include "FFTlog.h"


// =====================================================================================


namespace cbl {

  /**
   * @brief generic distribution function
   * @return distribution function
   */
  typedef std::function<double(double, std::shared_ptr<void>, std::vector<double>)> distribution_func;

  /**
   * @brief distribution function used for a n-dimensional
   * distribution
   * @return distribution function
   */
  typedef std::function<double(std::vector<double>, std::shared_ptr<void>)> nDim_distribution_func;

  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> random numbers </B>
   *  
   *  The \e random namespace contains all the functions and classes
   *  used to handle random numbers
   */
  namespace random {

    /**
     *  @class RandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class RandomNumbers
     *
     *  The base class to generate random numbers
     *
     */
    class RandomNumbers
    {
    protected:

      /// pseudo-random numbers generator
      std::mt19937_64 m_generator;

      /// seed
      int m_seed;

      /// minimum value to generate
      double m_MinVal;

      /// maximum value to generate
      double m_MaxVal;

    public:

      /**
       * @brief default constructor
       * 
       */
      RandomNumbers () = default;
      
      /**
       *  @brief constructor
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  
       */
      RandomNumbers (const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  
       */
      virtual ~RandomNumbers () = default;


      /**
       * @brief extract number from the distribution
       * @return random values
       */
      virtual double operator () () = 0;

      /**
       *  @brief set the random number generator seed
       *  @param seed the random number generator seed
       */
      void set_seed (const int seed);

      /**
       *  @brief set the range for the random number extraction
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       */
      void set_range (const double MinVal, const double MaxVal);

      /**
       *  @brief set the value for constant distribution
       *  @param value the value to be returned
       */
      void set_value (const double value) 
      { (void)value; ErrorCBL("error!", "set_value", "RandomNumbers.h"); }

      /**
       *  @brief set the mean for Poisson distribution
       *  @param mean the Poisson distribution mean
       */
      virtual void set_mean (const double mean)
      { (void)mean; ErrorCBL("error!", "set_mean", "RandomNumbers.h"); }

      /**
       *  @brief set parameters for Normal distribution
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       */
      virtual void set_mean_sigma (const double mean, const double sigma)
      { (void)mean; (void)sigma; ErrorCBL("error!", "set_mean_sigma", "RandomNumbers.h"); }

      /**
       *  @brief set parameters for Discrete distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       */
      virtual void set_discrete_values (const std::vector<double> values, const std::vector<double> weights)
      { (void)values; (void)weights; ErrorCBL("error!", "set_discrete_values", "RandomNumbers.h"); }

      /**
       *  @brief set the parameters for the interpolated distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       */
      virtual void set_interpolated_distribution (const std::vector<double> values, const std::vector<double> weights)
      { (void)values; (void)weights; ErrorCBL("error!", "set_interpolated_distribution", "RandomNumbers.h"); }

      /**
       *  @brief set parameters for interpolated distribution
       *
       *  @param func the probability distribution function
       *
       *  @param modelInput pointer to the data possibly used to
       *  construct the probability distribution function
       *
       *  @param parameter the parameters of the probability
       *  distribution function
       */
      virtual void set_custom_distribution (const distribution_func func, const std::shared_ptr<void> modelInput, const std::vector<double> parameter)
      { (void)func; (void)modelInput; (void)parameter; ErrorCBL("error!", "set_custom_distribution", "RandomNumbers.h"); }

    };

    /**
     *  @class ConstantRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class ConstantRandomNumbers
     *
     *  This class return always the same value, 
     *  written for consistency
     *
     */
    class ConstantRandomNumbers : public RandomNumbers
    {
    protected:

      /// returned value
      double m_value;

    public:

      /**
       *  @brief constructor
       *  @param value the value to be returned
       */
      ConstantRandomNumbers (const double value) : RandomNumbers(1)
      {
	m_value = value;
      }

      /**
       *  @brief default destructor
       */
      ~ConstantRandomNumbers () = default; 

      /**
       *  @brief set the value for constant distribution
       *  @param value the value to be returned
       */
      void set_value(const double value)
      { m_value = value; }

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ()
      {
	return m_value;
      }

    };

    
    /**
     *  @class UniformRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class UniformRandomNumbers
     *
     *  The base class to generate random numbers
     *  on a interval
     *
     */
    class UniformRandomNumbers : public RandomNumbers
    {
    protected:

      /// uniform distribution
      std::shared_ptr<std::uniform_real_distribution<double>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @param seed the random number generator seed
       */
      UniformRandomNumbers (double MinVal, const double MaxVal, const int seed);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~UniformRandomNumbers () = default; 

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    
    /**
     *  @class UniformRandomNumbers_Int RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class UniformRandomNumbers_Int
     *
     *  The base class to generate random integers
     *  on a interval
     *
     */
    class UniformRandomNumbers_Int : public RandomNumbers
    {
    protected:

      /// uniform distribution
      std::shared_ptr<std::uniform_int_distribution<int>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param MinVal lower limit of the random numbers range;
       *  the corresponding lower limit will be ceil(MinVal)
       *  @param MaxVal upper limit of the random numbers range;
       *  the corresponding upper limit will be floor(MaxVal)
       *  @param seed the random number generator seed
       */
      UniformRandomNumbers_Int (double MinVal, const double MaxVal, const int seed);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~UniformRandomNumbers_Int () = default; 

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    
    /**
     *  @class PoissonRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class PoissonRandomNumbers
     *
     *  The base class to generate random numbers
     *  following a Poisson distribution
     *
     */
    class PoissonRandomNumbers : public RandomNumbers
    {
      
    protected:

      /// mean
      double m_mean;

      /// Poisson distribution
      std::shared_ptr<std::poisson_distribution<int> > m_distribution;

    public:
      
      /**
       *  @brief constructor
       *  @param mean the Poisson distribution mean
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       */
      PoissonRandomNumbers (const double mean, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~PoissonRandomNumbers () = default;

      /**
       *  @brief set the mean for Poisson distribution
       *  @param mean the Poisson distribution mean
       *  
       */
      void set_mean (const double mean);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    
    /**
     *  @class NormalRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class NormalRandomNumbers
     *
     *  The base class to generate random numbers
     *  following a Normal distribution
     */
    class NormalRandomNumbers : public RandomNumbers
    {
      
    protected:

      /// mean
      double m_mean;

      /// standard deviation
      double m_sigma;

      /// normal distributionnormal distributionnormal distributionnormal distribution
      std::shared_ptr<std::normal_distribution<double>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       */
      NormalRandomNumbers (const double mean, const double sigma, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~NormalRandomNumbers () = default;

      /**
       *  @brief set parameters for Normal distribution
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       *  
       */
      void set_mean_sigma (const double mean, const double sigma);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    
    /**
     *  @class DiscreteRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class DiscreteRandomNumbers
     *
     *  The base class to generate random numbers
     *  from a discrete sample
     *
     */
    class DiscreteRandomNumbers : public RandomNumbers
    {
      
    protected:
      
      /// discrete values
      std::vector<double> m_values;

      /// weights for the values
      std::vector<double> m_weights;

      /// discrete distribution
      std::shared_ptr<std::discrete_distribution<int>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param values values of the discrete distribution
       *  @param weights weights of the discrete distribution
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  
       */
      DiscreteRandomNumbers (const std::vector<double> values, const std::vector<double> weights, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~DiscreteRandomNumbers () = default;

      /**
       *  @brief set parameters for Discrete distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       *  
       */
      void set_discrete_values (const std::vector<double> values, const std::vector<double> weights);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
      
    };

    
    /**
     *  @class DistributionRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class DistributionRandomNumbers
     *
     *  The base class to generate random numbers
     *  from a tabulated density distribution
     *
     */
    class DistributionRandomNumbers : public RandomNumbers
    {

    protected:

      /// Uniform random number generator
      std::shared_ptr<UniformRandomNumbers> m_uniform_generator;

      /// interpolated distribution
      std::shared_ptr<glob::FuncGrid> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param xx values of the distribution
       *  @param distribution_function density probability function at xx 
       *  @param interpolation_method the method to interpolate
       *  @param seed the random number generator seed
       *  
       */
      DistributionRandomNumbers (const std::vector<double> xx, const std::vector<double> distribution_function, const std::string interpolation_method, const int seed);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~DistributionRandomNumbers () = default;

      /**
       *  @brief set the random number generator seed
       *  @param seed the random number generator seed
       *  
       */
      void set_seed (const int seed);

      /**
       *  @brief set parameters for interpolated distribution
       *  @param xx vector containing the values with known distribution function 
       *  @param distribution_function vector containing the distribution function at xx 
       *  @param interpolation_method the method of interpolation
       *  
       */
      void set_interpolated_distribution (const std::vector<double> xx, const std::vector<double> distribution_function, const std::string interpolation_method);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    /**
     *  @class CustomDistributionRandomNumbers RandomNumbers.h
     *  "Headers/RandomNumbers.h"
     *
     *  @brief The class CustomDistributionRandomNumbers
     *
     *  The base class to generate random numbers
     *  from an user-defined function 
     *
     */
    class CustomDistributionRandomNumbers : public RandomNumbers
    {

    protected:

      /// Uniform random number generator
      std::shared_ptr<UniformRandomNumbers> m_uniform_generator;

      /// the probability distribution function
      distribution_func m_func;

      /// parameters of the probability distribution function
      std::vector<double> m_func_parameter;
      
      /// pointer to the data possibly used to construct the probability the distribution function
      std::shared_ptr<void> m_func_modelInput;

      /// the distribution normalization
      double m_normalization;

    public:

      /**
       *  @brief constructor
       *
       *  @param func the probability distribution function
       *
       *  @param modelInput pointer to the data possibly used to
       *  construct the probability distribution function
       *
       *  @param parameter the parameters of the probability
       *  distribution function
       *
       *  @param seed the random number generator seed
       *
       *  @param MinVal minimum value
       *
       *  @param MaxVal maximum value
       *
       *  
       */
      CustomDistributionRandomNumbers (const distribution_func func, const std::shared_ptr<void> modelInput, const std::vector<double> parameter, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  
       */
      ~CustomDistributionRandomNumbers () = default;

      /**
       *  @brief set the random number generator seed
       *  @param seed the random number generator seed
       *  
       */
      void set_seed (const int seed);

      /**
       *  @brief set parameters for interpolated distribution
       *
       *  @param func the probability distribution function
       *
       *  @param modelInput pointer to the data possibly used to
       *  construct the probability distribution function
       *
       *  @param parameter the parameters of the probability
       *  distribution function
       *
       *  
       */
      void set_custom_distribution (const distribution_func func, const std::shared_ptr<void> modelInput, const std::vector<double> parameter);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();

    };

  }
}

#endif
