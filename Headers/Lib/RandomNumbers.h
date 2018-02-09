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
 *  @file Headers/Lib/RandomNumbers.h
 *
 *  @brief Class functions used to generate random numbers
 * 
 *  This file contains the classes to generate random numbers from various 
 *  distributions
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __RANNUM__
#define __RANNUM__ 


// =====================================================================================


namespace cosmobl {

  /**
   * @brief generic distribution function
   * @return distribution function
   */
  typedef function<double(double, shared_ptr<void>, vector<double>)> distribution_func;

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
     *  "Headers/Lib/RandomNumbers.h"
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
      mt19937_64 m_generator;

      /// seed
      int m_seed;

      /// minimum value to generate
      double m_MinVal;

      /// maximum value to generate
      double m_MaxVal;

    public:

      /**
       * @brief default constructor
       * @return object of class RandomNumbers
       */
      RandomNumbers () = default;
      
      /**
       *  @brief constructor
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @return object of class RandomNumbers
       */
      RandomNumbers (const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  @return none
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
       *  @return none
       */
      void set_seed (const int seed);

      /**
       *  @brief set the range for the random number extraction
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @return none
       */
      void set_range (const double MinVal, const double MaxVal);

      /**
       *  @brief set the value for constant distribution
       *  @param value the value to be returned
       *  @return none
       */
      void set_value(const double value) 
      { (void)value; ErrorCBL("Error in set_value() of RandomNumbers.h"); }

      /**
       *  @brief set the mean for Poisson distribution
       *  @param mean the Poisson distribution mean
       *  @return none
       */
      virtual void set_mean (const double mean)
      { (void)mean; ErrorCBL("Error in set_mean() of RandomNumbers.h"); }

      /**
       *  @brief set parameters for Normal distribution
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       *  @return none
       */
      virtual void set_mean_sigma (const double mean, const double sigma)
      { (void)mean; (void)sigma; ErrorCBL("Error in set_mean_sigma() of RandomNumbers.h"); }

      /**
       *  @brief set parameters for Discrete distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       *  @return none
       */
      virtual void set_discrete_values (const vector<double> values, const vector<double> weights)
      { (void)values; (void)weights; ErrorCBL("Error in set_discrete_values() of RandomNumbers.h"); }

      /**
       *  @brief set the parameters for the interpolated distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       *  @return none
       */
      virtual void set_interpolated_distribution (const vector<double> values, const vector<double> weights)
      { (void)values; (void)weights; ErrorCBL("Error in set_interpolated_distribution() of RandomNumbers.h"); }

      /**
       *  @brief set parameters for interpolated distribution
       *  @param func the distribution func 
       *  @param fixed_pars fixed parameters
       *  @param pars distribution free parameters
       *  @return none
       */
      virtual void set_custom_distribution (const distribution_func func, const shared_ptr<void> fixed_pars, const vector<double> pars)
      { (void)func; (void)fixed_pars; (void)pars; ErrorCBL("Error in set_custom_distribution() of RandomNumbers.h"); }

    };

    /**
     *  @class ConstantRandomNumbers RandomNumbers.h
     *  "Headers/Lib/RandomNumbers.h"
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
       *  @return object of class ConstantRandomNumbers
       */
      ConstantRandomNumbers (const double value) : RandomNumbers(1)
      {
	m_value=value;
      }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~ConstantRandomNumbers () = default; 

      /**
       *  @brief set the value for constant distribution
       *  @param value the value to be returned
       *  @return none
       */
      void set_value(const double value) {m_value = value;}

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
     *  "Headers/Lib/RandomNumbers.h"
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
      shared_ptr<uniform_real_distribution<double>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @param seed the random number generator seed
       *  @return object of class UniformRandomNumbers
       */
      UniformRandomNumbers (double MinVal, const double MaxVal, const int seed);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~UniformRandomNumbers () = default; 

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    
    /**
     *  @class PoissonRandomNumbers RandomNumbers.h
     *  "Headers/Lib/RandomNumbers.h"
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
      shared_ptr<poisson_distribution<int> > m_distribution;

    public:
      
      /**
       *  @brief constructor
       *  @param mean the Poisson distribution mean
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @return object of class PoissonRandomNumbers
       */
      PoissonRandomNumbers (const double mean, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~PoissonRandomNumbers () = default;

      /**
       *  @brief set the mean for Poisson distribution
       *  @param mean the Poisson distribution mean
       *  @return none
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
     *  "Headers/Lib/RandomNumbers.h"
     *
     *  @brief The class NormalRandomNumbers
     *
     *  The base class to generate random numbers
     *  following a Normal distribution
     *
     */
    class NormalRandomNumbers : public RandomNumbers
    {
      
    protected:

      /// mean
      double m_mean;

      /// standard deviation
      double m_sigma;

      /// normal distributionnormal distributionnormal distributionnormal distribution
      shared_ptr<normal_distribution<double>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @return object of class NormalRandomNumbers
       */
      NormalRandomNumbers (const double mean, const double sigma, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~NormalRandomNumbers () = default;

      /**
       *  @brief set parameters for Normal distribution
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       *  @return none
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
     *  "Headers/Lib/RandomNumbers.h"
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
      vector<double> m_values;

      /// weights for the values
      vector<double> m_weights;

      /// discrete distribution
      shared_ptr<discrete_distribution<int>> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param values values of the discrete distribution
       *  @param weights weights of the discrete distribution
       *  @param seed the random number generator seed
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @return object of class RandomNumbers
       */
      DiscreteRandomNumbers (const vector<double> values, const vector<double> weights, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~DiscreteRandomNumbers () = default;

      /**
       *  @brief set parameters for Discrete distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       *  @return none
       */
      void set_discrete_values (const vector<double> values, const vector<double> weights);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
      
    };

    
    /**
     *  @class DistributionRandomNumbers RandomNumbers.h
     *  "Headers/Lib/RandomNumbers.h"
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
      shared_ptr<UniformRandomNumbers> m_uniform_generator;

      /// interpolated distribution
      shared_ptr<glob::FuncGrid> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param xx values of the distribution
       *  @param distribution_function density probability function at xx 
       *  @param interpolation_method the method to interpolate
       *  @param seed the random number generator seed
       *  @return object of class RandomNumbers
       */
      DistributionRandomNumbers (const vector<double> xx, const vector<double> distribution_function, const string interpolation_method, const int seed);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~DistributionRandomNumbers () = default;

      /**
       *  @brief set the random number generator seed
       *  @param seed the random number generator seed
       *  @return none
       */
      void set_seed (const int seed);

      /**
       *  @brief set parameters for interpolated distribution
       *  @param xx vector containing the values with known distribution function 
       *  @param distribution_function vector containing the distribution function at xx 
       *  @param interpolation_method the method of interpolation
       *  @return none
       */
      void set_interpolated_distribution (const vector<double> xx, const vector<double> distribution_function, const string interpolation_method);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();
    };

    /**
     *  @class CustomDistributionRandomNumbers RandomNumbers.h
     *  "Headers/Lib/RandomNumbers.h"
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
      shared_ptr<UniformRandomNumbers> m_uniform_generator;

      /// the probability distribution function
      distribution_func m_func;

      /// parameters of the distribution function
      vector<double> m_func_pars;
      
      /// void pointer for the distribution function
      shared_ptr<void> m_func_fixed_pars;

      /// the distribution normalization
      double m_normalization;

    public:

      /**
       *  @brief constructor
       *
       *  @param func function
       *  @param fixed_pars function fixed parameters
       *  @param pars function free parameters
       *  @param seed the random number generator seed
       *  @param MinVal minimum value
       *  @param MaxVal maximum value
       *
       *  @return object of class RandomNumbers
       */
      CustomDistributionRandomNumbers (const distribution_func func, const shared_ptr<void> fixed_pars, const vector<double> pars, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~CustomDistributionRandomNumbers () = default;

      /**
       *  @brief set the random number generator seed
       *  @param seed the random number generator seed
       *  @return none
       */
      void set_seed (const int seed);

      /**
       *  @brief set parameters for interpolated distribution
       *
       *  @param func function
       *  @param fixed_pars function fixed parameters
       *  @param pars function free parameters
       *
       *  @return none
       */
      void set_custom_distribution (const distribution_func func, const shared_ptr<void> fixed_pars, const vector<double> pars);

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ();

    };

  }
}

#endif
