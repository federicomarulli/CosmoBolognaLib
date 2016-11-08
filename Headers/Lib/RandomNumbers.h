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
      default_random_engine m_generator;

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
      RandomNumbers (const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble)
      {
	set_seed(seed);
	set_range(MinVal, MaxVal);
      }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~RandomNumbers () = default;


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
      void set_seed (const int seed)
      {
	m_seed = seed;
	m_generator.seed(m_seed);
      }

      /**
       *  @brief set the range for the random number extraction
       *  @param MinVal lower limit of the random numbers range
       *  @param MaxVal upper limit of the random numbers range
       *  @return none
       */
      void set_range (const double MinVal, const double MaxVal)
      {
	m_MinVal = MinVal;
	m_MaxVal = MaxVal;
      }

      /**
       *  @brief set the mean for Poisson distribution
       *  @param mean the Poisson distribution mean
       *  @return none
       */
      virtual void set_mean (const double mean)
      { (void)mean; ErrorCBL("Error in set_parameters() of RandomNumbers.h"); }

      /**
       *  @brief set parameters for Normal distribution
       *  @param mean the Normal distribution mean
       *  @param sigma the Normal distribution standard deviation
       *  @return none
       */
      virtual void set_mean_sigma (const double mean, const double sigma)
      { (void)mean; (void)sigma; ErrorCBL("Error in set_parameters() of RandomNumbers.h"); }

      /**
       *  @brief set parameters for Discrete distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       *  @return none
       */
      virtual void set_discrete_values (const vector<double> values, const vector<double> weights)
      { (void)values; (void)weights; ErrorCBL("Error in set_parameters() of RandomNumbers.h"); }

      /**
       *  @brief set the parameters for the interpolated distribution
       *  @param values the values to be extracted
       *  @param weights the values weights
       *  @return none
       */
      virtual void set_interpolated_distribution (const vector<double> values, const vector<double> weights)
      { (void)values; (void)weights; ErrorCBL("Error in set_parameters() of RandomNumbers.h"); }
      
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
      UniformRandomNumbers (double MinVal, const double MaxVal, const int seed) : RandomNumbers(seed, MinVal, MaxVal)
      {
	m_distribution = make_shared<uniform_real_distribution<double> >(uniform_real_distribution<double>(0,1));
      }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~UniformRandomNumbers () = default; 

      double operator () ()
      {
	return (m_MaxVal-m_MinVal)*m_distribution->operator()(m_generator)+m_MinVal;
      }

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
      PoissonRandomNumbers (const double mean, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble) : RandomNumbers(seed, MinVal, MaxVal)
      {
	set_mean(mean);
      }

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
      void set_mean (const double mean)
      {
	m_mean = mean;
	m_distribution = make_shared<poisson_distribution<int> >(poisson_distribution<int>(mean));
      }

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ()
      {
	double val = m_distribution->operator()(m_generator);

	while (val>=m_MaxVal || val<=m_MinVal)
	  val = m_distribution->operator()(m_generator);

	return val;
      }
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
      NormalRandomNumbers (const double mean, const double sigma, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble) : RandomNumbers(seed, MinVal, MaxVal)
	{
	  set_mean_sigma(mean, sigma);
	}

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
      void set_mean_sigma (const double mean, const double sigma)
      {
	m_mean = mean;
	m_sigma = sigma;
	m_distribution = make_shared<normal_distribution<double> >(normal_distribution<double>(mean, sigma));
      }

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ()
      {
	double val = m_distribution->operator()(m_generator);

	while (val>=m_MaxVal || val<=m_MinVal)
	  val = m_distribution->operator()(m_generator);

	return val;
      }
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
      DiscreteRandomNumbers (const vector<double> values, const vector<double> weights, const int seed, const double MinVal = par::defaultDouble, const double MaxVal = -par::defaultDouble) : RandomNumbers(seed, MinVal, MaxVal)
	{
	  set_discrete_values(values, weights);
	}

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
      void set_discrete_values (const vector<double> values, const vector<double> weights)
      {
	if (weights.size()==0) {
	  m_values = values;
	  m_weights.erase(m_weights.begin(), m_weights.end());
	  fill(m_weights.begin(), m_weights.end(), 1.);
	}
	else if (weights.size()!=values.size())
	  ErrorCBL("Error in set_parameters of DiscreteRandomNumbers.h: value and weight vectors have different sizes!");
	else {
	  m_values = values;
	  m_weights = weights;
	}
      
	for (size_t i=0; i<m_values.size(); i++)
	  if (m_values[i]>=m_MaxVal || m_values[i]<=m_MinVal)
	    m_weights[i]=0;
      
	m_distribution = make_shared<discrete_distribution<int> >(discrete_distribution<int>(m_weights.begin(), m_weights.end()));
      }

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ()
      {
	return m_values[m_distribution->operator()(m_generator)];
      }
      
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
      shared_ptr<classfunc::func_grid_GSL> m_distribution;

    public:

      /**
       *  @brief constructor
       *  @param xx values of the distribution
       *  @param distribution_function density probability function at xx 
       *  @param interpolation_method the method to interpolate
       *  @param seed the random number generator seed
       *  @return object of class RandomNumbers
       */
      DistributionRandomNumbers (const vector<double> xx, const vector<double> distribution_function, const string interpolation_method, const int seed) : RandomNumbers()
	{
	  set_interpolated_distribution(xx, distribution_function, interpolation_method);
	  m_uniform_generator = make_shared<UniformRandomNumbers>(0., 1., seed);
	}

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
      void set_seed (const int seed)
      {
	m_uniform_generator->set_seed(seed);
      }

      /**
       *  @brief set parameters for interpolated distribution
       *  @param xx vector containing the values with known distribution function 
       *  @param distribution_function vector containing the distribution function at xx 
       *  @param interpolation_method the method of interpolation
       *  @return none
       */
      void set_interpolated_distribution (const vector<double> xx, const vector<double> distribution_function, const string interpolation_method)
      {
	classfunc::func_grid_GSL ff(xx,distribution_function,interpolation_method);
	double norm = ff.integrate_qag(Min(xx), Max(xx));

	vector<double> FX;
	for (size_t i=0; i<xx.size(); i++)
	  FX.push_back(ff.integrate_qag(Min(xx), xx[i])/norm);

	m_distribution = make_shared<classfunc::func_grid_GSL>(FX, xx, interpolation_method);
      }

      /**
       * @brief extract number from the distribution
       * @return random values
       */
      double operator () ()
      {
	return m_distribution->operator()(m_uniform_generator->operator()());
      }
    };
    
  }
}

#endif
