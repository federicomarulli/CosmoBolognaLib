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
 *  @file Headers/Lib/Prior.h
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

#include "Chain.h"


// ============================================================================================


namespace cosmobl {

  namespace statistics {

    /**
     * @var typedef prior_func
     * @brief definition of a function for the priors 
     */
    typedef function<double(double, shared_ptr<void>, vector<double>)> prior_func;

    /**
     * @enum PriorType
     * @brief the two-point correlation function error type
     */
    enum PriorType {
      
      /// Identity function
      _UniformPrior_,

      ///Gaussian function
      _GaussianPrior_,

      /// Poisson function
      _PoissonPrior_,
      
      /// User defined Tabulated Function
      _InterpolatedPrior_,

      /// Discrete prior
      _DiscretePrior_
      
    };

    /**
     *  @class Prior Prior.h "Headers/Lib/Prior.h"
     *
     *  @brief The class Prior
     *
     *  This class is used to define the prior
     */
    class Prior {

    protected:

      /// the prior function
      prior_func m_func;

      /// the prior random generator
      shared_ptr<random::RandomNumbers> m_prior_random;

      /// the prior lower limit
      double m_xmin;

      /// the prior upper limit
      double m_xmax;

      /// parameters of the prior func
      vector<double> m_prior_func_pars;

      /// void pointer for the prior func
      shared_ptr<void> m_prior_func_fixed_pars;

      /// prior normalization
      double m_prior_normalization;

      /**
       * @brief set prior normalization 
       * @return none
       */
      void m_set_prior_normalization ();


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
      Prior () : Prior(statistics::PriorType::_UniformPrior_, 0., 0.) {}

      /**
       *  @brief constructor of a flat prior
       *
       *  @param priorType the type of prior to be created
       *
       *  @param xmin lower limit of the prior
       *
       *  @param xmax upper limit of the prior
       *
       *  @param seed the prior seed for random sampling
       *
       *  @return object of class Prior
       */
      Prior (const PriorType priorType, const double xmin, const double xmax, const int seed=1);

      /**
       *  @brief constructor
       *
       *  @param priorType the type of prior to be created
       *
       *  @param prior_params parameters of the prior function or discrete
       *  list of values for discrete prior
       *
       *  @param xmin lower limit of the prior
       *
       *  @param xmax upper limit of the prior
       *
       *  @param seed the prior seed for random sampling
       *
       *  @return object of class Prior
       */
      Prior (const PriorType priorType, const vector<double> prior_params, const double xmin, const double xmax, const int seed=1);

      /**
       *  @brief constructor
       *
       *  @param priorType the type of prior to be created
       *
       *  @param discrete_values list of discrete values 
       *
       *  @param weights list of weights for discrete values
       *
       *  @param seed the prior seed for random sampling
       *
       *  @return object of class Prior
       */
      Prior (const PriorType priorType, const vector<double> discrete_values, const vector<double> weights, const int seed=1);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Prior () = default;

      ///@}

      /**
       * @brief evaluate prior 
       *
       * @param xx the value for prior calculation
       *
       * @return the prior value
       */
      double operator() (double xx) 
      {
	if (xx<m_xmin || xx>m_xmax) return 0;
	else return m_func(xx, m_prior_func_fixed_pars, m_prior_func_pars)/m_prior_normalization;
      }

      /**
       * @brief set prior seed
       * @param seed the prior seed
       * return none
       */
      void set_seed (const int seed) {m_prior_random->set_seed(seed);}

      /**
       * @brief set the prior limits 
       *
       * @param xmin lower limit of the prior
       *
       * @param xmax upper limit of the prior
       *
       * @return none
       */
      void set_limits (const double xmin, const double xmax);

      /**
       * @brief set an uniform prior with input limits and seed
       *
       * @param xmin lower limit of the prior
       *
       * @param xmax upper limit of the prior
       *
       * @param seed the prior seed for random sampling
       *
       * @return none
       */
      void set_uniform_prior (const double xmin, const double xmax, const int seed=1);

      /**
       * @brief set normal prior 
       *
       * @param mean the normal distribution mean
       *
       * @param sigma the normal distribution standard deviation
       *
       * @param seed the prior seed for random sampling
       *
       * @return none
       */
      void set_gaussian_prior (const double mean, const double sigma, const int seed=1);

      /**
       * @brief set poisson prior  
       *
       * @param mean the poisson distribution mean
       *
       * @param seed the prior seed for random sampling
       *
       * @return none
       */       
      void set_poisson_prior (const double mean, const int seed=1);

      /**
       * @brief set discrete prior values and weights 
       *
       * @param discrete_values vector containing discrete values
       *
       * @param weights list of weights for discrete values
       *
       * @param seed the prior seed for random sampling
       *
       * @return none
       */       
      void set_discrete_values (const vector<double> discrete_values, const vector<double> weights, const int seed=1);

      /**
       * @brief return the private member m_xmin
       *
       * @return the prior lower limit
       */
      double xmin () const { return m_xmin; }

      /**
       * @brief return the private member m_xmax
       *
       * @return the prior upper limit
       */
      double xmax () const { return m_xmax; }

      /**
       * @brief check if a value is included in the prior limits
       * @param value the value to be checked
       * @return 0 &rarr; not included in prior range; 1 &rarr; included in prior range
       */
      bool isIncluded (const double value) const;

      /**
       * @brief sample a value from the prior
       *
       * @return value of the parameter
       */
      double sample () const;

      /**
       * @brief sample a value from the prior
       *
       * @param seed the seed for random number generation
       *
       * @return value of the parameter
       */
      double sample(const int seed);

      /**
       * @brief sample values from the prior
       *
       * @param nvalues the number of points to be generated
       *
       * @return values of the parameter
       */
      vector<double> sample_vector(const int nvalues);

    };
  }
}

#endif
