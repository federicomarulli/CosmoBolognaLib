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

  /**
   *  @brief The namespace of functions and classes used for statistical
   *  analysis
   *  
   *  The \e statistic namespace contains all the functions and classes
   *  used for statistical analyis
   */
  namespace statistics {

    typedef function<double(double, shared_ptr<void>, vector<double>)> prior_func;

    /**
     * @enum PriorType
     * @brief the two-point correlation function error type
     */
    enum PriorType {
      
      /// Identity function
      _IdentityPrior_,

      ///Gaussian function
      _GaussianPrior_,

      /// Poisson function
      _PoissonPrior_,
 
      /// User defined Function
      _FunctionPrior_,
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

      /// parameters of the prior func
      vector<double> m_prior_func_pars;

      /// the prior lower limit
      double m_xmin;

      /// the prior upper limit
      double m_xmax;

      /// discrete prior
      bool m_Discrete;

      /// discrete values
      vector<double> m_discrete_values;

      /// Normalization of the prior distribution
      double m_prior_normalization;

      /// Maximum probability of the prior distribution
      double m_prior_max_prob;

      /// Inverse cumulative distribution of the prior
      shared_ptr<classfunc::func_grid_GSL> m_prior_inverse_cumulative;

      /* function to set normalization and
       * maximum probability of the prior function, used
       * to sample values
       *
       * return none
       */
      void set_distribution_parameters();

      /* function to set the inverse comulative
       * distribution of the prior
       * return none
       */
      void set_inverse_cumulative_distribution();

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
      Prior (); 

      /**
       *  @brief constructor of a flat prior
       *
       *  @param priorType the type of prior to be created
       *
       *  @param pmin lower limit of the prior
       *
       *  @param pmax upper limit of the prior
       *
       *  @param discrete_values discrete values for the prior
       *
       *  @return object of class Prior
       */
      Prior (const PriorType priorType, const double pmin, const double pmax, const vector<double> discrete_values = {});

      /**
       *  @brief constructor
       *
       *  @param priorType the type of prior to be created
       *
       *  @param prior_params parameters of the prior function or discrete
       *  list of values for discrete prior
       *
       *  @param Limits limits of the prior
       *
       *  @param discrete_values discrete values for the prior
       *
       *  @return object of class Prior
       */
      Prior (const PriorType priorType, const vector<double> prior_params, const vector<double> Limits = {}, const vector<double> discrete_values = {});

      /**
       *  @brief constructor
       *
       *  @param priorType the type of prior to be created
       *  @param func user-defined function for the prior
       *  @param prior_params parameters of the prior function or discrete
       *  list of values for discrete prior
       *  @param Limits limits of the prior
       *
       *  @param discrete_values discrete values for the prior
       *
       *  @return object of class Prior
       */
      Prior (const PriorType priorType, const prior_func func, const vector<double> prior_params, const vector<double> Limits = {}, const vector<double> discrete_values = {});

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Prior () {} 

      ///@}

      /**
       * @brief evaluate prior 
       *
       * @param xx the value for prior calculation
       *
       * @return the prior value
       */
      double operator() (double xx) {
	shared_ptr<void> pp = NULL;
	if (xx<m_xmin || xx>m_xmax) return 1.e30;
	else return m_func(xx, pp, m_prior_func_pars)/m_prior_normalization;
      }

      /**
       * @brief set the prior limits 
       *
       * @param pmin lower limit of the prior
       *
       * @param pmax upper limit of the prior
       *
       * @return none
       */
      void set_limits (const double pmin, const double pmax);

      /**
       * @brief set parameters for gaussian prior 
       *
       * @param prior_params parameters of the prior function
       *
       * @return none
       */
      void set_parameters (const vector<double> prior_params);

      /**
       * @brief set prior function 
       *
       * @param func the prior function
       *
       * @param prior_params the prior function parameters
       *
       * @return none
       */       
      void set_func_parameters (const prior_func func, const vector<double> prior_params);

      /**
       * @brief set discrete values 
       *
       * @param discrete_values vector containing discrete values
       *
       * @return none
       */       
      void set_discrete_values (const vector<double> discrete_values);

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
       * @brief check if the prior is discrete
       *
       * return 0 &rarr; continuous; 1 &rarr; discrete
       */
      bool isDiscrete () const { return m_Discrete;} 

      /**
       * @brief check if a value is included in the prior limits
       *
       * @return 0 &rarr; not included in prior range; 1 &rarr; included in prior range
       */
      bool isIncluded (const double value) const;

      /**
       * @brief sample a value from the prior
       *
       * @return value of the parameter
       */
      double sample ();

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

      /**
       * @brief sample values from the prior
       *
       * @param prob the probability
       *
       * @return values of the parameter
       */
      double sample(const double prob);

      /**
       * @brief sample a value from the prior
       *
       * @param value value to be checked
       *
       * @return value of the parameter
       */
      double apply_discrete (const double value) const;
	  
    };
  }
}

#endif
