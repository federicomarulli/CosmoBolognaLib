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
 *  @file Headers/Lib/BaseParameter.h
 *
 *  @brief The class BaseParameter
 *
 *  This file defines the interface of the class BaseParameter, used
 *  to manage input model parameters in statistical analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __INPARAM__
#define __INPARAM__

#include "Parameter.h"


// ===================================================================================================


namespace cosmobl {

  namespace statistics {

    /**
     *  @class BaseParameter BaseParameter.h
     *  "Headers/Lib/BaseParameter.h"
     *
     *  @brief The class BaseParameter
     *
     *  This class is used to define the input parameters of models
     */
    class BaseParameter : public Parameter {

      protected:

	/// the parameter value
	double m_value = 0.;

	/// the parameter prior
	shared_ptr<statistics::Prior> m_prior;

	/// the parameter prior 
	shared_ptr<statistics::Prior> m_prior_free;

	/// false &rarr; the parameter is free; true &rarr; the parameter is fixed
	bool m_fixed = false;

	/**
	 * @brief set the prior 
	 * 
	 * @return none
	 */
	void m_set_prior ();

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class BaseParameter
	 */
	BaseParameter () = default;

	/**
	 *  @brief constructor for BaseParameter
	 *
	 *  @param value parameter value
	 *
	 *  @param name parameter name
	 *
	 *  @return object of class BaseParameter, set fixed
	 */
	BaseParameter (const double value, const string name);

	/**
	 *  @brief constructor for free parameters 	
	 *
	 *  @param prior the parameter prior
	 *
	 *  @param name parameter name
	 *
	 *  @return object of class BaseParameter
	 */
	BaseParameter (const statistics::Prior prior, const string name = "parameter");

	/**
	 *  @brief constructor for BaseParameter
         *
         *  @param xmin lower limit of the distribution
         *
         *  @param xmax upper limit of the distribution
         *
         *  @param seed the distribution seed for random sampling
         *
         *  @param name the parameter name
	 *
	 *  @return object of class BaseParameter with uniform prior
	 */
	BaseParameter (const double xmin, const double xmax, const int seed, const string name="parameter"); 

	/**
	 *  @brief constructor for BaseParameter, with Poisson or
	 *  Gaussiand distribution as a prior (see the class
	 *  Distribution)
	 *
	 *  @param priorType the type of distribution to be
	 *  created
	 *
	 *  @param prior_params parameters of the distribution
	 *  function or discrete list of values for discrete
         *  distribution
         *
         *  @param xmin lower limit of the distribution
         *
         *  @param xmax upper limit of the distribution
         *
         *  @param seed the distribution seed for random sampling
         *
         *  @param name the parameter name
         *
	 *  @return object of class BaseParameter
	 */
	 BaseParameter (const glob::DistributionType priorType, const vector<double> prior_params, const double xmin, const double xmax, const int seed, const string name="parameter");

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~BaseParameter () = default;

	///@}

	
	/**
	 *  @name Member functions used to set private/protected of the the Parameter
	 */
	///@{

	/**
	 *  @brief return private member m_fixed;
	 *
	 *  @return false &rarr; the parameter is free; true &rarr;
	 *  the parameter is fixed
	 */
	bool fixed () const;

	/**
	 * @brief set the parameter as fixed
	 *
	 * @param fix 1 &rarr; fix the parameter; 
	 * 0 &rarr; do not fix the parameter
	 *
	 * @return none
	 */
	void set_fixed (const bool fix);

	/**
	 * @brief free the parameter
	 *
	 * @return none
	 */
	void free ();

	/**
	 * @brief fix the parameter 
	 * at the m_value;
	 *
	 * @return none
	 */
	void fix ();

	/**
	 * @brief fix the parameter at the 
	 * input value;
	 * 
	 * @param value the input value
	 *
	 * @return none
	 */
	void fix (const double value);

	/**
	 * @brief fix the parameter at the 
	 * input value;
	 * 
	 * @return none
	 */
	void fix_at_bestfit ();
	
	/**
	 * @brief get the protected member m_value
	 *
	 * @return the value of the parameter
	 */
	double value () const;

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param value parameter value
	 *
	 *  @return none
	 */
	void set_value (const double value); 

	///@}
	
	/**
	 *  @name Member functions used to interact with prior distribution
	 */
	///@{
	
	/**
	 * @brief get the defined prior
	 * 
	 * @return a shared pointer to an object of class
	 * glob::Distribution
	 */
	shared_ptr<statistics::Prior> prior () const;
	
	/**
	 *  @brief set user defined prior
	 *
	 *  @param prior user defined prior
	 *
	 *  @return none
	 */
	void set_prior (const statistics::Prior prior);

	/**
	 *  @brief set an uniform prior
	 *
	 *  @param xmin lower limit of the distribution
	 *
	 *  @param xmax upper limit of the distribution
	 *
	 *  @param seed the distribution seed for random sampling
	 *
	 *  @return none
	 */
	void set_prior (const double xmin, const double xmax, const int seed);

	/**
	 *  @brief set Poisson or Gaussiand distribution prior: 
	 *  see the Documentation for the class Distribution
	 *
	 *  @param priorType the type of distribution to be created
	 *
	 *  @param prior_params parameters of the distribution function or discrete
	 *  list of values for discrete distribution
	 *
	 *  @param xmin lower limit of the distribution
	 *
	 *  @param xmax upper limit of the distribution
	 *
	 *  @param seed the distribution seed for random sampling
	 *
	 *  @return none
	 */
	void set_prior (const glob::DistributionType priorType, const vector<double> prior_params, const double xmin, const double xmax, const int seed);

	/**
	 * @brief value of the prior at a proposed value
	 *
	 * @param value proporsed parameter value
	 *
	 * @return prior value
	 */ 
	double PriorProbability (const double value) const;

	/**
	 * @brief value of the prior at a proposed value
	 *
	 * @param value proporsed parameter value
	 *
	 * @return prior value
	 */ 
	double LogPriorProbability (const double value) const;

	/**
	 *  @brief set prior seed
	 *
	 *  @param seed the prior seed
	 *
	 *  @return none
	 */
	void set_prior_seed (const int seed);

	/**
	 *  @brief extract a parameter value from the prior
	 *
	 *  @return value extracted from the prior distribution
	 */
	double prior_sample ();

	/**
	 *  @brief extract values from the prior 	 
	 *
	 *  @param sample_size the size of the extracted sample
	 *
	 *  @return a vector of values extracted from the 
	 *  prior distribution
	 */
	vector<double> prior_sample (const int sample_size);

	/**
	 *  @brief return the size of 
	 *  the prior range times epsilon 	 
	 *
	 *  @param epsilon the fraction of the interval size
	 *
	 *  @return the prior range
	 */
	double prior_range (const double epsilon=1.);

	///@}

	/**
	 *  @name Member functions used to interact with posterior distribution
	 */
	///@{ 

	/**
	 * @brief set the posterior distribution
	 * from the chain
	 *
	 * @param start the chain values
	 *
	 * @param thin the number of parallel walkers
	 *
	 * @param seed the distribution seed
	 *
	 * @return none
	 */
	void set_posterior (const int start, const int thin, const int seed);

	/**
	 * @brief set the posterior distribution
	 * from the chain
	 *
	 * @param chain_values the chain values
	 *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @param seed the distribution seed
	 *
	 * @return none
	 */
	void set_posterior (const vector<double> chain_values, const int nwalkers, const int seed);

	/**
	 * @brief set the posterior distribution
	 * from the chain
	 *
	 * @param chain_values the chain values
	 *
	 * @param seed the distribution seed
	 *
	 * @return none
	 */
	void set_posterior (const vector<vector<double>> chain_values, const int seed);

	///@}
    };
  }
}

#endif
