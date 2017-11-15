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
 *  @file Headers/Lib/LikelihoodParameters.h
 *
 *  @brief The class LikelihoodParameters
 *
 *  This file defines the interface of the class LikelihoodParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LPARAM__
#define __LPARAM__

#include "BaseParameter.h"
#include "DerivedParameter.h"

namespace cosmobl {

  namespace statistics {

    /**
     * @class LikelihoodParameters LikelihoodParameters.h
     * "Headers/Lib/LikelihoodParameters.h"
     *
     *  @brief The class LikelihoodParameters
     *
     *  This class is used to define the parameters of the likelihood
     */
    class LikelihoodParameters {

      protected:

	/// likelihood parameters
	vector<shared_ptr<Parameter>> m_parameters;

	/// number of parameters
	int m_nparameters = 0;

	/// number of free parameters
	int m_nparameters_free = 0;

	/// number of fixed parameters
	int m_nparameters_fixed = 0;

	/// number of output parameters
	int m_nparameters_output = 0;

	/// the chain size
	int m_chain_size = 0;

	/// the number of parallel walkers
	int m_nwalkers = 0;

	/// function parameter covariance matrix
	vector<vector<double>> m_parameter_covariance;

	/// indexes of fixed parameters
	vector<unsigned int> m_fixed_parameters;

	/// indexes of the free parameters
	vector<unsigned int> m_free_parameters;

	/// indexes of the output parameters
	vector<unsigned int> m_output_parameters;

	/**
	 * @brief private member to set the parameters
	 *
	 * @param parameters a vector containing pointers
	 * to the parameters
	 *
	 * @return none
	 */
	void m_set_parameter_type(vector<shared_ptr<Parameter>> parameters);

	/**
	 * @brief private member to set posterior
	 * from chains
	 *
         * @param start the starting position
         *
         * @param thin number of jumped indexes in the chain
         *
         * @param seed the posterior seed
         *
	 * @return none
	 */
	void m_set_parameter_posterior (const int start, const int thin=1, const int seed=3214);

	/**
	 * @brief private member to set internal variable 
	 * m_parameter_covariance
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @return none
	 */
	void m_set_parameter_covariance (const int start, const int thin=1);

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class LikelihoodParameters
	 */
	LikelihoodParameters () = default;

	/**
	 *  @brief constructor for LikelihoodParameters
	 *
	 * @param parameters a vector containing pointers
	 * to the parameters
	 *
	 *  @return object of class LikelihoodParameters
	 */
	LikelihoodParameters (const vector<shared_ptr<Parameter>> parameters);

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~LikelihoodParameters () = default;

	///@}
	
	/**
	 *  @name Member functions used to set private/protected of the LikelihoodParameters
	 */
	///@{

	/**
	 * @brief return the total number of
	 * parameters
	 *
	 * @return the total number of parameters
	 */
	int nparameters () const;

	/**
	 * @brief return the number of free
	 * parameters
	 *
	 * @return the number of free parameters
	 */
	int nparameters_free () const;

	/**
	 * @brief return all the likelihood parameters
	 * 
	 * @param free_parameter_values vector of free parameters
	 *
	 * @return all the likelihood parameters value
	 */
	vector<double> full_parameters (const vector<double> free_parameter_values) const;

	/**
	 * @brief return the number of fixed
	 * parameters
	 *
	 * @return the number of fixed parameters
	 */
	int nparameters_fixed() const;

	/**
	 * @brief return the number of output
	 * parameters
	 *
	 * @return the number of output parameters
	 */
	int nparameters_output() const;

	/**
	 * @brief return the i-th
	 * parameter
	 *
	 * @param i the parameter index
	 *
	 * @return the i-th parameters
	 */
	shared_ptr<Parameter> parameter (const unsigned int i) const;

	/**
	 * @brief return the parameters
	 *
	 * @return vector containing pointers to the
	 * parameters
	 */
	vector<shared_ptr<Parameter>> parameters () const;

	/**
	 *  @brief set the parameter
	 *
	 *  @param parameters vector containing pointers to the
	 * parameters
	 *
	 *  @return none
	 */
	void set_parameters (const vector<shared_ptr<Parameter>> parameters);

	/**
	 * @brief return the protected member m_parameter_covariance
	 *
	 * @param i index i
	 *
	 * @param j index j
	 *
	 * @return the value of the parameter covariance at i,j
	 */
	double parameter_covariance (const int i, const int j) const;

	/**
	 * @brief return the protected member m_parameter_covariance
	 *
	 * @return vector containing the parameter covariance
	 * matrix
	 */
	vector<vector<double>> parameter_covariance () const;

	/**
	 * @brief return the chain size
	 *
	 * @return the chain size
	 */
	int chain_size () const;

	/**
	 * @brief return the number
	 * of parallel walkers
	 *
	 * @return the number of parallel walkers
	 */
	int nwalkers () const;

	///@}
	
	/**
	 *  @name Member functions to manage fixed/free parameters
	 */
	///@{
		
	/**
	 * @brief set m_fixed to false;
	 * 
	 * @param p the p-th parameter
	 *
	 * @return none
	 */
	void free (const int p);

	/**
	 * @brief set m_fixed to true
	 *
	 * @param p the p-th parameter
	 * 
	 * @return none
	 */
	void fix (const int p);

	/**
	 * @brief fix the parameter at the input value;
	 * 
	 * @param p the p-th parameter
	 *
	 * @param value the input value
	 *
	 * @return none
	 */
	void fix (const int p, const double value);

	/**
	 * @brief fix the parameter at the bestfit value,
	 * contained in m_bestfit_value;
	 *
	 * @param p the p-th parameter
	 *
	 * @return none
	 */
	void fix_at_bestfit (const int p);

	///@}

	/**
	 *  @name Member functions used to interact with prior distribution
	 */
	///@{

	/**
	 * @brief value of the prior at a proposed value
	 *
	 * @param p the p-th parameter
	 *
	 * @param value proposed parameter value
	 *
	 * @return prior value
	 */ 
	double PriorProbability (const int p, const double value) const;

	/**
	 * @brief value of the prior at a proposed value
	 *
	 * @param value proposed parameter values
	 *
	 * @return prior value
	 */ 
	vector<double> PriorProbability (const vector<double> value) const;
	
	/**
	 * @brief value of the prior at a proposed value
	 *
	 * @param p the p-th parameter
	 *
	 * @param value proposed parameter value
	 *
	 * @return prior value
	 */ 
	double LogPriorProbability (const int p, const double value) const;

	/**
	 * @brief value of the prior at a proposed value
	 *
	 * @param value proposed parameter values
	 *
	 * @return prior value
	 */ 
	vector<double> LogPriorProbability (const vector<double> value) const;

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
         *  @param p the parameter index
	 *
	 *  @return a parameter value
	 */
	double prior_sample (const int p);

	/**
	 *  @brief extract a parameter value from the prior
	 *
	 *  @return a vector containing parameter values
	 */
	vector<double> prior_sample ();

	/**
	 *  @brief extract values from the prior 	 
	 *
	 *  @param p the p-th parameter
	 *
	 *  @param sample_size the size of the extracted sample
	 *
	 *  @return a vector containing parameter values
	 */
	vector<double> prior_sample (const int p, const int sample_size);

	/**
	 *  @brief return the size of 
	 *  the prior interval times epsilon 	 
	 *
	 *  @param p the p-th parameter
	 *
	 *  @param epsilon the fraction of the interval size
	 *
	 *  @return the scaled prior interval
	 */
	double prior_range (const int p, const double epsilon=1);

	/**
	 *  @brief return the size of 
	 *  the prior interval times epsilon 	 
	 *
	 *  @param epsilon the fraction of the interval size
	 *
	 *  @return  the scaled prior intervals
	 */
	vector<double> prior_range (const double epsilon=1);

	///@}

	/**
	 *  @name Member functions used to get/set best fit values
	 */
	///@{
	
	/**
	 * @brief get the protected member m_value
	 *
	 *  @param p the p-th parameter
	 *
	 * @return the bestfit value of the parameter
	 */
	double bestfit_value (const int p) const;
	
	/**
	 * @brief get the protected member m_value
	 *
	 * @return the parameter bestfit values
	 */
	vector<double> bestfit_values () const;

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param p the p-th parameter
	 *
	 *  @param bestfit_value parameter bestfit value
	 *
	 *  @return none
	 */
	void set_bestfit_value (const int p, const double bestfit_value);

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param bestfit_value parameter bestfit values
	 *
	 *  @return none
	 */
	void set_bestfit_value (const vector<double> bestfit_value);

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param bestfit_value parameter value
	 *
	 *  @return none
	 */
	void write_bestfit_info(const vector<double> bestfit_value);

	///@}
	
	/**
	 *  @name Member functions used to interact with posterior distribution
	 */
	///@{

	/**
	 * @brief value of the posterior at a proposed value
	 *
	 * @param p the p-th parameter
	 *
	 * @param value proporsed parameter value
	 *
	 * @return posterior value
	 */ 
	double PosteriorProbability (const int p, const double value) const;

	/**
	 * @brief value of the posterior at a proposed value
	 *
	 * @param value proporsed parameter value
	 *
	 * @return posterior value
	 */ 
	vector<double> PosteriorProbability (const vector<double> value) const;

	/**
	 * @brief get the posterior distribution mean
	 *
	 * @param p the p-th parameter
	 *
	 * @return the mean of the parameter
	 * posterior distribution
	 */
	double posterior_mean (const int p) const;

	/**
	 * @brief get the posterior distribution mean
	 *
	 * @return the mean of the parameter
	 * posterior distribution
	 */
	vector<double> posterior_mean () const;

	/**
	 * @brief get the posterior distribution median
	 *
	 * @param p the p-th parameter
	 * 
	 * @return the median value of the parameter
	 * posterior
	 */   
	double posterior_median (const int p) const;

	/**
	 * @brief get the posterior distribution median
	 *
	 * @return the median value of the parameter
	 * posterior
	 */   
	vector<double> posterior_median () const;

	/**
	 * @brief get the posterior distribution standard deviation
	 *
	 * @param p the p-th parameter
	 *
	 * @return the standard deviation value of the 
	 * parameter posterior
	 */
	double posterior_std (const int p) const;

	/**
	 * @brief get the posterior distribution standard deviation
	 *
	 * @return the standard deviation value of the 
	 * parameter posterior
	 */
	vector<double> posterior_std () const;

	/**
	 * @brief get the posterior percentile
	 *
	 * @param p the p-th parameter
	 *
	 * @param i the i-th percentile
	 *
	 * @return the posterior i-th percentile
	 */
	double posterior_percentile (const int p, const unsigned int i) const;

	/**
	 * @brief get the posterior percentile
	 *
	 * @param i the i-th percentile
	 *
	 * @return the posterior i-th percentile
	 */
	vector<double> posterior_percentile (const unsigned int i) const;

	/**
	 *  @brief extract a parameter value from the posterior
	 *
	 *  @param p the p-th parameter
	 *
	 *  @return a parameter value
	 */
	double posterior_sample (const int p);

	/**
	 * @brief extract values from the posterior 	 
	 *
	 * @param p the p-th parameter
	 *
	 * @param sample_size the size of the extracted sample
	 *
	 * @return a vector of values
	 */
	vector<double> posterior_sample (const int p, const int sample_size);

	///@}
		
	/**
	 *  @name Member functions used to interact with arameter chains
	 */
	///@{

	/**
	 * @brief return the chain values
	 *
	 * @param pos the position in the chain
	 *
	 * @param ww the walker
	 * 
	 * @param par the parameter index
	 *
	 * @return the chain values
	 */
	double chain_value (const int pos, const int ww, const int par);

	/**
	 * @brief return the chain values
	 *
	 * @param par the parameter index
	 *
	 * @param pos the position in the chain
	 *
	 * @return the chain values
	 */
	vector<double> chain_values (const int pos, const int par);
	
	/**
	 * @brief return the chain values
	 *
	 * @param pos the position in the chain
	 *
	 * @return the chain values
	 */
	vector<vector<double>> chain_values (const int pos);
	
	/**
	 * @brief set the chain values
	 *
	 * @param pos the position in the chain
	 *
	 * @param ww the walker
	 *
	 * @param par the parameter index
	 *
	 * @param value the input chain value
	 *
	 * @return none
	 */
	void set_chain_value (const int pos, const int ww, const int par, const double value);
		
	/**
	 * @brief set the chain values
	 *
	 * @param pos the position in the chain
	 *
	 * @param par the parameter index
	 *
	 * @param values the input chain values
	 *
	 * @return none
	 */
	void set_chain_values (const int pos, const int par, const vector<double> values);

	/**
	 * @brief set the chain values
	 *
	 * @param pos the position in the chain
	 *
	 * @param values the input chain values
	 *
	 * @return none
	 */
	void set_chain_values (const int pos, const vector<vector<double>> values);
	
	/**
	 * @brief set the chain
	 *
	 * @param chain_size the chain size
	 *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @return none
	 */
	void set_chains (const int chain_size, const int nwalkers);

	/**
	 * @brief set the chain
	 *
	 * @param chain_values the chain values
	 *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @return none
	 */
	void set_chains (const vector<vector<double>> chain_values, const int nwalkers);

	/**
	 * @brief set the chain
	 *
	 * @param chain_values the chain values
	 *
	 * @return none
	 */
	void set_chains (const vector<vector<vector<double>>> chain_values);

	/**
	 * @brief initialize chains setting the internal
	 * variable m_chain extracting
	 * randomly from the parameter priors
	 *
	 * @param seed base seed for prior random
	 * extraction
	 *
	 * @return none
	 */
	void initialize_chains_from_prior (const int seed);
  
	/**
	 * @brief initialize chains setting the internal
	 * variable m_chain extracting randomly close 
	 * to the bestfit values
	 *
	 * @param radius the stardand distribution for
	 * normal distribution random extraction
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void initialize_chains_around_bestfit_values (const double radius, const int seed);
  
	/**
	 * @brief initialize chains setting the internal
	 * variable m_chain extracting randomly close 
	 * to input values
	 *
	 * @param values the means  for
	 * normal distribution random extraction
	 *
	 * @param radius the stardand deviation for
	 * normal distribution random extraction
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void initialize_chains_around_values (const vector<double> values, const double radius, const int seed);

	/**
	 * @brief show the results on the standard output
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void show_results (const int start, const int thin, const int seed=34121);

	/**
	 * @brief write results on files
	 *
	 * @param dir name of the output folder
	 *
	 * @param file name of the output file
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void write_results (const string dir, const string file, const int start, const int thin, const int seed=34121);

	///@}

    };
  }
}

#endif
