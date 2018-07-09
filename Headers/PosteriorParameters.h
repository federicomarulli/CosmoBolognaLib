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
 *  @file Headers/PosteriorParameters.h
 *
 *  @brief The class PosteriorParameters
 *
 *  This file defines the interface of the class PosteriorParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PPARAM__
#define __PPARAM__

#include "ModelParameters.h"
#include "PosteriorDistribution.h"
#include "Prior.h"

namespace cbl {

  namespace statistics {

    /**
     * @class PosteriorParameters PosteriorParameters.h
     * "Headers/PosteriorParameters.h"
     *
     *  @brief The class PosteriorParameters
     *
     *  This class is used to define the model parameters in
     *  a bayesian framework
     */
    class PosteriorParameters : public ModelParameters {

      protected:
	
	/// model parameter bestfit values
	std::vector<double> m_parameter_bestfit_value;

	/// numbers of free parameters
	size_t m_nparameters_free = 0;

	/// number of fixed parameters
	size_t m_nparameters_fixed = 0;

	/// indexes of fixed parameters
	std::vector<unsigned int> m_fixed_parameters;

	/// indexes of the free parameters
	std::vector<unsigned int> m_free_parameters;

	/// parameter prior distributions
	std::vector<std::shared_ptr<PriorDistribution>> m_parameter_prior;

	/// parameter posterior distributions
	std::vector<std::shared_ptr<PosteriorDistribution>> m_parameter_posterior;

	/// function parameter covariance matrix
	std::vector<std::vector<double>> m_parameter_covariance;

	/// the lenght of the chain
	size_t m_chain_size;

	/// the number of parallel walkers
	size_t m_chain_nwalkers;

	/// content of the chain
	std::vector<std::vector<double>> m_chain_values;

	/**
	 * @brief private member to set the parameter
	 * types
	 *
	 * @return none
	 */
	void m_set_parameter_type() override;

	/**
	 *  @brief get the position in the vector
	 *  m_values from position index and walker index
	 *
	 *  @param pp the position in the chain
	 *
	 *  @param ww the walker
	 *
	 *  @return object of class Chain.  
	 */
	int m_inds_to_index(const int pp, const int ww) const { return pp*m_chain_nwalkers+ww;}

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class PosteriorParameters
	 */
	PosteriorParameters () = default;

	/**
	 *  @brief constructor for PosteriorParameters
	 *
	 *  @param nparameters the number of parameters
	 *
	 *  @param priorDistributions the priorDistributions
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names
	 *
	 *  @return object of class PosteriorParameters
	 */
	PosteriorParameters (const size_t nparameters, const std::vector<std::shared_ptr<PriorDistribution>> priorDistributions, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames);

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~PosteriorParameters () = default;

	///@}
	
	/**
	 * @brief return the number of free
	 * parameters
	 *
	 * @return the number of free parameters
	 */
	size_t nparameters_free () const override;

	/**
	 * @brief return the number of fixed
	 * parameters
	 *
	 * @return the number of fixed parameters
	 */
	size_t nparameters_fixed() const override;

	/**
	 * @brief return all the model parameters
	 * 
	 * @param parameter_values vector of free parameters
	 *
	 * @return all the parameter values
	 */
	std::vector<double> full_parameters (const std::vector<double> parameter_values) const override;

	/**
	 *  @brief set the parameter
	 *
	 *  @param nparameters the number of parameters
	 *
	 *  @param priorDistributions vector containing the parameter
	 *  priors
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names
	 *
	 *  @return none
	 */
	void set_parameters (const size_t nparameters, const std::vector<std::shared_ptr<PriorDistribution>> priorDistributions, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames) override;

	/**
	 *  @name Member functions used to set private/protected of the PosteriorParameters
	 */
	///@{
			
	/**
	 * @brief set the internal method
	 *  m_parameter_covariance
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @return none
	 */
	void set_parameter_covariance(const int start=0, const int thin=1);

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
	std::vector<std::vector<double>> parameter_covariance () const;

	///@}


	/**
	 *  @name Member functions used to interact with prior distribution
	 */
	///@{
	
	/**
	 * @brief set the prior distribution for the p-th 
	 * parameter
	 *
	 * @param p the p-th parameter
	 *
	 * @param priorDistribution the prior distribution
	 *
	 * @return none
	 */
	void set_prior_distribution (const int p, const std::shared_ptr<PriorDistribution> priorDistribution);

	/**
	 * @brief set the prior distributions for the 
	 * parameters
	 *
	 * @param priorDistributions the prior distributions
	 *
	 * @return none
	 */
	void set_prior_distribution (const std::vector<std::shared_ptr<PriorDistribution>> priorDistributions);

	/**
	 * @brief set the prior distributions for the 
	 * parameters
	 *
	 * @param ran_generator the random generator
	 *
	 * @return none
	 */
	void set_prior_distribution_seed (const std::shared_ptr<random::UniformRandomNumbers_Int> ran_generator);
	
	/**
	 * @brief get the prior distribution for the p-th 
	 * parameter
	 *
	 * @param p the p-th parameter
	 *
	 * @return pointer to the prior distribution of the p-th
	 * parameter
	 */
	std::shared_ptr<PriorDistribution> prior_distribution (const int p) const {return m_parameter_prior[p];}

	/**
	 * @brief get the prior distribution for the p-th 
	 * parameter
	 *
	 * @return vector containing pointers to the prior distributions
	 */
	std::vector<std::shared_ptr<PriorDistribution>> prior_distribution () const {return m_parameter_prior;}

	/**
	 * @brief get the prior function
	 *
	 * @return pointer to a class Prior
	 */
	std::shared_ptr<Prior> prior () const;

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
	std::vector<double> bestfit_values () const;

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param bestfit_value parameter bestfit values
	 *
	 *  @return none
	 */
	void set_bestfit_value (const std::vector<double> bestfit_value);

	/**
	 *  @brief write bestfit info
	 *
	 *  @return none
	 */
	void write_bestfit_info();

	///@}
	
	/**
	 *  @name Member functions used to interact with posterior distribution
	 */
	///@{
			
	/**
	 * @brief set the posterior distribution from
	 * the chains
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param nbins the number of bins
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void set_posterior_distribution (const int start, const int thin, const int nbins, const int seed=34121);

	/**
	 *  @brief get the posterior distribution for the 
	 *  chosen parameter
	 *
	 *  @param par the index of the parameter
	 *
	 *  @return the protected member m_parameter_posterior[param]
	 */
	std::shared_ptr<PosteriorDistribution> posterior_distribution (const int par) const { return m_parameter_posterior[par]; }

	///@}

	/**
	 * @brief show the results on the standard output
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param nbins the number of bins
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void show_results (const int start, const int thin, const int nbins, const int seed=34121);

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
	 * @param nbins the number of bins
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	void write_results (const std::string dir, const std::string file, const int start, const int thin, const int nbins, const int seed=34121);

	/**
	 * @brief return the private member m_chain_size
	 *
	 * @return the chain size
	 */
	size_t chain_size () const { return m_chain_size; }

	/**
	 * @brief return the private member m_chain_nwalkers
	 *
	 * @return the chain size
	 */
	size_t chain_nwalkers () const { return m_chain_nwalkers; }

	/**
	 * @brief set the chain
	 *
	 * @param size the chain lenght 
         *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @return none
	 */
	void set_chain (const size_t size, const size_t nwalkers);

	/**
	 * @brief reset the chain using m_size and m_nwalkers
	 *
	 * @return none
	 */
	void reset_chain ();

	/**
	 * @brief expand the already existing chain
	 *
	 * @param append the lenght of the empty chunk of the chain 
	 *
	 * @return none
	 */
	void expand_chain (const int append);

	/**
	 * @brief return the private member m_chain_values at the pos step
	 * for the ww-th walker, for the chosen parameter
	 *
	 * @param param the parameter index
	 * 
	 * @param pos the position in the chain
	 * 
	 * @param ww the walker index
	 *
	 * @return the chain value
	 */
	double chain_value (const int param, const int pos, const int ww) const { return m_chain_values[param][m_inds_to_index(pos, ww)]; }

	/**
	 * @brief return the private member m_values at the pp-th step
	 * for the ww-th step for all the parameters
	 * 
	 * @param pos the position in the chain
	 * 
	 * @param ww the walker index
	 *
	 * @return the chain value
	 */
	std::vector<double> chain_value_parameters (const int pos, const int ww) const;

	/**
	 * @brief return all the chain values for a parameter
	 * 
	 * @param param the parameter index
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @return the chain value
	 */
	std::vector<double> parameter_chain_values (const int param, const int start=0, const int thin = 1) const;

	/**
	 * @brief set the private member m_chain_values at the pp-th step
	 * for the ww-th step
	 *
	 * @param param the parameter index
	 * 
	 * @param pos the position in the chain
	 * 
	 * @param ww the walker index
	 *
	 * @param value the chain value
	 *
	 * @return none
	 */
	void set_chain_value (const int param, const int pos, const int ww, const double value) { m_chain_values[param][m_inds_to_index(pos, ww)] = value; }

	/**
	 *  @brief set the chain values
	 *
	 *  @param values the input chain values
	 *
	 *  @param nwalkers the number of parallel walkers
	 *
	 *  @return none
	 */
	void set_chain_values (const std::vector<std::vector<double>> values, const int nwalkers);

	/**
	 *  @brief set the chain values
	 *
	 *  @param values the input chain values
	 *
	 *  @return none
	 */
	void set_chain_values (const std::vector<std::vector<std::vector<double>>> values);

	/**
	 * @brief initialize the chain values
	 *
	 * @param values the starting values
	 *
	 * @return none
	 */
	void initialize_chain (const std::vector<std::vector<double>> values);

	/**
	 * @brief initialize the chain values
	 * random sampling the parameter priors
	 *
	 * @return none
	 */
	void initialize_chain_from_prior ();

	/**
	 * @brief initialize the chain values
	 *
	 * @param center the ball center
	 *
	 * @param radius the ball radius
	 *
	 * @param seed the random number generator seed
	 *
	 * @return none
	 */
	void initialize_chain_ball (const std::vector<double> center, const double radius, const double seed);

	/**
	 * @brief initialize the chain values around bestfit 
	 * values
	 *
	 * @param radius the ball radius
	 *
	 * @param seed the random number generator seed
	 *
	 * @return none
	 */
	void initialize_chain_ball_bestfit (const double radius, const double seed);
    };
  }
}

#endif
