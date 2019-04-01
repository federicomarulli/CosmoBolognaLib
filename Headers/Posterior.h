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
 *  @file Headers/Posterior.h
 *
 *  @brief The class Posterior 
 *
 *  This file defines the interface of the class Posterior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __POSTERIOR__
#define __POSTERIOR__

#include "RandomNumbers.h"
#include "Prior.h"
#include "Likelihood.h"


// ===================================================================================================


namespace cbl {

  namespace statistics {

    /**
     *  @class Posterior Posterior.h "Headers/Posterior.h"
     *
     *  @brief The class Posterior
     *
     *  This class is used to define the distribution
     */
    class Posterior : public Likelihood {

      protected:

	/// the prior distribution
	std::shared_ptr<cbl::statistics::Prior> m_prior;

	/// value of the log
	std::vector<double> m_logposterior_values;

	/// the MCMC acceptance rate
	std::vector<double> m_acceptance;

        /// general seed for prior/posterior distribution and sampler
	int m_seed;

	/// seed generator
	std::shared_ptr<cbl::random::UniformRandomNumbers_Int> m_seed_generator;

	/**
	 * @brief set the internal attribute m_seed
	 * and related quantities
	 *
	 * @param seed the seed
	 *
	 * @return none
	 */
	void m_set_seed (const int seed);

	/**
	 * @brief generate a seed
	 *
	 * @return a seed generated from m_seed_generator
	 */
	int m_generate_seed () { return m_seed_generator->operator()(); }

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class Posterior
	 */
	Posterior () {}

	/**
	 *  @brief constructor
	 *
	 *  @param prior_distributions vector containing the priors for the 
	 *  likelihood parameters
	 *
	 *  @param likelihood pointer to an object of type Likelihood
	 *
	 *  @param seed general seed for prior/posterior distribution and sampler
	 *
	 *  @return object of class Posterior
	 */
	Posterior (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const Likelihood &likelihood, const int seed=5341);

	/**
	 *  @brief constructor
	 *
	 *  @param prior_distributions vector containing the priors for the 
	 *  likelihood parameters
	 *
	 *  @param data object of type data
	 *
	 *  @param model object of type model
	 *
	 *  @param likelihood_type type of the likelihood
	 *
	 *  @param x_index index(s) of the extra info std::vector containing 
	 *  the point(s) where the model is evaluated
	 *
	 *  @param w_index std::vector containing the data point weight
	 *
	 *  @param seed general seed for prior/posterior distribution and sampler
	 *
	 *  @return object of class Posterior
	 */
	Posterior (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const int seed=5341);

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~Posterior () = default;

	///@}
	
	/**
	 * @brief return the posterior parameters
	 *
	 * @return pointer containing the posterior parameters
	 */
	std::shared_ptr<ModelParameters> parameters () const { return m_model_parameters; }

	/**
	 *  @brief evaluate the un-normalized posterior
	 *
	 *  \f[ P((\vec{theta} | \vec{d}) =
	 *  \mathcal{L}(\vec{d}|\vec{theta}) \cdot Pr(\vec{theta}) \f]
	 *
	 *  where \f$P\f$ is the
	 *  posterior,\f$\mathcal{L}(\vec{d}|\vec{theta})\f$ is the
	 *  likelihood and \f$Pr(\vec{theta})\f$ is the prior
	 *
	 *  @param pp the parameters
	 *
	 *  @return the value of the un-normalized posterior
	 */
	double operator () (std::vector<double> &pp) const;

	/**
	 *  @brief evaluate the logarithm of the un-normalized
	 *  posterior
	 *
	 *  \f[ P((\vec{theta} | \vec{d}) =
	 *  \mathcal{L}(\vec{d}|\vec{theta}) \cdot Pr(\vec{theta}) \f]
	 *
	 *  where \f$P\f$ is the
	 *  posterior,\f$\mathcal{L}(\vec{d}|\vec{theta})\f$ is the
	 *  likelihood and \f$Pr(\vec{theta})\f$ is the prior
	 *
	 *  @param pp the parameters
	 *
	 *  @return the logarithm of the un-normalized posterior
	 */
	double log (std::vector<double> &pp) const;

	/**
	 * @brief set the model for the likelihood analysis 
	 *
	 * @param model pointer to the model
	 *
	 * @param model_parameters parameters of the model
	 *
	 * @return none
	 */
	void set_model (std::shared_ptr<Model> model=NULL, const std::shared_ptr<ModelParameters> model_parameters=NULL);

	/**
	 * @brief set the posterior type using the LikelihoodType object 
	 *
	 * @param prior_distributions vector containing the priors for the 
	 * likelihood parameters
	 *
	 * @param data pointer to an object of type Data
	 *
	 * @param model pointer to an object of type model
	 *
	 * @param likelihood_type the likelihood type, specified with the 
	 * LikelihoodType object
	 *
	 * @param x_index index(s) of the extra info std::vector containing the point(s) where to evaluate the model
	 *
	 * @param w_index std::vector containing the data point weight 
	 *
	 * @param seed the seed
	 *
	 * @return none
	 */
	void set (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const int seed);

	/**
	 *  @brief function that maximize the likelihood, find the
	 *  best-fit parameters and store them in model
	 *
	 *  @param start std::vector containing initial values for
	 *  the likelihood maximization
	 *
	 *  @param parameter_limits limits for the parameters
	 *
	 *  @param max_iter the maximum number of iterations
	 *
	 *  @param tol the tolerance in finding convergence 
	 *
	 *  @param epsilon the relative fraction of the interval size
	 *
	 *  @return none
	 */
	void maximize (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3) {(void)start; (void)parameter_limits; (void)max_iter; (void)tol; (void)epsilon; ErrorCBL("Error in maximize() of Posterior.h, use the method without parameter_ranges!");}

	/**
	 *  @brief function that maximize the posterior, find the
	 *  best-fit parameters and store them in model
	 *
	 *  @param start std::vector containing initial values for
	 *  the posterior maximization
	 *
	 *  @param max_iter the maximum number of iterations
	 *
	 *  @param tol the tolerance in finding convergence 
	 *  
	 *  @param epsilon the relative size of the initial trial step
	 *
	 *  @return none
	 */
	void maximize (const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-4);

	/**
	 * @brief initialize the chains by drawing from the prior
	 * distributions
	 *
	 * the starting values of the chain are extracted from the
	 * (possibly different) distributions of the priors
	 *
	 * @param chain_size the chain lenght
	 *
	 * @param nwalkers the number of parallel
	 * chains
	 *
	 * @return none
	 */
	void initialize_chains (const int chain_size, const int nwalkers);

	/**
	 * @brief initialize the chains in a ball around the posterior
	 * best-fit parameter values
	 *
	 * the starting values of the chain are extracted from uniform
	 * distributions in the range [parameter-radius,
	 * parameter+radius] (for each likelihood parameter)
	 *
	 * this function first maximizes the posterior, starting the
	 * computation at the values of the input vector 'start', then
	 * it inizializes the chain
	 *
	 * @param chain_size the chain lenght
	 *
	 * @param nwalkers the number of parallel
	 * chains
	 *
	 * @param radius radius of the ball in parameter space
	 *
	 * @param start std::vector containing initial values for
	 * the posterior maximization
	 *
	 * @param max_iter the maximum number of iterations
	 *
	 * @param tol the tolerance in finding convergence 
	 *
	 * @return none
	 */
	void initialize_chains (const int chain_size, const int nwalkers, const double radius, const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6);

	/**
	 * @brief initialize the chains in a ball around the input
	 * parameter values
	 *
	 * the starting values of the chain are extracted from uniform
	 * distributions in the range [value[i]-radius,
	 * value[i]+radius] (for each i-th likelihood parameter)
	 *
	 * @param chain_size the chain lenght
	 *
	 * @param nwalkers the number of parallel
	 * chains
	 *
	 * @param value vector containing the input values, centres of
	 * the ball in the parameter space
	 *
	 * @param radius radius of the ball in the parameter space
	 *
	 * @return none
	 */
	void initialize_chains (const int chain_size, const int nwalkers, std::vector<double> &value, const double radius);

	/**
	 * @brief initialize the chains with input values
	 *
	 * the starting values of the chain are the elements of the
	 * input matrix 'chain_values'
	 *
	 * @param chain_size the chain lenght
	 *
	 * @param chain_value matrix of size (nwalkers, nparameters),
	 * starting values of the chain
	 *
	 * @return none
	 */
	void initialize_chains (const int chain_size, const std::vector<std::vector<double>> chain_value);

	/**
	 * @brief initialize the chains reading from an input file 
	 *
	 * the starting values of the chain are get from the last
	 * lines of an input chain file; it can be used to continue an
	 * MCMC sampling computation
	 *
	 * @param chain_size the chain lenght
	 *
	 * @param nwalkers the number of parallel
	 * chains
	 *
	 * @param input_dir input directory
	 *
	 * @param input_file the input file
	 *
	 * @return none
	 */
	void initialize_chains (const int chain_size, const int nwalkers, const std::string input_dir, const std::string input_file);

	/**
	 *  @brief sample the posterior using the stretch-move sampler
	 *  (Foreman-Mackey et al. 2012)
	 *
	 *  @param aa the parameter of the \f$g(z)\f$ distribution
	 *
	 *  @param parallel false \f$\rightarrow\f$ non-parallel
	 *  sampler; true \f$\rightarrow\f$ parallel sampler
	 *
	 *  @param outputFile output file where the chains are written
         * 	   during run-time. Leave it to default value to have
         * 	   no output.  WARNING: this option is intended for
         * 	   debug. It only works for the non-parallelized
         * 	   stretch-move algorithm.  The chain in output will
         * 	   be written in a different format with respect to
         * 	   the method
         * 	   cbl::statistics::Posterior::write_chain::ascii
         * 	   col1) chain step col2) walker index col3-npar)
         * 	   parameter values col npar+3) value of the
         * 	   posterior.  @return none
	 *	
	 *  @param start the minimum chain position used to compute
	 *  the median
	 *
	 *  @param thin the step used for chain dilution
	 *
	 *  @param nbins the number of bins to estimate the posterior
	 *  distribution, used to assess its properties 
	 *
	 *  @warning if parallel is set true, than pointers cannot be
	 *  used inside the posterior function
	 */
	void sample_stretch_move (const double aa=2, const bool parallel=true, const std::string outputFile=par::defaultString, const int start=0, const int thin=1, const int nbins=50);

	/**
	 * @brief write the chains obtained after 
	 * the MCMC sampling on an ascii file
	 *
	 * @param output_dir the output directory
	 *
	 * @param output_file the output file
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution
	 *
	 * @return none
	 */
	void write_chain_ascii (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1);

	/**
	 * @brief write the chains obtained after 
	 * the MCMC sampling on a fits file
	 *
	 * @param output_dir the output directory
	 *
	 * @param output_file the output file
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution
	 *
	 * @return none
	 */
	void write_chain_fits (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1);

	/**
	 * @brief write the chains obtained after 
	 * the MCMC sampling
	 *
	 * @param output_dir the output directory
	 *
	 * @param output_file the output file
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution
	 *
	 * @param fits false \f$\rightarrow\f$ ascii file; true
	 * \f$\rightarrow\f$ fits file
	 *
	 * @return none
	 */
	void write_chain (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1, const bool fits=false);

	/**
	 * @brief read the chains from an ascii file
	 *
	 * @param input_dir the intput directory
	 *
	 * @param input_file the intput file
	 *
	 * @param nwalkers the number of parallel chains
	 *
	 * @param skip_header the lines to be skipped in
	 * the chain file
	 *
	 * @return none
	 */
	void read_chain_ascii (const std::string input_dir, const std::string input_file, const int nwalkers, const int skip_header=1);

	/**
	 * @brief read the chains from an ascii file
	 *
	 * @param input_dir the input directory
	 *
	 * @param input_file the input file
	 *
	 * @param nwalkers the number of parallel chains
	 *
	 * @return none
	 */
	void read_chain_fits (const std::string input_dir, const std::string input_file, const int nwalkers);

	/**
	 * @brief read the chains
	 *
	 * @param input_dir the input directory
	 *
	 * @param input_file the input file
	 *
	 * @param nwalkers the number of parallel chains
	 *
	 * @param skip_header the lines to be skipped in
	 * the chain file
	 *
	 * @param fits false \f$\rightarrow\f$ ascii file; true
	 * \f$\rightarrow\f$ fits file
	 *
	 * @return none
	 */
	void read_chain (const std::string input_dir, const std::string input_file, const int nwalkers, const int skip_header=1, const bool fits=false);

	/**
	 * @brief write maximization results on a file
	 *
	 * @param output_dir the output directory 
	 *
	 * @param root_file the root of the output file to be written
	 *
	 * @return none
	 */
	void write_maximization_results (const std::string output_dir, const std::string root_file);

	/**
	 * @brief show the results of the MCMC sampling on screen
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution on screen
	 *
	 * @param nbins the number of bins to estimate the posterior
	 * distribution, used to assess its properties
	 *
	 * @param show_mode true \f$\rightarrow\f$ show the posterior
	 * mode; false \f$\rightarrow\f$ do not show the posterior
	 * mode
	 *
	 * @return none
	 */
	void show_results (const int start, const int thin, const int nbins=50, const bool show_mode=false);

	/**
	 * @brief store the results of the MCMC sampling to file
	 * 
	 * this function stores to file the posterior mean, standard
	 * deviation, median, 18th and 82th percentiles, and
	 * optionally the mode
	 *
	 * @param output_dir the output directory 
	 *
	 * @param root_file the root of the output file to be written
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution on screen
	 *
	 * @param nbins the number of bins to estimate the posterior
	 * distribution, used to assess its properties
	 *
	 * @param fits false \f$\rightarrow\f$ ascii file; true
	 * \f$\rightarrow\f$ fits file
	 *
	 * @param compute_mode true \f$\rightarrow\f$ compute the
	 * posterior mode; false \f$\rightarrow\f$ do not compute the
	 * posterior mode
	 *
	 * @return none
	 */
	void write_results (const std::string output_dir, const std::string root_file, const int start=0, const int thin=1, const int nbins=50, const bool fits=false, const bool compute_mode=false);

	/**
	 * @brief write the model at xx, yy computing 16th, 50th and
	 * 84th percentiles from the MCMC chains
	 *
	 * @param output_dir the output directory
	 *
	 * @param output_file the output file
	 *
	 * @param xx vector of points at which the model is computed
	 *
	 * @param yy vector of points at which the model is computed
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution on screen
	 *
	 * @return none
	 */
	void write_model_from_chain (const std::string output_dir, const std::string output_file, const std::vector<double> xx={}, const std::vector<double> yy={}, const int start=0, const int thin=1);
    };
  }
}

#endif
