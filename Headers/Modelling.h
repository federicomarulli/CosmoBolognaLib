/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Modelling.h
 *
 *  @brief The class Modelling
 *
 *  This file defines the interface of the class Modelling, used for
 *  modelling any kind of measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING__
#define __MODELLING__


#include "Posterior.h"


// ===================================================================================================


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for <B>
   *  modelling </B>
   *  
   * The \e modelling namespace contains all the functions and classes
   * used to model any kind of measurements
   */
  namespace modelling {

    /**
     *  @class Modelling Modelling.h
     *  "Headers/Modelling.h"
     *
     *  @brief The class Modelling
     *
     *  This file defines the interface of the base class Modelling,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling {

      protected:

        /// input data to be modelled
        std::shared_ptr<data::Data> m_data = NULL;

        /// check if fit range has been set
        bool m_fit_range = false;

        /// input data to be modelled
        std::shared_ptr<data::Data> m_data_fit;

        /// input model
        std::shared_ptr<statistics::Model> m_model = NULL;

        /// likelihood
        std::shared_ptr<statistics::Likelihood> m_likelihood = NULL;

        /// prior
        std::vector<std::shared_ptr<statistics::PriorDistribution>> m_parameter_priors;

        /// posterior
        std::shared_ptr<statistics::Posterior> m_posterior = NULL;

        /**
         *  @brief set the interal variable m_parameter_priors
	 *
	 *  @param prior_distribution vector containing the prior distributions
         *
	 *  @return none
         */
        void m_set_prior (std::vector<statistics::PriorDistribution> prior_distribution);

        /**
         *  @brief set the interal variable m_posterior
	 *
	 *  @param seed the seed
         *
	 *  @return none
         */
        void m_set_posterior (const int seed);

      public:

        /**
         *  @name Constructors/destructors
         */
        ///@{

        /**
         *  @brief default constuctor
         *  @return object of class Modelling
         */
        Modelling () = default;

        /**
         *  @brief default destructor
         *  @return none
         */
        virtual ~Modelling () = default;

        ///@}

        /**
         *  @name Member functions used to get the protected members
         */
        ///@{

        /**
         * @brief return the dataset
         * @return pointer to the current dataset
         */
        std::shared_ptr<data::Data> data () { return m_data; }

        /**
         * @brief return the dataset
         * @return pointer to the current dataset
         */
        std::shared_ptr<data::Data> data_fit () 
        { 
          if (!m_fit_range) 
            ErrorCBL("Error in data_fit of Modelling, no fit range has been set!");

          return m_data_fit; 
        }

        /**
         * @brief return the likelihood parameters
         * @return pointer to the likelihood parameters
         */    
        std::shared_ptr<statistics::Likelihood> likelihood ();

        /**
         * @brief return the posterior parameters
         * @return pointer to the posterior parameters
         */    
        std::shared_ptr<statistics::Posterior> posterior ();

        /**
         * @brief return the likelihood parameters
         * @return pointer to the likelihood parameters
         */    
        std::shared_ptr<statistics::ModelParameters> likelihood_parameters ();

        /**
         * @brief return the posterior parameters
         * @return pointer to the posterior parameters
         */    
        std::shared_ptr<statistics::ModelParameters> posterior_parameters ();

        ///@}


        /**
         *  @name Member functions used to set internal parameters
         */
        ///@{

        /**
         *  @brief reset the fit range 
         *
         *  set m_fit_range = false, that means that the fit range is
         *  unset
         *
         *  @return none
         */
        void reset_fit_range () { m_fit_range = false; }

        /**
         *  @brief set the fit range 
         *
         *  @param xmin minimum x value used for the fit
         *
         *  @param xmax maximum x value used for the fit
         *
         *  @return none
         */
        virtual void set_fit_range (const double xmin, const double xmax)
        { (void)xmin; (void)xmax; ErrorCBL("Error in set_fit_range of Modelling.h."); }

        /**
         *  @brief set the fit range 
         *
         *  @param xmin minimum x value used for the fit
         *
         *  @param xmax maximum x value used for the fit
         *
         *  @param ymin minimum y value used for the fit
         *
         *  @param ymax maximum y value used for the fit
         *
         *  @return none
         */
        virtual void set_fit_range (const double xmin, const double xmax, const double ymin, const double ymax)
        { (void)xmin; (void)xmax; (void)ymin; (void)ymax; ErrorCBL("Error in set_fit_range of Modelling.h."); }

        /**
         * @brief set the dataset
         *
         * @param dataset the dataset 
         *
         * @return none
         */
        void set_data (const std::shared_ptr<data::Data> dataset) { m_data = move(dataset); }

        /**
         *  @brief set the interal variable m_likelihood
         *
         * @param likelihood_type the likelihood type, specified with the 
         * LikelihoodType object
         *
         *  @param x_index index(s) of the extra info std::vector
         *  containing the point(s) where to evaluate the model
         *
         *  @param w_index index of the extra info std::vector
         *  containing the data point weight
         *
         *  @return none
         */
        void set_likelihood (const statistics::LikelihoodType likelihood_type, const std::vector<size_t> x_index={0,2}, const int w_index=-1);

        ///@}

        /**
         *  @brief function that maximizes the likelihood, finds the
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
        void maximize_likelihood (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3);

        /**
         *  @brief function that maximizes the posterior, finds the
         *  best-fit parameters and stores them in model
	 *
	 *  @param start vector containing initial values for
	 *  the posterior maximization
	 *
	 *  @param max_iter the maximum number of iterations
	 *
	 *  @param tol the tolerance to find convergence
	 *
	 *  @param epsilon the relative fraction of the interval size
	 *
	 *  @param seed the seed
	 *
	 *  @return none
         */
        void maximize_posterior (const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3, const int seed=34123);

        /**
         * @brief sample the posterior, initializing the chains by
         * drawing from the prior distributions
	 *
	 *  the starting values of the chain are extracted from the
	 *  (possibly different) distributions of the priors
	 * 
         * @param chain_size the chain lenght
         *
         * @param nwalkers the number of parallel
         * chains
	 *
         * @param seed the seed
         *
         * @param aa the parameter of the \f$g(z)\f$ distribution
         *
         * @param parallel false \f$\rightarrow\f$ non-parallel
         * sampler; true \f$\rightarrow\f$ parallel sampler
         *
         * @return none
         */
        void sample_posterior (const int chain_size, const int nwalkers, const int seed=34123, const double aa=2, const bool parallel=true);

        /**
         * @brief sample the posterior, initializing the chains in a
         * ball around the posterior best-fit parameters values
         *
	 *  the starting values of the chain are extracted from
	 *  uniform distributions in the range [parameter-radius,
	 *  parameter+radius] (for each likelihood parameter)
	 *
	 *  this function first maximizes the posterior, starting the
	 *  computation at the values of the input vector 'start',
	 *  then it inizializes the chain
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
         * @param seed the seed
         *
         * @param aa the parameter of the \f$g(z)\f$ distribution
         *
         * @param parallel false \f$\rightarrow\f$ non-parallel
         * sampler; true \f$\rightarrow\f$ parallel sampler
	 *
         * @return none
         */
        void sample_posterior (const int chain_size, const int nwalkers, const double radius, const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const int seed=34123, const double aa=2, const bool parallel=true);

        /**
         * @brief sample the posterior, initializing the chains in a
         * ball around the input parameter values
         *
	 * the starting values of the chain are extracted from uniform
	 * distributions in the range [value[i]-radius,
	 * value[i]+radius] (for each i-th likelihood parameter)
         * @param chain_size the chain lenght
	 *
         * @param nwalkers the number of parallel
         * chains
         *
         * @param value input values, center of the ball in parameter
         * space
         *
         * @param radius radius of the ball in parameter space
	 *
         * @param seed the seed
         *
         * @param aa the parameter of the \f$g(z)\f$ distribution
         *
         * @param parallel false \f$\rightarrow\f$ non-parallel
         * sampler; true \f$\rightarrow\f$ parallel sampler
         *
         * @return none
         */
        void sample_posterior (const int chain_size, const int nwalkers, std::vector<double> &value, const double radius, const int seed=34123, const double aa=2, const bool parallel=true);

        /**
         * @brief sample the posterior, initializing the chains with
         * input values
         *	 
	 * the starting values of the chain are the elements of the
	 * input matrix 'chain_values'
	 *
	 * @param chain_size the chain lenght
         * @param chain_value std::vector of size (nwalkers,
         * nparameters) starting values of the chain
	 *
         * @param seed the seed
         *
         * @param aa the parameter of the \f$g(z)\f$ distribution
         *
         * @param parallel false \f$\rightarrow\f$ non-parallel
         * sampler; true \f$\rightarrow\f$ parallel sampler
         *
         * @return none
         */
        void sample_posterior (const int chain_size, const std::vector<std::vector<double>> chain_value, const int seed=34123, const double aa=2, const bool parallel=true);

        /**
         * @brief sample the posterior, initializing the chains
         * reading the input values from an input file
	 *	 
	 * the starting values of the chain are get from the last
	 * lines of an input chain file; it can be used to continue an
	 * MCMC sampling computation
	 *
	 * @param chain_size the chain lenght
         * @param nwalkers the number of parallel
         * chains
         *
         * @param input_dir input directory
         *
         * @param input_file the input file
	 *
         * @param seed the seed
         *
         * @param aa the parameter of the \f$g(z)\f$ distribution
         *
         * @param parallel false \f$\rightarrow\f$ non-parallel
         * sampler; true \f$\rightarrow\f$ parallel sampler
         *
         * @return none
         */
        void sample_posterior (const int chain_size, const int nwalkers, const std::string input_dir, const std::string input_file, const int seed=34123, const double aa=2, const bool parallel=true);

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
         * @brief show the results of the MCMC sampling on screen
         *
         * @param start the minimum chain position to be written
         *
         * @param thin the step used for dilution on screen
         *
         * @param nbins the number of bins
         *
	 * @param show_mode true \f$\rightarrow\f$ show the posterior
	 * mode; false \f$\rightarrow\f$ do not show the posterior
	 * mode
	 *
         * @return none
         */
        void show_results (const int start=0, const int thin=1, const int nbins=50, const bool show_mode=false);



	/**
	 * @brief write the results of the MCMC sampling to file
	 * 
	 * this function stores to file the posterior mean, standard
	 * deviation, median, 18th and 82th percentiles, and
	 * optionally the mode
	 *
	 * @param output_dir the output director
	 *
	 * @param root_file the root of the output files: 
	 * - file_root_parameters.dat file containing the output of
	 * the MCMC sampling for each parameter
	 * - file_root_covariance.dat file containing the covariance
	 * of the parameters
	 * - file_root_chain file containing the chains: the extention
	 * can be .dat or .fits
	 *
	 * @param start the minimum chain position to be written
	 *
	 * @param thin the step used for dilution on screen
	 *
	 * @param nbins the number of bins
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
         *  @brief write the model at xx
         *  for given parameters
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  @param parameters vector containing the input parameters
         *  used to compute the model; if this vector is not provided,
         *  the model will be computed using the best-fit parameters
         *
         *  @return none
         */
        virtual void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters)
        { (void)output_dir; (void)output_file; (void)xx; (void)parameters; cbl::ErrorCBL("Error in write_model() of Modelling.h!"); }

        /**
         *  @brief write the model at xx, yy
         *  for given parameters
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  first axis
         *  @param yy vector of points at which the model is computed,
         *  second axis
         *  @param parameters vector containing the input parameters
         *  used to compute the model; if this vector is not provided,
         *  the model will be computed using the best-fit parameters
         *
         *  @return none
         */
        virtual void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> parameters)
        { (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)parameters; cbl::ErrorCBL("Error in write_model() of Modelling.h!"); }

        /**
         *  @brief write the model at xx 
         *  with best-fit parameters obtained from likelihood
         *  maximization
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *
         *  @return none
         */
        virtual void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx)
        { (void)output_dir; (void)output_file; (void)xx; cbl::ErrorCBL("Error in write_model_at_bestfit() of Modelling.h!"); }

        /**
         *  @brief write the model at xx, yy
         *  with best-fit parameters obtained from likelihood
         *  maximization       
         *  
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  first axis
         *  @param yy vector of points at which the model is computed,
         *  second axis
         *
         *  @return none
         */
        virtual void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy)
        { (void)output_dir; (void)output_file; (void)xx; (void)yy; cbl::ErrorCBL("Error in write_model_at_bestfit of Modelling.h!"); }

        /**
         *  @brief write the model at xx 
         *  computing 16th, 50th and 84th percentiles
         *  from the chains.
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  @param start the starting position for each chain
         *  @param thin the position step
         *
         *  @return none
         */
        virtual void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start=0, const int thin=1)
        { (void)output_dir; (void)output_file; (void)xx; (void)start; (void)thin; cbl::ErrorCBL("Error in write_model_from_chains of Modelling.h!"); }

        /**
         *  @brief write the model at xx, yy
         *  computing 16th, 50th and 84th percentiles
         *  from the chains.
         *
         *  @param output_dir the output directory
         *  @param output_file the output file
         *  @param xx vector of points at which the model is computed,
         *  first axis
         *  @param yy vector of points at which the model is computed,
         *  second axis
         *  @param start the starting position for each chain
         *  @param thin the position step
         *
         *  @return none
         */
        virtual void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const int start=0, const int thin=1)
        { (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)start; (void)thin; cbl::ErrorCBL("Error in write_model_from_chains of Modelling.h!"); }

        ///@}

  };
}
}

#endif
