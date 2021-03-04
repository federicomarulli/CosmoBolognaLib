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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
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

      /// the logarithm of the likelihood
      std::vector<double> m_log_likelihood;

      /// the MCMC acceptance rate
      std::vector<double> m_acceptance;

      /// general seed for prior/posterior distribution and sampler
      int m_seed;

      /// seed generator
      std::shared_ptr<cbl::random::UniformRandomNumbers_Int> m_seed_generator;

      /// the logarithm of the posterior
      std::vector<double> m_log_posterior;

      /// the number of parameters
      int m_Nparameters;

      /// the parameters
      std::vector<std::vector<double>> m_parameters;

      /// the chain weight
      std::vector<double> m_weight;

      /**
       * @brief set the internal attribute m_seed
       * and related quantities
       *
       * @param seed the seed
       *
       * 
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
       *  
       */
      Posterior () = default;

      /**
       *  @brief constructor
       *
       *  @param prior_distributions vector containing the priors
       *  for the likelihood parameters
       *
       *  @param likelihood pointer to an object of type Likelihood
       *
       *  @param seed general seed for prior/posterior distribution
       *  and sampler
       *
       *  
       */
      Posterior (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const Likelihood &likelihood, const int seed=5341);

      /**
       *  @brief constructor
       *
       *  @param prior_distributions vector containing the priors
       *  for the likelihood parameters
       *
       *  @param data object of type data
       *
       *  @param model object of type model
       *
       *  @param likelihood_type type of the likelihood
       *
       *  @param x_index index(s) of the extra info std::vector
       *  containing the point(s) where the model is evaluated
       *
       *  @param w_index std::vector containing the data point
       *  weight
       *
       *  @param seed general seed for prior/posterior distribution
       *  and sampler
       *
       *  
       */
      Posterior (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const int seed=5341);

      /**
       *  @brief constructor
       *
       *  @param file_name the name of the input file
       *
       *  @param path_file the path of the input file
       *
       *  @param parameter_names the vector containing the names of parameters
       *
       *  @param usecols the columns to read
       *
       *  @param skip_header the number of lines to skip
       *
       */
      Posterior (const std::string file_name, const std::string path_file, const std::vector<int> usecols,  const std::vector<std::string> parameter_names, const int skip_header);

      /**
       *  @brief default destructor
       *
       *  
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
       *  @brief the un-normalized posterior
       *
       *  \f[ P((\vec{\theta} | \vec{d}) =
       *  \mathcal{L}(\vec{d}|\vec{\theta}) \cdot Pr(\vec{\theta})
       *  \f]
       *
       *  where \f$P\f$ is the
       *  posterior,\f$\mathcal{L}(\vec{d}|\vec{\theta})\f$ is the
       *  likelihood and \f$Pr(\vec{\theta})\f$ is the prior
       *
       *  @param pp the parameters
       *
       *  @return the value of the un-normalized posterior
       */
      double operator () (std::vector<double> &pp) const;

      /**
       * @brief return the internal member m_weights
       *
       * @param start the minimum chain position to be written
       *
       * @param thin the step used for dilution
       *
       * @return vector containing the posterior weights
       */
      std::vector<double> weight (const int start=0, const int thin=1)  const;

      /**
       *  @brief the logarithm of the un-normalized posterior
       *
       *  \f[ P(\vec{\theta} | \vec{d}) =
       *  \mathcal{L}(\vec{d}|\vec{\theta}) \cdot Pr(\vec{\theta})
       *  \f]
       *
       *  where \f$P(\vec{\theta} | \vec{d})\f$ is the posterior,
       *  \f$\mathcal{L}(\vec{d}|\vec{\theta})\f$ is the likelihood,
       *  \f$Pr(\vec{\theta})\f$ is the prior, \f$\vec{\theta}\f$ is
       *  the set of model parameters and \f$\vec{d}\f$ is the
       *  dataset
       *
       *  @param pp the parameters
       *
       *  @return the logarithm of the un-normalized posterior
       */
      double log (std::vector<double> &pp) const;

      /**
       *  @brief the \f$\chi^2\f$
       *
       *  this function computes 
       *
       *  \f[ \chi^2 \equiv -2\log(\frac{P(\vec{\theta} |
       *  \vec{d})}{Pr(\vec{\theta})}) \f]
       *
       *  where \f$P(\vec{\theta} | \vec{d})\f$ is the posterior and
       *  \f$Pr(\vec{\theta})\f$ is the prior
       *
       *  @param parameter vector containing the model parameters; if
       *  not provided in input, the chi2 is evaluated at best-fit
       *  parameter values
       *
       *  @return \f$\chi^2\f$
       */
      double chi2 (const std::vector<double> parameter={}) const;

      /**
       * @brief set the model for the posterior analysis 
       *
       * @param model pointer to the model
       *
       * @param model_parameters parameters of the model
       *
       * 
       */
      void set_model (std::shared_ptr<Model> model=NULL, const std::shared_ptr<ModelParameters> model_parameters=NULL);

      /**
       * @brief set the posterior type using the LikelihoodType
       * object
       *
       * @param prior_distributions vector containing the priors for
       * the likelihood parameters
       *
       * @param data pointer to an object of type Data
       *
       * @param model pointer to an object of type model
       *
       * @param likelihood_type the likelihood type, specified with
       * the LikelihoodType object
       *
       * @param x_index index(s) of the extra info std::vector
       * containing the point(s) where to evaluate the model
       *
       * @param w_index std::vector containing the data point weight 
       *
       * @param seed the seed
       *
       * 
       */
      void set (const std::vector<std::shared_ptr<PriorDistribution>> prior_distributions, const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index, const int w_index, const int seed);

      /**
       * @brief get the number of parameters used
       *
       * @return the value of the internal variable m_Nparameters
       */
      int get_Nparameters () { return m_Nparameters; }

      /**
       * @brief get the value of the internal variable m_seed
       *
       * @return the value of the internal variable m_seed
       */
      int m_get_seed () { return m_seed; }

      /**
       * @brief get the value of the internal variable m_use_grid
       *
       * @return the boolean value of m_use_grid
       */
      bool get_m_use_grid() { return m_use_grid; }

      /**
       * @brief get the values of the internal variable m_parameters
       *
       * @return vector of vector containing parameters for each MCMC step
       */
      std::vector<std::vector<double>> get_parameters () { return m_parameters; }

      /**
       * @brief get the values of the internal variable m_log_posterior
       *
       * @return vector of double containing the logPosterior for each MCMC step
       */
      std::vector<double> get_log_posterior () { return m_log_posterior; }

      /**
       * @brief get the values of the internal variable m_likelihood_inputs
       *
       * @return pointer containing m_likelihood_inputs
       */
      std::shared_ptr<void> get_m_likelihood_inputs() { return m_likelihood_inputs; }

      /**
       * @brief get the values of the internal variable m_prior
       *
       * @return pointer containing the prior
       */
      std::shared_ptr<cbl::statistics::Prior> get_m_prior() { return m_prior; }

      /**
       * @brief get the values of the internal variable m_log_likelihood_function
       *
       * @return m_log_likelihood_function
       */
      LogLikelihood_function get_m_log_likelihood_function () { return m_log_likelihood_function; };

      /**
       * @brief get the values of the internal variable m_likelihood_function
       *
       * @return m_likelihood_function
       */
      Likelihood_function get_m_likelihood_function () { return m_likelihood_function; };

      /**
       * @brief get the values of the internal variable m_likelihood_function_grid
       *
       * @return m_likelihood_function_grid
       */
      Likelihood_function get_m_likelihood_function_grid () { return m_likelihood_function_grid; };

      /**
       * @brief get the values of the internal variable m_log_likelihood_function_grid
       *
       * @return m_log_likelihood_function_grid
       */
      LogLikelihood_function get_m_log_likelihood_function_grid () { return m_log_likelihood_function_grid; };


      /**
       *  @brief function that maximizes the posterior, finds the
       *  best-fit parameters and store them in the model
       *
       *  this function exploits the Nelder-Mead method
       *  https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
       *
       *  the algorithm defines a simplex (i.e a k-dimensional
       *  polytope which is the convex hull of its k+1 vertices) in
       *  the parameter space. At each step, it identifies the
       *  simplex vertex at which the function to be minimised
       *  (i.e. the negative posterior in this case) has the
       *  greatest value, and moves it, via reflections and scaling,
       *  to a new position in which the function has a lower
       *  value. This iteration stops when the simplex area becomes
       *  lower than the tolerance. For instance, in 2D, the
       *  starting vertices of the simplex (a triangle in 2D) are
       *  the following: (start[0], start[1]) ; (start[0]+epsilon,
       *  start[1]) ; (start[0], start[1]+epsilon)
       *
       *  @param start std::vector containing initial values for the
       *  posterior maximization
       *
       *  @param parameter_limits limits for the parameters
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the simplex side
       *
       *  
       */
      void maximize (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3)
      { (void)start; (void)parameter_limits; (void)max_iter; (void)tol; (void)epsilon; ErrorCBL("the method is used without parameter_ranges!", "maximize", "Posterior.h"); }

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
       *  @param epsilon the simplex side
       *
       *  
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
       * @param n_walkers the number of parallel
       * chains
       *
       * 
       */
      void initialize_chains (const int chain_size, const int n_walkers);

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
       * @param n_walkers the number of parallel
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
       * @param epsilon the simplex side
       *
       * 
       */
      void initialize_chains (const int chain_size, const int n_walkers, const double radius, const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3);

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
       * @param n_walkers the number of parallel
       * chains
       *
       * @param value vector containing the input values, centres of
       * the ball in the parameter space
       *
       * @param radius radius of the ball in the parameter space
       *
       * 
       */
      void initialize_chains (const int chain_size, const int n_walkers, std::vector<double> &value, const double radius);

      /**
       * @brief initialize the chains with input values
       *
       * the starting values of the chain are the elements of the
       * input matrix 'chain_values'
       *
       * @param chain_size the chain lenght
       *
       * @param chain_value matrix of size (n_walkers, n_parameters),
       * starting values of the chain
       */
      void initialize_chains (const int chain_size, const std::vector<std::vector<double>> chain_value);

      /**
       * @brief initialize the chains reading from an input file 
       *
       * the starting values of the chain are get from the last lines
       * of an input chain file; it can be used to continue an MCMC
       * sampling computation
       *
       * @param chain_size the chain lenght
       *
       * @param n_walkers the number of parallel
       * chains
       *
       * @param input_dir the input directory
       *
       * @param input_file the input file
       */
      void initialize_chains (const int chain_size, const int n_walkers, const std::string input_dir, const std::string input_file);

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
       *  during run-time. Leave it to default value to have no
       *  output.  WARNING: this option is intended for debug. It
       *  only works for the non-parallelized stretch-move
       *  algorithm.  The chain in output will be written in a
       *  different format with respect to the method
       *  cbl::statistics::Posterior::write_chain::ascii col1) chain
       *  step col2) walker index col3-npar) parameter values col
       *  npar+3) value of the posterior
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
       * @brief perform importance sampling
       *
       * Importance sampling is a convenient technique to join
       * independet datasets.
       *
       * This function takes in input a chain and computes the
       * posterior looping over all entries. It's possible to
       * specify the columns, in case the input chain has different
       * ordering, or larger number of parameters.
       *
       * It's possible either to sum or overwrite the log
       * likelihood.
       *
       * @param input_dir the input directory
       *
       * @param input_file the input file
       *
       * @param column the column of the input file to be read
       *
       * @param header_lines_to_skip the lines to be skipped in the
       * chain file
       * 
       * @param is_FITS_format true \f$\rightarrow\f$ the format of
       * the input file is FITS; false \f$\rightarrow\f$ the format
       * of the input file is ASCII
       *
       * @param apply_to_likelihood true \f$\rightarrow\f$ the
       * likelihood ratio is used as weight; false \f$\rightarrow\f$
       * the posterior ratio is used as weight
       *
       * @param n_walkers this parameter is used here only to
       * parallelize the computation
       *
       * @warning column is used only for ASCII chain files
       */
      void importance_sampling (const std::string input_dir, const std::string input_file, const std::vector<size_t> column={}, const int header_lines_to_skip=1, const bool is_FITS_format=false, const bool apply_to_likelihood=false, const int n_walkers=100);

      /**
       * @brief read the chains
       *
       * @param input_dir the input directory
       *
       * @param input_file the input file
       *
       * @param n_walkers the number of parallel chains
       *
       * @param column the columns of the input file to be read
       *
       * @param header_lines_to_skip the lines to be skipped in the
       * chain file
       *
       * @param is_FITS_format true \f$\rightarrow\f$ the format of
       * the input file is FITS; false \f$\rightarrow\f$ the format of
       * the input file is ASCII
       */
      void read_chain (const std::string input_dir, const std::string input_file, const int n_walkers, const std::vector<size_t> column={}, const int header_lines_to_skip=1, const bool is_FITS_format=false);
	
      /**
       * @brief read the chains from an ascii file
       *
       * @param input_dir the intput directory
       *
       * @param input_file the intput file
       *
       * @param n_walkers the number of parallel chains
       *
       * @param column the columns of the input file to be read.
       *
       * @param header_lines_to_skip the lines to be skipped in the
       * chain file
       */
      void read_chain_ascii (const std::string input_dir, const std::string input_file, const int n_walkers, const std::vector<size_t> column={}, const int header_lines_to_skip=1);

      /**
       * @brief read the chains from an ascii file
       *
       * @param input_dir the input directory
       *
       * @param input_file the input file
       *
       * @param n_walkers the number of parallel chains
       *
       * @param column the columns of the input file to be read
       */
      void read_chain_fits (const std::string input_dir, const std::string input_file, const int n_walkers, const std::vector<size_t> column);

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
       * @param is_FITS_format true \f$\rightarrow\f$ the format of
       * the input file is FITS; false \f$\rightarrow\f$ the format of
       * the input file is ASCII
       *
       * @param prec decimal precision
       *
       * @param ww number of characters to be used as field width
       *
       * @warning column only work for ascii chain file
       */
      void write_chain (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1, const bool is_FITS_format=false, const int prec=5, const int ww=14);
	
      /**
       * @brief write the chains obtained after the MCMC sampling on
       * an ascii file
       *
       * @param output_dir the output directory
       *
       * @param output_file the output file
       *
       * @param start the minimum chain position to be written
       *
       * @param thin the step used for dilution
       *
       * @param prec decimal precision
       *
       * @param ww number of characters to be used as field width
       */
      void write_chain_ascii (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1, const int prec=5, const int ww=14);

      /**
       * @brief write the chains obtained after the MCMC sampling on
       * a FITS file
       *
       * @param output_dir the output directory
       *
       * @param output_file the output file
       *
       * @param start the minimum chain position to be written
       *
       * @param thin the step used for dilution
       */
      void write_chain_fits (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1);

	
      /**
       * @brief write maximization results on a file
       *
       * @param output_dir the output directory 
       *
       * @param root_file the root of the output file to be written
       */
      void write_maximization_results (const std::string output_dir, const std::string root_file);

      /**
       * @brief show the results of the MCMC sampling on screen
       *
       * if the covariance matrix has been estimated from a set of
       * mock catalogues, and the input parameters ns (number of
       * samples used to estimate the covariance matrix) and nb
       * (number of data measurements, e.g. the bins of the dataset)
       * are provided (>0), then the parameter errors (\f$\sigma_p\f$)
       * will be corrected to take into account the uncertainities in
       * the covariance estimate (Percival et al. 2014):
       *
       * \f[ \sigma_p = \sqrt{\frac{1+B(n_b-n_p)}{1+A+B(n_p+1)}} \f]
       *
       * where 
       *
       * \f[ A = \frac{2}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       * \f[ B = \frac{(n_s-n_b-2)}{(n_s-n_b-1)(n_s-n_b-4)} \,. \f]
       *
       * this correction can be applied only if the likelihood is
       * Gaussian. Morever, the inverce covariance matrix estimator
       * has to be corrected to take into account the inverse
       * Wishart distribution (Hartlap, Simon and Schneider 2006).
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
       * @param ns number of samples used to estimate the covariance
       * matrix
       *
       * @param nb number of data measurements, e.g. the bins of the
       * dataset
       */
      void show_results (const int start, const int thin, const int nbins=50, const bool show_mode=false, const int ns=-1, const int nb=-1);

      /**
       * @brief store the results of the MCMC sampling to file
       * 
       * this function stores to file the posterior mean, the
       * posterior standard deviation, the posterior median, 18th
       * and 82th posterior percentiles, and, optionally, the
       * posterior mode.
       *
       * If the covariance matrix has been estimated from a set of
       * mock catalogues, and the input parameters ns (number of
       * samples used to estimate the covariance matrix) and nb
       * (number of data measurements, e.g. the bins of the dataset)
       * are provided (>0), then the parameter errors
       * (\f$\sigma_p\f$) will be corrected to take into account the
       * uncertainities in the covariance estimate (Percival et
       * al. 2014):
       *
       * \f[ \sigma_p = \sqrt{\frac{1+B(n_b-n_p)}{1+A+B(n_p+1)}} \f]
       *
       * where 
       *
       * \f[ A = \frac{2}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       * \f[ B = \frac{(n_s-n_b-2)}{(n_s-n_b-1)(n_s-n_b-4)} \,. \f]
       *
       * this correction can be applied only if the likelihood is
       * Gaussian. Morever, the inverce covariance matrix estimator
       * has to be corrected to take into account the inverse
       * Wishart distribution (Hartlap, Simon and Schneider 2006).
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
       * @param ns number of samples used to estimate the covariance
       * matrix
       *
       * @param nb number of data measurements, e.g. the bins of the
       * dataset
       */
      void write_results (const std::string output_dir, const std::string root_file, const int start=0, const int thin=1, const int nbins=50, const bool fits=false, const bool compute_mode=false, const int ns=-1, const int nb=-1);

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
       */
      void write_model_from_chain (const std::string output_dir, const std::string output_file, const std::vector<double> xx={}, const std::vector<double> yy={}, const int start=0, const int thin=1);
    };
  }
}

#endif
