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
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
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
   *  The \e modelling namespace contains all the functions and
   *  classes used to model any kind of measurements
   */
  namespace modelling {

    /**
     *  @class Modelling Modelling.h "Headers/Modelling.h"
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

      /// input data restricted to the range used for the fit
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
       *  
       */
      void m_set_prior (std::vector<statistics::PriorDistribution> prior_distribution);

      /**
       *  @brief set the interal variable m_posterior
       *
       *  @param seed the seed
       *
       *  
       */
      void m_set_posterior (const int seed);

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  
       */
      Modelling () = default;

      /**
       *  @brief default destructor
       *  
       */
      virtual ~Modelling () = default;

      ///@}

      /**
       *  @name Member functions used to get the protected members
       */
      ///@{

      /**
       *  @brief return the dataset
       *  @return pointer to the current dataset
       */
      std::shared_ptr<data::Data> data () { return m_data; }

      /**
       *  @brief return the dataset
       *  @return pointer to the current dataset
       */
      std::shared_ptr<data::Data> data_fit () 
      { 
	if (!m_fit_range) 
	  ErrorCBL("no fit range has been set!", "data_fit", "Modelling.h");

	return m_data_fit; 
      }

      /**
       *  @brief return the likelihood parameters
       *  @return pointer to the likelihood parameters
       */    
      std::shared_ptr<statistics::Likelihood> likelihood ();

      /**
       *  @brief return the posterior parameters
       *  @return pointer to the posterior parameters
       */    
      std::shared_ptr<statistics::Posterior> posterior ();

      /**
       *  @brief return the likelihood parameters
       *  @return pointer to the likelihood parameters
       */    
      std::shared_ptr<statistics::ModelParameters> likelihood_parameters ();

      /**
       *  @brief return the posterior parameters
       *  @return pointer to the posterior parameters
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
       */
      void reset_fit_range () { m_fit_range = false; }

      /**
       *  @brief set the fit range 
       *
       *  @param xmin minimum x value used for the fit
       *
       *  @param xmax maximum x value used for the fit
       * 
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_fit_range (const double xmin, const double xmax)
      { (void)xmin; (void)xmax; ErrorCBL("", "set_fit_range", "Modelling.h"); }
      
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
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void set_fit_range (const double xmin, const double xmax, const double ymin, const double ymax)
      { (void)xmin; (void)xmax; (void)ymin; (void)ymax; ErrorCBL("", "set_fit_range", "Modelling.h"); }

      /**
       *  @brief set the dataset
       *
       *  @param dataset the dataset 
       */
      void set_data (const std::shared_ptr<data::Data> dataset) { m_data = move(dataset); }

      /**
       *  @brief set the likelihood function
       *
       *  if the input parameter Nres is larger than 0, the inverted
       *  covariance matrix is corrected as follows (Hartlap, Simon
       *  and Schneider 2006):
       *
       *  \f[ \hat{C}^{-1} =
       *  \left(1-\frac{n_b+1}{N_{res}-1}\right)C^{-1} \f]
       *
       *  where \f$n_b\f$ is the number of bins and \f$N_{res}\f$ is
       *  the number of resamplings
       *
       *  @param likelihood_type the likelihood type, specified with
       *  the LikelihoodType object
       *
       *  @param x_index index(s) of the extra info std::vector
       *  containing the point(s) where to evaluate the model
       *
       *  @param w_index index of the extra info std::vector
       *  containing the data point weight
       *
       *  @param prec the precision required in the inversion of the
       *  covariance matrix
       *
       *  @param Nres \f$N_{res}\f$, the number of catalogue
       *  resamplings used to estimate the covariance matrix;
       *  \f$N_{res}=-1\f$ if the covariance matrix has not been
       *  estimated with resampling methods
       */
      void set_likelihood (const statistics::LikelihoodType likelihood_type, const std::vector<size_t> x_index={0, 2}, const int w_index=-1, const double prec=1.e-10, const int Nres=-1);

      ///@}
      

      /**
       *  @name Member functions used to manage likelihood/posterior
       *  distributions
       */
      ///@{
      
      /**
       *  @brief function that maximizes the posterior, finds the
       *  best-fit parameters and stores them in the model
       *
       *  this function exploits the Nelder-Mead method
       *  https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
       *
       *  the algorithm defines a simplex (i.e a k-dimensional
       *  polytope which is the convex hull of its k+1 vertices) in
       *  the parameter space. At each step, it identifies the
       *  simplex vertex at which the function to be minimised
       *  (i.e. the negative likelihood in this case) has the
       *  greatest value, and moves it, via reflections and scaling,
       *  to a new position in which the function has a lower
       *  value. This iteration stops when the simplex area becomes
       *  lower than the tolerance. For instance, in 2D, the
       *  starting vertices of the simplex (a triangle in 2D) are
       *  the following: (start[0], start[1]) ; (start[0]+epsilon,
       *  start[1]) ; (start[0], start[1]+epsilon)
       *
       *  @param start std::vector containing the initial values for
       *  the likelihood maximization
       *  
       *  @param parameter_limits limits of the parameter space in
       *  where the maximum of likelihood will be searched
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the simplex side
       *
       *  
       */
      void maximize_likelihood (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3);

      /**
       *  @brief function that maximizes the posterior, finds the
       *  best-fit parameters and stores them in the model
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
       *  @param start vector containing initial values for
       *  the posterior maximization
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance to find convergence
       *
       *  @param epsilon the relative fraction of the initial
       *  simplex size
       *
       *  @param seed the seed
       *
       *  
       */
      void maximize_posterior (const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3, const int seed=666);

      /**
       *  @brief sample the posterior, initializing the chains by
       *  drawing from the prior distributions
       *
       *  the starting values of the chain are extracted from the
       *  (possibly different) distributions of the priors
       * 
       *  @param chain_size the chain lenght
       *
       *  @param nwalkers the number of parallel chains
       *
       *  @param seed the seed
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       *
       *  
       */
      void sample_posterior (const int chain_size, const int nwalkers, const int seed=666, const double aa=2, const bool parallel=true);

      /**
       *  @brief sample the posterior, initializing the chains in a
       *  ball around the posterior best-fit parameters values
       *
       *  the starting values of the chain are extracted from
       *  uniform distributions in the range [parameter-radius,
       *  parameter+radius] (for each likelihood parameter)
       *
       *  this function first maximizes the posterior, starting the
       *  computation at the values of the input vector 'start',
       *  then it inizializes the chain
       *
       *  @param chain_size the chain lenght
       *
       *  @param nwalkers the number of parallel chains
       *
       *  @param radius radius of the ball in parameter space
       *
       *  @param start std::vector containing initial values for the
       *  posterior maximization
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the simplex side
       *
       *  @param seed the seed
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       *
       *  
       */
      void sample_posterior (const int chain_size, const int nwalkers, const double radius, const std::vector<double> start, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3, const int seed=666, const double aa=2, const bool parallel=true);

      /**
       *  @brief sample the posterior, initializing the chains by
       *  drawing from the prior distributions
       *
       *  the starting values of the chain are extracted from the
       *  (possibly different) distributions of the priors
       *
       *  @param chain_size the chain lenght
       *
       *  @param nwalkers the number of parallel chains
       *
       *  @param value input values, center of the ball in parameter
       *  space
       *
       *  @param radius radius of the ball in parameter space
       *
       *  @param seed the seed
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       *
       *  
       */
      void sample_posterior (const int chain_size, const int nwalkers, std::vector<double> &value, const double radius, const int seed=666, const double aa=2, const bool parallel=true);

      /**
       *  @brief sample the posterior, initializing the chains with
       *  input values
       *	 
       *  the starting values of the chain are the elements of the
       *  input matrix 'chain_values'
       *
       *  @param chain_size the chain lenght
       *  
       *  @param chain_value std::vector of size (nwalkers,
       *  nparameters) starting values of the chain
       *
       *  @param seed the seed
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       *
       *  
       */
      void sample_posterior (const int chain_size, const std::vector<std::vector<double>> chain_value, const int seed=666, const double aa=2, const bool parallel=true);

      /**
       *  @brief sample the posterior, initializing the chains
       *  reading the input values from an input file
       *	 
       *  the starting values of the chain are get from the last
       *  lines of an input chain file; it can be used to continue
       *  an MCMC sampling computation
       *
       *  @param chain_size the chain lenght
       *
       *  @param nwalkers the number of parallel chains
       *
       *  @param input_dir input directory
       *
       *  @param input_file the input file
       *
       *  @param seed the seed
       *
       *  @param aa the parameter of the \f$g(z)\f$ distribution
       *
       *  @param parallel false \f$\rightarrow\f$ non-parallel
       *  sampler; true \f$\rightarrow\f$ parallel sampler
       *
       *  
       */
      void sample_posterior (const int chain_size, const int nwalkers, const std::string input_dir, const std::string input_file, const int seed=666, const double aa=2, const bool parallel=true);

      /**
       * @brief perform importance sampling
       *
       * Importance sampling is a convenient technique to join
       * independet datasets.
       *
       * This function takes in input a chain and computes the
       * likelihood or posterior looping over all entries. It's
       * possible to specify the columns, in case the input chain has
       * different ordering, or larger number of parameters.
       *
       * @param input_dir input directory
       *
       * @param input_file the input file
       *
       * @param seed the seed
       *
       * @param column the columns of the input file to be read
       *
       * @param header_lines_to_skip the lines to be skipped in the
       * chain file
       *
       * @param is_FITS_format true \f$\rightarrow\f$ the format of
       * the input file is FITS; false \f$\rightarrow\f$ the format of
       * the input file is ASCII
       *
       * @param apply_to_likelihood true \f$\rightarrow\f$ the
       * likelihood ratio is used as weight; false \f$\rightarrow\f$
       * the posterior ratio is used as weight
       *
       * @warning column is used only for ASCII chain files
       *
       * 
       */
      void importance_sampling (const std::string input_dir, const std::string input_file, const int seed=666, const std::vector<size_t> column={}, const int header_lines_to_skip=1, const bool is_FITS_format=false, const bool apply_to_likelihood=false);

      ///@}

      
      /**
       *  @name Member functions used for Input/Output
       */
      ///@{
      
      /**
       *  @brief write the chains obtained after the MCMC sampling
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       *
       *  @param start the minimum chain position to be written
       *
       *  @param thin the step used for dilution
       *
       *  @param is_FITS_format true \f$\rightarrow\f$ the format of
       *  the input file is FITS; false \f$\rightarrow\f$ the format
       *  of the input file is ASCII
       *
       *  @param prec decimal precision
       *
       *  @param ww number of characters to be used as field width
       *
       *  
       */
      void write_chain (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1, const bool is_FITS_format=false, const int prec=5, const int ww=14);

      /**
       *  @brief read the chains
       *
       *  @param input_dir the input directory
       *
       *  @param input_file the input file
       *
       *  @param nwalkers the number of parallel chains
       *
       * @param columns the columns of the input file to be read.
       *
       *  @param skip_header the lines to be skipped in the chain
       *  file
       *
       *  @param fits false \f$\rightarrow\f$ ascii file; true
       *  \f$\rightarrow\f$ fits file
       *
       *  
       */
      void read_chain (const std::string input_dir, const std::string input_file, const int nwalkers, const std::vector<size_t> columns={}, const int skip_header=1, const bool fits=false);

      /**
       *  @brief show the results of the MCMC sampling on screen
       *
       *  if the covariance matrix has been estimated from a set of
       *  mock catalogues, and the input parameters ns (\f$n_s\f$, the
       *  number of samples used to estimate the covariance matrix) is
       *  provided (>0), then the parameter errors (\f$\sigma_p\f$)
       *  will be corrected to take into account the uncertainities in
       *  the covariance estimate (Percival et al. 2014):
       *
       *  \f[ \sigma_p = \sqrt{\frac{1+B(n_b-n_p)}{1+A+B(n_p+1)}} \f]
       *
       *  where 
       *
       *  \f[ A = \frac{2}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       *  \f[ B = \frac{(n_s-n_b-2)}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       *  where \f$n_b\f$ is number of data measurements (e.g. the
       *  bins of the dataset).
       *
       *  This correction can be applied only if the likelihood is
       *  Gaussian. Morever, the inverce covariance matrix estimator
       *  has to be corrected to take into account the inverse Wishart
       *  distribution (Hartlap, Simon and Schneider 2006).
       *
       *  @param start the minimum chain position to be written
       *
       *  @param thin the step used for dilution on screen
       *
       *  @param nbins the number of bins
       *
       *  @param show_mode true \f$\rightarrow\f$ show the posterior
       *  mode; false \f$\rightarrow\f$ do not show the posterior mode
       *
       *  @param ns number of samples used to estimate the covariance
       *  matrix
       *
       *  
       */
      void show_results (const int start=0, const int thin=1, const int nbins=50, const bool show_mode=false, const int ns=-1);

      /**
       *  @brief write the results of the MCMC sampling to file
       *  
       *  this function stores to file the posterior mean, the
       *  posterior standard deviation, the posterior median, 18th and
       *  82th posterior percentiles, and, optionally, the posterior
       *  mode.
       *
       *  If the covariance matrix has been estimated from a set of
       *  mock catalogues, and the input parameters ns (\f$n_s\f$, the
       *  number of samples used to estimate the covariance matrix) is
       *  provided (>0), then the parameter errors (\f$\sigma_p\f$)
       *  will be corrected to take into account the uncertainities in
       *  the covariance estimate (Percival et al. 2014):
       *
       *  \f[ \sigma_p = \sqrt{\frac{1+B(n_b-n_p)}{1+A+B(n_p+1)}} \f]
       *
       *  where 
       *
       *  \f[ A = \frac{2}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       *  \f[ B = \frac{(n_s-n_b-2)}{(n_s-n_b-1)(n_s-n_b-4)} \,, \f]
       *
       *  where \f$n_b\f$ is number of data measurements (e.g. the
       *  bins of the dataset).
       *
       *  This correction can be applied only
       *  if the likelihood is Gaussian. Morever, the inverce
       *  covariance matrix estimator has to be corrected to take into
       *  account the inverse Wishart distribution (Hartlap, Simon and
       *  Schneider 2006).
       *
       *  @param output_dir the output director
       *
       *  @param root_file the root of the output files: -
       *  file_root_parameters.dat file containing the output of the
       *  MCMC sampling for each parameter -
       *  file_root_covariance.dat file containing the covariance of
       *  the parameters - file_root_chain file containing the
       *  chains: the extention can be .dat or .fits
       *
       *  @param start the minimum chain position to be written
       *
       *  @param thin the step used for dilution on screen
       *
       *  @param nbins the number of bins
       *
       *  @param fits false \f$\rightarrow\f$ ascii file; true
       *  \f$\rightarrow\f$ fits file
       *
       *  @param compute_mode true \f$\rightarrow\f$ compute the
       *  posterior mode; false \f$\rightarrow\f$ do not compute the
       *  posterior mode
       *
       *  @param ns number of samples used to estimate the covariance
       *  matrix
       *
       *  
       */
      void write_results (const std::string output_dir, const std::string root_file, const int start=0, const int thin=1, const int nbins=50, const bool fits=false, const bool compute_mode=false, const int ns=-1);

      /**
       *  @brief write the model at xx for given parameters
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       *
       *  @param xx vector of points at which the model is computed,
       *  
       *  @param parameters vector containing the input parameters
       *  used to compute the model; if this vector is not provided,
       *  the model will be computed using the best-fit parameters
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters)
      { (void)output_dir; (void)output_file; (void)xx; (void)parameters; cbl::ErrorCBL("", "write_model", "Modelling.h"); }

      /**
       *  @brief write the model at xx, yy for given parameters
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *
       *  @param yy vector of points at which the model is computed,
       *  second axis
       *
       *  @param parameters vector containing the input parameters
       *  used to compute the model; if this vector is not provided,
       *  the model will be computed using the best-fit parameters
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> parameters)
      { (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)parameters; cbl::ErrorCBL("", "write_model", "Modelling.h"); }

      /**
       *  @brief write the model at xx with best-fit parameters
       *  obtained from posterior maximization
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       *
       *  @param xx vector of points at which the model is computed
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx)
      { (void)output_dir; (void)output_file; (void)xx; cbl::ErrorCBL("", "write_model_at_bestfit", "Modelling.h"); }

      /**
       *  @brief write the model at xx, yy with best-fit parameters
       *  obtained from likelihood maximization
       *  
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *
       *  @param yy vector of points at which the model is computed,
       *  second axis
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy)
      { (void)output_dir; (void)output_file; (void)xx; (void)yy; cbl::ErrorCBL("", "write_model_at_bestfit", "Modelling.h"); }

      /**
       *  @brief write the model at xx computing 16th, 50th and 84th
       *  percentiles from the chains
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file 
       *
       *  @param xx vector of points at which the model is computed
       *
       *  @param start the starting position for each chain 
       *
       *  @param thin the position step
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start=0, const int thin=1)
      { (void)output_dir; (void)output_file; (void)xx; (void)start; (void)thin; cbl::ErrorCBL("", "write_model_from_chains", "Modelling.h"); }

      /**
       *  @brief write the model at xx, yy computing 16th, 50th and
       *  84th percentiles from the chains
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *
       *  @param yy vector of points at which the model is computed,
       *  second axis
       *
       *  @param start the starting position for each chain 
       *
       *  @param thin the position step
       *
       *  @return none, or an error message if the derived object does
       *  not have this member
       */
      virtual void write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const int start=0, const int thin=1)
      { (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)start; (void)thin; cbl::ErrorCBL("", "write_model_from_chains", "Modelling.h"); }

      ///@}
	
    };
  }
}

#endif
