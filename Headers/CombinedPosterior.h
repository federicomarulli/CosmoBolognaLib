/********************************************************************
 *  Copyright (C) 2020 by Davide Pelliciari and Sofia Contarini     *
 *  davide.pelliciari@studio.unibo.it                               *
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
 *  @file Headers/CombinedPosterior.h
 *
 *  @brief The class CombinedPosterior
 *
 *  This file defines the interface of the class CombinedPosterior
 *
 *  @authors Davide Pelliciari, Sofia Contarini, Giorgio Lesci
 *
 *  @authors davide.pelliciari@studio.unibo.it
 *
 *  @authors sofia.contarini3@unibo.it
 *
 *  @authors giorgio.lesci2@unibo.it
 */

#ifndef __COMBINEDPOSTERIOR__
#define __COMBINEDPOSTERIOR__

#include "RandomNumbers.h"
#include "Prior.h"
#include "Likelihood.h"
#include "Posterior.h"
#include "CovarianceMatrix.h"
#include "SuperSampleCovariance.h"


// ===================================================================================================


namespace cbl {

  namespace statistics {
  
  
    /**
     *  @struct STR_DependentProbes_data_model
     *  @brief the structure STR_DependentProbes_data_model
     *
     *  This structure contains the data of statistically
     *  dependent probes
     */
    struct STR_DependentProbes_data_model {
      
      /// pointer to the CovarianceMatrix object 
      std::shared_ptr<data::CovarianceMatrix> covariance;
      
      /// vector of pointers to the Model objects
      std::vector<std::shared_ptr<Model>> models;
      
      /// vector of pointers to the response functions, used for super-sample covariance
      std::vector<std::shared_ptr<Model>> responses;
      
      /// pointer to cbl::cosmology::SuperSampleCovariance object
      std::shared_ptr<cbl::cosmology::SuperSampleCovariance> Sij;
	
      /// values where the models are computed
      std::vector<std::vector<double>> xx;
	
      /// flat xx vector
      std::vector<double> flat_xx;
	
      /// flat data vector
      std::vector<double> flat_data;
	
      /// parameter indexes for the dependent probes
      std::vector<std::vector<int>> par_indexes;
      
      /// indexes of the cosmological parameters
      std::vector<int> cosmoPar_indexes;

      /**
       *  @brief default constructor
       */
      STR_DependentProbes_data_model () = default;
    };

    /**
     *  @class CombinedPosterior CombinedPosterior.h "Headers/CombinedPosterior.h"
     *
     *  @brief The class CombinedPosterior
     *
     *  This class is used to define the distribution
     */
    class CombinedPosterior : public Posterior {

    private:
    
      /// the container of parameters for statistically dependent probes
      STR_DependentProbes_data_model m_data_model;

      /// shared pointers vector containing  input Posterior objects
      std::vector<std::shared_ptr<Posterior>> m_posteriors;

      /// number of posteriors
      int m_Nposteriors;
     
      /// indexes of the parameters, for each posterior object, in the total parameter vector
      std::vector<std::vector<int>> m_parameter_indexes;
      
      /// indexes of the parameters, for each posterior object, in the total parameter vector
      std::vector<std::vector<int>> m_parameter_indexes2;
      
      /// indexes of the cosmological parameters
      std::vector<int> m_cosmoPar_indexes;
      
      /// names of the parameters, for each posterior object
      std::vector<std::vector<std::string>> m_parameter_names;

      /// importance sampling mode
      bool impsampling = false;

      /// likelihood inputs
      std::vector<std::shared_ptr<void>> m_likelihood_inputs;

      /// log-likelihood functions
      std::vector<LogLikelihood_function> m_log_likelihood_functions;

      /// likelihood functions
      std::vector<Likelihood_function> m_likelihood_functions;

      /// likelihood functions on a grid
      std::vector<LogLikelihood_function> m_likelihood_functions_grid;

      /// log-likelihood functions on a grid
      std::vector<LogLikelihood_function> m_log_likelihood_functions_grid;

      /// data containers
      std::vector<std::shared_ptr<data::Data>> m_datasets;

      /// models to test
      std::vector<std::shared_ptr<Model>> m_models;

      /// use_grid vector
      std::vector<bool> m_use_grid;


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
      CombinedPosterior () = default;

      /**
       *  @brief Constructor used to set the modelling of
       *  statistically independent probes, or used to perform
       *  importance sampling.
       *
       *  @param posteriors pointers to Posterior objects
       *
       *  @param repeated_par parameters shared by different probes, for which
       *  the user wants different posteriors for each probe. For example, if
       *  the probes A and B depend on the same parameter \f$p\f$, whose identification
       *  string is "par", then "par" must be given in input if the user desires to have
       *  different posteriors of \f$p\f$ for the two probes A and B. Every probe depending
       *  on \f$p\f$ will provide a different posterior on \f$p\f$. This is useful
       *  in the case of astrophysical parameters, e.g. the parameters describing
       *  galaxy cluster profiles.
       *
       *  @param common_repeated_par for each argument in repeated_par, a vector of vectors of
       *  integers can be defined here. Such vectors define the sets of probes for which
       *  the same priors and posteriors are provided for the same parameters. 
       *  For example, let us consider the probes {A1, A2, A3, A4}, provided in the posteriors parameter.
       *  All the probes (A1, A2, A3, A4) depend on the parameter \f$p\f$, 
       *  whose identification string is "par", and we set repeated_par = {"par"}.
       *  If we want the pair of probes {A1, A2}, {A3, A4}, to provide
       *  a different posterior on \f$p\f$, we must set common_repeated_par = { { {0,1}, {2,3} } }.
       *  This is useful when different probes in the same bin provide constraints
       *  on the same parameters. 
       *  If common_repeated_par is not provided, every probe depending on \f$p\f$ 
       *  will provide a different posterior on \f$p\f$.
       *  If two parameters are provided in repeated_par, e.g. repeated_par = {"par1", "par2"}, and
       *  only "par2" must be shared by more than one probe, then leave blank the
       *  vector of vectors corresponding to "par1" in common_repeated_par.
       *
       */
      CombinedPosterior (const std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par={}, const std::vector<std::vector<std::vector<int>>> common_repeated_par={});
      
      /**
       *  @brief Constructor used to set the modelling of
       *  statistically dependent probes. For each vector of cbl::statistics::Posterior
       *  objects in the \e posteriors argument, a cbl::data::CovarianceMatrix
       *  object must be defined. The probes in each vector within
       *  \e posteriors are described by the same likelihood function.
       *  The final log-likelihood is given by the sum of the logarithms
       *  of each likelihood describing a set of dependent probes.
       *
       *  If sets of probes are described by user-defined likelihoods,
       *  then the cbl::data::CovarianceMatrix objects must not be provided
       *  for such sets.
       *
       *  @param posteriors vector of vectors of pointers 
       *  to Posterior objects. In each vector, a set of
       *  probes is contained, with a covariance matrix
       *  defined by the corresponding cbl::data::CovarianceMatrix object
       *  given as input in the second argument of this constructor
       *
       *  @param covariance objects defining the covariance
       *  matrices for the Posterior objects given in input
       *  
       *  @param likelihood_types likelihood types describing each set of probes
       *
       *  @param repeated_par parameters shared by different probes, for which
       *  the user wants different posteriors for each probe. For example, if
       *  the probes A and B depend on the same parameter \f$p\f$, whose identification
       *  string is "par", then "par" must be given in input if the user desires to have
       *  different posteriors of \f$p\f$ for the two probes A and B. This is useful
       *  in the case of astrophysical parameters, e.g. the parameters describing
       *  galaxy cluster profiles.
       *
       *  @param common_repeated_par for each argument in repeated_par, a vector of vectors of
       *  integers can be defined here. Such vectors define the sets of probes for which
       *  the same priors and posteriors are provided for the same parameters. 
       *  For example, let us consider the two sets of probes {A1, A2, A3} and {B1, B2, B3}, provided
       *  in the posteriors parameter.
       *  All the probes (A1, A2, A3, B1, B2, B3) depend on the parameter \f$p\f$, 
       *  whose identification string is "par", and we set repeated_par = {"par"}.
       *  If we want each pair of probes {A1, B1}, {A2, B2}, {A3, B3}, to provide
       *  a different posterior on \f$p\f$, we must set common_repeated_par = { { {0,3}, {1,4}, {2,5} } }.
       *  This is useful when different probes in the same bin provide constraints
       *  on the same parameters.
       *  If common_repeated_par is not provided, every probe depending on \f$p\f$ 
       *  will provide a different posterior on \f$p\f$.
       *  If two parameters are provided in repeated_par, e.g. repeated_par = {"par1", "par2"}, and
       *  only "par2" must be shared by more than one probe, then leave blank the
       *  vector of vectors corresponding to "par1" in common_repeated_par.
       *
       *  @param SSC vector of pointers to cbl::cosmology::SuperSampleCovariance
       *  objects, for the computation of the \f$S_{ij}\f$ matrices. If, for example,
       *  the sets of probes A, B, C, are considered, and only for A and C the 
       *  super-sample covariance must be computed, then set the second element of
       *  \e SSC equal to \e NULL.
       *
       */
      CombinedPosterior (const std::vector<std::vector<std::shared_ptr<Posterior>>> posteriors, const std::vector<std::shared_ptr<data::CovarianceMatrix>> covariance, const std::vector<cbl::statistics::LikelihoodType> likelihood_types, const std::vector<std::string> repeated_par={}, const std::vector<std::vector<std::vector<int>>> common_repeated_par={}, const std::vector<std::shared_ptr<cbl::cosmology::SuperSampleCovariance>> SSC={});

      /**
       *  @brief default destructor
       *
       *
       */
      ~CombinedPosterior () = default;

      ///@}


      /**
       * @brief set the internal values of m_log_posterior as the concatenation
       * of the logposterior vectors of two cosmological chains
       *
       * @param logpostA the first logposterior distribution
       *
       * @param logpostB the second logposterior distribution
       *
       */
      void set_log_posterior(const std::vector<double> logpostA, const std::vector<double> logpostB);

      /**
       * @brief set the internal values of m_parameters as the concatenation
       * of the parameters vectors of two cosmological chains
       *
       * @param parametersA the first parameters vector
       *
       * @param parametersB the second parameter vector
       *
       */
      void set_parameters(const std::vector<std::vector<double>> parametersA, const std::vector<std::vector<double>> parametersB);

      /**
       * @brief set the internal values of m_weight as the concatenation
       * of the weights vectors of two MCMC chains
       *
       * @param weightsA the first weights vector
       *
       * @param weightsB the second weights vector
       *
       */
      void set_weight(const std::vector<double> weightsA, const std::vector<double> weightsB);
      
      /**
       * @brief check the repeated parameters
       *
       * @param dummy_Nposteriors number of posterior objects
       *
       * @param posteriors pointers to Posterior objects
       *
       * @param repeated_par repeated parameters
       *
       */
      void m_check_repeated_par (int dummy_Nposteriors, std::vector<std::shared_ptr<Posterior>> posteriors, const std::vector<std::string> repeated_par);
      
      /**
       * @brief check the common repeated parameters
       *
       * @param dummy_Nposteriors number of posterior objects
       *
       * @param posteriors vector of input Posterior objects
       *
       * @param repeated_par repeated parameters
       *
       * @param common_repeated_par indices of the common 
       * repeated parameters
       *
       */
      void m_check_common_repeated_par (int dummy_Nposteriors, std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par, std::vector<std::vector<std::vector<int>>> common_repeated_par);
      
      /**
       * @brief add a prior
       *
       * @param par_is_repeated if true, the considered
       * parameter is repeated
       *
       * @param posteriors pointer to the Posterior objects
       *
       * @param N index of the Posterior object
       *
       * @param k index of the parameter of the N-th Posterior
       *
       * @param prior_distributions vector containing all the Posterior objects
       *
       */
      void m_add_prior (bool par_is_repeated, std::vector<std::shared_ptr<Posterior>> posteriors, const int N, const int k, std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> &prior_distributions);
      
      /**
       * @brief set a common repeated parameter
       *
       * @param posteriors pointer to the Posterior objects
       *
       * @param is_in_parnames if true, the parameter is already
       * in the vector of parameter names
       *
       * @param N index of the Posterior object
       *
       * @param k index of the parameter of the N-th Posterior
       *
       * @param repeated_par repeated parameter names
       *
       * @param common_repeated_par indices of the common 
       * repeated parameters
       *
       * @param prior_distributions vector containing all the Posterior objects
       *
       * @param parameter_names vector containing all the parameter names
       *
       * @param original_names vector of original parameter names
       *
       * @param parameter_types vector containing all the parameter types
       *
       */
      void m_set_common_repeated_par (std::vector<std::shared_ptr<Posterior>> posteriors, const bool is_in_parnames, const int N, const int k, const std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par, std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> &prior_distributions, std::vector<std::string> &parameter_names, std::vector<std::string> &original_names, std::vector<ParameterType> &parameter_types);
      
      /**
       * @brief set a repeated parameter
       *
       * @param posteriors pointer to the Posterior objects
       *
       * @param is_in_parnames if true, the parameter is already
       * in the vector of parameter names
       *
       * @param N index of the Posterior object
       *
       * @param k index of the parameter of the N-th Posterior
       *
       * @param repeated_par repeated parameter names
       *
       * @param common_repeated_par indices of the common 
       * repeated parameters
       *
       * @param prior_distributions vector containing all the Posterior objects
       *
       * @param parameter_names vector containing all the parameter names
       *
       * @param original_names vector of original parameter names
       *
       * @param parameter_types vector containing all the parameter types
       *
       */
      void m_set_repeated_par (std::vector<std::shared_ptr<Posterior>> posteriors, const bool is_in_parnames, const int N, const int k, const std::vector<std::string> repeated_par, const std::vector<std::vector<std::vector<int>>> common_repeated_par, std::vector<std::shared_ptr<cbl::statistics::PriorDistribution>> &prior_distributions, std::vector<std::string> &parameter_names, std::vector<std::string> &original_names, std::vector<ParameterType> &parameter_types);
      
      /**
       * @brief set the parameters and the priors
       *
       * @param posteriors pointers to Posterior objects
       *
       *  @param repeated_par parameters shared by different probes, for which
       *  the user wants different posteriors for each probe. For example, if
       *  the probes A and B depend on the same parameter \f$p\f$, whose identification
       *  string is "par", then "par" must be given in input if the user desires to have
       *  different posteriors of \f$p\f$ for the two probes A and B. This is useful
       *  in the case of astrophysical parameters, e.g. the parameters describing
       *  galaxy cluster profiles.
       *
       *  @param common_repeated_par for each argument in repeated_par, a vector of vectors of
       *  integers can be defined here. Such vectors define the sets of probes for which
       *  the same priors and posteriors are provided for the same parameters. 
       *  For example, let us consider the two sets of probes {A1, A2, A3} and {B1, B2, B3}, provided
       *  in the posteriors parameter.
       *  All the probes (A1, A2, A3, B1, B2, B3) depend on the parameter \f$p\f$, 
       *  whose identification string is "par", and we set repeated_par = {"par"}.
       *  If we want each pair of probes {A1, B1}, {A2, B2}, {A3, B3}, to provide
       *  a different posterior on \f$p\f$, we must set common_repeated_par = { { {0,1}, {2,3}, {4,5} } }.
       *  This is useful when different probes in the same bin provide constraints
       *  on the same parameters.
       *  If common_repeated_par is not provided, every probe depending on \f$p\f$ 
       *  will provide a different posterior on \f$p\f$.
       *  If two parameters are provided in repeated_par, e.g. repeated_par = {"par1", "par2"}, and
       *  only "par2" must be shared by more than one probe, then leave blank the
       *  vector of vectors corresponding to "par1" in common_repeated_par.
       *
       */
      void m_set_parameters_priors (std::vector<std::shared_ptr<Posterior>> posteriors, std::vector<std::string> repeated_par={}, const std::vector<std::vector<std::vector<int>>> common_repeated_par={});
      
      
      /**
       * @brief set all the internal variables needed for modelling
       * independent probes
       *
       */
      void m_set_independent_probes ();
      

      /**
       * @brief do the importance sampling for two Posterior objects, which has been read externally
       *
       * @param distNum number of points for the averaging
       *
       * @param cell_size cell dimension
       *
       * @param rMAX searching radius
       *
       * @param cut_sigma for the weight's distribution cut
       *
       */
      void importance_sampling(const int distNum, const double cell_size, const double rMAX, const double cut_sigma=-1);

      /**
       * @brief do the importance sampling for two Posterior objects, that are constructed
       * from dataset and models
       *
       * @param output_path the path where the chains will be stored
       *
       * @param model_nameA the name of the first model
       *
       * @param model_nameB the name of the second model
       *
       * @param start std::vector containing initial values for
       * the posterior maximization
       *
       * @param chain_size the size of the chains that has to be ran
       *
       * @param nwalkers the number of parallal walkers for the chain
       *
       * @param burn_in the minimum chain position to be written
       *
       * @param thin the step used for dilution on screen
       *
       */
      void importance_sampling(const std::string output_path, const std::string model_nameA, const std::string model_nameB, const std::vector<double> start, const int chain_size, const int nwalkers, const int burn_in, const int thin);

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
       */
      void initialize_chains (const int chain_size, const int n_walkers);
      
      /**
       * @brief initialize the chains reading the input values from an input file
       *
       * @param chain_size the chain lenght
       *
       * @param n_walkers the number of parallel
       * chains
       *
       * @param input_dir input directory
       *
       * @param input_file input file
       *
       * @param seed the seed
       *
       */
      void initialize_chains (const int chain_size, const int n_walkers, const std::string input_dir, const std::string input_file, const int seed);

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
       *
       *
       *  @warning if parallel is set true, than pointers cannot be
       *  used inside the posterior function
       */
      void sample_stretch_move (const double aa=2, const bool parallel=true, const std::string outputFile=par::defaultString, const int start=0, const int thin=1, const int nbins=50);


      /**
       *  @brief evaluate the un-normalized posterior as the product of the
       *  un-normalized entry posteriors. For N Posterior objects, at a given point \f$ \vec{\theta_i} \f$ in the parameter space:
       *
       *  \f[ P(\vec{\theta_i} | \vec{d}) =
       *  \prod_{j=1}^{N} \mathcal{L_j}(\vec{d}|\vec{\theta_i}) \cdot Pr(\vec{\theta_i})
       *  \f]
       *
       *  where \f$P\f$ is the
       *  combined posterior,\f$\mathcal{L_j}(\vec{d}|\vec{\theta})\f$ is the
       *  j-th likelihood function and \f$Pr(\vec{\theta_i})\f$ is the prior distribution.
       *
       *  @param pp the parameters
       *
       *  @return the value of the un-normalized posterior
       */
      double operator () (std::vector<double> &pp) const;

      /**
       *  @brief evaluate the logarithm of the un-normalized
       *  posterior as the sum of all the logarithm of the un-normalized
       *  posteriors of the N entry objects. Considering a given point \f$ \vec{\theta_i} \f$ in the parameter space:
       *
       *  \f[ P(\vec{\theta_i} | \vec{d}) =
       *  \sum_{j=1}^{N} \log{\mathcal{L_j}(\vec{d}|\vec{\theta_i})} + \log{Pr(\vec{\theta_i})}
       *  \f]
       *
       *  where \f$P\f$ is the
       *  combined posterior,\f$\mathcal{L_j}(\vec{d}|\vec{\theta})\f$ is the
       *  j-th likelihood function and \f$Pr(\vec{\theta})\f$ is the prior distribution.
       *
       *  @param pp the parameters
       *
       *  @return the logarithm of the un-normalized posterior
       */
      double log (std::vector<double> &pp) const;

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
       *
       *
       */
      void write_results (const std::string output_dir, const std::string root_file, const int start=0, const int thin=1, const int nbins=50, const bool fits=false, const bool compute_mode=false, const int ns=-1, const int nb=-1);

      /**
       * @brief show the results of the MCMC combination on the screen.
       *
       * In the case of the sum of logPosteriors:
       *
       * if the covariance matrix has been estimated from a set of
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
       * In the case of importance sampling:
       *
       * The weighted average of the chains parameters are shown.
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
       *
       *
       */
      void show_results (const int start, const int thin, const int nbins=50, const bool show_mode=false, const int ns=-1, const int nb=-1);
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
       * the input file is FITS; false \f$\rightarrow\f$ the format
       * of the input file is ASCII
       *
       * @param prec decimal precision
       *
       * @param ww number of characters to be used as field width
       *
       * @warning column only work for ascii chain file
       *
       *
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
       *
       *
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
       *
       *
       */
      void write_chain_fits (const std::string output_dir, const std::string output_file, const int start=0, const int thin=1);

      /**
       * @brief write the model computing 16th, 50th and
       * 84th percentiles from the MCMC
       *
       * @param output_dir the output directory
       *
       * @param output_file the output file
       *
       * @param start the minimum chain position to be written
       *
       * @param thin the step used for dilution on screen
       *
       *
       */
      void write_model_from_chain (const std::string output_dir, const std::string output_file, const int start, const int thin);

      /**
       * @brief write maximization results on a file
       *
       * @param dir_output the output directory
       *
       * @param file the name of the output file to be written
       *
       *
       */
      void write_maximization_results (const std::string dir_output, const std::string file);

      /**
       * @brief read an entire 2d table data (MCMC chain) of unknown dimension
       *
       * @param path the path in which the file is stored
       *
       * @param filename the name of the file that has to be read
       *
       * @return two-dimensional vector containing all the read data
       */
      std::vector<std::vector<double>> read (const std::string path, const std::string filename);


    };
    
    /**
     *  @name log-likelihood functions
     */
    ///@{

    /**
     * @brief Gaussian log-likelihood
     *
     * @param likelihood_parameter the parameter values
     * set in the modelling
     *
     * @param input fixed parameters
     *
     * @return the value of the log-likelihood
     *
     */
    double LogLikelihood_Gaussian_combined (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);
    
    /**
     * @brief Poissonian log-likelihood
     *
     * @param likelihood_parameter the parameter values
     * set in the modelling
     *
     * @param input fixed parameters
     *
     * @return the value of the log-likelihood
     *
     */
    double LogLikelihood_Poissonian_combined (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);
    
    /**
     * @brief Poissonian log-likelihood with super-sample covariance
     *
     * @param likelihood_parameter the parameter values
     * set in the modelling
     *
     * @param input fixed parameters
     *
     * @return the value of the log-likelihood
     *
     */
    double LogLikelihood_Poissonian_SSC_combined (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);
    
    ///@}
        
  }
}

#endif
