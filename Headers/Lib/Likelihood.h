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
 *  @file Headers/Lib/Likelihood.h
 *
 *  @brief The class Likelihood 
 *
 *  This file defines the interface of the class Likelihood, used for
 *  statistical analyses and Bayesian inference
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LIKEL__
#define __LIKEL__

#include "LikelihoodFunction.h"


// ============================================================================================


namespace cosmobl {

  namespace statistics {
    
    /**
     *  @class Likelihood Likelihood.h "Headers/Lib/Likelihood.h"
     *
     *  @brief The class Likelihood
     *
     *  This class is used to handle objects of type likelihood. It is
     *  used for all kind of likelihood analyses, sample and
     *  maximization
     */
    class Likelihood
    {
      
    protected:

      /// data containers 
      shared_ptr<data::Data> m_data;

      /// model to test
      shared_ptr<Model> m_model; 

      /// likelihood parameters
      shared_ptr<void> m_likelihood_parameters;

      /// likelihood parameters
      shared_ptr<LikelihoodParameters> m_parameters;

      /// type of the likelihood
      LikelihoodType m_likelihood_type = LikelihoodType::_Likelihood_NotSet_;

      /// likelihood function
      LogLikelihood_function m_log_likelihood_function;

      /// value of the likelihood
      vector<double> m_likelihood_values;

      /// acceptance values
      vector<double> m_acceptance;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Likelihood
       */
      Likelihood () : m_likelihood_type(LikelihoodType::_Likelihood_NotSet_) {}

      /**
       *  @brief constructor
       *  
       *  @param data pointers to the data container
       *  
       *  @param model pointers to the model 
       *
       *  @param parameters vector containing likelihood parameters 
       *
       *  @param likelihood_type type of likelihood
       *
       *  @return object of class Likelihood
       */
      Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const vector<shared_ptr<Parameter>> parameters, const LikelihoodType likelihood_type);

      /**
       *  @brief constructor
       *
       *  @param data pointers to the data container
       *
       *  @param model pointers to the model 
       *
       *  @param parameters vector containing likelihood parameters 
       *
       *  @param loglikelihood_function function of type loglikelihood_function
       *
       *  @return object of class Likelihood
       */
      Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const vector<shared_ptr<Parameter>> parameters, const LogLikelihood_function loglikelihood_function);

      /**
       *  @brief constructor
       *  
       *  @param data pointers to the data container
       *  
       *  @param model pointers to the model 
       *
       *  @param parameters likelihood parameters 
       *
       *  @param likelihood_type type of likelihood
       *
       *  @return object of class Likelihood
       */
      Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const shared_ptr<LikelihoodParameters> parameters, const LikelihoodType likelihood_type);

      /**
       *  @brief constructor
       *
       *  @param data pointers to the data container
       *
       *  @param model pointers to the model 
       *
       *  @param parameters likelihood parameters 
       *
       *  @param loglikelihood_function function of type loglikelihood_function
       *
       *  @return object of class Likelihood
       */
      Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const shared_ptr<LikelihoodParameters> parameters, const LogLikelihood_function loglikelihood_function);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Likelihood () = default;

      ///@}


      /**
       * @brief evaluate the prior
       *
       * @param pp the likelihood parameters
       *
       * @return none
       */
      double prior (vector<double> &pp) const;

      /**
       * @brief evaluate the prior
       *
       * @param pp the likelihood parameters
       *
       * @return none
       */
      double log_prior (vector<double> &pp) const;

      /**
       * @brief evaluate the likelihood
       *
       * @param pp the likelihood parameters
       *
       * @return the likelihood \f$ \mathcal{L} \f$
       */
      double likelihood (vector<double> &pp) const;

      /**
       * @brief evaluate the log-likelihood
       *
       * @param pp the likelihood parameters
       *
       * @return the log-likelihood \f$ \log(\mathcal{L}) \f$
       */
      double log_likelihood (vector<double> &pp) const;

      /**
       * @brief evaluate the likelihood
       *
       * @param pp the likelihood parameters
       *
       * @return  \f$ \mathcal{L} \times P(\vec{\theta}) \f$
       */
      double likelihood_and_priors (vector<double> &pp) const;

      /**
       * @brief evaluate the likelihood times the priors
       *
       * @param pp the likelihood parameters
       *
       * @return \f$ \log(\mathcal{L}) + log(P(\vec{\theta})) \f$
       */
      double log_likelihood_and_priors (vector<double> &pp) const;

      /**
       * @brief set the data for the likelihood analysis 
       *
       * @param data pointer to the dataset
       *
       * @return none
       */
      void set_data (shared_ptr<data::Data> data);

      /**
       * @brief set the model for the likelihood analysis 
       *
       * @param model pointer to the model
       *
       * @return none
       */
      void set_model (shared_ptr<Model> model);

      /**
       * @brief set the likelihood parameters  
       *
       * @param parameters vector containing pointers
       * to the likelihood parameters
       *
       * @return none
       */
      void set_parameters (vector<shared_ptr<Parameter>> parameters);

      /**
       * @brief set the likelihood type using the LikelihoodType object 
       *
       * @param likelihood_type the likelihood type, specified with the 
       * LikelihoodType object
       *
       * @return none
       */
      void set_likelihood_function (const LikelihoodType likelihood_type);

      /**
       * @brief set the likelihood function 
       *
       * @param loglikelihood_function the likelihood function
       *
       * @return none
       */
      void set_likelihood_function (const LogLikelihood_function loglikelihood_function);

      /**
       *  @brief function that maximize the likelihood, find the
       *  best-fit parameters and store them in model
       *
       *  @param guess vector containing initial guess values and
       *  result of the likelihood maximization
       *
       *  @param usePriors false &rarr; minimize \f$
       *  -\log(\mathcal{L}) \f$; true &rarr; minimize \f$
       *  -\log(\mathcal{L})- \log(\mathcal{L}) \f$:
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void maximize (vector<double> &guess, const bool usePriors=false, const unsigned int max_iter=100, const double tol=1.e-6, const double epsilon=1.e-3); 

      /**
       *  @brief function that maximize the likelihood, find the
       *  best-fit parameters and store them in model
       *
       *  @param guess vector containing the results of the likelihood
       *  maximization
       *
       *  @param ntry number of prior extractions to find the initial guess
       *
       *  @param prior_seed base prior seed
       *
       *  @param usePriors false &rarr; minimize \f$ -\log(\mathcal{L}) \f$; 
       *  true &rarr; minimize \f$ -\log(\mathcal{L})- \log(\mathcal{L}) \f$:
       *
       *  @param max_iter maximum number of iterations 
       *
       *  @param tol the tolerance in finding convergence
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void maximize (vector<double> &guess, const int ntry, const int prior_seed=431412, const bool usePriors=false, const unsigned int max_iter=100, const double tol=1.e-6, const double epsilon=1.e-3);

      /**
       *  @brief function that 
       *  samples likelihood, using stretch-move
       *  algorithm on n-dimensional parameter space, and stores chain
       *  parameters.
       *
       *  @param seed the seed for random number generator
       *
       *  @param aa parameter for the stretch-move distribution
       *
       *  @return averace acceptance ratio
       */
      void sample_stretch_move (const int seed, const double aa=2);

      /**
       *  @brief parallel version of the function that 
       *  samples likelihood, using stretch-move
       *  algorithm on n-dimensional parameter space, and stores chain
       *  parameters.
       *
       *  @param seed the seed for random number generator
       *
       *  @param aa parameter for the stretch-move distribution
       *
       *  @return averace acceptance ratio
       */
      void sample_stretch_move_parallel (const int seed, const double aa=2);

      /**
       *  @brief function that samples likelihood, using stretch-move
       *  algorithm on n-dimensional parameter space, and stores chain
       *  parameters.
       *  @param nstep_p1 number of step for the first parameter
       *  @param nstep_p2 number of step for the second parameter
       *  @param interpolation_method the interpolation method
       *  @param nchains number of chains to sample the parameter space 
       *  @param chain_size number of step in each chain 
       *  @param seed the seed for random number generator
       *  @param do_write_chain 0 &rarr; don't write chains at each step 
       *  1 &rarr; write chains at each step
       *  @param output_dir directory of output for chains 
       *  @param output_file file of output for chains 
       *  @return averace acceptance ratio
       */
      void sample_tabulated_likelihood (const int nstep_p1, const int nstep_p2, const string interpolation_method, const int nchains, const int chain_size, const int seed, bool do_write_chain = 0, const string output_dir=par::defaultString, const string output_file=par::defaultString);

      /**
       *  @brief function that write sampled parameters 
       *
       *  @param output_dir directory of output for chains 
       *
       *  @param output_file file of output for chains 
       *
       *  @param start starting position in the chains
       *
       *  @param thin take one chain step every \e thin steps
       *
       *  @return none
       */
      void write_chain (const string output_dir, const string output_file, const int start, const int thin=1);
      
      /**
       *  @brief function that write sampled parameters 
       *  on a FITS file
       *
       *  @param output_dir directory of output for chains 
       *
       *  @param output_file file of output for chains 
       *
       *  @param table_name table name
       *
       *  @param burn_in the number of the first chain steps to be
       *  discarded
       *
       *  @param thin take one chain step every \e thin steps
       *
       *  @return none
       */
      void write_chain_FITS (const string output_dir, const string output_file, const string table_name, const int burn_in=0, const int thin=1)
      { (void)output_dir; (void)output_file; (void)table_name; (void)burn_in; (void)thin; ErrorCBL("Error in write_chain_FITS! Work in progress...", glob::_workInProgress_); }

      /**
       *  @brief function that write sampled parameters 
       *
       *  @param output_dir directory of output for chains 
       *
       *  @param output_file file of output for chains 
       *
       *  @param nwalkers the number of parallel walkers
       *  
       *  @param skip_header header lines
       *
       *  @return none
       */
      void read_chain (const string output_dir, const string output_file, const int nwalkers, const int skip_header);
      
      /**
       *  @brief function that read sampled parameters 
       *  from a FITS file
       *
       *  @param output_dir directory of output for chains 
       *
       *  @param output_file file of output for chains 
       *
       *  @param table_name table name
       *
       *  @param burn_in the number of the first chain steps to be
       *  discarded
       *
       *  @param thin take one chain step every \e thin steps
       *
       *  @return none
       */
      void read_chain_FITS (const string output_dir, const string output_file, const string table_name, const int burn_in=0, const int thin=1)
      { (void)output_dir; (void)output_file; (void)table_name; (void)burn_in; (void)thin; ErrorCBL("Error in read_chain_FITS() of Likelihood.h! Work in progress...", glob::_workInProgress_); }

      /**
       *  @brief initialize the chains
       *
       *  @param chain_size the chain size
       *
       *  @param nwalkers the number of parallel
       * walkers
       *
       *  @param seed the base seed for initialization
       *
       *  @return none
       */
      void initialize_chains (const int chain_size, const int nwalkers, const int seed);

      /**
       *  @brief initialize the chains
       *
       *  @param chain_size the chain size
       *
       *  @param nwalkers the number of parallel
       * walkers
       *
       *  @param seed the base seed for initialization
       *
       *  @param radius the stardand distribution for
       * normal distribution random extraction
       *
       *  @param ntry number of prior extractions to find the initial guess
       *
       *  @param prior_seed base prior seed
       *
       *  @param max_iter maximum number of iteration 
       *
       *  @param tol the tolerance for minima finding convergence 
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void initialize_chains (const int chain_size, const int nwalkers, const int seed, const double radius, const int ntry, const int prior_seed=23442341, const unsigned int max_iter=100, const double tol=1.e-6, const double epsilon=1.e-3);

      /**
       *  @brief initialize the chains
       *
       *  @param chain_size the chain size
       *
       *  @param nwalkers the number of parallel
       * walkers
       *
       *  @param seed the base seed for initialization
       *
       *  @param radius the stardand distribution for
       * normal distribution random extraction   
       * 
       *  @param guess vector containing initial guess and result
       *  of the likelihood maximization
       *
       *  @param max_iter maximum number of iteration 
       *
       *  @param tol the tolerance for minima finding convergence 
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void initialize_chains (const int chain_size, const int nwalkers, const int seed, const double radius, vector<double> &guess, const unsigned int max_iter=100, const double tol=1.e-6, const double epsilon=1.e-3);

      /**
       *  @brief initialize the chains
       *
       *  @param chain_size the chain size
       *
       *  @param nwalkers the number of parallel
       * walkers
       *
       *  @param seed the base seed for initialization
       *
       *  @param values the means  for
       * normal distribution random extraction
       *
       *  @param radius the stardand distribution for
       * normal distribution random extraction
       *
       *  @return none
       */
      void initialize_chains (const int chain_size, const int nwalkers, const int seed, vector<double> &values, const double radius);

      /**
       *  @brief function that initialize chains
       *
       *  @param chain_size the chain size
       *
       *  @param chain_values the value of the chain 
       *
       *  @return none
       */
      void initialize_chains (const int chain_size, const vector<vector<double>> chain_values);

      /**
       *  @brief function that initialize chains
       *  reading from an input chain file 
       *
       *  @param chain_size the chain size
       *
       *  @param nwalkers the number of parallel
       * walkers
       *
       *  @param input_dir directory of input for chains 
       *
       *  @param input_file file of input for chains 
       *
       *  @return none
       */
      void initialize_chains (const int chain_size, const int nwalkers, const string input_dir, const string input_file);

      /**
       *  @brief show the results on the standard output
       *
       *  @param start starting position in the chains
       *
       *  @param thin take one chain step every \e thin steps
       *
       *  @param seed the base seed for initialization
       *
       *  @return none
       */
      void show_results (const int start, const int thin, const int seed);

      /**
       *  @brief write results on files
       *
       *  @param dir name of the output directory
       *
       *  @param file name of the output file
       *
       *  @param start starting position in the chains
       *
       *  @param thin take one chain step every \e thin steps
       *
       *  @param seed the base seed for initialization
       *
       *  @return none
       */
      void write_results (const string dir, const string output, const int start, const int thin, const int seed);

    };
  }
}

#endif
