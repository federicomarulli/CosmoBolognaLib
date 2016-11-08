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

#include "Chi2.h"


// ============================================================================================


namespace cosmobl {

  namespace statistics {
    
    /**
     * @enum LikelihoodType
     * @brief the type of likelihood function
     */
    enum LikelihoodType {

      /// not set
      _Likelihood_NotSet_,

      /// Gaussian likelihood model
      _GaussianLikelihood_Model_,

      /// Gaussian likelihood error
      _GaussianLikelihood_Error_,
      
      /// Gaussian likelihood covariance
      _GaussianLikelihood_Covariance_,
      
      /// Likelihood function defined by the user
      _UserDefinedLikelihood_
    };

    
    /**
     * @struct STR_likelihood_parameters
     * @brief the struct STR_likelihood_parameters
     *
     * This struct contains the data
     * and the model for the likelihood analysis
     */
    struct STR_likelihood_parameters
    {
      /// data containers
      shared_ptr<data::Data> data;

      /// model to test
      shared_ptr<Model> model;

      /**
       *  @brief constructor
       *  @param _data pointers to the data container
       *  @param _model pointers to the model 
       *  @return object of type STR_likelihood_parameters
       */ 
      STR_likelihood_parameters (const shared_ptr<data::Data> _data, const shared_ptr<Model> _model)
      : data(_data), model(_model) {}
    };


    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  @param model_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood; 
     */
    double LogLikelihood_Gaussian_1D_model (vector<double> model_parameters, const shared_ptr<void> inputs);

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  @param model_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_1D_error (vector<double> model_parameters, const shared_ptr<void> inputs);

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  @param model_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_1D_covariance (vector<double> model_parameters, const shared_ptr<void> inputs);

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  @param model_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_2D_model (vector<double> model_parameters, const shared_ptr<void> inputs);

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  model with one parameter  &chi;&sup2; 
     *  @param model_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_2D_error (vector<double> model_parameters, const shared_ptr<void> inputs);

    /**
     * @var typedef LogLikelihood_function
     * @brief definition of a function for computation of 
     * the LogLikelihood
     */
    typedef function<double (const vector<double>, const shared_ptr<void>)> LogLikelihood_function;
    
    /**
     *  @class Likelihood Likelihood.h "Headers/Lib/Likelihood.h"
     *
     *  @brief The class Likelihood
     *
     *  This class is used to handle objects of type likelihood. It is
     *  used for all kind of likelihood analyses, sample and minimization
     */
    class Likelihood
    {
      
    protected:

      /// data containers 
      shared_ptr<data::Data> m_data;

      /// model to test
      shared_ptr<Model> m_model; 

      /// type of the likelihood
      LikelihoodType m_likelihood_type;

      /// likelihood function
      LogLikelihood_function m_log_likelihood_function;

      /// number of chains
      int m_nchains;

      /// size of the chains
      int m_chain_size;

      /// number of parameters
      bool m_npar;

      /// use data covariance matrix; 0 &rarr; don't use data covariance matrix 
      /// 1 &rarr; use data covariance matrix; 
      bool m_cov;

      /**
       *  @brief function that write sampled parameters 
       *  @param output_dir directory of output for chains 
       *  @param output_file file of output for chains 
       *  @param start starting position in the chains
       *  @param stop final position in the chains
       *  @param thin interval of parameter in output
       *  @return none
       */
      void m_write_chain (const string output_dir, const string output_file, const int start=-1, const int stop=-1, const int thin=1);

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
       *  @param data pointers to the data container
       *  @param model pointers to the model 
       *  @param likelihood_type type of likelihood
       *  @return object of class Likelihood
       */
      Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const LikelihoodType likelihood_type);

      /**
       *  @brief constructor
       *  @param data pointers to the data container
       *  @param model pointers to the model 
       *  @param likelihood_type type of likelihood
       *  @param loglikelihood_function function of type loglikelihood_function
       *  @param cov if cov=0 &rarr do not invert covariance matrix; else &rarr invert the covariance matrix;
       *  @param cov 
       *  @return object of class Likelihood
       */
      Likelihood (const shared_ptr<data::Data> data, const shared_ptr<Model> model, const LikelihoodType likelihood_type, const LogLikelihood_function loglikelihood_function, const bool cov);

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Likelihood () = default;

      ///@}

      /**
       * @brief return the best-fit value of the i-th parameter
       *
       * @param i the parameter index
       *
       * @return the best-fit value of the i-th parameter
       */
      double best_fit (const int i) { return m_model->parameter(i)->value(); }

	
      /**
       * @brief evaluate the likelihood
       *
       * @param pp the model parameters
       *
       * @return none
       */
      double operator () (vector<double> pp) 
      {	
	shared_ptr<void> pars = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data, m_model));
	return m_log_likelihood_function(pp, pars);
      }

      /**
       * @brief set the likelihood type using the LikelihoodType object 
       *
       * @param likelihood_type the likelihood type, specified with the 
       * LikelihoodType object
       *
       * @param cov use data covariance matrix; 0 &rarr; don't use data covariance
       * matrix, 1 &rarr; use data covariance matrix; 
       *
       * @return none
       */
      void set_likelihood_type (const LikelihoodType likelihood_type, const bool cov);

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
       * @brief set the likelihood function 
       * according among pre-defined likelihood functions 
       *
       * @return none
       */
      void set_likelihood_function ();

      /**
       * @brief set the likelihood function 
       *
       * @param loglikelihood_function the likelihood function
       *
       * @param cov use data covariance matrix; 0 &rarr; don't use data covariance
       * matrix, 1 &rarr; use data covariance matrix; 
       *
       * @return none
       */
      void set_likelihood_function (const LogLikelihood_function loglikelihood_function, const bool cov);

      /**
       *  @brief return the Log prior probability
       *  @return value of the Log prior probability
       */
      double LogPriorProbability ();

      /**
       *  @brief return the Log prior probability
       *  @param parameter_values values of the parameters
       *  @return value of the Log Prior probability
       */
      double LogPriorProbability (const vector<double> parameter_values);

      /**
       *  @brief function that maximize Likelihood, find best fit
       *  parameters and store them in model       
       *  @param usePriors minimize the Log likelihood using the Priors
       *  @param max_iter maximum number of iteration 
       *  @param tol the tolerance for minima finding convergence 
       *  @return none
       */
      void minimize_LogLikelihood (const bool usePriors=true, const unsigned int max_iter=100, const double tol=1.e-6); 

      /**
       *  @brief function that samples likelihood, using stretch-move
       *  algorithm on n-dimensional parameter space, and stores chain
       *  parameters.
       *  @param nchains number of chains to sample the parameter space 
       *  @param chain_size number of step in each chain 
       *  @param seed the seed for random number generator
       *  @param do_write_chain 0 &rarr; don't write chains at each step 
       *  1 &rarr; write chains at each step
       *  @param output_dir directory of output for chains 
       *  @param output_file file of output for chains 
       *  @return averace acceptance ratio
       */
      double sample_stretch_move (const int nchains, const int chain_size, const int seed, bool do_write_chain = 0, const string output_dir=par::defaultString, const string output_file=par::defaultString);

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
      double sample_tabulated_likelihood (const int nstep_p1, const int nstep_p2, const string interpolation_method, const int nchains, const int chain_size, const int seed, bool do_write_chain = 0, const string output_dir=par::defaultString, const string output_file=par::defaultString);
      
      /**
       *  @brief function that write sampled parameters 
       *  @param output_dir directory of output for chains 
       *  @param output_file file of output for chains 
       *  @param start starting position in the chains
       *  @param stop final position in the chains
       *  @param thin interval of parameter in output
       *  @return none
       */
      void write_chain (const string output_dir, const string output_file, const double start=0.5, const double stop=1, const int thin=1);
      
    };
  }
}

#endif
