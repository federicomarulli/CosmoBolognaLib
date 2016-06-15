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

  /**
   *  @brief The namespace of functions and classes used for statistical
   *  analysis
   *  
   * The \e statistic namespace contains all the functions and classes
   * used for statistical analyis
   */
  namespace statistics {
    
    /// type of the likelihood to be used
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

    struct STR_likelihood_parameters
    {
      shared_ptr<Data> data;
      shared_ptr<Model> model;
      
      STR_likelihood_parameters (shared_ptr<Data> _data, shared_ptr<Model> _model) :
      data(_data), model(_model) {}
    };


    double LogLikelihood_Gaussian_1D_model (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double LogLikelihood_Gaussian_1D_error (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double LogLikelihood_Gaussian_1D_covariance (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);

    double LogLikelihood_Gaussian_2D_model (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double LogLikelihood_Gaussian_2D_error (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);

    typedef function<double(const vector<double>, const shared_ptr<void>)> LogLikelihood_function;

    class Likelihood
    {
      
    protected:

      shared_ptr<Data> m_data;
      shared_ptr<Model> m_model; 

      LikelihoodType m_likelihood_type;

      LogLikelihood_function m_log_likelihood_function;

      int m_nchains;
      int m_chain_size;
	
      bool m_npar;

      bool m_cov;

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
      Likelihood (const shared_ptr<Data> data, const shared_ptr<Model> model, const LikelihoodType likelihood_type);

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
      Likelihood (const shared_ptr<Data> data, const shared_ptr<Model> model, const LikelihoodType likelihood_type, const LogLikelihood_function loglikelihood_function, const bool cov);

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Likelihood () {}

      ///@}

      double operator () (vector<double> pp) 
      {	
	shared_ptr<void> pars = make_shared<STR_likelihood_parameters>(STR_likelihood_parameters(m_data, m_model));
	return m_log_likelihood_function(pp, pars);
      }

      void set_likelihood_type (const LikelihoodType likelihood_type, const bool cov);

      void set_data (shared_ptr<Data> data);

      void set_model (shared_ptr<Model> model);

      void set_likelihood_function ();

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
       *  @return averace acceptance ratio
       */
      double sample_stretch_move (const int nchains, const int chain_size, const int seed);

      /*
       *  @brief function that samples likelihood, using
       *  Metropolis-Hastings algorithm on uni-dimensional parameter
       *  space, and stores chain parameters.       
       *  @param nchains number of chains to sample the parameter space 
       *  @param chain_size number of step in each chain 
       *  @param shift  percentage shift for proposed parameter 
       *  @param nsubstep number of substeps in each iteration 
       *  @return averace acceptance ratio
       */
      // double sample (const int nchains, const int chain_size, const double shift, const int nsubstep = 100); 

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
