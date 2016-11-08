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
 *  @file Headers/Lib/Modelling.h
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


#include "Likelihood.h"


// ===================================================================================================


namespace cosmobl {
  
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
     *  "Headers/Lib/Modelling.h"
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
      shared_ptr<data::Data> m_data;

      /// input model
      shared_ptr<statistics::Model> m_model;

      /// likelihood
      shared_ptr<statistics::Likelihood> m_likelihood;

      
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

      /**
       * @brief return the dataset
       * @return pointer to the current dataset
       */
      shared_ptr<data::Data> data () { return m_data; }

      /**
       * @brief return the model
       * @return pointer to the current model
       */    
      shared_ptr<statistics::Model> model () { return m_model; }

      /**
       * @brief return the likelihood
       * @return pointer to the current likelihood
       */    
      shared_ptr<statistics::Likelihood> likelihood () { return m_likelihood; }

      ///@}

      
      /**
       * @brief set the dataset
       *
       * @param dataset the dataset 
       *
       * @return none
       */
      void set_data (const shared_ptr<data::Data> dataset) { m_data = move(dataset); }

      
      /**
       * @brief set the model
       *
       * @param model the model
       *
       * @return none
       */    
      void set_model (const shared_ptr<statistics::Model> model) { m_model = move(model); }

      
      /**
       * @brief set the likelihood
       *
       * @param xmin lower bound for the fit
       * @param xmax upper bound for the fit
       * @param likelihood_type type of likelihood
       *
       * @return none
       */
      void set_likelihood (const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type);

      
      /**
       * @brief set the likelihood
       *
       * @param xmin lower bound for the fit, first dimension
       * @param xmax upper bound for the fit, first dimension
       * @param ymin lower bound for the fit, second dimension
       * @param ymax upper bound for the fit, second dimension
       * @param likelihood_type type of likelihood
       *
       * @return none
       */
      void set_likelihood (const double xmin, const double xmax, const double ymin, const double ymax, const statistics::LikelihoodType likelihood_type);

      
      /**
       * @brief sample the likelihood and write
       * chains on a file
       *
       * @param xmin lower bound for the fit
       * @param xmax upper bound for the fit
       * @param likelihood_type type of likelihood
       * @param n_chains number of parallel chains
       * @param chain_size size of the chains
       * @param seed the seed for random number generation
       * @param dir_output output directory
       * @param chain_file output file
       * @param do_write_chain 0 &rarr; do not write the chains during sampling
       * 1 &rarr; write the chains during sampling
       * @param start starting position in the chains
       * @param stop final position in the chains
       * @param thin interval of parameter in output
       *
       * @return none
       */
      void sample_likelihood (const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain=0, const double start=0.5, const double stop=1, const int thin=1);

      
      /**
       * @brief sample the likelihood and write
       * chains on a file
       *
       * @param xmin lower bound for the fit, first dimension
       * @param xmax upper bound for the fit, first dimension
       * @param ymin lower bound for the fit, second dimension
       * @param ymax upper bound for the fit, second dimension
       * @param likelihood_type type of likelihood
       * @param n_chains number of parallel chains
       * @param chain_size size of the chains
       * @param seed the seed for random number generation
       * @param dir_output output directory
       * @param chain_file output file
       * @param do_write_chain 0 &rarr; do not write the chains during sampling
       * 1 &rarr; write the chains during sampling
       * @param start starting position in the chains
       * @param stop final position in the chains
       * @param thin interval of parameter in output
       *
       * @return none
       */
      void sample_likelihood (const double xmin, const double xmax, const double ymin, const double ymax, const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain=0, const double start=0.5, const double stop=1, const int thin=1);

      
      /**
       * @brief sample the likelihood and write
       * chains on a file
       *
       * @param likelihood_type type of likelihood
       * @param n_chains number of parallel chains
       * @param chain_size size of the chains
       * @param seed the seed for random number generation
       * @param dir_output output directory
       * @param chain_file output file
       * @param do_write_chain 0 &rarr; do not write the chains during sampling
       * 1 &rarr; write the chains during sampling
       * @param start starting position in the chains
       * @param stop final position in the chains
       * @param thin interval of parameter in output
       *
       * @return none
       */
      void sample_likelihood (const statistics::LikelihoodType likelihood_type, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain=0, const double start=0.5, const double stop=1, const int thin=1);

      
      /**
       * @brief sample the likelihood and write the chains on a file
       *
       * @param xmin lower bound for the fit
       * @param xmax upper bound for the fit
       * @param likelihood_type type of likelihood
       * @param loglikelihood_function the user-defined loglikelihood
       * function
       * @param cov 0 &rarr; do not use data covariance, 1 &rarr; use data covariance 
       * @param n_chains number of parallel chains
       * @param chain_size size of the chains
       * @param seed the seed for random number generation
       * @param do_write_chain 0 &rarr; do not write the chains during sampling
       * 1 &rarr; write the chains during sampling
       * @param dir_output output directory
       * @param chain_file output file
       * @param start starting position in the chains
       * @param stop final position in the chains
       * @param thin interval of parameter in output
       *
       * @return none
       */
      void sample_likelihood (const double xmin, const double xmax, const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain=0, const double start=0.5, const double stop=1, const int thin=1);

      
      /**
       * @brief sample the likelihood and write
       * chains on a file
       *
       * @param xmin lower bound for the fit, first dimension
       * @param xmax upper bound for the fit, first dimension
       * @param ymin lower bound for the fit, second dimension
       * @param ymax upper bound for the fit, second dimension
       * @param likelihood_type type of likelihood
       * @param loglikelihood_function the user-defined loglikelihood
       * function
       * @param cov 0 &rarr; do not use data covariance, 1 &rarr; use data covariance 
       * @param n_chains number of parallel chains
       * @param chain_size size of the chains
       * @param seed the seed for random number generation
       * @param do_write_chain 0 &rarr; do not write the chains during sampling
       * 1 &rarr; write the chains during sampling
       * @param dir_output output directory
       * @param chain_file output file
       * @param start starting position in the chains
       * @param stop final position in the chains
       * @param thin interval of parameter in output
       *
       * @return none
       */
      void sample_likelihood (const double xmin, const double xmax, const double ymin, const double ymax, const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain=0, const double start=0.5, const double stop=1, const int thin=1);

      
      /**
       * @brief sample the likelihood and write
       * chains on a file
       *
       * @param likelihood_type type of likelihood
       * @param loglikelihood_function the user-defined loglikelihood
       * function
       * @param cov 0 &rarr; do not use data covariance, 1 &rarr; use
       * data covariance
       * @param n_chains number of parallel chains
       * @param chain_size size of the chains
       * @param seed the seed for random number generation
       * @param dir_output output directory
       * @param chain_file output file
       * @param do_write_chain 0 &rarr; do not write the chains during
       * sampling 1 &rarr; write the chains during sampling
       * @param start starting position in the chains
       * @param stop final position in the chains
       * @param thin interval of parameter in output
       *
       * @return none
       */
      void sample_likelihood (const statistics::LikelihoodType likelihood_type, const statistics::LogLikelihood_function loglikelihood_function, const bool cov, const int n_chains, const int chain_size, const int seed, const string dir_output, const string chain_file, const bool do_write_chain=0, const double start=0.5, const double stop=1, const int thin=1);

      
      /**
       *  @name Member functions used to write the outputs
       */
      ///@{
      
      /**
       *  @brief store the best-fit values of the model parameters
       *
       *  @param dir the output directory
       *  @param file the output file
       * 
       *  @return none
       */
      void write_parameters (const string dir, const string file) const;

      /**
       *  @brief compute and write the model
       *
       *  @param xx vector of points at which the model is computed
       *  @param dir the output directory
       *  @param file the output file
       *
       *  @param parameter vector containing input parameters used to
       *  compute the model; if it is not provided, internal
       *  parameters will be used
       *
       *  @return none
       */
      virtual void write_model (const vector<double> xx, const string dir, const string file, const vector<double> parameter={}) const
      { (void)xx; (void)dir; (void)file; (void)parameter; cosmobl::ErrorCBL("Error in write_model of Modelling.h!"); }

      /**
       *  @brief compute and write the model 
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *  @param yy vector of points at which the model is computed,
       *  second axis
       *  @param dir the output directory
       *  @param file the output file
       *
       *  @param parameter vector containing input parameters used to
       *  compute the model; if it is not provided, internal
       *  parameters will be used
       *
       *  @return none
       */
      virtual void write_model (const vector<double> xx, const vector<double> yy, const string dir, const string file, const vector<double> parameter={}) const
      { (void)xx; (void)yy; (void)dir; (void)file; (void)parameter; cosmobl::ErrorCBL("Error in write_model of Modelling.h!"); }

      ///@}
      
    };
  }
}

#endif
