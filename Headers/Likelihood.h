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
 *  @file Headers/Likelihood.h
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


namespace cbl {

  namespace statistics {
    
    /**
     *  @class Likelihood Likelihood.h "Headers/Likelihood.h"
     *
     *  @brief The class Likelihood
     *
     *  This class is used to handle objects of type Likelihood. It is
     *  used for all kind of likelihood maximization
     */
    class Likelihood
    {
      
    protected:

      /// data containers 
      std::shared_ptr<data::Data> m_data;

      /// model to test
      std::shared_ptr<Model> m_model; 

      /// likelihood inputs
      std::shared_ptr<void> m_likelihood_inputs;

      /// likelihood parameters
      std::shared_ptr<ModelParameters> m_model_parameters;

      /// type of the likelihood
      LikelihoodType m_likelihood_type = LikelihoodType::_NotSet_;

      /// log-likelihood function
      LogLikelihood_function m_log_likelihood_function;

      /// likelihood function
      Likelihood_function m_likelihood_function;

      /// likelihood function on a grid
      LogLikelihood_function m_likelihood_function_grid;

      /// log-likelihood function on a grid
      LogLikelihood_function m_log_likelihood_function_grid;

      /// the index in extra info for the coordinates in which
      /// the model is computed
      std::vector<size_t> m_x_index;

      /// the index in extra info where the bin weight is stored
      int m_w_index;

      /// use grid
      bool m_use_grid = false;

      /**
       * @brief set the likelihood grid
       * and write the grid on a file
       *
       * @param npoints the number of points
       * to bin the parameter space
       *
       * @param parameter_limits the limits for
       * the parameter space
       *
       * @param output_file the input file
       *
       * @return none
       */
      void m_set_grid_likelihood_1D (const int npoints, const std::vector<std::vector<double>> parameter_limits, const std::string output_file);
      
      /**
       * @brief set the likelihood grid
       * from file
       *
       * @param input_file the input file
       *
       * @return none
       */
      void m_set_grid_likelihood_1D (const std::string input_file);

      /**
       * @brief set the likelihood grid
       * and write the grid on a file
       *
       * @param npoints the number of points
       * to bin the parameter space
       *
       * @param parameter_limits the limits for
       * the parameter space
       *
       * @param output_file the input file
       *
       * @return none
       */
      void m_set_grid_likelihood_2D (const int npoints, const std::vector<std::vector<double>> parameter_limits, const std::string output_file);

      /**
       * @brief set the likelihood grid
       * from file
       *
       * @param input_file the input file
       *
       * @return none
       */
      void m_set_grid_likelihood_2D (const std::string input_file);

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Likelihood
       */
      Likelihood () : m_likelihood_type(LikelihoodType::_NotSet_) {}

      /**
       *  @brief constructor
       *  
       *  @param data pointers to the data container
       *  
       *  @param model pointers to the model 
       *
       *  @param likelihood_type type of likelihood
       *
       *  @param x_index index(s) of the extra info std::vector containing the point(s) where to evaluate the model
       *
       *  @param w_index std::vector containing the data point weight 
       *
       * @param model_parameters parameters of the likelihood
       *
       *  @return object of class Likelihood
       */
      Likelihood (const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const LikelihoodType likelihood_type, const std::vector<size_t> x_index={0,2}, const int w_index=-1, const std::shared_ptr<ModelParameters> model_parameters=NULL);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      virtual ~Likelihood () = default;

      ///@}

      /**
       * @brief return the likelihood parameters
       *
       * @return pointer containing the likelihood parameters
       */
      std::shared_ptr<ModelParameters> parameters () const { return m_model_parameters; }

      /**
       * @brief evaluate the likelihood
       *
       * @param pp the likelihood parameters
       *
       * @return the likelihood \f$ \mathcal{L} \f$
       */
      double operator() (std::vector<double> &pp) const;

      /**
       * @brief evaluate the log-likelihood
       *
       * @param pp the likelihood parameters
       *
       * @return the log-likelihood \f$ \log(\mathcal{L}) \f$
       */
      double log (std::vector<double> &pp) const;

      /**
       * @brief set the data for the likelihood analysis 
       *
       * @param data pointer to the dataset
       *
       * @return none
       */
      void set_data (std::shared_ptr<data::Data> data);

      /**
       * @brief set the model for the likelihood analysis 
       *
       * @param model pointer to the model
       *
       * @param model_parameters parameters of the likelihood
       *
       * @return none
       */
      void set_model (std::shared_ptr<Model> model=NULL, const std::shared_ptr<ModelParameters> model_parameters=NULL);

      /**
       * @brief set the likelihood function 
       * with internal values of LikelihoodType 
       *
       * @return none
       */
      void unset_grid ();

      /**
       * @brief set the likelihood function 
       * as a grid, to speed up computation: this only
       * works for one or two free parameters
       *
       * @param npoints the number of grid points
       * @param parameter_limits
       * @param file the file to read/write the likelihood
       * computed on a grid
       * @param read if true, read the likelihood grid
       * from a list
       *
       * @return none
       */
      void set_grid (const int npoints, const std::vector<std::vector<double>> parameter_limits, const std::string file, const bool read=false);

      /**
       * @brief set the likelihood type using the LikelihoodType object 
       *
       * @param likelihood_type the likelihood type, specified with the 
       * LikelihoodType object
       *
       *  @param x_index index(s) of the extra info std::vector containing the point(s) where to evaluate the model
       *
       *  @param w_index index of the extra info std::vector containing the data point weight 
       *
       * @return none
       */
      void set_function (const LikelihoodType likelihood_type, const std::vector<size_t> x_index={0,2}, const int w_index=-1);

      /**
       * @brief set the likelihood function 
       *
       * @param loglikelihood_function the loglikelihood function
       *
       * @return none
       */
      void set_function (const LogLikelihood_function loglikelihood_function);

      /**
       *  @brief write best-fit results on a file
       *
       *  @param dir_output output directory
       *
       *  @param file output file
       *
       *  @return none
       */
      void write_results (const std::string dir_output, const std::string file); 

      /**
       *  @brief function that maximizes the likelihood, finds the
       *  best-fit parameters and stores them in the model
       *
       *  this function exploits the Nelder-Mead method
       *  https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
       *
       *  the algorithm defines a simplex (i.e a k-dimensional
       *  polytope which is the convex hull of its k+1 vertices) in
       *  the parameter space. At each step, it identifies the simplex
       *  vertex at which the function to be minimised (i.e. the
       *  negative likelihood in this case) has the greatest value,
       *  and moves it, via reflections and scaling, to a new position
       *  in which the function has a lower value. This iteration
       *  stops when the simplex area becomes lower than the
       *  tolerance. For instance, in 2D, the starting vertices of the
       *  simplex (a triangle in 2D) are the following: (start[0],
       *  start[1]) ; (start[0]+epsilon, start[1]) ; (start[0],
       *  start[1]+epsilon)
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
       *  @param epsilon the simplex side
       *
       *  @return none
       */
      void maximize (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3); 

      /**
       *  @brief write the model at xx at bestfit
       *
       *  @param output_dir the output directory
       *  @param output_file the output file
       *  @param parameters the model parameters
       *  @param xx vector of points at which the model is computed
       *  @param yy vector of points at which the model is computed
       *
       *  @return none
       */
      void write_model (const std::string output_dir, const std::string output_file, const std::vector<double> parameters, const std::vector<double> xx={}, const std::vector<double> yy={});

      /**
       *  @brief write the model at xx at bestfit
       *
       *  @param output_dir the output directory
       *  @param output_file the output file
       *  @param xx vector of points at which the model is computed
       *  @param yy vector of points at which the model is computed
       *
       *  @return none
       */
      void write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx={}, const std::vector<double> yy={})
      { write_model(output_dir, output_file, m_model_parameters->bestfit_values(), xx, yy); }

    };
  }
}

#endif
