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
 *  @file Headers/Model.h
 *
 *  @brief The class Model
 *
 *  This file defines the interface of the class Model
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODEL__
#define __MODEL__

#include "ModelParameters.h"


// ======================================================================================


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for
   *  <B> statistical analyses </B>
   *  
   *  The \e statistic namespace contains all the functions and
   *  classes used for statistical analyes
   */
  namespace statistics {

    /**
     *  @brief 1D function: the inputs are a vector of values at which
     *  the function is computed, a pointer to a set of data used to
     *  construct the function and a vector of parameters
     */
    typedef std::function<std::vector<double>(std::vector<double>, std::shared_ptr<void>, std::vector<double> &)> model_function_1D;
    
    /**
     *  @brief 2D function: the inputs are the values at which the
     *  function is computed, a pointer to a set of data used to
     *  construct the function and a vector of parameters
     */
    typedef std::function<std::vector<std::vector<double>>(std::vector<double>, std::vector<double>, std::shared_ptr<void>, std::vector<double> &)> model_function_2D;

    /**
     *  @brief generic function: the inputs are the values at which the
     *  function is computed, a pointer to a set of data used to
     *  construct the function and a vector of parameters
     */
    typedef std::function<std::vector<std::vector<double>>(std::vector<std::vector<double>>, std::shared_ptr<void>, std::vector<double> &)> model_function_generic;

    /**
     *  @class Model Model.h "Headers/Model.h"
     *
     *  @brief The class Model
     *
     *  This class is used to define models
     */
    class Model {

    protected:
      
      /// inputs of the model
      std::shared_ptr<void> m_inputs;

      /// parameters of the model
      std::shared_ptr<ModelParameters> m_parameters;

      /// the model function
      model_function_generic m_function;

      /// the model dimension
      Dim m_dimension;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{ 
      
      /**
       *  @brief default constructor
       */
      Model () = default;

      /**
       * @brief constructor
       *
       * @param nparameters the number of parameters
       * 
       * @param parameterTypes the parameter types
       *
       * @param parameterNames the parameter names
       *
       * @param inputs inputs of the model
       */
      Model (const int nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames, const std::shared_ptr<void> inputs)
	{ set_inputs(inputs); set_parameters(nparameters, parameterTypes, parameterNames); }

      /**
       *  @brief default destructor
       *
       *  
       */
      virtual ~Model () = default;

      ///@}

      /**
       * @brief return the model dimension
       *
       * @return the model dimension
       */
      Dim dimension() const { return m_dimension; }

      /**
       * @brief set the model function, for Model1D
       *
       * @param function pointer to a model_function_1D
       * 
       *  
       */
      void set_function (const model_function_1D function)
      { (void)function; ErrorCBL("", "set_function", "Model.h"); } 

      /**
       * @brief set the model function, for Model1D
       *
       * @param function pointer to a model_function_1D
       *
       * @warning function passed in this way doesn't allow for
       * input parameters 
       * 
       *  
       */
      void set_function (const std::vector<double>(*function)(const std::vector<double> xx, std::vector<double> &val))
      { (void)function; ErrorCBL("", "set_function", "Model.h"); } 

      /**
       * @brief set the model function, for Model2D
       *
       * @param function pointer to a model_function_2D
       */
      void set_function (const model_function_2D function)
      { (void)function; ErrorCBL("", "set_function", "Model.h"); } 

      /**
       *  @brief set the model inputs
       *
       *  @param inputs the inputs of the model
       */
      void set_inputs (std::shared_ptr<void> inputs) { m_inputs = inputs; }

      /**
       *  @brief set the model parameters
       * 
       *  @param parameters pointer to parameters
       */
      void set_parameters (const std::shared_ptr<ModelParameters> parameters) {m_parameters = parameters;}

      /**
       *  @brief set the model parameters
       * 
       *  @param nparameters the number of parameters
       *
       *  @param parameterTypes the parameter types
       *
       *  @param parameterNames the parameter names
       */
      void set_parameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames);

      /**
       *  @brief return the model inputs
       *
       *  @return pointer to the model inputs
       */
      std::shared_ptr<void> inputs () { return m_inputs; }

      /**
       *  @brief return the model parameters
       *
       *  @return pointer to an object of class 
       *  ModelParameters
       */
      std::shared_ptr<ModelParameters> parameters () { return m_parameters; }

      /**
       * @brief evalueate the model function at xx, for Model1D
       *
       * @param xx vector containing the positions at which the model 
       * is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      std::vector<std::vector<double>> operator () (const std::vector<std::vector<double>> xx, std::vector<double> &parameters) const;

      /**
       * @brief evaluate the model function at xx, for Model1D
       *
       * @param xx the positions at which the model 
       * is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      virtual double operator () (const double xx, std::vector<double> &parameters) const
      { (void)xx; (void)parameters; ErrorCBL("", "operator ()", "Model.h"); return 0; } 

      /**
       * @brief evaluate the model function at xx, for Model1D
       *
       * @param xx vector containing the positions at which the
       * model is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      virtual std::vector<double> operator () (const std::vector<double> xx, std::vector<double> &parameters) const
      { (void)xx; (void)parameters; ErrorCBL("", "operator ()", "Model.h"); std::vector<double> vv; return vv; } 

      /**
       * @brief evaluate the model function at xx, yy, for Model2D
       *
       * @param xx the positions at which the model is computed
       *
       * @param yy the positions at which the model is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      virtual double operator () (const double xx, const double yy, std::vector<double> &parameters) const
      { (void)xx; (void)yy; (void)parameters; ErrorCBL("", "operator ()", "Model.h"); return 0; } 

      /**
       *  @brief evaluate the model function at xx, yy, for Model2D
       *
       *  @param xx vector containing the positions at which the
       *  model is computed
       *
       *  @param yy vector containing the positions at which the
       *  model is computed
       *
       *  @param parameters vector containing the model parameters
       *
       *  @return vector containing model values
       */
      virtual std::vector<std::vector<double>> operator () (const std::vector<double> xx, const std::vector<double> yy, std::vector<double> &parameters) const
      { (void)xx; (void)yy; (void)parameters; ErrorCBL("", "operator ()", "Model.h"); std::vector<std::vector<double>> vv; return vv; } 

      /**
       *  @brief compute the median and percentiles of the model
       *  from MCMC chains
       *
       *  @param xx vector of points at which the model is computed,
       *  @param median_model the median model
       *  @param low_model the model at 16th percentile
       *  @param up_model the model at 84th percentile
       *  @param start the starting position for each chain
       *  @param thin the position step
       *  
       */
      void stats_from_chains (const std::vector<std::vector<double>> xx, std::vector<std::vector<double>> &median_model, std::vector<std::vector<double>> &low_model, std::vector<std::vector<double>> &up_model, const int start=0, const int thin=1);

      /**
       *  @brief compute the median and percentiles of the model
       *  from MCMC chains
       *
       *  @param xx vector of points at which the model is computed,
       *  @param median_model the median model
       *  @param low_model the model at 16th percentile
       *  @param up_model the model at 84th percentile
       *  @param start the starting position for each chain
       *  @param thin the position step
       */
      virtual void stats_from_chains (const std::vector<double> xx, std::vector<double> &median_model, std::vector<double> &low_model, std::vector<double> &up_model, const int start=0, const int thin=1)
      { (void)xx; (void)median_model; (void)low_model; (void)up_model; (void)start; (void)thin; ErrorCBL("", "stats_from_chains", "Model.h"); } 


      /**
       *  @brief compute the median and percentiles of the model
       *  from MCMC chains
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *  @param yy vector of points at which the model is computed,
       *  second axis
       *  @param median_model the median model
       *  @param low_model the model at 16th percentile
       *  @param up_model the model at 84th percentile
       *  @param start the starting position for each chain
       *  @param thin the position step
       */
      virtual void stats_from_chains (const std::vector<double> xx, const std::vector<double> yy, std::vector<std::vector<double>> &median_model, std::vector<std::vector<double>> &low_model, std::vector<std::vector<double>> &up_model, const int start=0, const int thin=1)
      { (void)xx; (void)yy; (void)median_model; (void)low_model; (void)up_model; (void)start; (void)thin; ErrorCBL("", "stats_from_chains", "Model.h"); } 

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
       */
      virtual void write (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters)
      { (void)output_dir; (void)output_file; (void)xx; (void)parameters; cbl::ErrorCBL("", "write", "Model.h"); }

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
       */
      virtual void write (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> parameters)
      { (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)parameters; cbl::ErrorCBL("", "write", "Model.h"); }

      /**
       *  @brief write the model at xx 
       *  with best-fit parameters obtained from likelihood
       *  maximization
       *
       *  @param output_dir the output directory
       *  @param output_file the output file
       *  @param xx vector of points at which the model is computed
       */
      virtual void write_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx)
      { (void)output_dir; (void)output_file; (void)xx; cbl::ErrorCBL("", "write_at_bestfit", "Model.h"); }

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
       */
      virtual void write_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy)
      { (void)output_dir; (void)output_file; (void)xx; (void)yy; cbl::ErrorCBL("", "write_at_bestfit", "Model.h"); }

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
       */
      virtual void write_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start=0, const int thin=1)
      { (void)output_dir; (void)output_file; (void)xx; (void)start; (void)thin; cbl::ErrorCBL("", "write_from_chains", "Model.h"); }

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
       */
      virtual void write_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const int start=0, const int thin=1)
      { (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)start; (void)thin; cbl::ErrorCBL("", "write_from_chains", "Model.h"); }

    };
  }
}

#endif
