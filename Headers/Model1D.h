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
 *  @file Headers/Model1D.h
 *
 *  @brief The class Model1D
 *
 *  This file defines the interface of the class Model1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODEL1D__
#define __MODEL1D__

#include "Model.h"


// ===================================================================================================


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
     *  @class Model1D Model1D.h "Headers/Model1D.h"
     *
     *  @brief The class Model1D
     *
     *  This class is used to define 1D models
     */
    class Model1D : public Model {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{ 

      /**
       *  @brief default constructor
       *
       *  @return object of class Model1D
       */
      Model1D () : Model() {m_dimension=Dim::_1D_;}

      /**
       *  @brief constructor
       *
       *  @param function model function
       *
       *  @param nparameters the number of parameters
       *
       *  @param parameterTypes the parameter types
       *
       *  @param parameterNames the parameter names
       *
       *  @param inputs inputs of the model
       *
       *  @return object of class Model1D
       */
      Model1D (const model_function_1D function, const size_t nparameters, std::vector<ParameterType> parameterTypes={}, std::vector<std::string> parameterNames={}, const std::shared_ptr<void> inputs=NULL) : Model(nparameters, parameterTypes, parameterNames, inputs)
      { set_function(function); m_dimension=Dim::_1D_; }

      /**
       *  @brief constructor
       *
       *  @param function model function
       *
       *  @param nparameters the number of parameters
       *
       *  @param parameterTypes the parameter types
       *
       *  @param parameterNames the parameter names
       *
       *  @param inputs inputs of the model
       *
       *  @return object of class Model1D
       */
      Model1D (const std::vector<double> (*function)(const std::vector<double> xx, std::vector<double> &val), const size_t nparameters, std::vector<ParameterType> parameterTypes={}, std::vector<std::string> parameterNames={}, const std::shared_ptr<void> inputs=NULL) : Model(nparameters, parameterTypes, parameterNames, inputs)
      { set_function(function); m_dimension=Dim::_1D_; }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Model1D () = default;

      ///@}

      /**
       * @brief set the model function
       *
       * @param function pointer to a model_function_1D
       * 
       * @return none 
       */
      void set_function (const model_function_1D function);

      /**
       * @brief set the model function
       *
       * @param function pointer to a model_function_1D
       *
       * @warning function passed in this way doesn't allow for
       * input parameters
       * 
       * @return none 
       */
      void set_function (const std::vector<double> (*function)(const std::vector<double> xx, std::vector<double> &val));

      /**
       * @brief evaluate the model function at xx
       *
       * @param xx the positions at which the model is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      double operator () (const double xx, std::vector<double> &parameters) const override
      {
	std::vector<std::vector<double>> xvec(1, std::vector<double>(1, xx));
	return Model::operator()(xvec, parameters)[0][0];
      }

      /**
       * @brief evaluate the model function 
       * at xx 
       *
       * @param xx vector containing the positions at which the model 
       * is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      std::vector<double> operator () (const std::vector<double> xx, std::vector<double> &parameters) const override
      { 
	std::vector<std::vector<double>> xvec(1, xx);
	return Model::operator()(xvec, parameters)[0];
      }

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
       *  @return none
       */
      void stats_from_chains (const std::vector<double> xx, std::vector<double> &median_model, std::vector<double> &low_model, std::vector<double> &up_model, const int start=0, const int thin=1) override;

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
      void write (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters) override;

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
      void write_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx) override;

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
      void write_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start=0, const int thin=1) override;
	
    };
  }
}

#endif
