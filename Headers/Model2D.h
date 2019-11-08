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
 *  @file Headers/Model2D.h
 *
 *  @brief The class Model2D
 *
 *  This file defines the interface of the class Model2D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODEL2D__
#define __MODEL2D__

#include "Model1D.h"


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
     *  @class Model2D Model2D.h "Headers/Model2D.h"
     *
     *  @brief The class Model2D
     *
     *  This class is used to define 2D models
     */
    class Model2D : public Model {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{ 

      /**
       *  @brief default constructor
       *
       *  @return object of class Model2D
       */
      Model2D () : Model() {m_dimension=Dim::_2D_;} 

      /**
       *  @brief constructor
       *
       *  @param function the model function
       *
       *  @param nparameters the number of parameters
       *
       *  @param parameterTypes the parameter types
       *
       *  @param parameterNames the parameter names
       *
       *  @param inputs list of model inputs
       *
       *  @return object of class Model2D
       */
      Model2D (model_function_2D function, const size_t nparameters, std::vector<ParameterType> parameterTypes={}, std::vector<std::string> parameterNames={}, std::shared_ptr<void> inputs=NULL)
	: Model(nparameters, parameterTypes, parameterNames, inputs) { set_function(function); m_dimension=Dim::_2D_; }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Model2D () = default;

      ///@}

      /**
       * @brief set the model function
       *
       * @param function the model function
       * 
       * @return none 
       */
      void set_function (const model_function_2D function);

      /**
       * @brief evaluate the model function at xx, yy
       *
       * @param xx the positions at which the model is computed
       *
       * @param yy the positions at which the model is computed
       *
       * @param parameters vector containing the model parameters
       *
       * @return vector containing model values
       */
      double operator () (const double xx, const double yy, std::vector<double> &parameters) const override
      {
	std::vector<std::vector<double>> xvec(2); xvec[0] = {xx}; xvec[1] = {yy};
	return Model::operator()(xvec, parameters)[0][0];
      }

      /**
       *  @brief evaluate the model function at xx, yy
       *
       *  @param xx vector containing the positions at which the model
       *  is computed
       *
       *  @param yy vector containing the positions at which the model 
       *  is computed
       *
       *  @param parameters vector containing the model parameters
       *
       *  @return vector containing model values
       */
      std::vector<std::vector<double>> operator () (const std::vector<double> xx, const std::vector<double> yy, std::vector<double> &parameters) const override
      {
	std::vector<std::vector<double>> xvec = {xx, yy};
	return Model::operator()(xvec, parameters);
      }

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
       *  @return none
       */
      void stats_from_chains (const std::vector<double> xx, const std::vector<double> yy, std::vector<std::vector<double>> &median_model, std::vector<std::vector<double>> &low_model, std::vector<std::vector<double>> &up_model, const int start=0, const int thin=1) override;

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
       *
       *  @return none
       */
      void write (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> parameters) override;

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
       *
       *  @return none
       */
      void write_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy) override;

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
       *
       *  @return none
       */
      void write_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> yy, const int start=0, const int thin=1) override;

    };
  }
}

#endif
