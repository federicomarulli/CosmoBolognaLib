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
 *  @file Headers/Lib/Model1D.h
 *
 *  @brief The class Model1D
 *
 *  This file defines the interface of the class Model1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODEL1D__
#define __MODEL1D__

#include "Model.h"


// ===================================================================================================


namespace cosmobl {
  
  /**
   *  @brief The namespace of the functions and classes used for
   *  <B> statistical analyses </B>
   *  
   *  The \e statistic namespace contains all the functions and
   *  classes used for statistical analyes
   */
  namespace statistics {

    /**
     *  @class Model1D Model1D.h "Headers/Lib/Model1D.h"
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
      Model1D () = default;

      /**
       *  @brief constructor
       *
       *  @param model model function
       *
       *  @param fixed_parameters list of fixed model parameters
       *
       *  @return object of class Model1D
       */
      Model1D (const model_function_1D model, const shared_ptr<void> fixed_parameters=NULL) : Model(fixed_parameters)
      { set_model(model); }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Model1D () = default;

      ///@}

      /**
       * @brief set model
       *
       * @param model pointer to a model_function_1D
       * 
       * @return none 
       */
      void set_model (const model_function_1D model);

      /**
       * @brief evaluate the model function at xx
       *
       * @param xx the positions at which the model is computed
       *
       * @param parameter vector containing the model parameters
       *
       * @return vector containing model values
       */
      double operator () (const double xx, vector<double> &parameter) const override
      {
	vector<vector<double>> xvec(1, vector<double>(1, xx));
	return m_function(xvec, m_fixed_parameters, parameter)[0][0];
      }

      /**
       * @brief evaluate the model function 
       * at xx 
       *
       * @param xx vector containing the positions at which the model 
       * is computed
       *
       * @param parameter vector containing the model parameters
       *
       * @return vector containing model values
       */
      vector<double> operator () (const vector<double> xx, vector<double> &parameter) const override
        { 
          vector<vector<double>> xvec(1, xx);
          return m_function({xx}, m_fixed_parameters, parameter)[0];
        }

      /**
       *  @brief write the model at xx
       *
       *  @param output_dir the output directory
       *
       *  @param output_file the output file
       * 
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *
       *  @param parameter vector containing the input parameters
       *  used to compute the model; if this vector is not provided,
       *  the model will be computed using the best-fit parameters
       *
       *  @return none
       */
      void write_model (const string output_dir, const string output_file, const vector<double> xx, vector<double> &parameter) const override;
	
    };
  }
}

#endif
