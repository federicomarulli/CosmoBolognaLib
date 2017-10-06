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
 *  @file Headers/Lib/Model.h
 *
 *  @brief The class Model
 *
 *  This file defines the interface of the class Model
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODEL__
#define __MODEL__

#include "Data2D_extra.h"


// ======================================================================================


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
     *  @brief 1D function: the inputs are a vector of values at which
     *  the function is computed, a pointer to a set of data used to
     *  construct the function and a vector of free parameters
     */
    typedef function< vector<double>(vector<double>, shared_ptr<void>, vector<double> &)> model_function_1D;

    /**
     *  @brief 2D function: the inputs are the values at which the
     *  function is computed, a pointer to a set of data used to
     *  construct the function and a vector of free parameters
     */
    typedef function<vector<vector<double>>(vector<double>, vector<double>, shared_ptr<void>, vector<double> &)> model_function_2D;

    /**
     *  @brief generic function: the inputs are the values at which the
     *  function is computed, a pointer to a set of data used to
     *  construct the function and a vector of free parameters
     */
    typedef function<vector<vector<double>>(vector<vector<double>>, shared_ptr<void>, vector<double> &)> model_function_generic;

    /**
     *  @class Model Model.h "Headers/Lib/Model.h"
     *
     *  @brief The class Model
     *
     *  This class is used to define models
     */
    class Model {

      protected:

	/// fixed parameters of the model
	shared_ptr<void> m_fixed_parameters;

	/// the model function
	model_function_generic m_function;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{ 

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class Model
	 */
	Model () = default;

	/**
	 *  @brief constructor
	 *
	 *  @param fixed_parameters fixed parameters of the model
	 *
	 *  @return object of class Model
	 */
	Model (const shared_ptr<void> fixed_parameters)
	  { set_fixed_parameters(fixed_parameters); }

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	virtual ~Model () = default;

	///@}


	/**
	 * @brief set model, for Model1D
	 *
	 * @param model pointer to a model_function_1D
	 * 
	 * @return none 
	 */
	void set_model (const model_function_1D model)
	{ (void)model; ErrorCBL("Error in set_model of Model.h!"); } 

	/**
	 * @brief set model, for Model2D
	 *
	 * @param model pointer to a model_function_2D
	 * 
	 * @return none 
	 */
	void set_model (const model_function_2D model)
	{ (void)model; ErrorCBL("Error in set_model of Model.h!"); } 

	/**
	 *  @brief set the model parameters
	 *
	 *  @param fixed_parameters the fixed parameters of the model
	 *
	 *  @return none
	 */
	void set_fixed_parameters (shared_ptr<void> fixed_parameters) { m_fixed_parameters = fixed_parameters; }

	/**
	 *  @brief return the model parameters
	 *
	 *  @return pointer to the model fixed parameters
	 */
	shared_ptr<void> fixed_parameters () { return m_fixed_parameters; }

	/**
	 * @brief evalueate the model function at xx, for Model1D
	 *
	 * @param xx vector containing the positions at which the model 
	 * is computed
	 *
	 * @param parameter vector containing the model parameters
	 *
	 * @return vector containing model values
	 */
	vector<vector<double>> operator () (const vector<vector<double>> xx, vector<double> &parameter) const
	{
	  vector<double> par = parameter;
	  return m_function(xx, m_fixed_parameters, par);
	}

  	/**
	 * @brief evaluate the model function at xx, for Model1D
	 *
	 * @param xx the positions at which the model 
	 * is computed
	 *
	 * @param parameter vector containing the model parameters
	 *
	 * @return vector containing model values
	 */
	virtual double operator () (const double xx, vector<double> &parameter) const
	{ (void)xx; (void)parameter; ErrorCBL("Error in operator() of Model.h!"); return 0; } 

   	/**
	 * @brief evaluate the model function at xx, for Model1D
	 *
	 * @param xx vector containing the positions at which the
	 * model is computed
	 *
	 * @param parameter vector containing the model parameters
	 *
	 * @return vector containing model values
	 */
	virtual vector<double> operator () (const vector<double> xx, vector<double> &parameter) const
	{ (void)xx; (void)parameter; ErrorCBL("Error in operator() of Model.h!"); vector<double> vv; return vv; } 

  	/**
	 * @brief evaluate the model function at xx, yy, for Model2D
	 *
	 * @param xx the positions at which the model is computed
	 *
	 * @param yy the positions at which the model is computed
	 *
	 * @param parameter vector containing the model parameters
	 *
	 * @return vector containing model values
	 */
	virtual double operator () (const double xx, const double yy, vector<double> &parameter) const
	{ (void)xx; (void)yy; (void)parameter; ErrorCBL("Error in operator() of Model.h!"); return 0; } 

   	/**
	 *  @brief evaluate the model function at xx, yy, for Model2D
	 *
	 *  @param xx vector containing the positions at which the
	 *  model is computed
	 *
	 *  @param yy vector containing the positions at which the
	 *  model is computed
	 *
	 *  @param parameter vector containing the model parameters
	 *
	 *  @return vector containing model values
	 */
	virtual vector<vector<double>> operator () (const vector<double> xx, const vector<double> yy, vector<double> &parameter) const
	{ (void)xx; (void)yy; (void)parameter; ErrorCBL("Error in operator_() of Model.h!"); vector<vector<double>> vv; return vv; } 

	/**
	 *  @brief write the model at xx (implemented by Model1D)
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
	virtual void write_model (const string output_dir, const string output_file, const vector<double> xx, vector<double> &parameter) const 
	{ (void)output_dir; (void)output_file; (void)xx; (void)parameter; cosmobl::ErrorCBL("Error in write_model() of Model.h!"); }

	/**
	 *  @brief write the model at xx, yy (implemented by Model2D)
	 *
	 *  @param output_dir the output directory
	 *
	 *  @param output_file the output file
	 *
	 *  @param xx vector of points at which the model is computed,
	 *  first axis
	 *
	 *  @param yy vector of points at which the model is computed,
	 *  second axis
	 *
	 *  @param parameter vector containing the input parameters
	 *  used to compute the model; if this vector is not provided,
	 *  the model will be computed using the best-fit parameters
	 *
	 *  @return none
	 */
	virtual void write_model (const string output_dir, const string output_file, const vector<double> xx, const vector<double> yy, vector<double> &parameter) const
	{ (void)output_dir; (void)output_file; (void)xx; (void)yy; (void)parameter; cosmobl::ErrorCBL("Error in write_model of Model.h!"); }

    };
  }
}

#endif
