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

#include "Parameter.h"

namespace cosmobl {

  /**
   *  @brief The namespace of functions and classes used for statistical
   *  analysis
   *  
   * The \e statistic namespace contains all the functions and classes
   * used for statistical analyis
   */

  namespace statistics {

    typedef function<double() > random_numbers;

    typedef function<double(double,  shared_ptr<void> , vector<double>)> model_function_1D;

    typedef function< vector<double>(vector<double> ,  shared_ptr<void> , vector<double>)> model_function_1D_vector;

    typedef function<double(double, double, shared_ptr<void> , vector<double>)> model_function_2D;

    /**
     *  @class Model Model.h "Headers/Lib/Model.h"
     *
     *  @brief The class Model
     *
     *  This class is used to define the model
     */
    class Model {

    protected:

      /// number of model parameters
      unsigned int m_npar;

      /// the effective number of model parameters
      unsigned int m_npar_eff;

      /// function parameters
      vector<shared_ptr<Parameter> > m_parameters;

      /// fixed parameters of the model
      shared_ptr<void> m_model_parameters;

      /// 0 &rarr; 1D functon; 1 &rarr; 2D function
      bool m_2d;

      /// Random number generator
      default_random_engine m_generator;

      /// Random number distribution;
      uniform_real_distribution<double> distribution;

      random_numbers m_random_numbers;



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
      Model () {}

      /**
       *  @brief constructor
       *
       *  @param parameters list of parameters for the model
       *  @param model_parameters fixed parameters of the model
       *
       *  @return object of class Model
       */
      Model (const vector<Parameter> parameters, const shared_ptr<void> model_parameters);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      virtual ~Model() {}

      ///@}

      virtual double operator () (const double xx) = 0;

      virtual double operator () (const double xx, const vector<double> pars) = 0;

      virtual vector<double> operator () (const vector<double> xx) = 0;

      virtual vector<double> operator () (const vector<double> xx, const vector<double> pars) = 0;

      virtual double operator () (const double xx, const double yy) = 0;

      virtual double operator () (const double xx, const double yy, const vector<double> pars) = 0;

      /**
       *  @brief Return the private member m_npar
       *
       *  @return the number of parameters
       */
      unsigned int npar() { return m_npar; }

      /**
       *  @brief Return the private member m_npar_eff
       *
       *  @return the number of effective parameters (freeze == 0)
       */
      unsigned int npar_eff ()
      {
	m_npar_eff = 0;
	for(auto &&pp : m_parameters)
	  m_npar_eff = (pp->isFreezed()) ? m_npar_eff : m_npar_eff+1;

	return m_npar_eff;
      }

      /**
       *  @brief Return the provate variable m_2d;
       *
       *  @return the number of model variables 
       */
      bool dimension () { return m_2d; }

      /**
       *  @brief Return the i-th element of the method m_parameters
       *
       *  @param i: the i-th parameter of the model
       *
       *  @return the i-th parameter of the model
       */
      shared_ptr<Parameter> parameter(const int i) { return m_parameters[i]; }

      /**
       *  @brief Return the private method m_parameters
       *
       *  @return the parameters of the model
       */
      vector<shared_ptr<Parameter> > parameters () { return m_parameters; }

      /**
       *  @brief Update one parameter and return the current parameter
       *  values
       *  @param new_parameter the new parameter
       *  @return the value of model parameter
       */
      double update_parameter (const double new_parameter); 

      /**
       *  @brief Update parameters and return a vector containing the
       *  current parameter values
       *  @param new_parameters the new parameters
       *  @return the values of model parameters
       */
      vector<double> update_parameters (const vector<double> new_parameters); 

      /**
       *  @brief set the chains for likelihood analysis
       *  @param nchains the number of chains
       *  @param chain_size the chain size
       *  @return none
       */
      void set_chains (const int nchains, const int chain_size);
	
      /**
       *  @brief set the random_number generator
       *  @param generator the random numbers generator
       *  @return none
       */
      void set_random_number_generator (const default_random_engine generator);
    };


    /**
     *  @class Model1D Model.h "Headers/Lib/Model.h"
     *
     *  @brief The class Model1D
     *
     *  This class is used to define the 1D model
     */
    class Model1D : public Model {

    protected:

      /// model: functional form for the 1D model
      model_function_1D m_model;

      model_function_1D_vector m_model_vector;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{ 

      /**
       *  @brief default constructor
       *
       * return object of class Model1D
       */
      Model1D () : Model() { m_2d=0; }

      /**
       *  @brief default constructor
       *
       *  @param parameters: list of model parameters
       *  @param model_parameters: list of fixed model parameters
       *  @param model: model function
       *
       * return object of class Model1D
       */
      Model1D (const vector<Parameter> parameters, const shared_ptr<void> model_parameters, const model_function_1D model) :
      Model(parameters, model_parameters), m_model(model) { m_2d = 0; }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Model1D () {}

      ///@}

      /**
       *  @brief evaluate the model for the value xx
       *
       *  @param xx: the variable position
       * 
       *  @return the value of the model at xx
       */
      double operator () (const double xx)
      {
	vector<double> parameters;
	for (auto &&p : m_parameters)
	  parameters.push_back(p->value());
	return m_model(xx, m_model_parameters, parameters);
      }

      /**
       *  @brief evaluate the model for the value xx
       *
       *  @param xx the variable position
       *  @param parameters vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */
      double operator () (const double xx, const vector<double> parameters)
      {
	return m_model(xx, m_model_parameters, parameters);
      }
	
      /**
       *  @brief evaluate the model for the value xx
       *
       *  @param xx: the variable position vector
       * 
       *  @return the value of the model at xx
       */
      vector<double> operator () (const vector<double> xx)
	{
	  vector<double> parameters;
	  for (auto &&p : m_parameters)
	    parameters.push_back(p->value());
	  return m_model_vector(xx, m_model_parameters, parameters);
	}

      /**
       *  @brief evaluate the model for the value xx
       *
       *  @param xx the variable position vector
       *
       *  @param parameters vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */
      vector<double> operator () (const vector<double> xx, const vector<double> parameters)
	{
	  return m_model_vector(xx, m_model_parameters, parameters);
	}

      double operator () (const double xx, const double yy)
      { cosmobl::ErrorMsg("Error in Model::operator() of Model1D.h!"); return 0; }

      double operator () (const double xx, const double yy, const vector<double> pars)
      { cosmobl::ErrorMsg("Error in Model::operator() of Model1D.h!"); return 0; }
    };


    /**
     *  @class Model2D Model.h "Headers/Lib/Model.h"
     *
     *  @brief The class Model2D
     *
     *  This class is used to define the 2D model
     */
    class Model2D : public Model{

    protected:

      /// the functional form for the 2D model
      model_function_2D m_model;

    public:
      /**
       *  @name Constructors/destructors
       */
      ///@{ 

      /** 
       *  @brief default constructor
       *
       * return object of class Model2D
       */
      Model2D () : Model() { m_2d = 1; } 

      /** 
       *  @brief default constructor
       *
       *  @param parameters: list of model parameters
       *  @param model_parameters: list of fixed model parameters
       *  @param model: model function
       *
       * return object of class Model2D
       */
      Model2D (vector<Parameter> parameters, shared_ptr<void> model_parameters, model_function_2D model)
	: Model(parameters,model_parameters), m_model(model) { m_2d = 1; }

      /**
       *  @brief default destructor
       *
       * return none
       */
      ~Model2D () {}

      ///@}


      /**
       *  @brief evaluate the model at the values x and y
       *
       *  @param xx the x variable
       *  @param yy the y variable
       * 
       *  @return the value of the model at x,y
       */
      double operator () (const double xx, const double yy)
      {
	vector<double> parameters;
	for (auto &&p : m_parameters)
	  parameters.push_back(p->value());

	return m_model(xx, yy, m_model_parameters, parameters);
      }

      /**
       *  @brief evaluate the model for the values x and y, using given
       *  parameters
       *
       *  @param xx the x variable
       *  @param yy the y variable
       *  @param parameters the model parameters
       * 
       *  @return the value of the model at xx
       */  
      double operator () (const double xx, const double yy, const vector<double> parameters)
      {
	return m_model(xx, yy, m_model_parameters, parameters);
      }

      virtual double operator () (const double xx)
      { cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return 0; }

      virtual double operator () (const double xx, const vector<double> pars)
      { cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return 0; }

      virtual vector<double> operator () (const vector<double> xx)
      { cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return {}; }

      virtual vector<double> operator () (const vector<double> xx, const vector<double> pars)
      { cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return {}; }

    };
  }
}

#endif
