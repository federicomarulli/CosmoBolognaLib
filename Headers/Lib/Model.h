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
   *  @brief The namespace of the functions and classes used for
   *  <B> statistical analyses </B>
   *  
   *  The \e statistic namespace contains all the functions and classes
   *  used for statistical analyes
   */
  namespace statistics {

    /**
     * @var typedef model_function_1D
     * @brief 1D functin that takes as an input the position as which it's computed,
     * a pointer to fixed parameters and a vector of free parameters
     */
    typedef function<double(double, shared_ptr<void>, vector<double>)> model_function_1D;

    /**
     * @var typedef model_function_1D_vector
     * @brief 1D function that takes as an input the positions as which it's computed,
     * a pointer to fixed parameters and a vector of free parameters
     */
    typedef function< vector<double>(vector<double>, shared_ptr<void>, vector<double>)> model_function_1D_vector;

    /**
     * @var typedef model_function_2D
     * @brief 2D function that takes as an input the positions as which it's computed,
     * a pointer to fixed parameters and a vector of free parameters
     */
    typedef function<double(double, double, shared_ptr<void>, vector<double>)> model_function_2D;

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
      shared_ptr<void> m_fixed_parameters;

      /// 0 &rarr; 1D functon; 1 &rarr; 2D function
      bool m_2d;

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
       *  @param parameters list of parameters for the model
       *  @param fixed_parameters fixed parameters of the model
       *
       *  @return object of class Model
       */
      Model (const vector<Parameter> parameters, const shared_ptr<void> fixed_parameters);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      virtual ~Model() {}

      ///@}


      /**
       *  @brief evaluate the model for the value xx
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx: the variable position
       * 
       *  @return the value of the model at xx
       */
      virtual double operator () (const double xx) = 0;

      /**
       *  @brief evaluate the model for the value xx
       *  the model uses input parameters
       *
       *  @param xx: the variable position
       *  @param pars vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */    
      virtual double operator () (const double xx, const vector<double> pars) = 0;

      /**
       *  @brief evaluate the model for the values xx
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx: the variable position vector
       * 
       *  @return the value of the model at xx
       */
      virtual vector<double> operator () (const vector<double> xx) = 0;

      /**
       *  @brief evaluate the model for the values xx
       *  the model uses input parameters
       *
       *  @param xx: the variable position vector
       *  @param pars vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */   
      virtual vector<double> operator () (const vector<double> xx, const vector<double> pars) = 0;

      /**
       *  @brief evaluate the model for the value (xx,yy)
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx: the X variable position
       *  @param yy: the Y variable position
       * 
       *  @return the value of the model at (xx,yy)
       */ 
      virtual double operator () (const double xx, const double yy) = 0;

      /**
       *  @brief evaluate the model for the value (xx,yy)
       *  the model uses input parameters
       *
       *  @param xx: the X variable position
       *  @param yy: the Y variable position
       *  @param pars vector containing the model parameters
       * 
       *  @return the value of the model at (xx,yy)
       */ 
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
      shared_ptr<Parameter> parameter (const int i) { return m_parameters[i]; }

      /**
       *  @brief Return the private method m_parameters
       *
       *  @return the parameters of the model
       */
      vector<shared_ptr<Parameter> > parameters () { return m_parameters; }

      /**
       *  @brief Return the stored values of model parameters
       *  @return vector containing the stored values of model parameters
       */
      vector<double> parameter_values ();

      /**
       *  @brief Return the values of model parameters in the i-th position of the n-th chain
       *  @param chain the n-th chain
       *  @param position the i-th position in the chain
       *  @return vector containing values of model parameters
       */
      vector<double> parameter_values_from_chain (const int chain, const int position);

      /**
       *  @brief Update parameters  
       *  @param parameter_values the new parameters
       *  @return none
       */
      void set_parameter_values (const vector<double> parameter_values); 

      /**
       *  @brief set the chains for likelihood analysis
       *  @param nchains the number of chains
       *  @param chain_size the chain size
       *  @return none
       */
      void set_chains (const int nchains, const int chain_size);

      /**
       *  @brief return the model parameters
       *  @return pointer to the model fixed parameters
       */
      shared_ptr<void> fixed_parameters () { return m_fixed_parameters; }

      /**
       * @brief compute and write the model using the stored 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model(const vector<double> xx, const string dir_model, const string file_model)
      { (void)xx; (void)dir_model; (void)file_model; cosmobl::ErrorMsg("Error in write_model of Model.h!"); }

      /**
       * @brief compute and write the model using input 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed
       * @param parameters vector of parameters values
       * at which the model is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model_parameters(const vector<double> xx, const vector<double> parameters, const string dir_model, const string file_model)
      { (void)xx; (void)parameters; (void)dir_model; (void)file_model; cosmobl::ErrorMsg("Error in write_model of Model.h!"); }

      /**
       * @brief compute and write the model using the stored 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed, first axis
       * @param yy vector of point at which the model
       * is computed, second axis
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model(const vector<double> xx, const vector<double> yy, const string dir_model, const string file_model)
      { (void)xx; (void)yy; (void)dir_model; (void)file_model; cosmobl::ErrorMsg("Error in write_model of Model.h!"); }

      /**
       * @brief compute and write the model using input 
       * parameter values
       *
       * @param xx vector of point at which the model is computed
       * @param yy vector of point at which the model is computed,
       * second axis
       * @param parameters vector of parameters values at which the
       * model is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model_parameters(const vector<double> xx, const vector<double> yy, const vector<double> parameters, const string dir_model, const string file_model)
      { (void)xx; (void)yy; (void)parameters; (void)dir_model; (void)file_model; cosmobl::ErrorMsg("Error in write_model of Model.h!"); }

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

      /// model: functional form for the 1D model
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
       *  @brief constructor
       *
       *  @param parameters list of model parameters
       *  @param fixed_parameters list of fixed model parameters
       *  @param model model function
       *
       * return object of class Model1D
       */
      Model1D (const vector<Parameter> parameters, const shared_ptr<void> fixed_parameters, const model_function_1D model) :
      Model(parameters, fixed_parameters), m_model(model) { m_2d = 0; }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Model1D () = default;

      ///@}


      /**
       *  @brief evaluate the model for the value xx
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx the variable position
       * 
       *  @return the value of the model at xx
       */
      double operator () (const double xx)
      {
	return m_model(xx, m_fixed_parameters, parameter_values());
      }

      /**
       *  @brief evaluate the model for the value xx
       *  the model uses input parameters
       *
       *  @param xx the variable position
       *  @param parameters vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */
      double operator () (const double xx, const vector<double> parameters)
      {
	return m_model(xx, m_fixed_parameters, parameters);
      }

      /**
       *  @brief evaluate the model for the values xx
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx the variable position vector
       * 
       *  @return the value of the model at xx
       */
      vector<double> operator () (const vector<double> xx)
	{
	  return m_model_vector(xx, m_fixed_parameters, parameter_values());
	}

      /**
       *  @brief evaluate the model for the values xx
       *  the model uses input parameters
       *
       *  @param xx the variable position vector
       *  @param parameters vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */
      vector<double> operator () (const vector<double> xx, const vector<double> parameters)
	{
	  return m_model_vector(xx, m_fixed_parameters, parameters);
	}

      /**
       *  @brief evaluate the model for the value (xx,yy)
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx the X variable position
       *  @param yy the Y variable position
       * 
       *  @return the value of the model at (xx,yy)
       */ 
      double operator () (const double xx, const double yy)
      { (void)xx; (void)yy; cosmobl::ErrorMsg("Error in Model::operator() of Model1D.h!"); return 0; }

      /**
       *  @brief evaluate the model for the value (xx,yy)
       *  the model uses input parameters
       *
       *  @param xx the X variable position
       *  @param yy the Y variable position
       *  @param pars vector containing the model parameters
       * 
       *  @return the value of the model at (xx,yy)
       */
      double operator () (const double xx, const double yy, const vector<double> pars)
      { (void)xx; (void)yy; (void)pars; cosmobl::ErrorMsg("Error in Model::operator() of Model1D.h!"); return 0; }

      /**
       * @brief compute and write the model using the stored 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      void write_model(const vector<double> xx, const string dir_model, const string file_model);

      /**
       * @brief compute and write the model using input 
       * parameter values
       *
       * @param xx vector of point at which the model
       * is computed
       * @param parameters vector of parameters values
       * at which the model is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model_parameters(const vector<double> xx, const vector<double> parameters, const string dir_model, const string file_model);
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
       *  @return object of class Model2D
       */
      Model2D () : Model() { m_2d = 1; } 

      /**
       *  @brief constructor
       *
       *  @param parameters list of model parameters
       *  @param fixed_parameters list of fixed model parameters
       *  @param model model function
       *
       *  @return object of class Model2D
       */
      Model2D (vector<Parameter> parameters, shared_ptr<void> fixed_parameters, model_function_2D model)
	: Model(parameters, fixed_parameters), m_model(model) { m_2d = 1; }

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Model2D () = default;

      ///@}


      /**
       *  @brief evaluate the model for the value (xx,yy)
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx the X variable position
       *  @param yy the Y variable position
       * 
       *  @return the value of the model at (xx,yy)
       */ 
      double operator () (const double xx, const double yy)
      {
	return m_model(xx, yy, m_fixed_parameters, parameter_values());
      }

      /**
       *  @brief evaluate the model for the value (xx,yy)
       *  the model uses input parameters
       *
       *  @param xx the X variable position
       *  @param yy the Y variable position
       *  @param parameters vector containing the model parameters
       * 
       *  @return the value of the model at (xx,yy)
       */
      double operator () (const double xx, const double yy, const vector<double> parameters)
      {
	return m_model(xx, yy, m_fixed_parameters, parameters);
      }

      /**
       *  @brief evaluate the model for the value xx
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx the variable position
       * 
       *  @return the value of the model at xx
       */
      virtual double operator () (const double xx)
      { (void)xx; cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return 0; }

      /**
       *  @brief evaluate the model for the value xx
       *  the model uses input parameters
       *
       *  @param xx the variable position
       *  @param pars vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */
      virtual double operator () (const double xx, const vector<double> pars)
      { (void)xx; (void)pars; cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return 0; }

      /**
       *  @brief evaluate the model for the values xx
       *  the model uses the values of parameters stored 
       *  in the private variable vector<Parameter>
       *
       *  @param xx the variable position vector
       * 
       *  @return the value of the model at xx
       */
      virtual vector<double> operator () (const vector<double> xx)
      { (void)xx; cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return {}; }

      /**
       *  @brief evaluate the model for the values xx
       *  the model uses input parameters
       *
       *  @param xx the variable position vector
       *  @param pars vector containing the model parameters
       * 
       *  @return the value of the model at xx
       */
      virtual vector<double> operator () (const vector<double> xx, const vector<double> pars)
      { (void)xx; (void)pars; cosmobl::ErrorMsg("Error in Model::operator() of Model2D.h!"); return {}; }

      /**
       * @brief compute and write the model using the stored 
       * parameter values
       *
       * @param xx vector of point at which the model is computed,
       * first axis
       * @param yy vector of point at which the model is computed,
       * second axis
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model (const vector<double> xx, const vector<double> yy, const string dir_model, const string file_model);

      /**
       * @brief compute and write the model using input 
       * parameter values
       *
       * @param xx vector of point at which the model is computed
       * @param yy vector of point at which the model is computed,
       * second axis
       * @param parameters vector of parameters values at which the
       * model is computed
       * @param dir_model the output directory of the model
       * @param file_model the name of the file
       *
       * @return none
       */
      virtual void write_model_parameters (const vector<double> xx, const vector<double> yy, const vector<double> parameters, const string dir_model, const string file_model);
    };
  }
}

#endif
