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
 *  @file Headers/Models/ModelCosmology.h
 *
 *  @brief The class ModelCosmology
 *
 *  This file defines the interface of the class ModelCosmology
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELCOSM__
#define __MODELCOSM__

#include "Model.h"

namespace cosmobl {
  
  /**
   *  @class ModelCosmology ModelCosmology.h "Headers/Lib/ModelCosmology.h"
   *
   *  @brief The class ModelCosmology
   *
   *  This class is used to model cosmological parameters
   */
  class ModelCosmology : public Model
  {

  protected:
    
    /// number of model parameters
    vector<cosmobl::CosmoPar> m_cosmological_parameters;

    
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
    ModelCosmology () {}

    /**
     *  @brief default destructor
     *
     *  @return none
     */
    virtual ~Model() {}

    ///@}
    
    virtual double operator () (const double xx)
    { cosmobl::ErrorMsg("Error in Model::operator() of Model.h!"); return 0; }

    virtual double operator () (const double xx, const vector<double> pars)
    { cosmobl::ErrorMsg("Error in Model::operator() of Model.h!"); return 0; }

    virtual double operator () (const double xx, const double yy)
    { cosmobl::ErrorMsg("Error in Model::operator() of Model.h!"); return 0; }

    virtual double operator () (const double xx, const double yy, const vector<double> pars)
    { cosmobl::ErrorMsg("Error in Model::operator() of Model.h!"); return 0; }

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
	m_npar_eff = (pp->freeze()) ? m_npar_eff : m_npar_eff+1;

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
     *  @brief Update parameter and return the current parameter values
     *
     *  @return the values of model parameters
     */
    double update_parameters(const double); 

    /**
     *  @brief Update parameters and return a vector containing the
     *  current parameter values
     *
     *  @return the values of model parameters
     */
    vector<double> update_parameters (const vector<double>); 
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
    Model1D (const vector<shared_ptr<Parameter> > parameters, const shared_ptr<void> model_parameters, const model_function_1D model) :
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
    Model2D (vector<shared_ptr<Parameter> > parameters, shared_ptr<void> model_parameters, model_function_2D model)
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
  };

}

#endif
