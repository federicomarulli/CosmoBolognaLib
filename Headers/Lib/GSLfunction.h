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
 *  @file Headers/Lib/GSLfunction.h
 *
 *  @brief The class GSLfunction 
 *
 *  This file defines the interface of the class GSLfunction, used to
 *  wrap gsl minimization function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __GSLfunc__
#define __GSLfunc__

#include "Func.h"

namespace cosmobl {

  /**
   * @var typedef func_1par_1
   * @brief definition of 1D function
   * that takes as an input the value at which the function
   * is computed and a pointer to fixed parameters
   */
  typedef function<double(double, shared_ptr<void>)> func_1par_1;

  /**
   * @var typedef func_1par_2
   * @brief definition of a one parameter function
   * that takes as an input the value at which the function
   * is computed a pointer to fixed parameters and a vector of free parameters 
   */
  typedef function<double(double, shared_ptr<void>, vector<double> )> func_1par_2;

  /**
   * @var typedef func_npar_1
   * @brief definition of a one parameter function
   * that takes as an input the values at which the function
   * is computed and a pointer to fixed parameters
   */
  typedef function<double(vector<double>, shared_ptr<void>)> func_npar_1;

   /**
   * @var typedef func_npar_2
   * @brief definition of a one parameter function
   * that takes as an input the values at which the function
   * is computed a pointer to fixed parameters and a vector of free parameters 
   */ 
  typedef function<double(vector<double>, shared_ptr<void>, vector<double> )> func_npar_2;

  /**
   *  @class GSLfunction GSLfunction.h "Headers/Lib/GSLfunction.h"
   *
   *  @brief The class GSLfunction
   *
   *  This class is used to wrap object of type GSL_function, ad use them in
   *  minimizing procedures
   */
  class GSLfunction{

    public:

      /**
       *  @brief static factory used to construct GSLfunction of type
       *  GSLfunction_1D_1
       *
       *  @param function object of type func_1par_1
       *
       *  @param model_parameters pointer to a container of parameters 
       *  for the function
       *
       *  @return a pointer to an object of class GSLfunction of
       *  a given type
       */
      static unique_ptr<GSLfunction> make_GSLfunction(func_1par_1 function , shared_ptr<void> model_parameters=NULL );

      /**
       *  @brief static factory used to construct GSLfunction of type
       *  GSLfunction_nD_1
       *
       *  @param npar number of function parameters
       *
       *  @param function object of type func_npar_1
       *
       *  @param model_parameters pointer to a container of parameters 
       *  for the function
       *
       *  @return a pointer to an object of class GSLfunction of
       *  a given type
       */    
      static unique_ptr<GSLfunction> make_GSLfunction(int npar, func_npar_1 function, shared_ptr<void> model_parameters=NULL);

      /**
       *  @brief static factory used to construct GSLfunction of type
       *  GSLfunction_1D_2
       *
       *  @param function object of type func_1par_2
       *
       *  @param params vector containing free parameters for the function
       *
       *  @param model_parameters pointer to a container of parameters 
       *  for the function
       *
       *  @return a pointer to an object of class GSLfunction of
       *  a given type
       */
      static unique_ptr<GSLfunction> make_GSLfunction(func_1par_2 function, vector<double> params, shared_ptr<void> model_parameters=NULL);

      /**
       *  @brief static factory used to construct GSLfunction of type
       *  GSLfunction_nD_2
       *
       *  @param npar number of function parameters
       *
       *  @param function object of type func_npar_2
       *
       *  @param params vector containing free parameters for the function
       *
       *  @param model_parameters pointer to a container of parameters 
       *  for the function
       *
       *  @return a pointer to an object of class GSLfunction of
       *  a given type
       */    
      static unique_ptr<GSLfunction> make_GSLfunction(int npar, func_npar_2 function, vector<double> params, shared_ptr<void> model_parameters=NULL);

      /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the value of the variable at which the function
       * has its minimum
       *
       * @param max_iter maximum number of iteration
       *
       * @param min lower bound for the minimization procedure
       *
       * @param max upper bound for the minimization procedure
       *
       * @return none
       */
      virtual void minimize(double &result, int max_iter=100, double min=-1.e30, double max=1.e30)
      {cosmobl::ErrorMsg ("Error in minimize of GSLfunction!");};
      
       /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the values of the variables at which the function
       * has its minimum
       *
       * @param max_iter maximum number of iteration
       *
       * @param tol tolerance of the minimization
       *
       * @return none
       */    
      virtual void minimize(vector<double> &result, unsigned int max_iter=100, double tol=1.e-6) 
      {cosmobl::ErrorMsg ("Error in minimize of GSLfunction!");};

       /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the values of the variables at which the function
       * has its minimum
       *
       * @param step_size size of the step for minima search
       *
       * @param max_iter maximum number of iteration
       *
       * @param tol tolerance of the minimization
       *
       * @return none
       */  
      virtual void minimize(vector<double> &result, vector<double> &step_size, unsigned int max_iter=100, double tol=1.e-6) 
      {cosmobl::ErrorMsg ("Error in minimize of GSLfunction!");};
  };

  /**
   *  @class GSLfunction_1D_1 GSLfunction.h "Headers/Lib/GSLfunction.h"
   *
   *  @brief The class GSLfunction_1D_1
   *
   *  This class is used to wrap object of type GSL_function, ad use them in
   *  minimizing procedures. It works with function of type func_1par_1
   */
  class GSLfunction_1D_1 : public GSLfunction
  {
    private:

      /// provided function of type func_1par_1 
      func_1par_1 m_function;

      /// fixed parameters of the function
      shared_ptr<void> m_function_parameters;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class GSLfunction_1D_1
       */
      GSLfunction_1D_1 () {}

      /**
       *  @brief default constructor of the class GSLfunction_1D_1
       *  
       *  @param function function of type func_1par_1
       *  
       *  @param function_parameters fixed parameters of the function
       *
       *  @return object of class GSLfunction_1D_1
       */
      GSLfunction_1D_1 (func_1par_1 function, shared_ptr<void> function_parameters = NULL) :
	m_function(function), m_function_parameters(function_parameters) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~GSLfunction_1D_1 () {}

      ///@}

      /**
       * @brief set the function 
       * @param function function of type func_1par_1
       * return none
       */
      void set_function(func_1par_1 function) {m_function = function;}

      /**
       * @brief set the function parameters
       * @param function_parameters pointer to the function parameters
       * return none
       */
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the value of the variable at which the function
       * has its minimum
       *
       * @param max_iter maximum number of iteration
       *
       * @param min lower bound for the minimization procedure
       *
       * @param max upper bound for the minimization procedure
       *
       * @return none
       */
      void minimize(double &result, int max_iter = 100, double min=-1.e30, double max=1.e30);
  };

  /**
   *  @class GSLfunction_1D_2 GSLfunction.h "Headers/Lib/GSLfunction.h"
   *
   *  @brief The class GSLfunction_1D_2
   *
   *  This class is used to wrap object of type GSL_function, ad use them in
   *  minimizing procedures. It works with function of type func_1par_2
   */
  class GSLfunction_1D_2 : public GSLfunction
  {
    private:

      /// provided function of type func_1par_2 
      func_1par_2 m_function;

      /// free parameters of the function
      vector<double> m_parameters;

      /// fixed parameters of the function
      shared_ptr<void> m_function_parameters;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class GSLfunction_1D_2
       */
      GSLfunction_1D_2 () {}

      /**
       *  @brief default constructor of the class GSLfunction_1D_2
       *  
       *  @param function function of type func_1par_2
       *  
       *  @param parameters vector containing function parameters
       *
       *  @param function_parameters fixed parameters of the function
       *
       *  @return object of class GSLfunction_1D_2
       */
      GSLfunction_1D_2 (func_1par_2 function, vector<double> parameters, shared_ptr<void> function_parameters = NULL) :
	m_function(function), m_parameters(parameters), m_function_parameters(function_parameters) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~GSLfunction_1D_2 () {}

      ///@}

      /**
       * @brief set the function 
       * @param function function of type func_1par_2
       * return none
       */
      void set_function(func_1par_2 function) {m_function = function;}

      /**
       * @brief set the function parameters
       * @param parameters vector containing the function parameters
       * return none
       */
      void set_parameters(vector<double> parameters) {m_parameters = parameters;}

      /**
       * @brief set the function parameters
       * @param function_parameters pointer to the function parameters
       * return none
       */
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the value of the parameter at which the function
       * has its minimum
       *
       * @param max_iter maximum number of iteration
       *
       * @param min lower bound for the minimization procedure
       *
       * @param max upper bound for the minimization procedure
       *
       * @return none
       */
      void minimize(double &result, int max_iter = 100, double min=-1.e30, double max=1.e30);
  };

  /**
   *  @class GSLfunction_nD_1 GSLfunction.h "Headers/Lib/GSLfunction.h"
   *
   *  @brief The class GSLfunction_nD_1
   *
   *  This class is used to wrap object of type GSL_function, ad use them in
   *  minimizing procedures. It works with function of type func_npar_1
   */
  class GSLfunction_nD_1 : public GSLfunction
  {
    private:

      /// number of free parameters of the function
      int m_npar;

      /// provided function of type func_1par_2 
      func_npar_1 m_function;
      
      /// fixed parameters of the function
      shared_ptr<void> m_function_parameters;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class GSLfunction_nD_1
       */
      GSLfunction_nD_1 () {}


      /**
       *  @brief default constructor of the class GSLfunction_nD_1
       *  
       *  @param npar the number of variables of the function
       *
       *  @param function function of type func_npar_1
       *
       *  @param function_parameters fixed parameters of the function
       *
       *  @return object of class GSLfunction_nD_1
       */
      GSLfunction_nD_1 (int npar, func_npar_1 function, shared_ptr<void> function_parameters = NULL) :
	m_npar(npar), m_function(function), m_function_parameters(function_parameters) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~GSLfunction_nD_1 () {}

      ///@}

      /**
       * @brief set the function 
       * @param function function of type func_npar_1
       * return none
       */
      void set_function(func_npar_1 function) {m_function = function;}

      /**
       * @brief set the function parameters
       * @param function_parameters pointer to the function parameters
       * return none
       */
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

       /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the values of the parameters at which the function
       * has its minimum
       *
       * @param max_iter maximum number of iteration
       *
       * @param tol tolerance of the minimization
       *
       * @return none
       */ 
      void minimize(vector<double> &result, unsigned int max_iter = 100, double tol=1.e-6); 

      /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the values of the parameters at which the function
       * has its minimum
       *
       * @param step_size size of the step for minima search
       *
       * @param max_iter maximum number of iteration
       *
       * @param tol tolerance of the minimization
       *
       * @return none
       */  
      void minimize(vector<double> &result, vector<double> step_size, unsigned int max_iter = 100, double tol=1.e-6); 
  };

  /**
   *  @class GSLfunction_nD_2 GSLfunction.h "Headers/Lib/GSLfunction.h"
   *
   *  @brief The class GSLfunction_nD_2
   *
   *  This class is used to wrap object of type GSL_function, ad use them in
   *  minimizing procedures. It works with function of type func_npar_2
   */
  class GSLfunction_nD_2 : public GSLfunction
  {
    private:

      /// number of free parameters of the function
      int m_npar;

      /// provided function of type func_1par_2 
      func_npar_2 m_function;

      /// free parameters of the function
      vector<double> m_parameters;

      /// fixed parameters of the function
      shared_ptr<void> m_function_parameters;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{
      
      /**
       *  @brief default constructor
       *  @return object of class GSLfunction_nD_1
       */
      GSLfunction_nD_2 () {}

      /**
       *  @brief default constructor of the class GSLfunction_nD_2
       *  
       *  @param npar the number of variables of the function
       *
       *  @param function function of type func_npar_2
       *
       *  @param parameters vetor containing parameters of the function
       *  
       *  @param function_parameters fixed parameters of the function
       *
       *  @return object of class GSLfunction_nD_2
       */
      GSLfunction_nD_2 (int npar, func_npar_2 function, vector<double> parameters, shared_ptr<void> function_parameters = NULL) :
	m_npar(npar), m_function(function), m_parameters(parameters), m_function_parameters(function_parameters) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~GSLfunction_nD_2 () {}

      ///@}
      
      /**
       * @brief set the function 
       * @param function function of type func_npar_2
       * return none
       */
      void set_function(func_npar_2 function)  {m_function = function;}

      /**
       * @brief set the function parameters
       * @param parameters vector containing the function parameters
       * return none
       */
      void set_parameters(vector<double> parameters) {m_parameters = parameters;}

      /**
       * @brief set the function parameters
       * @param function_parameters pointer to the function parameters
       * return none
       */
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the values of the parameters at which the function
       * has its minimum
       *
       * @param max_iter maximum number of iteration
       *
       * @param tol tolerance of the minimization
       *
       * @return none
       */ 
      void minimize(vector<double> &result, unsigned int max_iter = 100, double tol=1.e-6); 

      /**
       * @brief minimize the provided function using GSL procedure
       *
       * @param result store the values of the parameters at which the function
       * has its minimum
       *
       * @param step_size size of the step for minima search
       *
       * @param max_iter maximum number of iteration
       *
       * @param tol tolerance of the minimization
       *
       * @return none
       */  
      void minimize(vector<double> &result, vector<double> step_size, unsigned int max_iter = 100, double tol=1.e-6); 
  };
}

#endif
