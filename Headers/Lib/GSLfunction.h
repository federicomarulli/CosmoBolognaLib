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

  typedef function<double(double, shared_ptr<void>)> func_1par_1;
  typedef function<double(double, shared_ptr<void>, vector<double> )> func_1par_2;
  typedef function<double(vector<double>, shared_ptr<void>)> func_npar_1;
  typedef function<double(vector<double>, shared_ptr<void>, vector<double> )> func_npar_2;

  class GSLfunction{
    public:
      static unique_ptr<GSLfunction> make_GSLfunction(func_1par_1 , shared_ptr<void> model_parameters=NULL );
      static unique_ptr<GSLfunction> make_GSLfunction(int, func_npar_1 , shared_ptr<void> model_parameters=NULL);
      static unique_ptr<GSLfunction> make_GSLfunction(func_1par_2 , vector<double>, shared_ptr<void> model_parameters=NULL);
      static unique_ptr<GSLfunction> make_GSLfunction(int, func_npar_2 ,vector<double> ,shared_ptr<void> model_parameters=NULL);

      virtual void minimize(double &result, int max_iter=100, double min=-1.e30, double max=1.e30) {cosmobl::ErrorMsg ("Error in minimize of GSLfunction!");};
      virtual void minimize(vector<double> &result, unsigned int max_iter=100, double tol=1.e-6) {cosmobl::ErrorMsg ("Error in minimize of GSLfunction!");};
      virtual void minimize(vector<double> &result,vector<double> &step_size,unsigned int max_iter=100, double tol=1.e-6) {cosmobl::ErrorMsg ("Error in minimize of GSLfunction!");};
  };


  class GSLfunction_1D_1 : public GSLfunction
  {
    private:
      func_1par_1 m_function;
      shared_ptr<void> m_function_parameters;

    public:
      GSLfunction_1D_1 () {}
      ~GSLfunction_1D_1 () {}

      GSLfunction_1D_1 (func_1par_1 function, shared_ptr<void> function_parameters = NULL) :
	m_function(function), m_function_parameters(function_parameters) {}

      void set_function(func_1par_1 function) {m_function = function;}
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      void minimize(double &result, int max_iter = 100, double min=-1.e30, double max=1.e30);
  };

  class GSLfunction_1D_2 : public GSLfunction
  {
    private:
      func_1par_2 m_function;
      vector<double> m_parameters;
      shared_ptr<void> m_function_parameters;

    public:
      GSLfunction_1D_2 () {}
      ~GSLfunction_1D_2 () {}

      GSLfunction_1D_2 (func_1par_2 function, vector<double> parameters, shared_ptr<void> function_parameters = NULL) :
	m_function(function), m_parameters(parameters), m_function_parameters(function_parameters) {}

      void set_function(func_1par_2 function) {m_function = function;}
      void set_parameters(vector<double> parameters) {m_parameters = parameters;}
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      void minimize(double &result, int max_iter = 100, double min=-1.e30, double max=1.e30);
  };

  class GSLfunction_nD_1 : public GSLfunction
  {
    private:
      int m_npar;
      func_npar_1 m_function;
      shared_ptr<void> m_function_parameters;

    public:
      GSLfunction_nD_1 () {}
      ~GSLfunction_nD_1 () {}

      GSLfunction_nD_1 (int npar, func_npar_1 function, shared_ptr<void> function_parameters = NULL) :
	m_npar(npar), m_function(function), m_function_parameters(function_parameters) {}

      void set_function(func_npar_1 function) {m_function = function;}
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      void minimize(vector<double> &result, unsigned int max_iter = 100, double tol=1.e-6); 
      void minimize(vector<double> &result, vector<double> step_size, unsigned int max_iter = 100, double tol=1.e-6); 
  };

  class GSLfunction_nD_2 : public GSLfunction
  {
    private:
      int m_npar;
      func_npar_2 m_function;
      vector<double> m_parameters;
      shared_ptr<void> m_function_parameters;

    public:
      GSLfunction_nD_2 () {}
      ~GSLfunction_nD_2 () {}

      GSLfunction_nD_2 (int npar, func_npar_2 function, vector<double> parameters, shared_ptr<void> function_parameters = NULL) :
	m_npar(npar), m_function(function), m_parameters(parameters), m_function_parameters(function_parameters) {}

      void set_function(func_npar_2 function) {m_function = function;}
      void set_parameters(vector<double> parameters) {m_parameters = parameters;}
      void set_function_parameters(shared_ptr<void> function_parameters) {m_function_parameters = function_parameters;}

      void minimize(vector<double> &result, unsigned int max_iter = 100, double tol=1.e-6); 
      void minimize(vector<double> &result, vector<double> step_size, unsigned int max_iter = 100, double tol=1.e-6); 
  };
}

#endif
