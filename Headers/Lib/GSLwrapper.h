/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/Lib/GSLwrapper.h
 *
 *  @brief functions that wrap GSL routines for integration,
 *  root finding and minimization
 *
 *  This file contains the wrappers of GSL routines for 
 *  integration, root finding and minimization
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __GSLwrap__
#define __GSLwrap__


namespace cosmobl{

  /**
   *  @brief The namespace of the <B> GSL wrappers </B>
   *  
   *  The \e gsl namespace contains all the wrapper functions of the
   *  GSL routines
   */
  namespace gsl {

    struct STR_generic_func_GSL
    {
      function<double(double)> f;
      double xx0;
      function<double(vector<double>)> fmin;
    };
    
    /**
     *  @brief function used to integrate interpolated function 
     *  @param xx the point in which function is defined
     *  @param params the parameters of the function 
     *  @return the value of the function at xx
     */
    double generic_function (const double xx, void *params);

    /**
     *  @brief generic roots
     *  @param xx the point in which function is defined
     *  @param params the parameters of the function 
     *  @return the value of the generic root
     */
    double generic_roots (double xx, void *params);

    /**
     *  @brief generic roots
     *  @param xx the point in which function is defined
     *  @param params the parameters of the function 
     *  @return the value of the generic minimizer
     */
    double generic_minimizer (const gsl_vector * xx, void * params);

    /**
     *  @brief integral, computed using the GSL qag method 
     *  @param Func the GSL function to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @param rule the rule of integration
     *  @return the definite integral of the function
     */
    double GSL_integrate_qag (gsl_function Func, const double a, const double b, const double prec=1.e-2, const int limit_size=1000, const int rule=6);

    /**
     *  @brief integral, computed using the GSL qaws method 
     *  @param Func the GSL function to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param alpha &alpha;
     *  @param beta &beta;
     *  @param mu &mu;
     *  @param nu &nu;
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qaws (gsl_function Func, const double a, const double b, const double alpha=0, const double beta=0, const int mu=1, const int nu =0, const double prec=1.e-2, const int limit_size=1000);

    /**
     *  @brief integral, computed using the GSL qagiu method 
     *  @param Func the GSL function to be integrated
     *  @param a the lower limit of the integral
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @return the integral of the function
     */
    double GSL_integrate_qagiu (gsl_function Func, const double a, const double prec=1.e-2, const int limit_size=1000);

    /**
     *  @brief integral, using the GSL qag method
     *
     *  it only works with function defined as function<double(double)>
     *  that doesn't use fixed parameters (useful for class members,
     *  when the external parameters could be attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @param rule the rule of integration
     *  @return the definite integral of the function
     */
    double GSL_integrate_qag (function<double(double)> func, const double a, const double b, const double prec=1.e-2, const int limit_size=1000, const int rule=6);

    /**
     *  @brief integral, using the GSL qagiu method
     *
     *  it only works with a function defined as
     *  function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param a the lower limit of the integral
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qagiu (function<double(double)> func, const double a, const double prec=1.e-2, const int limit_size=1000);

    /**
     *  @brief integral, using the GSL qaws method 
     *  @param func the function to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param alpha &alpha;
     *  @param beta &beta;
     *  @param mu &mu;
     *  @param nu &nu;
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qaws (function<double(double)> func, const double a, const double b, const double alpha=0, const double beta=0, const int mu=1, const int nu =0, const double prec=1.e-2, const int limit_size=1000);

    /**
     *  @brief integral, using the GSL qag method
     *
     *  it only works with a function defined as
     *  function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param pp a void pointer 
     *  @param par a vector containing the coefficients
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @param rule the rule of integration
     *  @return the definite integral of the function
     */
    double GSL_integrate_qag (function<double(double, shared_ptr<void>, vector<double>)> func, shared_ptr<void> pp, vector<double> par, const double a, const double b, const double prec=1.e-2, const int limit_size=1000, const int rule=6);

    /**
     *  @brief integral, using the GSL qagiu method
     *
     *  it only works with a function defined as
     *  function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param pp a void pointer 
     *  @param par a vector containing the coefficients
     *  @param a the lower limit of the integral
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qagiu (function<double(double, shared_ptr<void>, vector<double>)> func, shared_ptr<void> pp, vector<double> par, const double a, const double prec=1.e-2, const int limit_size=1000);

    /**
     *  @brief integral, using the GSL qag method 
     *  @param func the function to be integrated
     *  @param pp a void pointer 
     *  @param par a vector containing the coefficients
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param alpha &alpha;
     *  @param beta &beta;
     *  @param mu &mu;
     *  @param nu &nu;
     *  @param prec the relative error tolerance
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qaws (function<double(double, shared_ptr<void>, vector<double>)> func, shared_ptr<void> pp, vector<double> par, const double a, const double b, const double alpha=0, const double beta=0, const int mu=1, const int nu =0, const double prec=1.e-2, const int limit_size=1000);

    /**
     *  @brief function to find roots using GSL qag method 
     *  @param Func the GSL function to be integrated
     *  @param low_guess the lower limit 
     *  @param up_guess the upper limit
     *  @param prec the relative error tolerance
     *  @return the function root
     */
    double GSL_root_brent (gsl_function Func, const double low_guess, const double up_guess, const double prec=1.e-3);

    /**
     *  @brief function to find roots using GSL brent method 
     *  @param func the function to be integrated
     *  @param xx0 the value of the zero
     *  @param low_guess the lower limit 
     *  @param up_guess the upper limit
     *  @param prec the relative error tolerance
     *  @return the function root
     */
    double GSL_root_brent (function<double(double)> func, double xx0, const double low_guess, const double up_guess, const double prec=1.e-3);

    /**
     * @brief minimize the provided function using GSL procedure
     *
     * @param func the function to minimize
     * @param start the starting point
     * @param step_size size of the step for minima search
     * @param max_iter maximum number of iteration
     * @param tol tolerance of the minimization
     *
     * @return vector containing the point that minimize the function
     */
    vector<double> GSL_minimize_nD (function<double(vector<double>)> func, const vector<double> start, const vector<double> step_size={}, const unsigned int max_iter=1000, const double tol=1.e-6);

    /**
     * @brief minimize the provided function using GSL procedure
     *
     * @param func the function to minimize
     * @param start the starting point
     * @param min the minimum of the search interval
     * @param max the maximum of the search interval
     * @param max_iter maximum number of iteration
     *
     * @return the point that minimize the function
     */
    double GSL_minimize_1D (function<double(double)> func, const double start, double min=par::defaultDouble, double max=-par::defaultDouble, const int max_iter=1000);

    double GSL_polynomial_eval(const double x, shared_ptr<void> fixed_parameters, const vector<double> coeff);

    void GSL_polynomial_root(const vector<double> coeff, vector<vector<double>> &root);
  }
}

#endif
