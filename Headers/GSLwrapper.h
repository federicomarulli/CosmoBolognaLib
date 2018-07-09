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
 *  @file Headers/GSLwrapper.h
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

#include "Kernel.h"

namespace cbl {

  /**
   *  @brief The namespace of the <B> GSL wrappers </B>
   *  
   *  The \e gsl namespace contains all the wrapper functions of the
   *  GSL routines
   */
  namespace gsl {

    struct STR_generic_func_GSL
    {
      FunctionDoubleDouble f;
      double xx0;

      FunctionDoubleVector fmin;

      FunctionDoubleVectorRef fmin_return;
      std::vector<double> parameters_return;
    };

    /**
     *  @brief Function used to check output of the wrapped
     *  GSL routines. 
     *  @param status the GSL routine status
     *  @param exit true &rarr; if the routine raise an error, exit the code; 
     *  false &rarr; show the status of the GSL routine without exiting
     *  @param CBLfunction the name of the function using the GSL routine 
     *  @param GSLroutine the name of the gsl routine
     *  @return none
     */
    void check_GSL_fail (const int status, const bool exit, const std::string CBLfunction, const std::string GSLroutine);
    
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
     *  @brief generic roots
     *  @param xx the point in which function is defined
     *  @param params the parameters of the function 
     *  @return the value of the generic minimizer
     */
    double generic_minimizer_return (const gsl_vector * xx, void * params);
    
    /**
     *  @brief the derivative of a function
     *
     *  This function computes the numerical derivative of the
     *  function Func at the point xx using an adaptive central
     *  difference algorithm with a step size of hh.
     *
     *  The initial value of hh is used to estimate an optimal
     *  step-size, based on the scaling of the truncation error and
     *  round-off error in the derivative calculation. The derivative
     *  is computed using a 5-point rule for equally spaced abscissae
     *  at xx-hh, xx-hh/2, xx, xx+hh/2, xx+hh, with an error estimate
     *  taken from the difference between the 5-point rule and the
     *  corresponding 3-point rule xx-hh, xx, xx+hh. Note that the
     *  value of the function at xx does not contribute to the
     *  derivative calculation, so only 4-points are actually used.
     *  (from the GSL documentation)
     *
     *  @param Func the GSL function to be derived 
     *
     *  @param xx point at which the derivative is computed
     *
     *  @param hh the initial value of the step size
     *
     *  @param prec the relative error tolerance
     *
     *  @return the definite integral of the function
     */
    double GSL_derivative (gsl_function Func, const double xx, const double hh, const double prec=1.e-2);
    
    /**
     *  @brief integral, using the GSL cquad method
     *
     *  it only works with a function defined as
     *  std::function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param Func the fuction to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param nevals the number of intervals
     *  @return the definite integral of the function
     */
    double GSL_integrate_cquad (gsl_function Func, const double a, const double b, const double rel_err=1.e-2, const double abs_err=1.e-6, const int nevals=100);

    /**
     *  @brief integral, computed using the GSL qag method 
     *  @param Func the GSL function to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @param rule the rule of integration
     *  @return the definite integral of the function
     */
    double GSL_integrate_qag (gsl_function Func, const double a, const double b, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000, const int rule=6);

    /**
     *  @brief integral, computed using the GSL qaws method 
     *  @param Func the GSL function to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param alpha &alpha;
     *  @param beta &beta;
     *  @param mu &mu;
     *  @param nu &nu;
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qaws (gsl_function Func, const double a, const double b, const double alpha=0, const double beta=0, const int mu=1, const int nu =0, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000);

    /**
     *  @brief integral, computed using the GSL qagiu method 
     *  @param Func the GSL function to be integrated
     *  @param a the lower limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @return the integral of the function
     */
    double GSL_integrate_qagiu (gsl_function Func, const double a, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000);

    /**
     *  @brief the derivative of a function
     *
     *  This function computes the numerical derivative of the
     *  function Func at the point xx using an adaptive central
     *  difference algorithm with a step size of hh.
     *
     *  The initial value of hh is used to estimate an optimal
     *  step-size, based on the scaling of the truncation error and
     *  round-off error in the derivative calculation. The derivative
     *  is computed using a 5-point rule for equally spaced abscissae
     *  at xx-hh, xx-hh/2, xx, xx+hh/2, xx+hh, with an error estimate
     *  taken from the difference between the 5-point rule and the
     *  corresponding 3-point rule xx-hh, xx, xx+hh. Note that the
     *  value of the function at xx does not contribute to the
     *  derivative calculation, so only 4-points are actually used.
     *  (from the GSL documentation)
     *
     *  It only works with function defined as std::function<double(double)>
     *  that doesn't use fixed parameters (useful for class members,
     *  when the external parameters could be attributes of the class)
     *
     *  @param func the function to be derived 
     *
     *  @param xx point at which the derivative is computed
     *
     *  @param hh the initial value of the step size
     *
     *  @param prec the relative error tolerance
     *
     *  @return the definite integral of the function
     */
    double GSL_derivative (FunctionDoubleDouble func, const double xx, const double hh, const double prec=1.e-2);

    /**
     *  @brief integral, using the GSL cquad method
     *
     *  it only works with function defined as std::function<double(double)>
     *  that doesn't use fixed parameters (useful for class members,
     *  when the external parameters could be attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param nevals the number of intervals
     *  @return the definite integral of the function
     */
    double GSL_integrate_cquad (FunctionDoubleDouble func, const double a, const double b, const double rel_err=1.e-2, const double abs_err=1.e-6, const int nevals=100);

    /**
     *  @brief integral, using the GSL qag method
     *
     *  it only works with function defined as std::function<double(double)>
     *  that doesn't use fixed parameters (useful for class members,
     *  when the external parameters could be attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @param rule the rule of integration
     *  @return the definite integral of the function
     */
    double GSL_integrate_qag (FunctionDoubleDouble func, const double a, const double b, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000, const int rule=6);

    /**
     *  @brief integral, using the GSL qagiu method
     *
     *  it only works with a function defined as
     *  std::function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param a the lower limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qagiu (FunctionDoubleDouble func, const double a, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000);

    /**
     *  @brief integral, using the GSL qaws method 
     *  @param func the function to be integrated
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param alpha &alpha;
     *  @param beta &beta;
     *  @param mu &mu;
     *  @param nu &nu;
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qaws (FunctionDoubleDouble func, const double a, const double b, const double alpha=0, const double beta=0, const int mu=1, const int nu =0, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000);

    /**
     *  @brief integral, using the GSL cquad method
     *
     *  it only works with a function defined as
     *  std::function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param pp a void pointer 
     *  @param par a vector containing the coefficients
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param nevals the number of intervals
     *  @return the definite integral of the function
     */
    double GSL_integrate_cquad (FunctionDoubleDoublePtrVectorRef func, std::shared_ptr<void> pp, std::vector<double> par, const double a, const double b, const double rel_err=1.e-2, const double abs_err=1.e-6, const int nevals=100);

    /**
     *  @brief integral, using the GSL qag method
     *
     *  it only works with a function defined as
     *  std::function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param pp a void pointer 
     *  @param par a vector containing the coefficients
     *  @param a the lower limit of the integral
     *  @param b the upper limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @param rule the rule of integration
     *  @return the definite integral of the function
     */
    double GSL_integrate_qag (FunctionDoubleDoublePtrVectorRef func, std::shared_ptr<void> pp, std::vector<double> par, const double a, const double b, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000, const int rule=6);

    /**
     *  @brief integral, using the GSL qagiu method
     *
     *  it only works with a function defined as
     *  std::function<double(double)> that doesn't use fixed parameters
     *  (useful for class members, when the external parameters can be
     *  attributes of the class)
     *
     *  @param func the fuction to be integrated
     *  @param pp a void pointer 
     *  @param par a vector containing the coefficients
     *  @param a the lower limit of the integral
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qagiu (FunctionDoubleDoublePtrVectorRef func, std::shared_ptr<void> pp, std::vector<double> par, const double a, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000);

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
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @param limit_size the maximum size of workspace
     *  @return the definite integral of the function
     */
    double GSL_integrate_qaws (FunctionDoubleDoublePtrVectorRef func, std::shared_ptr<void> pp, std::vector<double> par, const double a, const double b, const double alpha=0, const double beta=0, const int mu=1, const int nu =0, const double rel_err=1.e-2, const double abs_err=1.e-6, const int limit_size=1000);

    /**
     *  @brief function to find roots using GSL qag method 
     *  @param Func the GSL function to be integrated
     *  @param low_guess the lower limit 
     *  @param up_guess the upper limit
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @return the function root
     */
    double GSL_root_brent (gsl_function Func, const double low_guess, const double up_guess, const double rel_err=1.e-3, const double abs_err=1.e-6);

    /**
     *  @brief function to find roots using GSL brent method 
     *  @param func the function to be integrated
     *  @param xx0 the value of the zero
     *  @param low_guess the lower limit 
     *  @param up_guess the upper limit
     *  @param rel_err the relative error
     *  @param abs_err the absolute error
     *  @return the function root
     */
    double GSL_root_brent (FunctionDoubleDouble func, double xx0, const double low_guess, const double up_guess, const double rel_err=1.e-3, const double abs_err=1.e-6);

    /**
     * @brief minimize the provided function using GSL procedure
     *
     * @param func the function to minimize
     * @param start the starting point
     * @param ranges limits for the parameters
     * @param max_iter maximum number of iteration
     * @param tol tolerance of the minimization
     * @param epsilon percentage of the range to be used
     * to define the step size
     *
     * @return vector containing the point that minimize the function
     */
    std::vector<double> GSL_minimize_nD (FunctionDoubleVector func, const std::vector<double> start, const std::vector<std::vector<double>> ranges, const unsigned int max_iter=1000, const double tol=1.e-6, const double epsilon=0.1);

    /**
     * @brief minimize the provided function using GSL procedure
     *
     * @param func the function to minimize
     * @param start the starting point
     * @param ranges limits for the parameters
     * @param max_iter maximum number of iteration
     * @param tol tolerance of the minimization
     * @param epsilon percentage of the range to be used
     * to define the step size
     *
     * @return vector containing the point that minimize the function
     */
    std::vector<double> GSL_minimize_nD (FunctionDoubleVectorRef func, const std::vector<double> start, const std::vector<std::vector<double>> ranges, const unsigned int max_iter=1000, const double tol=1.e-6, const double epsilon=0.1);

    /**
     * @brief minimize the provided function using GSL procedure
     *
     * @param func the function to minimize
     * @param start the starting point
     * @param min the minimum of the search interval
     * @param max the maximum of the search interval
     * @param max_iter maximum number of iteration
     * @param verbose show output
     *
     * @return the point that minimize the function
     */
    double GSL_minimize_1D (FunctionDoubleDouble func, const double start, double min=par::defaultDouble, double max=-par::defaultDouble, const int max_iter=1000, const bool verbose=false);

    /**
     * @brief evaluate a polynomial 
     *
     * @param x the independent variable
     * @param fixed_parameters fixed parameters of 
     * the polynomial
     * @param coeff polynomial coefficients
     *
     * @return the value of the polynomial at x
     */
    double GSL_polynomial_eval (const double x, const std::shared_ptr<void> fixed_parameters, const std::vector<double> coeff);

    /**
     * @brief find polynomial roots
     *
     * @param coeff polynomial coefficients 
     * @param root polynomial roots
     *
     * @return none
     */
    void GSL_polynomial_root (const std::vector<double> coeff, std::vector<std::vector<double>> &root);
  }
}

#endif
