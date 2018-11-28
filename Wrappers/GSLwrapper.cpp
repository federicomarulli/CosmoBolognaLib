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
 *  @file Wrappers/GSLwrapper.cpp
 *
 *  @brief functions that wrap GSL routines for integration,
 *  root finding and minimization
 *
 *  This file contains the implementation of 
 *  wrappers of  GSL routines for integration, 
 *  root finding and minimization
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "GSLwrapper.h"

using namespace std;

using namespace cbl;
using namespace gsl;

// ============================================================================


void cbl::gsl::check_GSL_fail (const int status, const bool exit, const std::string CBLfunction, const std::string GSLroutine) 
{
  if (exit) {
    if (status) {
      ErrorCBL((GSLroutine != par::defaultString) ? "Error in the gsl routine "+GSLroutine+" used in "+CBLfunction+":"+string(gsl_strerror(status)): "Error in "+CBLfunction+":"+string(gsl_strerror(status)));
    }
  }
  else 
    WarningMsg("The gsl routine "+GSLroutine+" used in "+CBLfunction+" exited with status "+string(gsl_strerror(status)));
}


// ============================================================================


double cbl::gsl::generic_function (const double xx, void *params)
{
  gsl::STR_generic_func_GSL *pp = (gsl::STR_generic_func_GSL *) params;
  return pp->f(xx);
}


// ============================================================================


double cbl::gsl::generic_roots (double xx, void *params)
{
  gsl::STR_generic_func_GSL *pp = (gsl::STR_generic_func_GSL *) params;
  return pp->f(xx)-pp->xx0;
}


// ============================================================================


double cbl::gsl::generic_minimizer (const gsl_vector * xx, void * params)
{
  vector<double> _xx;
  for (size_t i=0; i<xx->size; i++)
    _xx.push_back(gsl_vector_get(xx, i));

  gsl::STR_generic_func_GSL *pp = (gsl::STR_generic_func_GSL *) params;
  pp->parameters_return = _xx;

  return pp->fmin(_xx);
}


// ============================================================================


double cbl::gsl::generic_minimizer_return (const gsl_vector * xx, void * params)
{
  vector<double> _xx;
  for (size_t i=0; i<xx->size; i++)
    _xx.push_back(gsl_vector_get(xx, i));

  gsl::STR_generic_func_GSL *pp = (gsl::STR_generic_func_GSL *) params;

  double val = pp->fmin_return(_xx);
  pp->parameters_return = _xx;
  
  return val;
}


// ============================================================================


double cbl::gsl::GSL_derivative (gsl_function Func, const double xx, const double hh, const double prec)
{
  gsl_set_error_handler_off();

  double Deriv, error;

  int status = gsl_deriv_central(&Func, xx, hh, &Deriv, &error);
  check_GSL_fail(status, true, "GSL_derivative", "gsl_deriv_central");

  return (Deriv!=0 && error/Deriv<prec) ? Deriv : ErrorCBL("Error in cosmobl::gsl::GSL_derivative of GSLwrapper! error/Deriv = "+conv(error/Deriv, par::fDP6)+" > prec = "+conv(prec, par::fDP3));
}


// ============================================================================


double cbl::gsl::GSL_integrate_cquad (gsl_function Func, const double a, const double b, const double rel_err, const double abs_err, const int nevals)
{
  gsl_set_error_handler_off();
  double Int, error;
  size_t nn = nevals;

  gsl_integration_cquad_workspace *ww = gsl_integration_cquad_workspace_alloc(nevals);
  
  int status = gsl_integration_cquad(&Func, a, b, abs_err, rel_err, ww, &Int, &error, &nn); 

  gsl_integration_cquad_workspace_free(ww);

  check_GSL_fail(status, true, "GSL_integrate_cquad", "gsl_integration_cquad");
  
  return Int;
}


// ============================================================================


double cbl::gsl::GSL_integrate_qag (gsl_function Func, const double a, const double b, const double rel_err, const double abs_err, const int limit_size, const int rule)
{
  gsl_set_error_handler_off();

  double Int, error;
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);

  int status = gsl_integration_qag(&Func, a, b, abs_err, rel_err, limit_size, rule, ww, &Int, &error); 

  check_GSL_fail(status, true, "GSL_integrate_qag", "gsl_integrate_qag");
  
  gsl_integration_workspace_free(ww);
  
  return Int;
}


// ============================================================================


double cbl::gsl::GSL_integrate_qagiu (gsl_function Func, const double a, const double rel_err, const double abs_err, const int limit_size)
{
  gsl_set_error_handler_off();

  double Int, error;
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
  
  int status = gsl_integration_qagiu(&Func, a, abs_err, rel_err, limit_size, ww, &Int, &error); 

  check_GSL_fail(status, true, "GSL_integrate_qagiu", "gsl_integrate_qagiu");

  gsl_integration_workspace_free(ww);
  
  return Int;
}


// ============================================================================


double cbl::gsl::GSL_integrate_qaws (gsl_function Func, const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double rel_err, const double abs_err, const int limit_size)
{
  gsl_set_error_handler_off();

  gsl_integration_qaws_table *T = gsl_integration_qaws_table_alloc(alpha, beta, mu, nu);

  double Int, error;
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);

  int status = gsl_integration_qaws (&Func, a, b, T, abs_err, rel_err, limit_size, ww, &Int, &error); 

  check_GSL_fail(status, true, "GSL_integrate_qaws", "gsl_integrate_qaws");

  gsl_integration_workspace_free(ww);
  gsl_integration_qaws_table_free(T);
  
  return Int;
}


// ============================================================================


double cbl::gsl::GSL_derivative (FunctionDoubleDouble func, const double xx, const double hh, const double prec)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;
  Func.function = generic_function;
  Func.params = &params;

  return GSL_derivative(Func, xx, hh, prec);
}


// ============================================================================


double cbl::gsl::GSL_integrate_cquad (FunctionDoubleDouble func, const double a, const double b, const double rel_err, const double abs_err, const int nevals)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;
  
  return GSL_integrate_cquad(Func, a, b, rel_err, abs_err, nevals);
}


// ============================================================================


double cbl::gsl::GSL_integrate_qag (FunctionDoubleDouble func, const double a, const double b, const double rel_err, const double abs_err, const int limit_size, const int rule)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;
  
  return GSL_integrate_qag(Func, a, b, rel_err, abs_err, limit_size, rule);
}


// ============================================================================


double cbl::gsl::GSL_integrate_qagiu (FunctionDoubleDouble func, const double a,  const double rel_err, const double abs_err, const int limit_size)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;

  return GSL_integrate_qagiu(Func, a, rel_err, abs_err, limit_size);
}


// ============================================================================


double cbl::gsl::GSL_integrate_qaws (FunctionDoubleDouble func, const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double rel_err, const double abs_err, const int limit_size)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;

  return GSL_integrate_qaws(Func, a, b, alpha, beta, mu, nu, abs_err, rel_err, limit_size);
}


// ============================================================================


double cbl::gsl::GSL_integrate_cquad (FunctionDoubleDoublePtrVectorRef func, const std::shared_ptr<void> pp, const std::vector<double> par, const double a, const double b, const double rel_err, const double abs_err, const int nevals)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_cquad(func_bind, a, b, abs_err, rel_err, nevals);
}


// ============================================================================


double cbl::gsl::GSL_integrate_qag (FunctionDoubleDoublePtrVectorRef func, const std::shared_ptr<void> pp, const std::vector<double> par, const double a, const double b, const double rel_err, const double abs_err, const int limit_size, const int rule)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_qag(func_bind, a, b, abs_err, rel_err, limit_size, rule);
}


// ============================================================================


double cbl::gsl::GSL_integrate_qagiu (FunctionDoubleDoublePtrVectorRef func, const std::shared_ptr<void> pp, const std::vector<double> par, const double a, const double rel_err, const double abs_err, const int limit_size)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_qagiu(func_bind, a, rel_err, abs_err, limit_size);
}


// ============================================================================


double cbl::gsl::GSL_integrate_qaws (FunctionDoubleDoublePtrVectorRef func, const std::shared_ptr<void> pp, const std::vector<double> par, const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double rel_err, const double abs_err, const int limit_size)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_qaws(func_bind, a, b, alpha, beta, mu, nu, rel_err, abs_err, limit_size);
}


// ============================================================================


double cbl::gsl::GSL_root_brent (gsl_function Func, const double low_guess, const double up_guess, const double rel_err, const double abs_err)
{
  gsl_set_error_handler_off();

  int status = 0, iter = 0, max_iter = 10000;
  const gsl_root_fsolver_type *T;
  double r;

  gsl_root_fsolver *s;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  
  double x_lo = low_guess;
  double x_hi = up_guess;

  gsl_root_fsolver_set(s, &Func, x_lo, x_hi); 

  do
  {
    iter ++;
    status = gsl_root_fsolver_iterate(s);

    if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE))
      check_GSL_fail(status, true, "GSL_root_brent", "gsl_root_fsolver_iterate");
    
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    
    status = gsl_root_test_interval(x_lo, x_hi, abs_err, rel_err);

    if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE))
      check_GSL_fail(status, true, "GSL_root_brent", "gsl_root_test_interval");
  }
  
  while (status==GSL_CONTINUE && iter<max_iter);

  gsl_root_fsolver_free(s);

  check_GSL_fail(status, true, "GSL_minimize_nD", par::defaultString);

  return r;
}


// ============================================================================


double cbl::gsl::GSL_root_brent (FunctionDoubleDouble func, const double xx0, const double low_guess, const double up_guess, const double rel_err, const double abs_err)
{
  gsl_set_error_handler_off();

  STR_generic_func_GSL params;
  params.f = func;
  params.xx0 = xx0;

  gsl_function Func;
  Func.function = &generic_roots;
  Func.params = &params;

  return GSL_root_brent(Func, low_guess, up_guess, rel_err, abs_err);
}


// ============================================================================


vector<double> cbl::gsl::GSL_minimize_nD (FunctionDoubleVector func, const std::vector<double> start, const std::vector<std::vector<double>> ranges, const unsigned int max_iter, const double tol, const double epsilon)
{
  if (ranges.size() != start.size() && ranges.size() != 0)
    ErrorCBL ("Error in GSL_minimize_nD of GSLwrapper.cpp, vector of ranges must have the same size of start vector.");
  gsl_set_error_handler_off();

  size_t npar = start.size();

  STR_generic_func_GSL params;
  params.fmin_return = (ranges.size() == start.size()) ? [&] (vector<double> par) {return (inRange(par, ranges)) ? func(par) : par::defaultDouble; } : func;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status = 0;
  double size;

  // Starting point
  
  x = gsl_vector_alloc(npar);
  ss = gsl_vector_alloc(npar);

  for (size_t i=0; i<npar; i++) {
    gsl_vector_set(x, i, start[i]);
    gsl_vector_set(ss, i, (epsilon>0 && ranges.size()>0) ? (ranges[i][1]-ranges[i][0])*epsilon : 1);
  }

  // Initialize the method and iterate 

  minex_func.n = npar;
  minex_func.f = &generic_minimizer_return;
  minex_func.params = &params;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter ++;
      status = gsl_multimin_fminimizer_iterate(s);

      if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE))
	check_GSL_fail(status, true, "GSL_minimize_nD", "gsl_multimin_fminimizer_iterate");

      size = gsl_multimin_fminimizer_size(s);
    
      status = gsl_multimin_test_size(size, tol);

      if ((status!=GSL_SUCCESS) && (status!=GSL_CONTINUE))
	check_GSL_fail(status, true, "GSL_minimize_nD", "gsl_multimin_fminimizer_iterate");
    }
  
  while (status == GSL_CONTINUE && iter < max_iter);

  check_GSL_fail(status, true, "GSL_minimize_nD", par::defaultString);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return params.parameters_return;
}


// ============================================================================


vector<double> cbl::gsl::GSL_minimize_nD (FunctionDoubleVectorRef func, const std::vector<double> start, const std::vector<std::vector<double>> ranges, const unsigned int max_iter, const double tol, const double epsilon)
{
  if (ranges.size() != start.size() && ranges.size() != 0)
    ErrorCBL ("Error in GSL_minimize_nD of GSLwrapper.cpp, vector of ranges must have the same size of start vector.");
  gsl_set_error_handler_off();

  size_t npar = start.size();

  STR_generic_func_GSL params;
  params.fmin_return = (ranges.size() == start.size()) ? [&] (vector<double> &par) { return (inRange(par, ranges)) ? func(par) : -par::defaultDouble; } : func;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status = 0;
  double size;

  // Starting point
  
  x = gsl_vector_alloc(npar);
  ss = gsl_vector_alloc(npar);

  for (size_t i=0; i<npar; i++) {
    gsl_vector_set(x, i, start[i]);
    double val = (epsilon>0 && ranges.size()>0) ? (ranges[i][1]-ranges[i][0])*epsilon : 1;
    gsl_vector_set(ss, i, val);
  }

  // Initialize the method and iterate 

  minex_func.n = npar;
  minex_func.f = &generic_minimizer_return;
  minex_func.params = &params;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter ++;

      if ((status != GSL_SUCCESS) && (status!=GSL_CONTINUE))
	check_GSL_fail(status, true, "GSL_minimize_nD", "gsl_multimin_test_size");

      status = gsl_multimin_fminimizer_iterate(s);

      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, tol);

      if ((status != GSL_SUCCESS) && (status!=GSL_CONTINUE))
	check_GSL_fail(status, true, "GSL_minimize_nD", "gsl_multimin_test_size");
    }
  
  while (status == GSL_CONTINUE && iter < max_iter);

  check_GSL_fail(status, true, "GSL_minimize_nD", par::defaultString);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return params.parameters_return;
}


// ============================================================================


double cbl::gsl::GSL_minimize_1D (FunctionDoubleDouble func, const double start, double min, double max, const int max_iter, const bool verbose)
{
  gsl_set_error_handler_off();

  const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);

  int status = 0;
  int iter = 0;
  double m = start;

  STR_generic_func_GSL params;
  params.f = func;

  gsl_function F;
  F.function = generic_function;
  F.params = &func;

  gsl_min_fminimizer_set(s, &F, m, min, max);
  
  do
    {
      iter ++;
      status = gsl_min_fminimizer_iterate(s);

      if ((status != GSL_SUCCESS) && (status != GSL_CONTINUE))
	check_GSL_fail(status, true, "GSL_minimize_1D", "gsl_min_fminimizer_iterate");

      m = gsl_min_fminimizer_x_minimum(s);
      min = gsl_min_fminimizer_x_lower(s);
      max = gsl_min_fminimizer_x_upper(s);

      status = gsl_min_test_interval(min, max, 0.001, 0.0);

      if ((status != GSL_SUCCESS) && (status != GSL_CONTINUE))
	check_GSL_fail(status, true, "GSL_minimize_1D", "gsl_min_test_interval");
    }
  
  while (status==GSL_CONTINUE && iter<max_iter);

  check_GSL_fail(status, true, "GSL_minimize_1D", par::defaultString);

  if ((status==GSL_SUCCESS) && verbose)
    coutCBL << "Converged to a minimum." << endl;

  gsl_min_fminimizer_free(s);
  
  return m;
}


// ============================================================================


double cbl::gsl::GSL_polynomial_eval (const double x, const std::shared_ptr<void> fixed_parameters, const std::vector<double> coeff)
{
  (void)fixed_parameters;
  return gsl_poly_eval(coeff.data(), coeff.size(), x);
}


// ============================================================================


void cbl::gsl::GSL_polynomial_root (const std::vector<double> coeff, std::vector<std::vector<double>> &root)
{
  gsl_set_error_handler_off();

  size_t size = coeff.size();
  size_t root_size = 2*(size-1);
  double *z = new double[root_size];
 
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(size);
  
  int status = gsl_poly_complex_solve(coeff.data(), size, w, z);

  check_GSL_fail(status, true, "GSL_polynomial", "gsl_poly_complex_solve");

  gsl_poly_complex_workspace_free(w);
  root.resize(size-1,vector<double>(2,0));

  for (size_t i=0; i<size-1; i++) {
    root[i][0] = z[2*i]; 
    root[i][1] = z[2*i+1];
  }

  delete[] z;
}
