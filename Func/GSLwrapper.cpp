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
 *  @file Func/GSLwrapper.cpp
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

#include "Func.h"


// ============================================================================


double cosmobl::gsl::generic_function (const double xx, void *params)
{
  cosmobl::gsl::STR_generic_func_GSL *pp = (cosmobl::gsl::STR_generic_func_GSL *) params;
  return pp->f(xx);
}


// ============================================================================


double cosmobl::gsl::generic_roots (double xx, void *params)
{
  cosmobl::gsl::STR_generic_func_GSL *pp = (cosmobl::gsl::STR_generic_func_GSL *) params;
  return pp->f(xx)-pp->xx0;
}


// ============================================================================


double cosmobl::gsl::generic_minimizer (const gsl_vector * xx, void * params)
{
  vector<double> _xx;
  for (size_t i=0; i<xx->size; i++)
    _xx.push_back(gsl_vector_get(xx, i));

  cosmobl::gsl::STR_generic_func_GSL *pp = (cosmobl::gsl::STR_generic_func_GSL *) params;
  pp->parameters_return = _xx;

  return pp->fmin(_xx);
}


// ============================================================================


double cosmobl::gsl::generic_minimizer_return (const gsl_vector * xx, void * params)
{
  vector<double> _xx;
  for (size_t i=0; i<xx->size; i++)
    _xx.push_back(gsl_vector_get(xx, i));

  cosmobl::gsl::STR_generic_func_GSL *pp = (cosmobl::gsl::STR_generic_func_GSL *) params;

  double val = pp->fmin_return(_xx);
  pp->parameters_return = _xx;
  
  return val;
}


// ============================================================================


double cosmobl::gsl::GSL_derivative (gsl_function Func, const double xx, const double hh, const double prec)
{
  double Deriv, error;

  gsl_deriv_central(&Func, xx, hh, &Deriv, &error);

  return (Deriv!=0 && error/Deriv<prec) ? Deriv : ErrorCBL("Error in cosmobl::gsl::GSL_derivative of CSLwrapper! error/Deriv = "+conv(error/Deriv, par::fDP6)+" > prec = "+conv(prec, par::fDP3));
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_cquad (gsl_function Func, const double a, const double b, const double prec, const int nevals)
{
  gsl_set_error_handler_off();

  double Int, error;
  size_t nn = nevals;

  gsl_integration_cquad_workspace *ww = gsl_integration_cquad_workspace_alloc(nevals);
  gsl_integration_cquad(&Func, a, b, 0., prec, ww, &Int, &error, &nn); 
  gsl_integration_cquad_workspace_free(ww);

  if (error/fabs(Int)>prec) {
    string Msg = "Warning in GSL_integrate_qag: unable to reach the requested precision: "+conv(error/fabs(Int), par::fDP6)+" "+conv(prec, par::fDP6);
    WarningMsg(Msg);
  }
  
  return Int;
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qag (gsl_function Func, const double a, const double b, const double prec, const int limit_size, const int rule)
{
  gsl_set_error_handler_off();

  double Int, error;
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
  gsl_integration_qag(&Func, a, b, 0., prec, limit_size, rule, ww, &Int, &error); 
  gsl_integration_workspace_free(ww);

  if (error/fabs(Int)>prec) {
    string Msg = "Warning in GSL_integrate_qag: unable to reach the requested precision: "+conv(error/fabs(Int), par::fDP6)+" "+conv(prec, par::fDP6);
    WarningMsg(Msg);
  }
  
  return Int;
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qagiu (gsl_function Func, const double a, const double prec, const int limit_size)
{
  gsl_set_error_handler_off();

  double Int, error;
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
  gsl_integration_qagiu(&Func, a, 0., prec, limit_size, ww, &Int, &error); 
  gsl_integration_workspace_free(ww);

  if (error/fabs(Int)>prec) {
     string Msg = "Warning in GSL_integrate_qagiu: unable to reach the requested precision: "+conv(error/fabs(Int), par::fDP6)+" "+conv(prec, par::fDP6);
     WarningMsg(Msg);
  }
  
  return Int;
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qaws (gsl_function Func, const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double prec, const int limit_size)
{
  gsl_set_error_handler_off();

  gsl_integration_qaws_table *T = gsl_integration_qaws_table_alloc(alpha, beta, mu, nu);

  double Int, error;
  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(limit_size);
  gsl_integration_qaws (&Func, a, b, T, 0., prec, limit_size, ww, &Int, &error); 

  gsl_integration_workspace_free(ww);
  gsl_integration_qaws_table_free(T);

  if (error/fabs(Int)>prec) {
     string Msg = "Warning in GSL_integrate_qagiu: unable to reach the requested precision: "+conv(error/fabs(Int), par::fDP6)+" "+conv(prec, par::fDP6);
     WarningMsg(Msg);
  }
  
  return Int;
}


// ============================================================================


double cosmobl::gsl::GSL_derivative (function<double(double)> func, const double xx, const double hh, const double prec)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;
  Func.function = generic_function;
  Func.params = &params;

  return GSL_derivative(Func, xx, hh, prec);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_cquad (function<double(double)> func, const double a, const double b, const double prec, const int nevals)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;
  
  return GSL_integrate_cquad(Func, a, b, prec, nevals);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qag (function<double(double)> func, const double a, const double b, const double prec, const int limit_size, const int rule)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;
  
  return GSL_integrate_qag(Func, a, b, prec, limit_size, rule);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qagiu (function<double(double)> func, const double a, const double prec, const int limit_size)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;

  return GSL_integrate_qagiu(Func, a, prec, limit_size);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qaws (function<double(double)> func, const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double prec, const int limit_size)
{
  STR_generic_func_GSL params;
  params.f = func;

  gsl_function Func;  
  Func.function = generic_function;
  Func.params = &params;

  return GSL_integrate_qaws(Func, a, b, alpha, beta, mu, nu, prec, limit_size);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_cquad (function<double(double, shared_ptr<void>, vector<double>)> func, const shared_ptr<void> pp, const vector<double> par, const double a, const double b, const double prec, const int nevals)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_cquad(func_bind, a, b, prec, nevals);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qag (function<double(double, shared_ptr<void>, vector<double>)> func, const shared_ptr<void> pp, const vector<double> par, const double a, const double b, const double prec, const int limit_size, const int rule)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_qag(func_bind, a, b, prec, limit_size, rule);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qagiu (function<double(double, shared_ptr<void>, vector<double>)> func, const shared_ptr<void> pp, const vector<double> par, const double a, const double prec, const int limit_size)
{

  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_qagiu(func_bind, a, prec, limit_size);
}


// ============================================================================


double cosmobl::gsl::GSL_integrate_qaws (function<double(double, shared_ptr<void>, vector<double>)> func, const shared_ptr<void> pp, const vector<double> par, const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double prec, const int limit_size)
{
  function<double(double)> func_bind = bind(func, placeholders::_1, pp, par);

  return GSL_integrate_qaws(func_bind, a, b, alpha, beta, mu, nu, prec, limit_size);
}


// ============================================================================


double cosmobl::gsl::GSL_root_brent (gsl_function Func, const double low_guess, const double up_guess, const double prec)
{
  int status;
  int iter = 0, max_iter = 10000;
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
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, prec);
  }
  
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  return r;
}


// ============================================================================


double cosmobl::gsl::GSL_root_brent (function<double(double)> func, double xx0, const double low_guess, const double up_guess, const double prec)
{
  STR_generic_func_GSL params;
  params.f = func;
  params.xx0 = xx0;

  gsl_function Func;
  Func.function = &generic_roots;
  Func.params = &params;

  return GSL_root_brent(Func, low_guess, up_guess, prec);
}


// ============================================================================


vector<double> cosmobl::gsl::GSL_minimize_nD (function<double(vector<double > &)> func, const vector<double> start, const vector<double> step_size, const unsigned int max_iter, const double tol)
{
  size_t npar = start.size();

  STR_generic_func_GSL params;
  params.fmin_return = func;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point
  
  x = gsl_vector_alloc(npar);
  ss = gsl_vector_alloc(npar);

  for (size_t i=0; i<npar; i++) {
    gsl_vector_set(x, i, start[i]);
    gsl_vector_set(ss, i, ( (step_size.size()==npar) ? step_size[i] : 1) );
  }

  
  // Initialize the method and iterate 

  minex_func.n = npar;
  minex_func.f = &generic_minimizer_return;
  minex_func.params = &params;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status) 
	break;

      size = gsl_multimin_fminimizer_size(s);
    
      status = gsl_multimin_test_size(size, tol);

      if (status == GSL_SUCCESS)
	coutCBL << "-----> converged to a minimum" << endl;
    }
  
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return params.parameters_return;
}


// ============================================================================


vector<double> cosmobl::gsl::GSL_minimize_nD (function<double(vector<double>)> func, const vector<double> start, const vector<double> step_size, const unsigned int max_iter, const double tol)
{
  size_t npar = start.size();

  STR_generic_func_GSL params;
  params.fmin = func;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  
  // Starting point
  
  x = gsl_vector_alloc(npar);
  ss = gsl_vector_alloc(npar);

  for (size_t i=0; i<npar; i++) {
    gsl_vector_set(x, i, start[i]);
    gsl_vector_set(ss, i, ( (step_size.size()==npar) ? step_size[i] : 1) );
  }

  
  // Initialize the method and iterate 

  minex_func.n = npar;
  minex_func.f = &generic_minimizer;
  minex_func.params = &params;

  s = gsl_multimin_fminimizer_alloc (T, npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status) 
	break;

      size = gsl_multimin_fminimizer_size(s);
    
      status = gsl_multimin_test_size(size, tol);

      if (status == GSL_SUCCESS)
	coutCBL << "-----> converged to a minimum" << endl;
    }
  
  while (status == GSL_CONTINUE && iter < max_iter);


  vector<double> result(npar);

  for (size_t i=0; i<npar; i++)
    result[i] = gsl_vector_get(s->x, i);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return result;
}


// ============================================================================


double cosmobl::gsl::GSL_minimize_1D (function<double(double)> func, const double start, double min, double max, const int max_iter, const bool verbose)
{
  int status;
  int iter = 0;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = start;

  STR_generic_func_GSL params;
  params.f = func;

  gsl_function F;
  F.function = generic_function;
  F.params = &func;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);
  gsl_min_fminimizer_set(s, &F, m, min, max);
  
  do
    {
      iter ++;
      status = gsl_min_fminimizer_iterate(s);

      m = gsl_min_fminimizer_x_minimum(s);
      min = gsl_min_fminimizer_x_lower(s);
      max = gsl_min_fminimizer_x_upper(s);

      status = gsl_min_test_interval(min, max, 0.001, 0.0);

      if (status==GSL_SUCCESS && verbose)
	coutCBL << "-----> Converged to minimum" << endl;
    }
  
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free(s);
  
  return m;
}


// ============================================================================


double cosmobl::gsl::GSL_polynomial_eval (const double x, shared_ptr<void> fixed_parameters, const vector<double> coeff)
{
  (void)fixed_parameters;
  return gsl_poly_eval(coeff.data(), coeff.size(), x);
}


// ============================================================================


void cosmobl::gsl::GSL_polynomial_root (const vector<double> coeff, vector<vector<double>> &root)
{
  size_t size = coeff.size();
  size_t root_size = 2*(size-1);
  double *z = new double[root_size];
 
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(size);
  
  gsl_poly_complex_solve(coeff.data(), size, w, z);

  gsl_poly_complex_workspace_free(w);
  root.resize(size-1,vector<double>(2,0));

  for (size_t i=0; i<size-1; i++) {
    root[i][0] = z[2*i]; 
    root[i][1] = z[2*i+1];
  }

  delete[] z;
}
