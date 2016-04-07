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

/** @file Func/GSLfunction.cpp
 *
 *  @brief Implementation of members of the class GSLfunction 
 *
 *  This file contains the implementation the class GSLfunction, used to
 *  wrap gsl minimization function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "GSLfunction.h"
using namespace cosmobl;


// ============================================================================================

unique_ptr<GSLfunction> cosmobl::GSLfunction::make_GSLfunction(func_1par_1 function, shared_ptr<void> function_parameters){
    return unique_ptr<GSLfunction_1D_1>(new GSLfunction_1D_1(function,function_parameters));
}


// ============================================================================================


unique_ptr<cosmobl::GSLfunction> cosmobl::GSLfunction::make_GSLfunction(int npar, func_npar_1 function, shared_ptr<void> function_parameters){
    return unique_ptr<GSLfunction_nD_1>(new GSLfunction_nD_1(npar,function,function_parameters));
}


// ============================================================================================


unique_ptr<GSLfunction> cosmobl::GSLfunction::make_GSLfunction(func_1par_2 function, vector<double> params, shared_ptr<void> function_parameters){
    return unique_ptr<GSLfunction_1D_2>(new GSLfunction_1D_2(function,params,function_parameters));
}


// ============================================================================================


unique_ptr<cosmobl::GSLfunction> cosmobl::GSLfunction::make_GSLfunction(int npar, func_npar_2 function, vector<double> params, shared_ptr<void> function_parameters){
    return unique_ptr<GSLfunction_nD_2>(new GSLfunction_nD_2(npar,function,params,function_parameters));
}


// ============================================================================================


void cosmobl::GSLfunction_1D_1::minimize(double &result, int max_iter, double min, double max)
{

  int status;
  int iter = 0;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = result;

  typedef function<double(double)> fun_type;
  fun_type func = [&](double x){return m_function(x,m_function_parameters);};

  gsl_function F;
  F.function = [](double x, void* p)
  {
    fun_type *f = (fun_type *)p; 
    return f->operator()(x);
  };

  F.params = &func;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, min, max);
  
  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate (s);

    m = gsl_min_fminimizer_x_minimum (s);
    min = gsl_min_fminimizer_x_lower (s);
    max = gsl_min_fminimizer_x_upper (s);

    status 
      = gsl_min_test_interval (min, max, 0.001, 0.0);

    if (status == GSL_SUCCESS)
    {
      printf ("-----> Converged to minimum \n");
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);
  
  result = m;
  
}


// ============================================================================================


void cosmobl::GSLfunction_1D_2::minimize(double &result, int max_iter, double min, double max)
{

  int status;
  int iter = 0;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = result;

  typedef function<double(double)> fun_type;
  fun_type func = [&](double x){return m_function(x,m_function_parameters,m_parameters);};

  gsl_function F;
  F.function = [](double x, void* p)
  {
    fun_type *f = (fun_type *)p; 
    return f->operator()(x);
  };

  F.params = &func;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, min, max);
  
  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate (s);

    m = gsl_min_fminimizer_x_minimum (s);
    min = gsl_min_fminimizer_x_lower (s);
    max = gsl_min_fminimizer_x_upper (s);

    status 
      = gsl_min_test_interval (min, max, 0.001, 0.0);

    if (status == GSL_SUCCESS)
    {
      printf ("-----> Converged to minimum \n");
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free (s);
  
  result = m;
  
}


// ============================================================================================


void cosmobl::GSLfunction_nD_1::minimize(vector<double> &result, unsigned int max_iter, double tol)
{

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  x = gsl_vector_alloc (m_npar);
  for(int i=0;i<m_npar;i++)
    gsl_vector_set(x,i,result[i]);
  
    

  // Set initial step sizes to 1 
  ss = gsl_vector_alloc (m_npar);
  gsl_vector_set_all (ss, 1);
  
  // Create the function to minimize 
  typedef function<double(vector<double> )> fun_type;
  fun_type func = [&](vector<double> x){return m_function(x,m_function_parameters);};

  // Initialize method and iterate 

  minex_func.n = m_npar;
  minex_func.f = [](const gsl_vector *gsl_x, void *p)
  {
    fun_type *f = (fun_type *)p;
    vector<double> xx;
    for(unsigned int i=0;i<gsl_x->size;i++)
      xx.push_back(gsl_vector_get(gsl_x,i));
    return f->operator()(xx);
  };
  minex_func.params = &func;

  s = gsl_multimin_fminimizer_alloc (T, m_npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (s);
    
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("-----> Converged to minimum \n");
    }

  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

}

void cosmobl::GSLfunction_nD_1::minimize(vector<double> &result, vector<double> step_size, unsigned int max_iter, double tol)
{

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  x = gsl_vector_alloc (m_npar);
  ss = gsl_vector_alloc (m_npar);
  for(int i=0;i<m_npar;i++){
    gsl_vector_set(x,i,result[i]);
    gsl_vector_set(ss,i,step_size[i]);
  }
    

  // Create the function to minimize 
  typedef function<double(vector<double> )> fun_type;
  fun_type func = [&](vector<double> x){return m_function(x,m_function_parameters);};

  // Initialize method and iterate 

  minex_func.n = m_npar;
  minex_func.f = [](const gsl_vector *gsl_x, void *p)
  {
    fun_type *f = (fun_type *)p;
    vector<double> xx;
    for(unsigned int i=0;i<gsl_x->size;i++)
      xx.push_back(gsl_vector_get(gsl_x,i));
    return f->operator()(xx);
  };
  minex_func.params = &func;

  s = gsl_multimin_fminimizer_alloc (T, m_npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (s);
    
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("-----> Converged to minimum \n");
    }

  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

}


// ============================================================================================


void cosmobl::GSLfunction_nD_2::minimize(vector<double> &result, unsigned int max_iter, double tol)
{

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  x = gsl_vector_alloc (m_npar);
  for(int i=0;i<m_npar;i++)
    gsl_vector_set(x,i,result[i]);
  
    

  // Set initial step sizes to 1 
  ss = gsl_vector_alloc (m_npar);
  gsl_vector_set_all (ss, 1);
  
  // Create the function to minimize 
  typedef function<double(vector<double> )> fun_type;
  fun_type func = [&](vector<double> x){return m_function(x,m_function_parameters,m_parameters);};

  // Initialize method and iterate 

  minex_func.n = m_npar;
  minex_func.f = [](const gsl_vector *gsl_x, void *p)
  {
    fun_type *f = (fun_type *)p;
    vector<double> xx;
    for(unsigned int i=0;i<gsl_x->size;i++)
      xx.push_back(gsl_vector_get(gsl_x,i));
    return f->operator()(xx);
  };
  minex_func.params = &func;

  s = gsl_multimin_fminimizer_alloc (T, m_npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (s);
    
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("-----> Converged to minimum \n");
    }

  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

}

void cosmobl::GSLfunction_nD_2::minimize(vector<double> &result, vector<double> step_size, unsigned int max_iter, double tol)
{

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  // Starting point 
  x = gsl_vector_alloc (m_npar);
  ss = gsl_vector_alloc (m_npar);
  for(int i=0;i<m_npar;i++){
    gsl_vector_set(x,i,result[i]);
    gsl_vector_set(ss,i,step_size[i]);
  }
    

  // Create the function to minimize 
  typedef function<double(vector<double> )> fun_type;
  fun_type func = [&](vector<double> x){return m_function(x,m_function_parameters,m_parameters);};

  // Initialize method and iterate 

  minex_func.n = m_npar;
  minex_func.f = [](const gsl_vector *gsl_x, void *p)
  {
    fun_type *f = (fun_type *)p;
    vector<double> xx;
    for(unsigned int i=0;i<gsl_x->size;i++)
      xx.push_back(gsl_vector_get(gsl_x,i));
    return f->operator()(xx);
  };
  minex_func.params = &func;

  s = gsl_multimin_fminimizer_alloc (T, m_npar);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (s);
    
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("-----> Converged to minimum \n");
    }

  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

}


