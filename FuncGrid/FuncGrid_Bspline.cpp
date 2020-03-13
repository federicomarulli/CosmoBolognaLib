/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file FuncGrid/FuncGrid_Bspline.cpp
 *
 *  @brief Methods of the class FuncGrid_Bspline
 *
 *  This file contains the implementation of the methods of the class
 *  FuncGrid_Bspline
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "FuncGrid_Bspline.h"

using namespace std;

using namespace cbl;
using namespace glob;


// ============================================================================


cbl::glob::FuncGrid_Bspline::FuncGrid_Bspline (const std::vector<double> x, const std::vector<double> fx, const int nbreakpoints, const int order, const double frac, const double xmin, const double xmax)
{
  this->set(x, fx, nbreakpoints, order, frac, xmin, xmax);
}


// ============================================================================


cbl::glob::FuncGrid_Bspline::FuncGrid_Bspline (const std::vector<double> x, const std::vector<double> fx, const std::vector<double> breakpoints, const int order, const double frac)
{
  this->set(x, fx, breakpoints, order, frac);
}


// ============================================================================


void cbl::glob::FuncGrid_Bspline::m_set_bspline(const vector<double> x, const vector<double> fx, const int nbreakpoints, const int order)
{
  m_x = x;
  m_fx = fx;
  m_order = order;
  m_nbreakpoints = nbreakpoints;
  m_ncoefficients = m_nbreakpoints-2+m_order;
 
  // Define gsl bspline workspace and vectors
  shared_ptr<gsl_bspline_workspace> bspline(gsl_bspline_alloc(m_order, m_nbreakpoints), gsl_bspline_free);
  shared_ptr<gsl_vector> Bcoeff(gsl_vector_alloc(m_ncoefficients), gsl_vector_free); 
  shared_ptr<gsl_vector> weights(gsl_vector_alloc(m_ncoefficients), gsl_vector_free); 
  shared_ptr<gsl_matrix> covariance(gsl_matrix_alloc(m_ncoefficients, m_ncoefficients), gsl_matrix_free);

  m_bspline = bspline;
  m_Bcoeff = Bcoeff;
  m_weights = weights;
  m_covariance = covariance;
}

// ============================================================================


void cbl::glob::FuncGrid_Bspline::m_set_knots (const double xmin, const double xmax)
{
  gsl_bspline_knots_uniform((xmin!=par::defaultDouble) ? xmin : Min(m_x),
      			    (xmax!=par::defaultDouble) ? xmax : Max(m_x),
      			    m_bspline.get());
}


// ============================================================================


void cbl::glob::FuncGrid_Bspline::m_set_knots (const std::vector<double> breakpoints)
{
  shared_ptr<gsl_vector> bp(gsl_vector_alloc(m_nbreakpoints), gsl_vector_free); 

  for (int i=0; i<m_nbreakpoints; i++)
    gsl_vector_set(bp.get(), i, breakpoints[i]);

  gsl_bspline_knots(bp.get(), m_bspline.get());
}


// ============================================================================


void cbl::glob::FuncGrid_Bspline::m_compute_func_integral ()
{
  auto func = [&] (const double &xx) {return this->operator()(xx, par::defaultDouble);};

  m_integral = wrapper::gsl::GSL_integrate_qag(func, Min(m_x), Max(m_x));
}

// ============================================================================


void cbl::glob::FuncGrid_Bspline::m_linear_fit(const double frac)
{
  cout << "Computing basis spline coefficients" << endl;

  size_t n = m_x.size();

  // Put input points in gsl vectors
  shared_ptr<gsl_vector> y(gsl_vector_alloc(n), gsl_vector_free); 
  shared_ptr<gsl_vector> w(gsl_vector_alloc(n), gsl_vector_free); 

  for (size_t i=0; i<m_x.size(); i++) {
    gsl_vector_set(y.get(), i, m_fx[i]);
    gsl_vector_set(w.get(), i, (frac>0) ? pow(frac*m_fx[i], -2) : 1.); //check!!!
  }

  // Construct the fit matrix X
  shared_ptr<gsl_matrix> XX(gsl_matrix_alloc(n, m_ncoefficients), gsl_matrix_free);
  for (size_t i = 0; i < n; ++i)
    {
      /* compute B_j(xi) for all j */
      gsl_bspline_eval(m_x[i], m_Bcoeff.get(), m_bspline.get());
      /* fill in row i of X */
      for (int j = 0; j < m_ncoefficients; ++j)
	gsl_matrix_set(XX.get(), i, j, gsl_vector_get(m_Bcoeff.get(), j));
    }

  // Compute coefficients
  double chisq;

  shared_ptr<gsl_multifit_linear_workspace> mw(gsl_multifit_linear_alloc(n, m_ncoefficients), gsl_multifit_linear_free);
  
  gsl_multifit_wlinear(XX.get(), w.get(), y.get(), m_weights.get(), m_covariance.get(), &chisq, mw.get());

  cout << "Done! Chi2/d.o.f. = " << chisq/(n-m_ncoefficients) << endl;
}


// ============================================================================


void cbl::glob::FuncGrid_Bspline::set (const std::vector<double> x, const std::vector<double> fx, const int nbreakpoints, const int order, const double frac, const double xmin, const double xmax)
{
  this->m_set_bspline(x, fx, nbreakpoints, order);
  this->m_set_knots(xmin, xmax);
  this->m_linear_fit(frac);
  this->m_compute_func_integral();
}


// ============================================================================


void cbl::glob::FuncGrid_Bspline::set (const std::vector<double> x, const std::vector<double> fx, const std::vector<double> breakpoints, const int order, const double frac)
{
  this->m_set_bspline(x, fx, static_cast<int>(breakpoints.size()), order);
  this->m_set_knots(breakpoints);
  this->m_linear_fit(frac);
  this->m_compute_func_integral();
}


// ============================================================================


double cbl::glob::FuncGrid_Bspline::operator () (const double xx, const double integral) const
{
  double fx, fx_err;
  
  gsl_bspline_eval(xx, m_Bcoeff.get(), m_bspline.get());
  gsl_multifit_linear_est(m_Bcoeff.get(), m_weights.get(), m_covariance.get(), &fx, &fx_err);

  return (integral!=par::defaultDouble) ? integral*fx/m_integral : fx;
}


// ============================================================================


vector<double> cbl::glob::FuncGrid_Bspline::eval_func (const vector<double> xx, const double integral) const
{
  vector<double> vv(xx.size(), 0), yerr;

  for (size_t i=0; i<xx.size(); i++)
    vv[i] = this->operator()(xx[i], integral);

  return vv;
}
