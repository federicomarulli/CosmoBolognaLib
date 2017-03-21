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
 *  @file Func/FuncGrid.cpp
 *
 *  @brief Methods of the class FuncGrid
 *
 *  This file contains the implementation of the methods of the class
 *  FuncGrid
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Func.h"

using namespace cosmobl;
using namespace glob;


// ============================================================================


cosmobl::glob::FuncGrid::FuncGrid (const vector<double> x, const vector<double> y, const string interpType)
  : m_x(x), m_y(y), m_size(x.size()), m_interpType(interpType)
{
  m_acc = gsl_interp_accel_alloc();
	  
  if (m_size<5 && interpType!="Linear") {
    WarningMsg("Warning in the constructor of FuncGrid in FuncGrid.h: the array size is less than 5 -> setting interpolation method to Linear");
    m_interpType = "Linear";
  }

  if (m_interpType=="Linear") 
    m_type = gsl_interp_linear;

  else if (m_interpType=="Poly") 
    m_type = gsl_interp_polynomial;

  else if (m_interpType=="Spline") 
    m_type = gsl_interp_cspline;

  else if (m_interpType=="Spline_periodic") 
    m_type = gsl_interp_cspline_periodic;

  else if (m_interpType=="Akima") 
    m_type = gsl_interp_akima;

  else if (m_interpType=="Akima_periodic") 
    m_type = gsl_interp_akima_periodic;

  else if (m_interpType=="Steffen") 
    m_type = gsl_interp_steffen;

  else 
    ErrorCBL("Error in the constructor of FuncGrid in FuncGrid.h: the value of m_interpType is not permitted!");

  m_xmin = Min(m_x);
  m_xmax = Max(m_x);
  m_spline = gsl_spline_alloc(m_type, m_size); 
  gsl_spline_init(m_spline, m_x.data(), m_y.data(), m_size);
}


// =====================================================================================

				 
void cosmobl::glob::FuncGrid::free ()
{
  gsl_spline_free(m_spline);
  gsl_interp_accel_free(m_acc);
}


// =====================================================================================
     

double cosmobl::glob::FuncGrid::operator () (const double xx) const
{	    
  if (xx<m_xmin) // perform a linear extrapolation
    return m_spline->y[0]+(xx-m_xmin)/(m_spline->x[1]-m_xmin)*(m_spline->y[1]-m_spline->y[0]);
  
  else if (xx>m_xmax) // perform a linear extrapolation
    return m_spline->y[m_size-2]+(xx-m_spline->x[m_size-2])/(m_xmax-m_spline->x[m_size-2])*(m_spline->y[m_size-1]-m_spline->y[m_size-2]);

  // performe an interpolation
  return gsl_spline_eval(m_spline, xx, m_acc);
}


// =====================================================================================


vector<double> cosmobl::glob::FuncGrid::eval_func (const vector<double> xx) const
{
  vector<double> yy;
  
  for (size_t i=0; i<xx.size(); i++)
    yy.push_back(gsl_spline_eval(m_spline, xx[i], m_acc));

  return yy;
}


// =====================================================================================


double cosmobl::glob::FuncGrid::D1v (const double xx) const
{
  return gsl_spline_eval_deriv(m_spline, xx, m_acc);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::D2v (const double xx) const
{
  return gsl_spline_eval_deriv2(m_spline, xx, m_acc);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::integrate_qag (const double a, const double b, const double prec, const int limit_size, const int rule)
{
  function<double(double)> f = bind(&FuncGrid::operator(), this, std::placeholders::_1);

  return gsl::GSL_integrate_qag(f, a, b, prec, limit_size, rule);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::integrate_qaws (const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double prec, const int limit_size)
{
  function<double(double)> f = bind(&FuncGrid::operator(), this, std::placeholders::_1);

  return gsl::GSL_integrate_qaws(f, a, b, alpha, beta, mu, nu, prec, limit_size);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::root (const double x_low, const double x_up, const double fx0, const double prec)
{
  function<double(double)> f = bind(&FuncGrid::operator(), this, std::placeholders::_1);

  return gsl::GSL_root_brent(f, fx0,  x_low, x_up, prec);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::root_D1v (const double x_low, const double x_up, const double fx0, const double prec)
{
  function<double(double)> f = bind(&FuncGrid::D1v, this, std::placeholders::_1);

  return gsl::GSL_root_brent(f, fx0,  x_low, x_up, prec);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::root_D2v (const double x_low, const double x_up, const double fx0, const double prec)
{
  function<double(double)> f = bind(&FuncGrid::D2v, this, std::placeholders::_1);

  return gsl::GSL_root_brent(f, fx0,  x_low, x_up, prec);
}
