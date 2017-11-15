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
 *  FuncGrid and FuncGrid_2D
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Func.h"

using namespace cosmobl;
using namespace glob;


// ============================================================================


cosmobl::glob::FuncGrid::FuncGrid (const vector<double> x, const vector<double> y, const string interpType, const cosmobl::binType bin_type)
  : m_x(x), m_y(y), m_size(x.size()), m_interpType(interpType), m_binType(bin_type)
{
  shared_ptr<gsl_interp_accel> acc(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  m_acc = acc;
	  
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
    ErrorCBL("Error in the constructor of cosmobl::glob::FuncGrid::FuncGrid() in FuncGrid.cpp: the value of m_interpType is not permitted!");

  vector<double> xx, yy;

  if (m_binType==cosmobl::binType::_linear_) {
    xx = m_x;
    yy = m_y;
  }
  else if (m_binType==cosmobl::binType::_logarithmic_) {
    for (size_t i=0; i<m_size; i++) {
      xx.emplace_back(log10(m_x[i]));
      yy.emplace_back(log10(m_y[i]));
    }
  }
  else
    ErrorCBL("Error in the constructor of cosmobl::glob::FuncGrid::FuncGrid() in FuncGrid.cpp: the value of m_binType is not permitted!");

  m_xmin = Min(m_x);
  m_xmax = Max(m_x);
  
  std::shared_ptr<gsl_spline> spline(gsl_spline_alloc(m_type, m_size), gsl_spline_free);
  gsl_spline_init(spline.get(), xx.data(), yy.data(), m_size);
  m_spline = spline;
}


// =====================================================================================

				 
void cosmobl::glob::FuncGrid::free ()
{
  //gsl_spline_free(m_spline);
  //gsl_interp_accel_free(m_acc);
}


// =====================================================================================
     

double cosmobl::glob::FuncGrid::operator () (const double xx) const
{
  const double _xx = (m_binType==cosmobl::binType::_logarithmic_) ? log10(xx) : xx;

  double val = -1.;
  
  if (xx<m_xmin) { // perform a linear extrapolation
    val = m_spline.get()->y[0]+(_xx-m_spline.get()->x[0])/(m_spline.get()->x[1]-m_spline.get()->x[0])*(m_spline.get()->y[1]-m_spline.get()->y[0]);
    val = (m_binType==cosmobl::binType::_logarithmic_) ? pow(10, val) : val;
    if (val!=val) return ErrorCBL("Error in cosmobl::glob::FuncGrid::operator () of FuncGrid.cpp: inside the xx<m_xmin condition, the return value is nan!");
    else return val;
  }
  
  else if (xx>m_xmax) { // perform a linear extrapolation
    val = m_spline.get()->y[m_size-2]+(_xx-m_spline.get()->x[m_size-2])/(m_spline.get()->x[m_size-1]-m_spline.get()->x[m_size-2])*(m_spline.get()->y[m_size-1]-m_spline.get()->y[m_size-2]);
    val = (m_binType==cosmobl::binType::_logarithmic_) ? pow(10, val) : val;
    if (val!=val) return ErrorCBL("Error in cosmobl::glob::FuncGrid::operator () of FuncGrid.cpp: inside the xx>m_xmax condition, the return value is nan!");
    else return val;
  }
  
  // performe an interpolation
  else {
    val = (m_binType==cosmobl::binType::_logarithmic_) ? pow(10., gsl_spline_eval(m_spline.get(), _xx, m_acc.get())) : gsl_spline_eval(m_spline.get(), _xx, m_acc.get());
    if (val!=val) return ErrorCBL("Error in cosmobl::glob::FuncGrid::operator () of FuncGrid.cpp: the return value is nan!");
    else return val;
  }
}


// =====================================================================================


vector<double> cosmobl::glob::FuncGrid::eval_func (const vector<double> xx) const
{
  vector<double> yy;
  
  for (size_t i=0; i<xx.size(); i++)
    yy.push_back(this->operator()(xx[i]));

  return yy;
}


// =====================================================================================


double cosmobl::glob::FuncGrid::D1v (const double xx) const
{
  double D1 = gsl_spline_eval_deriv(m_spline.get(), xx, m_acc.get());

  return ((m_binType==cosmobl::binType::_logarithmic_) ? xx*D1/this->operator()(xx): D1);
}


// =====================================================================================


double cosmobl::glob::FuncGrid::D2v (const double xx) const
{
  return gsl_spline_eval_deriv2(m_spline.get(), xx, m_acc.get());
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


// =====================================================================================


cosmobl::glob::FuncGrid2D::FuncGrid2D (const vector<double> x, const vector<double> y, const vector<vector<double>> fxy, const string interpType)
{
  // set internal variables
  m_x = x;

  m_size_x = m_x.size();
  m_xmin = Min(m_x);
  m_xmax = Max(m_x);

  m_y = y;
  m_size_y = m_y.size();
  m_ymin = Min(m_y);
  m_ymax = Max(m_y);

  std::shared_ptr<double> grid( new double[m_size_x*m_size_y], []( double *p ) { delete[] p; });

  for (size_t i=0; i<m_size_x; i++)
    for (size_t j=0; j<m_size_y; j++)
      grid.get()[i+j*m_size_x] = fxy[i][j];

  m_fxy = grid;

  m_interpType = interpType;

  if (m_interpType=="Linear") 
    m_type = gsl_interp2d_bilinear;

  else if (m_interpType=="Cubic") 
    m_type = gsl_interp2d_bicubic;

  shared_ptr<gsl_interp_accel> acc_x(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  shared_ptr<gsl_interp_accel> acc_y(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  m_acc_x = acc_x;
  m_acc_y = acc_x;

  shared_ptr<gsl_interp2d> interp(gsl_interp2d_alloc(m_type, m_size_x, m_size_y), gsl_interp2d_free);
  gsl_interp2d_init(interp.get(), m_x.data(), m_y.data(), m_fxy.get(), m_size_x, m_size_y);
  m_interp = interp;
}


// =====================================================================================

				 
void cosmobl::glob::FuncGrid2D::free ()
{
  /*
  gsl_interp2d_free(m_interp);
  gsl_interp_accel_free(m_acc_x);
  gsl_interp_accel_free(m_acc_y);
  std::free(m_fxy);
  */
}


// =====================================================================================
     

double cosmobl::glob::FuncGrid2D::operator () (const double xx, const double yy) const
{
  bool extr = false;
  if ((xx>m_xmax || xx<m_xmin) ||(yy>m_ymax || yy<m_ymin))
    extr = true;

  double val;
  if (extr)
    val= gsl_interp2d_eval_extrap(m_interp.get(), m_x.data(), m_y.data() , m_fxy.get(), xx, yy, m_acc_x.get(), m_acc_y.get());
  else
    val= gsl_interp2d_eval(m_interp.get(), m_x.data(), m_y.data() , m_fxy.get(), xx, yy, m_acc_x.get(), m_acc_y.get());

  return val;

}


// =====================================================================================


vector<double> cosmobl::glob::FuncGrid2D::eval_func (const vector<vector<double>> xx) const
{
  vector<double> yy;
  
  for (size_t i=0; i<xx.size(); i++)
    yy.push_back(this->operator()(xx[i][0], xx[i][1]));

  return yy;
}

