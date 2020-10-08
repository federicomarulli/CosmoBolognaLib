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
 *  @file FuncGrid/FuncGrid.cpp
 *
 *  @brief Methods of the class FuncGrid
 *
 *  This file contains the implementation of the methods of the class
 *  FuncGrid and FuncGrid_2D
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "FuncGrid.h"

using namespace std;

using namespace cbl;
using namespace glob;


// ============================================================================


cbl::glob::FuncGrid::FuncGrid (const std::vector<double> x, const std::vector<double> y, const std::string interpType, const cbl::BinType bin_type)
  : m_x(x), m_y(y), m_size(x.size()), m_interpType(interpType), m_binType(bin_type)
{
  if (!is_sorted(x.begin(), x.end()))
    ErrorCBL("the x array is not sorted!", "FuncGrid", "FuncGrid.cpp");

  shared_ptr<gsl_interp_accel> acc(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  m_acc = acc;
	  
  if (m_size<5 && interpType!="Linear") {
    WarningMsgCBL("the array size is less than 5 -> setting interpolation method to Linear!", "FuncGrid", "FuncGrid.cpp");
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
    ErrorCBL("the value of m_interpType is not permitted!", "FuncGrid", "FuncGrid.cpp");

  vector<double> xx, yy;

  if (m_binType==cbl::BinType::_linear_) {
    xx = m_x;
    yy = m_y;
  }
  else if (m_binType==cbl::BinType::_logarithmic_) {
    for (size_t i=0; i<m_size; i++) {
      xx.emplace_back(log10(m_x[i]));
      yy.emplace_back(log10(m_y[i]));
    }
  }
  else
    ErrorCBL("the value of m_binType is not permitted!", "FuncGrid", "FuncGrid.cpp");

  m_xmin = Min(m_x);
  m_xmax = Max(m_x);
  
  std::shared_ptr<gsl_spline> spline(gsl_spline_alloc(m_type, m_size), gsl_spline_free);
  gsl_spline_init(spline.get(), xx.data(), yy.data(), m_size);
  m_spline = spline;
}


// =====================================================================================

				 
void cbl::glob::FuncGrid::free ()
{
  //gsl_spline_free(m_spline);
  //gsl_interp_accel_free(m_acc);
}


// =====================================================================================
     

double cbl::glob::FuncGrid::operator () (const double xx) const
{
  const double _xx = (m_binType==cbl::BinType::_logarithmic_) ? log10(xx) : xx;

  double val = -1.;
  
  if (xx<m_xmin) { // perform a linear extrapolation
    val = m_spline.get()->y[0]+(_xx-m_spline.get()->x[0])/(m_spline.get()->x[1]-m_spline.get()->x[0])*(m_spline.get()->y[1]-m_spline.get()->y[0]);
    val = (m_binType==cbl::BinType::_logarithmic_) ? pow(10, val) : val;
    if (val!=val) return ErrorCBL("inside the xx<m_xmin condition, the return value is nan!", "operator ()", "FuncGrid.cpp");
    else return val;
  }
  
  else if (xx>m_xmax) { // perform a linear extrapolation
    val = m_spline.get()->y[m_size-2]+(_xx-m_spline.get()->x[m_size-2])/(m_spline.get()->x[m_size-1]-m_spline.get()->x[m_size-2])*(m_spline.get()->y[m_size-1]-m_spline.get()->y[m_size-2]);
    val = (m_binType==cbl::BinType::_logarithmic_) ? pow(10, val) : val;
    if (val!=val) return ErrorCBL("inside the xx>m_xmax condition, the return value is nan!", "operator ()", "FuncGrid.cpp");
    else return val;
  }
  
  // performe an interpolation
  else {
    val = (m_binType==cbl::BinType::_logarithmic_) ? pow(10., gsl_spline_eval(m_spline.get(), _xx, m_acc.get())) : gsl_spline_eval(m_spline.get(), _xx, m_acc.get());
    if (val!=val) coutCBL << "xx = " << conv(_xx, par::fDP3) << endl;
    if (val!=val) return ErrorCBL("the return value is nan!", "operator ()", "FuncGrid.cpp");
    else return val;
  }
}


// =====================================================================================


std::vector<double> cbl::glob::FuncGrid::eval_func (const std::vector<double> xx) const
{
  vector<double> yy;
  
  for (size_t i=0; i<xx.size(); i++)
    yy.push_back(this->operator()(xx[i]));

  return yy;
}


// =====================================================================================


double cbl::glob::FuncGrid::D1v (const double xx) const
{
  double D1 = gsl_spline_eval_deriv(m_spline.get(), xx, m_acc.get());

  return ((m_binType==cbl::BinType::_logarithmic_) ? xx*D1/this->operator()(xx): D1);
}


// =====================================================================================


double cbl::glob::FuncGrid::D2v (const double xx) const
{
  return gsl_spline_eval_deriv2(m_spline.get(), xx, m_acc.get());
}


// =====================================================================================


double cbl::glob::FuncGrid::integrate_qag (const double a, const double b, const double rel_err, const double abs_err, const int limit_size, const int rule)
{
  function<double(double)> f = bind(&FuncGrid::operator(), this, std::placeholders::_1);

  return wrapper::gsl::GSL_integrate_qag(f, a, b, rel_err, abs_err, limit_size, rule);
}


// =====================================================================================


double cbl::glob::FuncGrid::integrate_qaws (const double a, const double b, const double alpha, const double beta, const int mu, const int nu, const double rel_err, const double abs_err, const int limit_size)
{
  function<double(double)> f = bind(&FuncGrid::operator(), this, std::placeholders::_1);

  return wrapper::gsl::GSL_integrate_qaws(f, a, b, alpha, beta, mu, nu, rel_err, abs_err, limit_size);
}


// =====================================================================================


double cbl::glob::FuncGrid::root (const double x_low, const double x_up, const double fx0, const double rel_err, const double abs_err)
{
  function<double(double)> f = bind(&FuncGrid::operator(), this, std::placeholders::_1);

  return wrapper::gsl::GSL_root_brent(f, fx0,  x_low, x_up, rel_err, abs_err);
}


// =====================================================================================


double cbl::glob::FuncGrid::root_D1v (const double x_low, const double x_up, const double fx0, const double rel_err, const double abs_err)
{
  function<double(double)> f = bind(&FuncGrid::D1v, this, std::placeholders::_1);

  return wrapper::gsl::GSL_root_brent(f, fx0,  x_low, x_up, rel_err, abs_err);
}


// =====================================================================================


double cbl::glob::FuncGrid::root_D2v (const double x_low, const double x_up, const double fx0, const double rel_err, const double abs_err)
{
  function<double(double)> f = bind(&FuncGrid::D2v, this, std::placeholders::_1);

  return wrapper::gsl::GSL_root_brent(f, fx0,  x_low, x_up, rel_err, abs_err);
}


// =====================================================================================


cbl::glob::FuncGrid2D::FuncGrid2D (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> fxy, const std::string interpType)
{
  if(!is_sorted(x.begin(), x.end()))
    ErrorCBL("the x array is not sorted!", "FuncGrid2D", "FuncGrid.cpp");
  if(!is_sorted(y.begin(), y.end()))
    ErrorCBL("the y array is not sorted!", "FuncGrid2D", "FuncGrid.cpp");

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

  else 
    ErrorCBL("the value of m_interpType is not permitted!", "FuncGrid2D", "FuncGrid.cpp");

  shared_ptr<gsl_interp_accel> acc_x(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  shared_ptr<gsl_interp_accel> acc_y(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  m_acc_x = acc_x;
  m_acc_y = acc_y;

  shared_ptr<gsl_spline2d> spline(gsl_spline2d_alloc(m_type, m_size_x, m_size_y), gsl_spline2d_free);
  gsl_spline2d_init(spline.get(), m_x.data(), m_y.data(), m_fxy.get(), m_size_x, m_size_y);
  m_spline = spline;
}


// =====================================================================================

				 
void cbl::glob::FuncGrid2D::free ()
{
  /*
  gsl_interp2d_free(m_interp);
  gsl_interp_accel_free(m_acc_x);
  gsl_interp_accel_free(m_acc_y);
  std::free(m_fxy);
  */
}


// =====================================================================================
     

double cbl::glob::FuncGrid2D::operator () (const double xx, const double yy) const
{
  bool extr = false;
  if ((xx>m_xmax || xx<m_xmin) ||(yy>m_ymax || yy<m_ymin))
    extr = true;

  double val;
  
  if (extr)
    ErrorCBL("Point ("+conv(xx, par::fDP2)+", "+conv(yy, par::fDP2)+") is outside the interpolation range...", "operator ()", "FuncGrid.cpp", glob::ExitCode::_workInProgress_);
  
  else
    val = gsl_spline2d_eval(m_spline.get(), xx, yy, m_acc_x.get(), m_acc_y.get());

  return val;

}


// =====================================================================================


std::vector<double> cbl::glob::FuncGrid2D::eval_func (const std::vector<std::vector<double>> xx) const
{
  vector<double> yy;
  
  for (size_t i=0; i<xx.size(); i++)
    yy.push_back(this->operator()(xx[i][0], xx[i][1]));

  return yy;
}


// =====================================================================================


double cbl::glob::FuncGrid2D::IntegrateVegas (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  auto integrand = [&] (const vector<double> var)
  {
    return this->operator()(var[0], var[1]);
  };

  wrapper::cuba::CUBAwrapper cuba_integral(integrand, 2);

  return cuba_integral.IntegrateVegas({{xmin, xmax}, {ymin, ymax}});
}


// =====================================================================================


double cbl::glob::FuncGrid2D::IntegrateSuave (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  auto integrand = [&] (const vector<double> var)
  {
    return this->operator()(var[0], var[1]);
  };

  wrapper::cuba::CUBAwrapper cuba_integral(integrand, 2);

  return cuba_integral.IntegrateSuave({{xmin, xmax}, {ymin, ymax}});
}


// =====================================================================================


double cbl::glob::FuncGrid2D::IntegrateDivonne (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  auto integrand = [&] (const vector<double> var)
  {
    return this->operator()(var[0], var[1]);
  };

  wrapper::cuba::CUBAwrapper cuba_integral(integrand, 2);

  return cuba_integral.IntegrateDivonne({{xmin, xmax}, {ymin, ymax}});
}


// =====================================================================================


double cbl::glob::FuncGrid2D::IntegrateCuhre (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  auto integrand = [&] (const vector<double> var)
  {
    return this->operator()(var[0], var[1]);
  };

  wrapper::cuba::CUBAwrapper cuba_integral(integrand, 2);

  return cuba_integral.IntegrateCuhre ({{xmin, xmax}, {ymin, ymax}});
}
