/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file CosmoBolognaLib/Func/LegendrePolynomials.cpp
 *
 *  @brief Methods of the class LegendrePolynomials used to compute
 *  the Legendre Polynomials and their integrals
 *
 *  This file contains the implementation of the methods of the class
 *  LegendrePolynomials used to measure the Legendre polynomials at a
 *  given value
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@uniroma3.it
 */

#include "LegendrePolynomials.h"

#include <boost/math/special_functions/beta.hpp>

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::glob::LegendrePolynomials::m_set_coefficients (const int lMax)
{

  auto old_handler = gsl_set_error_handler_off();
  (void)old_handler;

  m_nOrders = lMax+1;

  m_coefficients.resize(m_nOrders, m_nOrders);

  auto binomial_coeff  = [&] (const double n, const double k) {

    gsl_sf_result lng_n, lng_k, lng_nmk;  
    double lng_n_sgn, lng_k_sgn, lng_nmk_sgn;

    int status_n   = gsl_sf_lngamma_sgn_e(n+1, &lng_n, &lng_n_sgn);
    int status_k   = gsl_sf_lngamma_sgn_e(k+1, &lng_k, &lng_k_sgn);
    int status_nmk = gsl_sf_lngamma_sgn_e(n-k+1, &lng_nmk, &lng_nmk_sgn);

    if(status_n != 0 or status_k != 0 or status_nmk != 0) {
      return 0.;
    }

    return lng_n_sgn/(lng_k_sgn*lng_nmk_sgn)*exp(lng_n.val-lng_k.val-lng_nmk.val);

  };

  for (int i=0; i<static_cast<int>(m_nOrders); i++) {
    double fact = pow(2, i);
    for (int j=0; j<static_cast<int>(m_nOrders); j++) 
      m_coefficients(i,j) = fact*binomial_coeff(i, j)*binomial_coeff( double(i+j-1)*0.5, i);
  }

  //gsl_set_error_handler(old_handler);
}


// ============================================================================================


cbl::glob::LegendrePolynomials::LegendrePolynomials (const int lMax, const bool safe)
{
  set(lMax);
  (void)safe;
}


// ============================================================================================


void cbl::glob::LegendrePolynomials::set (const int lMax)
{
  m_set_coefficients (lMax);
}


// ============================================================================================


double cbl::glob::LegendrePolynomials::operator() (const double x, const int ell)
{
  Eigen::VectorXd m_powers(m_nOrders);

  m_powers[0] = 1.;
  m_powers[1] = x;

  for (size_t i=2; i<m_nOrders; i++)
    m_powers[i] = m_powers[i-1]*x;
  
  return m_coefficients.row(ell).dot(m_powers);
}


// ============================================================================================


vector<double> cbl::glob::LegendrePolynomials::operator() (const double x)
{
  Eigen::VectorXd m_powers(m_nOrders);

  m_powers[0] = 1.;
  m_powers[1] = x;

  for (size_t i=2; i<m_nOrders; i++)
    m_powers[i] = m_powers[i-1]*x;

  Eigen::VectorXd leg_pols = m_coefficients*m_powers;

  vector<double> vv(leg_pols.data(), leg_pols.data()+leg_pols.size());
  return vv;
}


// ============================================================================================


double cbl::glob::LegendrePolynomials::integral (const double x_min, const double x_max, const int ell)
{
  Eigen::VectorXd min(m_nOrders);
  Eigen::VectorXd max(m_nOrders);
  Eigen::VectorXd norm(m_nOrders);

  min[0] = x_min;
  max[0] = x_max;
  norm[0] = 1.;

  for (size_t i=1; i<m_nOrders; i++) {
    min[i] = min[i-1]*x_min;
    max[i] = max[i-1]*x_max;
    norm[i] = norm[i-1]+1;
  }

  return m_coefficients.row(ell).dot( (max-min).cwiseQuotient(norm) )/(x_max-x_min);
}


// ============================================================================================


vector<double> cbl::glob::LegendrePolynomials::integral (const double x_min, const double x_max)
{
  Eigen::VectorXd min(m_nOrders);
  Eigen::VectorXd max(m_nOrders);
  Eigen::VectorXd norm(m_nOrders);

  min[0] = x_min;
  max[0] = x_max;
  norm[0] = 1.;

  for (size_t i=1; i<m_nOrders; i++) {
    min[i] = min[i-1]*x_min;
    max[i] = max[i-1]*x_max;
    norm[i] = norm[i-1]+1;
  }

  Eigen::VectorXd leg_pols = m_coefficients*( (max-min).cwiseQuotient(norm) )/(x_max-x_min);

  vector<double> vv(leg_pols.data(), leg_pols.data()+leg_pols.size());
  return vv;
}


// ============================================================================================

vector<vector<double>> cbl::glob::LegendrePolynomials::triangle (const double r12, const double r13, const double r23)
{
  double mu_23 = (r12*r12+r13*r13-r23*r23)/(2*r12*r13);
  double mu_13 = (r12*r12+r23*r23-r13*r13)/(2*r12*r23);
  double mu_12 = (r23*r23+r13*r13-r12*r12)/(2*r23*r13);

  vector<vector<double>> leg_pols;
  leg_pols.push_back(this->operator()(mu_23));
  leg_pols.push_back(this->operator()(mu_13));
  leg_pols.push_back(this->operator()(mu_12));

  return leg_pols;
}


// ============================================================================================


vector<double> cbl::glob::LegendrePolynomials::triangle_integral (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const double r23_min, const double r23_max, const double rel_err, const int nevals)
{
  double norm = 2*(pow(r12_max, 3)-pow(r12_min,3))*(pow(r13_max, 3)-pow(r13_min, 3))/9;

  Eigen::VectorXd m_powers(m_nOrders);

  for (int i=0; i<static_cast<int>(m_nOrders); i++) {

    size_t exponent = i+1;

    auto r12_integrand = [&] (const double r12) {
      auto r13_integrand = [&] (const double r13) {
	double x_min = (r12*r12+r13*r13-r23_max*r23_max)/(2*r12*r13);
	double x_max = (r12*r12+r13*r13-r23_min*r23_min)/(2*r12*r13);

	if ( x_min>1 or x_max < -1) 
	  return 0.;

	//double low = pow(max(x_min, -1.), exponent);
	//double up = pow(min(1., x_max), exponent);
	double low = gsl_pow_uint(max(x_min, -1.), exponent);
	double up = gsl_pow_uint(min(x_max, 1.), exponent);
	return r13*r13*(up-low)/exponent;
      };
      return r12*r12*cbl::wrapper::gsl::GSL_integrate_cquad(r13_integrand, r13_min, r13_max, rel_err, 0., nevals);
    };

    m_powers[i] = cbl::wrapper::gsl::GSL_integrate_cquad(r12_integrand, r12_min, r12_max, rel_err, 0., nevals)/norm;
  }

  Eigen::VectorXd leg_pols = m_coefficients*m_powers;

  vector<double> vv(leg_pols.data(), leg_pols.data()+leg_pols.size());
  return vv;
}
